import pandas as pd
from enum import Enum
import yaml
import numpy as np
import os
from itertools import product
from HTPolyNet.topocoord import TopoCoord, BTRC
from HTPolyNet.plot import trace
from HTPolyNet.gromacs import gromacs_distance, mdp_modify
from HTPolyNet.configuration import ReactionList
from HTPolyNet.molecule import MoleculeDict
from HTPolyNet.reaction import reaction_stage
from multiprocessing import Pool
from functools import partial
import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.stringthings import my_logger
import HTPolyNet.checkpoint as cp

import logging

logger=logging.getLogger(__name__)

class state(Enum):
    """Enumerated CURE state
    """
    cure_bondsearch=0
    cure_drag=1
    cure_update=2
    cure_relax=3
    cure_equilibrate=4
    cap_bondsearch=5
    cap_update=6  # dragging is not needed/allowed in caping since these are likely intramolecular bonds...
    cap_relax=7
    cap_equilibrate=8
    finished=9
    unknown=99
    def __str__(self):
        return self.name
    def basename(self):
        name=str(self)
        if 'cap_' in name:
            return name[len('cap_'):]
        if 'cure_' in name:
            return name[len('cure_'):]

class CureController:
    default_equilibration_sequence = [ { 'ensemble': 'min' }, 
        { 'ensemble': 'nvt', 'temperature': 600, 'nsteps': 1000 },
        { 'ensemble': 'npt', 'temperature': 600, 'pressure': 1, 'nsteps': 2000 }
    ]
    curedict_defaults={
        'output': {
            'bonds_file': 'bonds.csv'
        },
        'controls': {
            'max_conversion_per_iteration': 1.0,
            'search_radius': 0.5,
            'radial_increment': 0.05,
            'late_threshold': 0.85,
            'max_iterations': 100,
            'desired_conversion': 0.5,
            'ncpu' : os.cpu_count()
        },
        'drag': {
            'limit': 0.0,
            'kb': 300000.0, # kJ/mol/nm
            'trigger_distance': 0.0,
            'nstages': 0,
            'increment': 0.0,
            'cutoff_pad': 0.2,
            'equilibration': default_equilibration_sequence
        },
        'relax': {
            'nstages': 6,
            'increment': 0.0,
            'cutoff_pad': 0.2,
            'equilibration': default_equilibration_sequence
        },
        'equilibrate': {
            'temperature': 300,
            'pressure': 1,
            'nsteps': 50000,
            'ensemble': 'npt'
        },
        'gromacs': {
            'rdefault': 0.9
        }
    }

    def __init__(self,curedict={}):
        self.iter=0
        self.max_nxlinkbonds=0
        self.cum_nxlinkbonds=0
        self.desired_nxlinkbonds=0
        self.max_search_radius=0.0
        self.max_radidx=0
        self.state=state.unknown
        self.curr_nxlinkbonds=0

        self.current_stage={}
        self.current_stage['drag']=0
        self.current_stage['relax']=0

        self.current_radidx=0
        self.current_radius=0.0

        self.bonds_df:pd.DataFrame=None
        self.bonds_are='nonexistent!'

        self.search_failed=False
        self.dicts={}
        for k,v in self.curedict_defaults.items():
            self.dicts[k]=curedict.get(k,v)  # loads entire default if dict is just missing
            if type(v)==dict: # assign subitem defaults if their keys are missing
                for kk,vv in v.items():
                    if not kk in self.dicts[k]:
                        self.dicts[k][kk]=vv

        self.dragging_enabled=False
        d=self.dicts['drag']
        if (d['nstages']>0 or d['increment']>0.0) and d['limit']>0.0:
            self.dragging_enabled=True

    def _to_yaml(self,filename='cure_controller_state.yaml'):
        with open(filename,'w') as f:
            f.write(yaml.dump(self))

    @classmethod
    def from_yaml(cls,filename='cure_controller_state.yaml'):
        with open(filename,'r') as f:
            yaml_string=f.read()
        return yaml.load(yaml_string,Loader=yaml.Loader)

    def setup(self,max_nxlinkbonds=0,desired_nxlinkbonds=0,max_search_radius=0.0):
        self.max_nxlinkbonds=max_nxlinkbonds
        self.desired_nxlinkbonds=desired_nxlinkbonds
        self.max_search_radius=max_search_radius
        d=self.dicts['controls']
        self.max_radidx=int((self.max_search_radius-d['search_radius'])/d['radial_increment'])

    @cp.enableCheckpoint
    def do_iter(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,gromacs_dict={}):
        my_logger(f'Iteration {self.iter} begins',logger.info,fill='~')
        TC.grab_files() # copy files locally
        self._do_bondsearch(TC,RL,MD)
        self._do_preupdate_dragging(TC,gromacs_dict)
        self._do_topology_update(TC,MD)
        self._do_relax(TC,gromacs_dict)
        self._do_equilibrate(TC,gromacs_dict)
        logger.info(f'Iteration {self.iter} current conversion {self._curr_conversion():.3f} or {self.cum_nxlinkbonds} bonds')
        self._to_yaml()
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    @cp.enableCheckpoint
    def do_capping(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,gromacs_dict={}):
        self._do_cap_bondsearch(TC,RL,MD)
        self._do_topology_update(TC,MD)
        self._do_relax(TC)
        self._do_equilibrate(TC,gromacs_dict)
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    def is_cured(self):
        logger.debug(f'search_failed? {self.search_failed}')
        logger.debug(f'cumxlinks {self.cum_nxlinkbonds} maxxlinks {self.max_nxlinkbonds}: {(self.cum_nxlinkbonds>=self.max_nxlinkbonds)}')
        finished=self.search_failed or (self.cum_nxlinkbonds>=self.desired_nxlinkbonds)
        if finished:
            self.state=state.cap_bondsearch
        return finished

    def _curr_conversion(self):
        if not self.max_nxlinkbonds: return 0
        return float(self.cum_nxlinkbonds)/self.max_nxlinkbonds

    def reset(self):
        self.iter=1
        self.state=state.cure_bondsearch
        self.curr_nxlinkbonds=0
        self.current_stage['drag']=0
        self.current_stage['relax']=0
        self.current_radidx=0
        self.current_radius=0.0
        self.bonds_df=pd.DataFrame()
        self.bonds_are='nonexistent!'
        self.search_failed=False

    def next_iter(self):
        i=self.iter+1
        self.reset()
        self.iter=i
        return self.iter>=self.dicts['controls']['max_iterations']

    def _do_bondsearch(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,reentry=False):
        if self.state!=state.cure_bondsearch: return
        opfx=self._pfx()
        d=self.dicts['controls']
        self.current_radius=d['search_radius']+self.current_radidx*d['radial_increment']
        logger.info(f'Bond search using radius {self.current_radius} nm initiated')
        curr_conversion=self._curr_conversion()
        apply_probabilities=curr_conversion<d['late_threshold']
        bond_limit=int(d['max_conversion_per_iteration']*self.max_nxlinkbonds)
        bond_target=int((d['desired_conversion']-curr_conversion)*self.max_nxlinkbonds)
        bond_limit=min([bond_limit,bond_target])
        logger.debug(f'Iteration limited to at most {bond_limit} new bonds')
        # TODO: new branch: multisearch:
        # calc and store distances once outside this loop
        # have a controller that runs loop until a select number
        # of bonds has been formed, or we've maxed out the radius
        # raw_bdf=self.make_candidates(TC,RL,MD,...)
        nbonds=0
        while nbonds==0 and self.current_radidx<self.max_radidx:
            # test_bdf=raw_bdf[raw_bdf['r']<self.current_radius]
            # result_bdf=self.apply_filters(TC,RL,MD,...)
            nbdf=self._searchbonds(TC,RL,MD,stage=reaction_stage.cure,abs_max=bond_limit,apply_probabilities=apply_probabilities,reentry=reentry)
            nbonds=nbdf.shape[0]
            # nbonds+=result_bdf.shape[0]
            # if nbonds<PARAMETER:
            if nbonds==0:
                self.current_radidx+=1
                self.current_radius+=d['radial_increment']
                logger.info(f'Radius increased to {self.current_radius} nm')
        if nbonds>0:
            ess='' if nbonds==1 else 's'
            logger.info(f'Iteration {self.iter} will generate {nbdf.shape[0]} new bond{ess}')
            pairs=pd.DataFrame() # empty placeholder
            TC.add_length_attribute(nbdf,attr_name='initial_distance')
            self._register_bonds(nbdf,pairs,f'{opfx}-bonds.csv',bonds_are='identified')
            self.state=state.cure_drag if self.dragging_enabled else state.cure_update
        else:
            self.search_failed=True
            self.state=state.cap_bondsearch # proceed to cap
        self.cum_nxlinkbonds+=nbonds
        self._to_yaml()
        logger.debug(f'next: {self.state}')

    def _do_preupdate_dragging(self,TC:TopoCoord,gromacs_dict={}):
        if self.state!=state.cure_drag: return
        nbdf=self.bonds_df
        assert nbdf.shape[0]>0
        d=self.dicts['drag']
        nogos=[not self.dragging_enabled,nbdf['initial_distance'].max()<d['trigger_distance']]
        logger.debug(f'{nogos} {any(nogos)}')
        if any(nogos):
            logger.debug(f'No dragging is necessary')
        else:
            self._distance_attenuation(TC,mode='drag',gromacs_dict=gromacs_dict)
        self.state=state.cure_update
        self._to_yaml()
        logger.debug(f'next: {self.state}')

    def _do_relax(self,TC:TopoCoord,gromacs_dict={}):
        if self.state!=state.cure_relax and self.state!=state.cap_relax: return
        nbdf=self.bonds_df
        if nbdf.shape[0]==0: # no bonds identified
            if self.state==state.cure_relax:
                self.state=state.cap_bondsearch # drop out of CURE loop
            else:
                self.state=state.finished
            return
        self._distance_attenuation(TC,mode='relax',gromacs_dict=gromacs_dict)
        if self.state==state.cure_relax:
            self.state=state.cure_equilibrate
        else:
            self.state=state.cap_equilibrate
        self._to_yaml()
        logger.debug(f'next: {self.state}')

    def _do_topology_update(self,TC:TopoCoord,MD:MoleculeDict):
        if self.state!=state.cure_update and self.state!=state.cap_update: return
        opfx=self._pfx()
        logger.debug(f'Topology update')
        bonds_df,pairs_df=TC.update_topology_and_coordinates(self.bonds_df,template_dict=MD,write_mapper_to=f'{opfx}-idx-mapper.csv')
        TC.add_length_attribute(bonds_df,attr_name='initial_distance')
        TC.add_length_attribute(pairs_df,attr_name='initial_distance')
        self._register_bonds(bonds_df,pairs_df,f'{opfx}-bonds.csv',bonds_are='unrelaxed')
        TC.write_gro(f'{opfx}.gro')
        TC.write_top(f'{opfx}.top')
        TC.write_grx_attributes(f'{opfx}.grx')
        if self.state==state.cure_update:
            self.state=state.cure_relax
        else:
            self.state=state.cap_relax
        self._to_yaml()

    def _do_equilibrate(self,TC:TopoCoord,gromacs_dict={}):
        if self.state!=state.cure_equilibrate and self.state!=state.cap_equilibrate: return
        d=self.dicts['equilibrate']
        opfx=self._pfx()
        TC.equilibrate(deffnm=f'{opfx}',edict=d,gromacs_dict=gromacs_dict,plot_pfx=f'iter-{self.iter}-{str(self.state)}')
        if self.state==state.cure_equilibrate:
            # go to next iteration -- this whole method is skipped if nbonds==0 in relax
            self.state=state.cure_bondsearch # if not self.search_failed else state.cap_bondsearch
        elif self.state==state.cap_equilibrate:
            self.state=state.finished
        self._to_yaml()

    def _do_cap_bondsearch(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict):
        if self.state!=state.cap_bondsearch: return
        opfx=self._pfx()
        # multi: nbdf=self.make_cadidates(...,stage='cap')
        nbdf=self._searchbonds(TC,RL,MD,stage=reaction_stage.cap)
        nbonds=nbdf.shape[0]
        ess='' if nbonds==1 else 's'
        logger.info(f'Capping will generate {nbdf.shape[0]} new bond{ess}')
        if nbonds>0:
            cwd=pfs.go_to(f'systems/capping')
            pairs=pd.DataFrame() # empty placeholder
            TC.add_length_attribute(nbdf,attr_name='initial_distance')
            self._register_bonds(nbdf,pairs,f'{opfx}-bonds.csv',bonds_are='identified')
            self.state=state.cap_update
        else:
            self.state=state.finished
        self._to_yaml()

    def _write_bonds_df(self,bondsfile='bonds.csv'):
        self.bonds_df.to_csv(bondsfile,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def _read_bonds_df(self,bonds_file_override=''):
        infile=self.dicts['output']['bonds_file'] if not bonds_file_override else bonds_file_override
        assert os.path.exists(infile),f'Error: {infile} not found'
        self.bonds_df=pd.read_csv(infile,sep='\s+',header=0)
        self.dicts['output']['bonds_file']=os.path.abspath(infile)

    def _register_bonds(self,bonds,pairs,bonds_file,bonds_are='unrelaxed'):
        self.bonds_df=bonds
        self.pairs_df=pairs
        self.bonds_are=bonds_are
        self._write_bonds_df(bonds_file)
        self.dicts['output']['bonds_file']=os.path.abspath(bonds_file)

    def _pfx(self):
        return f'{self.state.value}-{self.state}'

    def _distance_attenuation(self,TC:TopoCoord,mode='drag',gromacs_dict={}):
        assert mode in ['drag','relax']
        statename=self.state.basename()
        opfx=self._pfx()
        nbdf=self.bonds_df.copy()
        pdf=self.pairs_df.copy()
        d=self.dicts[mode]
        nbdf['current_lengths']=nbdf['initial_distance'].copy()
        maxL,minL,meanL=nbdf['current_lengths'].max(),nbdf['current_lengths'].min(),nbdf['current_lengths'].mean()
        logger.debug(f'Lengths avg/min/max: {meanL:.3f}/{minL:.3f}/{maxL:.3f}')
        ess='' if nbdf.shape[0]==1 else 's'
        logger.info(f'{str(self.state).capitalize()} initiated on {nbdf.shape[0]} distance{ess} (max {maxL:.3f} nm)')
        roptions=[self.dicts['gromacs']['rdefault'],maxL]
        if mode=='drag':
            TC.add_restraints(nbdf,typ=6)
            logger.info('     Stage  Max-distance (nm)')
        else:
            pdf['current_lengths']=pdf['initial_distance'].copy()
            pmaxL,pminL,pmeanL=pdf['current_lengths'].max(),pdf['current_lengths'].min(),pdf['current_lengths'].mean()
            logger.debug(f'1-4 distances lengths avg/min/max: {pmeanL:.3f}/{pminL:.3f}/{pmaxL:.3f}')
            roptions.append(pmaxL)
            logger.info('     Stage  Max-distance (nm)  Max-1-4-distance (nm)')
        rcommon=max(roptions)
        for stg_dict in d['equilibration']:
            ensemble=stg_dict['ensemble']
            impfx=f'{statename}-{ensemble}' # e.g., drag-min, drag-nvt, drag-npt
            pfs.checkout(f'mdp/{impfx}.mdp')
            mdp_mods_dict={'rvdw':rcommon,'rcoulomb':rcommon,'rlist':rcommon}
            if ensemble!='min':
                mdp_mods_dict.update({'gen-temp':stg_dict['temperature'],'ref_t':stg_dict['temperature'],'gen-vel':'yes','nsteps':stg_dict['nsteps']})
            if ensemble=='npt':
                mdp_mods_dict['ref_p']=stg_dict['pressure']
            mdp_modify(f'{impfx}.mdp',mdp_mods_dict,add_if_missing=(ensemble!='min'))
        this_nstages=int(maxL/d['increment'])
        this_firststage=self.current_stage[mode]
        logger.debug(f'{self.state} {this_nstages} {this_firststage}')
        saveT=TC.copy_bond_parameters(self.bonds_df)
        for i in range(this_firststage,this_nstages):
            self.current_stage[mode]=i
            if mode=='drag':
                TC.attenuate_bond_parameters(self.bonds_df,i,this_nstages,minimum_distance=d['limit'],init_colname='initial_distance')
            else:
                TC.attenuate_bond_parameters(self.bonds_df,i,this_nstages,init_colname='initial_distance')
            stagepfx=f'{opfx}-stage-{i+1}'
            TC.write_top(f'{stagepfx}.top')
            for stg_dict in d['equilibration']:
                ensemble=stg_dict['ensemble']
                impfx=f'{statename}-{ensemble}' # e.g., drag-min, drag-nvt, drag-npt
                TC.grompp_and_mdrun(out=f'{stagepfx}-{ensemble}',mdp=f'{impfx}',**gromacs_dict)
                # logger.debug(f'{TC.files["gro"]}')
            TC.restore_bond_parameters(saveT)
            TC.add_length_attribute(nbdf,attr_name='current_lengths')
            maxL,minL,meanL=nbdf['current_lengths'].max(),nbdf['current_lengths'].min(),nbdf['current_lengths'].mean()
            logger.debug(f'Distances avg/min/max: {meanL:.3f}/{minL:.3f}/{maxL:.3f}')
            if mode=='drag':
                logger.info(f'{i+1:>10d}  {maxL:>17.3f}')
            else:
                TC.add_length_attribute(pdf,attr_name='current_lengths')
                pmaxL,pminL,pmeanL=pdf['current_lengths'].max(),pdf['current_lengths'].min(),pdf['current_lengths'].mean()
                logger.debug(f'1-4 distances lengths avg/min/max: {pmeanL:.3f}/{pminL:.3f}/{pmaxL:.3f}')
                logger.info(f'{i+1:>10d}  {maxL:>17.3f}  {pmaxL:>21.3f}')
            self._to_yaml()
        if mode=='drag':
            TC.remove_restraints(self.bonds_df)
        TC.write_top(f'{opfx}-complete.top')
        self._register_bonds(nbdf,pdf,f'{opfx}-{mode}-bonds.csv',bonds_are=('relaxed' if mode=='relax' else 'dragged'))
        self.current_stage[mode]=0
        self._to_yaml()

    # TODO: move to searchbonds.by; split into make_candidates() and apply_filters()
    def _searchbonds(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,stage=reaction_stage.cure,abs_max=0,apply_probabilities=False,reentry=False):

        adf=TC.gro_DataFrame('atoms')
        gro=TC.files['gro']
        ncpu=self.dicts['controls']['ncpu']
        if stage==reaction_stage.cure:
            TC.linkcell_initialize(self.current_radius,ncpu=ncpu,force_repopulate=reentry)
        raset=adf[adf['z']>0]  # this view will be used for downselecting to potential A-B partners
        bdf=pd.DataFrame()
        Rlist=[x for x in RL if (x.stage==stage and x.probability>0.0)]
        logger.debug(f'reactioncount {len(Rlist)} atomscount {raset.shape[0]}')
        for R in Rlist:
            logger.debug(f'Reaction {R.name} with {len(R.bonds)} bond(s)')
            prob=R.probability
            for bond in R.bonds:
                A=R.atoms[bond['atoms'][0]]
                B=R.atoms[bond['atoms'][1]]
                order=bond['order']
                aname=A['atom']
                areactantname_template=R.reactants[A['reactant']]
                aresid_template=A['resid']
                aresname=MD[areactantname_template].get_resname(aresid_template)
                az=A['z']
                bname=B['atom']
                breactantname_template=R.reactants[B['reactant']]
                bresid_template=B['resid']
                bresname=MD[breactantname_template].get_resname(bresid_template)
                bz=B['z']
                if stage==reaction_stage.cap:
                    assert areactantname_template==breactantname_template,f'Error: capping reaction {R.name} lists a bond whose atoms are in different reactants'
                    assert aresname==bresname,f'Error: capping reaction {R.name} lists a bond whose atoms are in different residues'

                Aset=raset[(raset['atomName']==aname)&(raset['resName']==aresname)&(raset['z']==az)&(raset['reactantName']==areactantname_template)]
                Bset=raset[(raset['atomName']==bname)&(raset['resName']==bresname)&(raset['z']==bz)&(raset['reactantName']==breactantname_template)]
                logger.debug(f'Aset {Aset.shape[0]} atoms')
                logger.debug(f'Bset {Bset.shape[0]} atoms')
                alist=list(zip(Aset['globalIdx'].to_list(),Aset['resNum'].to_list(),Aset['molecule'].to_list()))
                blist=list(zip(Bset['globalIdx'].to_list(),Bset['resNum'].to_list(),Bset['molecule'].to_list()))
                all_possible_pairs=list(product(alist,blist))

                idf=pd.DataFrame({'ai':           [int(x[0][0]) for x in all_possible_pairs],
                                  'ri':           [int(x[0][1]) for x in all_possible_pairs],
                                  'mi':           [int(x[0][2]) for x in all_possible_pairs],
                                  'aj':           [int(x[1][0]) for x in all_possible_pairs],
                                  'rj':           [int(x[1][1]) for x in all_possible_pairs],
                                  'mj':           [int(x[1][2]) for x in all_possible_pairs],
                                  'prob':         [prob for _ in all_possible_pairs],
                                  'reactantName': [R.product for _ in  all_possible_pairs],
                                  'order':        [order for _ in all_possible_pairs]})
                if stage==reaction_stage.cure:
                    # exclude atom pairs that have same resid or molid
                    idf=idf[(idf['ri']!=idf['rj'])&(idf['mi']!=idf['mj'])].copy()
                    logger.debug(f'Examining {idf.shape[0]} bond-candidates of order {order}')
                    if idf.shape[0]>0:
                        ess='' if idf.shape[0]!=1 else 's'
                        bondtestoutcomes={k:0 for k in BTRC}
                        p=Pool(processes=ncpu)
                        idf_split=np.array_split(idf,ncpu)
                        packets=[(i,idf_split[i]) for i in range(ncpu)]
                        logger.debug(f'Decomposed dataframe lengths: {", ".join([str(x.shape[0]) for x in idf_split])}')
                        results=p.map(partial(gromacs_distance,gro=gro,new_column_name='r'),packets)
                        p.close()
                        p.join()
                        # reassemble final dataframe:
                        # logger.debug(f'Checking dataframe lengths: {", ".join([str(x.shape[0]) for x in results])}')
                        idf=pd.concat(results,ignore_index=True)
                        # gromacs_distance(idf,gro,new_column_name='r') # use "gmx distance" to very quickly get all lengths
                        logger.debug(f'{idf.shape[0]} bond-candidate length{ess} avg/min/max: {idf["r"].mean():0.3f}/{idf["r"].min():0.3f}/{idf["r"].max():0.3f}')
                        idf=idf[idf['r']<self.current_radius].copy().reset_index(drop=True)
                        ess='' if idf.shape[0]!=1 else 's'
                        logger.debug(f'{idf.shape[0]} bond-candidate{ess} with lengths below {self.current_radius} nm')
                        p=Pool(processes=ncpu)
                        idf_split=np.array_split(idf,ncpu)
                        # logger.debug(f'Decomposed dataframe lengths: {", ".join([str(x.shape[0]) for x in idf_split])}')
                        results=p.map(partial(TC.bondtest_df),idf_split)
                        p.close()
                        p.join()
                        # reassemble final dataframe:
                        # logger.debug(f'Checking dataframe lengths: {", ".join([str(x.shape[0]) for x in results])}')
                        idf=pd.concat(results,ignore_index=True)
                        if not idf.empty:
                            logger.debug(f'Bond-candidate test outcomes:')
                            for k in bondtestoutcomes:
                                bondtestoutcomes[k]=idf[idf['result']==k].shape[0]
                                logger.debug(f'   {str(k)}: {bondtestoutcomes[k]}')
                            idf=idf[idf['result']==BTRC.passed].copy()
                elif stage==reaction_stage.cap:
                    idf=idf[idf['ri']==idf['rj']].copy()
                    logger.debug(f'Examining {idf.shape[0]} bond-candidates of order {order}')
                if not idf.empty:
                    bdf=pd.concat((bdf,idf),ignore_index=True)
        # logger.debug('Filtered (by single-bond tests) bonds:')
        # for ln in bdf.to_string().split('\n'):
        #     logger.debug(ln)

        if stage==reaction_stage.cure and bdf.shape[0]>0:
            bdf=bdf.sort_values('r',axis=0,ignore_index=True).reset_index(drop=True)
            bdf['allowed']=[True for x in range(bdf.shape[0])]
            unique_atomidx=set(bdf.ai.to_list()).union(set(bdf.aj.to_list()))
            unique_resids=set(bdf.ri.to_list()).union(set(bdf.rj.to_list()))
            for i,r in bdf.iterrows():
                if r.ai in unique_atomidx:
                    unique_atomidx.remove(r.ai)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated atom index {r.ai}')
                    bdf.loc[i,'allowed']=False
                if r.aj in unique_atomidx:
                    unique_atomidx.remove(r.aj)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated atom index {r.aj}')
                    bdf.loc[i,'allowed']=False
                if r.ri in unique_resids:
                    unique_resids.remove(r.ri)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated residue index {r.ri}')
                    bdf.loc[i,'allowed']=False
                if r.rj in unique_resids:
                    unique_resids.remove(r.rj)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated residue index {r.rj}')
                    bdf.loc[i,'allowed']=False

            logger.debug(f'{bdf[bdf["allowed"]==False].shape[0]} out of {bdf.shape[0]} bonds disallowed due to repeated atom indexes or residue indexes')

            bdf=bdf[bdf['allowed']==True].copy().reset_index(drop=True)
            # logger.debug('Allowed bonds:')
            # for ln in bdf.to_string().split('\n'):
            #     logger.debug(ln)

            bdf=TC.cycle_collective(bdf)
            logger.debug(f'{bdf[bdf["remove-to-uncyclize"]==True].shape[0]} out of {bdf.shape[0]} bonds removed to break nascent cycles')
            bdf=bdf[bdf['remove-to-uncyclize']==False].copy().reset_index(drop=True)
            # logger.debug('Non-cyclizing bonds:')
            # for ln in bdf.to_string().split('\n'):
            #     logger.debug(ln)

            ''' roll the dice '''
            bdf['lucky']=[True for _ in range(bdf.shape[0])]
            if apply_probabilities:
                for i,r in bdf.iterrows():
                    x=np.random.random()
                    if x>r.prob:
                        bdf.loc[i,'lucky']=False

            logger.debug(f'{bdf[bdf["lucky"]==True].shape[0]} bonds survive probability application')
            bdf=bdf[bdf['lucky']==True].copy().reset_index(drop=True)
            # logger.debug('Lucky bonds:')
            # for ln in bdf.to_string().split('\n'):
            #     logger.debug(ln)

            ''' apply the stated limit '''
            if abs_max>-1:
                if abs_max<bdf.shape[0]:
                    bdf=bdf.loc[:abs_max].copy().reset_index(drop=True)
                    logger.debug(f'Limiting to {bdf.shape[0]} allowed bonds')
            logger.debug('Final bonds:')
            for ln in bdf.to_string().split('\n'):
                logger.debug(ln)

        if stage==reaction_stage.cure:
            TC.linkcell_cleanup()

        return bdf