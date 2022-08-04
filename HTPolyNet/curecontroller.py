import pandas as pd
from enum import Enum
import numpy as np
import os
from itertools import product
from HTPolyNet.topocoord import TopoCoord, BTRC
from HTPolyNet.gromacs import gromacs_distance, mdp_modify
from HTPolyNet.configuration import ReactionList
from HTPolyNet.molecule import MoleculeDict
from multiprocessing import Pool
from functools import partial
import HTPolyNet.projectfilesystem as pfs

import logging

logger=logging.getLogger(__name__)

class state(Enum):
    """Enumerated CURE state
    """
    bondsearch=0
    drag=1
    update=2
    relax=3
    equilibrate=4
    postcure_bondsearch=5
    postcure_update=6
    postcure_relax=7
    postcure_equilibrate=8
    finished=9
    unknown=99
    def __str__(self):
        return self.name

class CureController:
    default_equilibration = [ { 'ensemble': 'min' }, 
        { 'ensemble': 'nvt', 'temperature': 600, 'nsteps': 1000 },
        { 'ensemble': 'npt', 'temperature': 600, 'pressure': 1, 'nsteps': 2000 }
    ]
    # internal name : name in cfg
    parameter_dicts={'cure':'CURE','drag':'drag','relax':'relax',
                    'equil':'equilibration','gromacs':'gromacs','pc_equil':'postcure_equilibration',
                    'anneal':'postcure_anneal','pa_equil':'postanneal_equilibration'}
    parameter_defaults={
        'cure': { 
            'max_conversion_per_iteration': 1.0,
            'search_radius': 0.5,
            'radial_increment': 0.05,
            'late_threshold': 0.85,
            'max_iterations': 100,
            'desired_conversion': 0.5
        },
        'drag': {
            'limit': 0.0,
            'kb': 300000.0,
            'trigger_distance': 0.0,
            'nstages': 0,
            'increment': 0.0,
            'cutoff_pad': 0.2,
            'equilibration': default_equilibration,
        },
        'relax': {
            'nstages': 6,
            'increment': 0.0,
            'cutoff_pad': 0.2,
            'equilibration': default_equilibration
        },
        'equil': {
            'temperature': 300,
            'pressure': 1,
            'nsteps': 50000,
            'ensemble': 'npt'
        },
        'gromacs': {
            'rdefault': 0.9
        },
        'postcure_equilibration': {
            'temperature': 300,
            'pressure': 1,
            'nsteps': 50000,
            'ensemble': 'npt'
        },
        'anneal': {
            'ncycles': 0 # default behavior is no annealing
        },
        'pa_equil': {
            'nsteps': 0 # default behavior is no post-anneal equilibration
        }
    }
    def __init__(self,basedict):
        self.iter=0
        self.max_nxlinkbonds=0
        self.cum_nxlinkbonds=0
        self.max_search_radius=0.0
        self.max_radidx
        self.state=state.unknown
        self.iter=0
        self.max_search_radius=0.0
        self.max_radidx=0
        self.curr_nxlinkbonds=0
        self.max_nxlinkbonds=0
        self.current_stage={}
        self.current_stage['drag']=0
        self.current_stage['relax']=0
        self.current_radidx=0
        self.radius=0.0
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'
        self.dicts={}
        for k,v in self.parameter_dicts.items():
            self.dicts[k]=basedict.get(v,{})
            for p,val in self.parameter_defaults[k].items():
                if not p in self.dicts[k]:
                    self.dicts[k][p]=basedict.get(f'{v}_{p}',val)
        self.dragging_enabled=False
        d=self.dicts['drag']
        if (d['nstages']>0 or d['increment']>0.0) and d['limit']>0.0:
            self.dragging_enabled=True
        self.ncpu=basedict.get('ncpu',os.cpu_count())
        self.bonds_file=basedict.get('bonds_file','bonds.csv')

    def setup(self,max_nxlinkbonds=0,max_search_radius=0.0):
        self.max_nxlinkbonds=max_nxlinkbonds
        self.max_search_radius=max_search_radius
        self.max_radidx=int((self.max_search_radius-self.dicts['cure']['search_radius'])/self.dicts['cure']['radial_increment'])

    def is_cured(self):
        return self.cum_nxlinkbonds==self.max_nxlinkbonds

    def curr_conversion(self):
        if not self.max_nxlinkbonds: return 0
        return self.cum_nxlinkbonds/self.max_nxlinkbonds

    def reset(self):
        self.iter=1
        self.state=state.bondsearch
        self.max_search_radius=0.0
        self.max_radidx=0
        self.curr_nxlinkbonds=0
        self.current_stage['drag']=0
        self.current_stage['relax']=0
        self.current_radidx=0
        self.radius=0.0
        self.bonds_file=None
        self.bonds_df=pd.DataFrame()
        self.bonds_are='nonexistent!'

    def next_iter(self):
        i=self.iter+1
        self.reset()
        self.iter=i
        return self.iter<self.dicts['cure']['max_iterations']

    def _write_bonds_df(self):
        self.bonds_df.to_csv(self.bonds_file,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def _read_bonds_df(self):
        assert os.path.exists(self.bonds_file),f'Error: {self.bonds_file} not found.'
        self.bonds_df=pd.read_csv(self.bonds_file,sep='\s+',header=0)

    def register_bonds(self,bonds,pairs,bonds_are='unrelaxed'):
        self.bonds_df=bonds
        self.pairs_df=pairs
        self.bonds_are=bonds_are

    def pfx(self):
        return f'{self.state.value}-{self.state}'

    def do_bondsearch(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict):
        if self.state!=state.bondsearch: return
        self.radius=self.dicts['cure']['search_radius']+self.current_radidx*self.dicts['cure']['radial_increment']
        logger.info(f'Bondsearch using radius {self.radius} nm initiated.')
        apply_probabilities=self.curr_conversion()<self.dicts['cure']['late_threshold']
        bond_limit=int(self.dicts['cure']['max_conversion_per_iteration']*self.max_nxlinkbonds)
        bond_target=int((self.dicts['cure']['desired_conversion']-self.curr_conversion)*self.max_nxlinkbonds)
        bond_limit=min([bond_limit,bond_target])
        logger.debug(f'Iteration limited to at most {bond_limit} new bonds')
        nbonds=0
        while nbonds==0 and self.current_radidx<self.max_radidx:
            nbdf=self.searchbonds(TC,RL,MD,stage='cure',abs_max=bond_limit,apply_probabilities=apply_probabilities)
            nbonds=nbdf.shape[0]
            if nbonds==0:
                self.current_radidx+=1
                self.radius+=self.dicts['cure']['radial_increment']
                logger.info(f'Increasing cutoff radius to {self.radius} nm')
        if nbonds>0:
            ess='' if nbonds==1 else 's'
            logger.info(f'CURE iteration {self.iter} will generate {nbdf.shape[0]} new bond{ess}.')
            pairs=pd.DataFrame() # empty placeholder
            self.register_bonds(nbdf,pairs,bonds_are='unrelaxed')
            self.bonds_df['initial_distance']=TC.return_bond_lengths(self.bonds_df)
            self.state=state.drag if self.dragging_enabled else state.update
        else:
            self.state=state.equilibrate

    def do_prebond_dragging(self,TC:TopoCoord):
        if self.state!=state.drag: return
        nbdf=self.bonds_df
        d=self.dicts['drag']
        nogos=[nbdf.shape[0]==0,not self.dragging_enabled,nbdf['initial_distance'].max()<d['trigger_distance']]
        if any(nogos):
            self.state=state.update
            return
        self.distance_attenuation(TC,mode='drag')
        self.state=state.update

    def do_relax(self,TC:TopoCoord):
        if self.state!=state.relax: return
        nbdf=self.bonds_df
        if nbdf.shape[0]==0:
            self.state=state.equilibrate
            return
        self.distance_attenuation(TC,mode='relax')
        self.state=state.equilibrate

    def distance_attenuation(self,TC:TopoCoord,mode='drag'):
        assert mode in ['drag','relax']
        opfx=self.pfx()
        nbdf=self.bonds_df
        pdf=self.pairs_df
        d=self.dicts[mode]
        nbdf['current_lengths']=nbdf['initial_distance'].copy()
        maxL,minL,meanL=nbdf['current_lengths'].max(),nbdf['current_lengths'].min(),nbdf['current_lengths'].mean()
        logger.debug(f'Lengths avg/min/max: {meanL:.3f}/{minL:.3f}/{maxL:.3f}')
        ess='' if nbdf.shape[0]==1 else 's'
        logger.info(f'{self.state} initiated on {nbdf.shape[0]} distance{ess} (max {maxL:.3f} nm).')
        roptions=[self.dicts['gromacs']['rdefault'],maxL]
        if mode=='drag':
            TC.add_restraints(nbdf,typ=6)
            logger.info('     Stage  Max-distance (nm)')
        else:
            pmaxL,pminL,pmeanL=pdf['current_lengths'].max(),pdf['current_lengths'].min(),pdf['current_lengths'].mean()
            logger.debug(f'1-4 distances lengths avg/min/max: {pmeanL:.3f}/{pminL:.3f}/{pmaxL:.3f}')
            roptions.append(pmaxL)
            logger.info('     Stage  Max-distance (nm)  Max-1-4-distance (nm)')
        rcommon=max(roptions)
        mdp_mods_dict={'rvdw':rcommon,'rcoulomb':rcommon,'rlist':rcommon,'gen-temp':stg_dict['temperature'],'ref_t':stg_dict['temperature'],'gen-vel':'yes','nsteps':stg_dict['nsteps']}
        if ensemble=='npt':
            mdp_mods_dict['ref_p']=stg_dict['pressure']
        for stg_dict in d['equilibration']:
            ensemble=stg_dict['ensemble']
            impfx=f'{self.state}-{ensemble}' # e.g., drag-min, drag-nvt, drag-npt
            pfs.checkout(f'mdp/{impfx}.mdp')
            mdp_modify(f'{impfx}.mdp',mdp_mods_dict,add_if_missing=(ensemble!='min'))
        this_nstages=int(maxL/d['increment'])
        this_firststage=self.current_stage[mode]
        saveT=TC.copy_bond_parameters(self.bonds_df)
        for i in range(this_firststage,this_nstages):
            TC.attenuate_bond_parameters(self.bonds_df,i,this_nstages,minimum_distance=d['limit'],init_colname='initial_distance')
            stagepfx=f'{opfx}-stage-{i+1}'
            for stg_dict in d['equilibration']:
                ensemble=stg_dict['ensemble']
                TC.grompp_and_mdrun(out=f'{stagepfx}-{ensemble}',mdp=f'{self.state}-{ensemble}')
                logger.debug(f'{TC.files["gro"]}')
            TC.restore_bond_parameters(saveT)
            self.bonds_df['current_lengths']=np.array(TC.return_bond_lengths(self.bonds_df))
            maxL,minL,meanL=self.bonds_df['current_lengths'].max(),self.bonds_df['current_lengths'].min(),self.bonds_df['current_lengths'].mean()
            logger.debug(f'Distances avg/min/max: {meanL:.3f}/{minL:.3f}/{maxL:.3f}')
            if mode=='drag':
                logger.info(f'{i+1:>10d}  {maxL:>17.3f}')
            else:
                pmaxL,pminL,pmeanL=pdf['current_lengths'].max(),pdf['current_lengths'].min(),pdf['current_lengths'].mean()
                logger.debug(f'1-4 distances lengths avg/min/max: {pmeanL:.3f}/{pminL:.3f}/{pmaxL:.3f}')
                logger.info(f'{i+1:>10d}  {maxL:>17.3f}  {pmaxL:>21.3f}')
        if mode=='drag':
            TC.remove_restraints(self.bonds_df)

    def do_topology_update(self,TC:TopoCoord,MD:MoleculeDict):
        if self.state!=state.update or self.state!=state.postcure_update: return
        opfx=self.pfx()
        logger.info(f'Topology update')
        self.bonds_df,self.pairs_df=TC.update_topology_and_coordinates(self.bonds_df,template_dict=MD,write_mapper_to=f'{opfx}-idx-mapper.csv')
        self.bonds_df['initial_distance']=TC.return_bond_lengths(self.bonds_df)
        if self.state==state.update:
            self.state=state.relax
        else:
            self.state=state.postcure_relax

    def do_equilibrate(self,TC:TopoCoord):
        if self.state!=state.equilibrate or self.state!=state.postcure_equilibrate: return
        d=self.dicts['equil']
        opfx=self.pfx()
        logger.info(f'Equilibration for {d["nsteps"]} steps at {d["temperature"]} K and {d["pressure"]} bar')
        pfx=f'{self.state}-d["ensemble"]'
        pfs.checkout(f'mdp/{pfx}.mdp')
        mod_dict={'gen-temp':d['temperature'],'gen-vel':'yes','ref_t':d['temperature'],'ref_p':d['pressure'],'nsteps':d['nsteps']}
        mdp_modify(f'{pfx}.mdp',mod_dict,new_filename=f'{opfx}.mdp')
        TC.grompp_and_mdrun(out=f'{opfx}-post',mdp=f'{opfx}',quiet=False)
        average_density=trace('Density',[f'{opfx}-post'],outfile='density.png')
        logger.info(f'  -> average density {average_density:.3f} kg/m^3')
        if self.state==state.equilibrate:
            self.next_iter()
        elif self.state==state.postcure_equilibrate:
            self.state=state.finished

    def do_postcure_bondsearch(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict):
        if self.state!=state.postcure_bondsearch: return
        nbdf=self.searchbonds(TC,RL,MD,stage='post-cure')
        nbonds=nbdf.shape[0]
        if nbonds>0:
            ess='' if nbonds==1 else 's'
            logger.info(f'CURE iteration {self.iter} will generate {nbdf.shape[0]} new bond{ess}.')
            pairs=pd.DataFrame() # empty placeholder
            self.register_bonds(nbdf,pairs,bonds_are='unrelaxed')
            self.bonds_df['initial_distance']=TC.return_bond_lengths(self.bonds_df)
            self.state=state.postcure_update
        else:
            self.state=state.postcure_equilibrate

    def searchbonds(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,stage='cure',abs_max=0,apply_probabilities=False):
        adf=TC.gro_DataFrame('atoms')
        gro=TC.files['gro']
        if stage=='cure':
            TC.linkcell_initialize(self.radius,ncpu=self.ncpu)
        raset=adf[adf['z']>0]  # this view will be used for downselecting to potential A-B partners
        bdf=pd.DataFrame()
        Rlist=[x for x in RL if (x.stage==stage and x.probability>0.0)]
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
                if stage=='post-cure':
                    assert areactantname_template==breactantname_template,f'Error: post-cure reaction {R.name} lists a bond whose atoms are in different reactants'
                    assert aresname==bresname,f'Error: post-cure reaction {R.name} lists a bond whose atoms are in different residues'

                Aset=raset[(raset['atomName']==aname)&(raset['resName']==aresname)&(raset['z']==az)&(raset['reactantName']==areactantname_template)]
                Bset=raset[(raset['atomName']==bname)&(raset['resName']==bresname)&(raset['z']==bz)&(raset['reactantName']==breactantname_template)]
                alist=list(zip(Aset['globalIdx'].to_list(),Aset['resNum'].to_list()))
                blist=list(zip(Bset['globalIdx'].to_list(),Bset['resNum'].to_list()))
                all_possible_pairs=list(product(alist,blist))

                idf=pd.DataFrame({'ai':           [int(x[0][0]) for x in all_possible_pairs],
                                  'ri':           [int(x[0][1]) for x in all_possible_pairs],
                                  'aj':           [int(x[1][0]) for x in all_possible_pairs],
                                  'rj':           [int(x[1][1]) for x in all_possible_pairs],
                                  'prob':         [prob for _ in all_possible_pairs],
                                  'reactantName': [R.product for _ in  all_possible_pairs],
                                  'order':        [order for _ in all_possible_pairs]})
                if stage=='cure':
                    # exclude atom pairs that have same resid
                    idf=idf[idf['ri']!=idf['rj']].copy()
                    logger.debug(f'Examining {idf.shape[0]} bond-candidates of order {order}')
                    if idf.shape[0]>0:
                        ess='' if idf.shape[0]!=1 else 's'
                        bondtestoutcomes={k:0 for k in BTRC}
                        gromacs_distance(idf,gro) # use "gmx distance" to very quickly get all lengths
                        # idf.info(buf=buffer)
                        # logger.debug(f'after gromacs distance:\n{buffer.getvalue()}')
                        logger.debug(f'{idf.shape[0]} bond-candidate length{ess} avg/min/max: {idf["r"].mean():0.3f}/{idf["r"].min():0.3f}/{idf["r"].max():0.3f}')
                        idf=idf[idf['r']<self.radius].copy().reset_index(drop=True)
                        # idf.info(buf=buffer)
                        # logger.debug(f'after distance filter:\n{buffer.getvalue()}')
                        
                        # logger.debug('Filtered (by gromacs distance) bonds:')
                        # for ln in bdf.to_string().split('\n'):
                        #     logger.debug(ln)
                        ess='' if idf.shape[0]!=1 else 's'
                        logger.debug(f'{idf.shape[0]} bond-candidate{ess} with lengths below {self.radius} nm')
                        p=Pool(processes=self.ncpu)
                        idf_split=np.array_split(idf,self.ncpu)
                        logger.debug(f'Decomposed dataframe lengths: {", ".join([str(x.shape[0]) for x in idf_split])}')
                        results=p.map(partial(TC.bondtest_df),idf_split)
                        p.close()
                        p.join()
                        # reassemble final dataframe:
                        logger.debug(f'Checking dataframe lengths: {", ".join([str(x.shape[0]) for x in results])}')
                        idf=pd.concat(results,ignore_index=True)
                        if not idf.empty:
                            logger.debug(f'Bond-candidate test outcomes:')
                            for k in bondtestoutcomes:
                                bondtestoutcomes[k]=idf[idf['result']==k].shape[0]
                                logger.debug(f'   {str(k)}: {bondtestoutcomes[k]}')
                            idf=idf[idf['result']==BTRC.passed].copy()
                elif stage=='post-cure':
                    idf=idf[idf['ri']==idf['rj']].copy()
                    logger.debug(f'Examining {idf.shape[0]} bond-candidates of order {order}')
                if not idf.empty:
                    bdf=pd.concat((bdf,idf),ignore_index=True)
        # logger.debug('Filtered (by single-bond tests) bonds:')
        # for ln in bdf.to_string().split('\n'):
        #     logger.debug(ln)

        if stage=='cure' and bdf.shape[0]>0:
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

            logger.debug(f'{bdf[bdf["lucky"]==True].shape[0]} bonds survive probability application.')
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

        return bdf