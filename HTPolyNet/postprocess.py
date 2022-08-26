import logging
import shutil
import numpy as np
from HTPolyNet.projectfilesystem import ProjectFileSystem
from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.gromacs import mdp_get, mdp_modify, gmx_energy_trace
import HTPolyNet.software as software
from HTPolyNet.configuration import Configuration
from HTPolyNet.plot import scatter

logger=logging.getLogger(__name__)

class tladder:
    def __init__(self,initstr):
        try:
            tok=initstr.split(',')
            self.T=[float(tok[0]),float(tok[1])]
            self.N=int(tok[2])
            self.ps_per_run=float(tok[3])
            self.ps_per_rise=float(tok[4])
            self.warmup_ps=float(tok[5])
            self.empty=False
        except:
            self.empty=True

    def __str__(self):
        if self.empty: return f'Tladder: empty'
        restr=f'Tladder: {self.T[0]}K -> {self.T[1]} K in {self.N} rung at {self.ps_per_run} ps per run and {self.ps_per_rise} ps per rise; warmup for {self.warmup_ps} ps'
        return restr

    def t_ladder(self,TC:TopoCoord,mdp_pfx='npt',**gromacs_dict):
        timestep=float(mdp_get(f'{mdp_pfx}.mdp','dt'))
        Tladder=np.linspace(self.T[0],self.T[1],self.N)
        timepoints=[0.0,self.warmup_ps]
        temppoints=[self.T[0],self.T[0]]
        for i in range(self.N):
            cT=Tladder[i]
            timepoints.append(timepoints[-1]+self.ps_per_rise)
            temppoints.append(cT)
            timepoints.append(timepoints[-1]+self.ps_per_run)
            temppoints.append(cT)
        duration=timepoints[-1]
        nsteps=int(duration/timestep)
        mod_dict={
            'ref_t':self.T[0],
            'gen-temp':self.T[0],
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'single',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(f'{mdp_pfx}.mdp',mod_dict)
        deffnm='t-ladder'
        msg=TC.grompp_and_mdrun(out=deffnm,mdp=mdp_pfx,quiet=False,mylogger=logger.info,**gromacs_dict)
        df=gmx_energy_trace(deffnm,['Temperature','Density','Volume'])
        scatter(df,'Temperature',['Density'],'rho_v_T.png')
        logger.info(f'Final coordinates in {deffnm}.gro')

def postprocess(args):
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(format='%(levelname)s> %(message)s',level=loglevel_numeric)
    ess='y' if len(args.proj)==0 else 'ies'
    cfg=Configuration.read(args.cfg,parse=False)
    # logger.info(f'{cfg.basedict}')
    logger.info(f'project director{ess}: {args.proj}')
    logger.info(f'{str(args.Tladder)}')
    software.sw_setup()
    software.set_gmx_preferences(cfg.basedict.get('gromacs',{}))
    for d in args.proj:
        pfs=ProjectFileSystem(projdir=d,topdirs=['molecules','systems','plots','postprocess'])
        pfs.go_to('postprocess')
        shutil.copy('../systems/final-results/final.gro','.')
        shutil.copy('../systems/final-results/final.top','.')
        shutil.copy('../systems/final-results/final.grx','.')
        TC=TopoCoord(topfilename='final.top',grofilename='final.gro',grxfilename='final.grx')
        logger.info(f'{TC.Coordinates.A.shape[0]} atoms {TC.total_mass(units="gromacs"):.2f} amu')
        if not args.Tladder.empty:
            pfs.library.checkout('mdp/npt.mdp')
            args.Tladder.t_ladder(TC)

