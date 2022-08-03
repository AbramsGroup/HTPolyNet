
class CureController:
    def __init__(self,basedict):
        self.cure_dict=basedict.get('CURE',{})
        self.drag_dict=basedict.get('drag',{})
        self.relax_dict=basedict.get('relax',{})
        self.equil_dict=basedict.get('postcure_equilibration',{})
        self.gromacs_dict=basedict.get('gromacs',{})
        self.pc_equil_dict=basedict.get('postcure_equilibration',{})
        if not self.cure_dict:
            self.cure_dict['max_conversion_per_iteration']=basedict.get('CURE_max_conversion_per_iteration',1.0)
            self.cure_dict['search_radius']=basedict.get('CURE_initial_search_radius',0.5)
            self.cure_dict['radial_increment']=basedict.get('CURE_radial_increment',0.25)
            self.cure_dict['late_threshold']=basedict.get('CURE_late_threshold',1.0)
            self.cure_dict['max_iterations']=basedict.get('CURE_max_iterations',100)
            self.cure_dict['desired_conversion']=basedict.get('CURE_desired_conversion',0.5)
        self.dragging_enabled=False
        if not self.drag_dict:
            self.drag_dict['limit']=basedict.get('drag_limit',0.0)
            self.drag_dict['trigger_distance']=basedict.get('drag_trigger_distance',0.0)
            self.drag_dict['nstages']=basedict.get('drag_nstages',0)
            self.drag_dict['increment']=basedict.get('drag_increment',0.0)
            self.drag_dict['temperature']=basedict.get('drag_temperature',300.0)
            self.drag_dict['pressure']=basedict.get('drag_pressure',1.0)
            self.drag_dict['cutoff_pad']=basedict.get('drag_cutoff_pad',0.2)
            if (self.drag_dict['nstages']>0 or self.drag_dict['increment']>0.0) and self.drag_dict['limit']>0.0:
                self.dragging_enabled=True
        if not self.relax_dict:
            self.relax_dict['nstages']=basedict.get('relax_nstages',6)
            self.relax_dict['increment']=basedict.get('relax_increment',0)
            self.relax_dict['temperature']=basedict.get('relax_temperature',300.0)
            self.relax_dict['pressure']=basedict.get('relax_pressure',1.0)
            self.relax_dict['cutoff_pad']=basedict.get('relax_cutoff_pad',0.2)
        if not self.equil_dict:
            self.equil_dict['temperature']=basedict.get('equilibration_temperature',300.0)
            self.equil_dict['pressure']=basedict.get('equilibration_pressure',1.0)
            self.equil_dict['nsteps']=basedict.get('equilibration_steps',50000)
        if not self.gromacs_dict:
            self.gromacs_dict['rdefault']=basedict.get('gromacs_rdefault',0.9)
        if not self.pc_equil_dict:
            self.pc_equil_dict['temperature']=basedict.get('postcure_equilibration_temperature',300.0)
            self.pc_equil_dict['pressure']=basedict.get('postcure_equilibration_pressure',1.0)
            self.pc_equil_dict['nsteps']=basedict.get('postcure_equilibration_steps',50000)
