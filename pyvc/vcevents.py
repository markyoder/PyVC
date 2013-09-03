#!/usr/bin/env python
from . import VCSys

class VCEvents(VCSys):
    def __init__(self, sim_file = ''):
        #print os.path.isfile(sim_file)
        
        super(VCEvents, self).__init__(sim_file = sim_file)
        
        #self.sys = sys
        #self._years = None
        #self._event_ids = None
    #self._mags = None
    #self.simple_event_list = None
    
    '''
        @property
        def events(self):
        f = h5py.File(self.sys.data_file, 'r')
        events = f['event_table']
        f.close()
        return events
        
        @property
        def sweeps(self):
        f = h5py.File(self.sys.data_file, 'r')
        sweeps = f['event_sweep_table']
        f.close()
        return sweeps
        '''
    
    @property
    def years(self):
        f = self.openSimFile()
        #f = h5py.File(self.sim_file, 'r')
        events = f['event_table']
        years = events['event_year']
        f.close()
        return years
    
    @property
    def mags(self):
        f = h5py.File(self.sys.data_file, 'r')
        events = f['event_table']
        mags = events['event_magnitude']
        f.close()
        return mags
    
    @property
    def event_ids(self):
        f = h5py.File(self.sys.data_file, 'r')
        events = f['event_table']
        event_ids = events['event_number']
        f.close()
        return event_ids
