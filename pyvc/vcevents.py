#!/usr/bin/env python
from . import VCSys

class VCEvents(VCSys):
    def __init__(self, sim_data):
        #print os.path.isfile(sim_file)
        
        super(VCEvents, self).__init__(sim_data)
    
        self.event_data = self.sim_data.file.root.event_table
        self.sweep_data = self.sim_data.file.root.event_sweep_table
        #self.sim_data.a += 'a has a events. '
        
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
    def unpack_event_range(self, event_range):
        #TODO: Better error checking here
        if event_range is not None:
            return event_range['type'], event_range['filter'][0], event_range['filter'][1]
        else:
            return 'ids', 0, self.event_data.nrows
            
    
    def get_event_elements(self, evnum):
        #print self.event_data[evnum]['start_sweep_rec'], self.event_data[evnum]['end_sweep_rec']
        return set([x['block_id'] for x in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]])
        #return [x['block_id'] for x in self.sweep_data.where('event_number == {}'.format(evnum))]
    
    def get_event_slip_area(self, evnum):
        areas = {}
        total_slip = 0.0
        slip_records = 0
        for sweep in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]:
            areas[sweep['block_id']] = sweep['area']
            total_slip += sweep['slip']
            slip_records += 1
        #print slip_records
        return (total_slip, sum(areas.values()), slip_records)
    
    def get_event_data(self, requested_data, event_range=None, magnitude_filter=None):
        type, start, stop = self.unpack_event_range(event_range)
        
        event_data = {}
        
        for col_name in requested_data:
            event_data[col_name] = []
            
        get_involved_elements = False
        if 'involved_elements' in requested_data:
            requested_data.remove('involved_elements')
            get_involved_elements = True
        
        query_str = ''
        if type == 'years':
            query_str += '(event_year >= %f) & (event_year < %f)'%(start, stop)
        elif type == 'ids':
            query_str += '(event_number >= %f) & (event_number < %f)'%(start, stop)
        
        if magnitude_filter is not None:
            query_str += ' & (event_magnitude %s)'%magnitude_filter
        '''
        if section_filter is not None:
            query_str += ' & ('
            for n, item in enumerate(section_filter['filter']):
                query_str += '(%s == %s)'%(section_filter['type'],item)
                if n != len(section_filter['filter']) - 1:
                    query_str += ' | '
            query_str += ')'
        '''
        
        for event in self.event_data.where(query_str):
            #if get_involved_elements:
                #involved_elements = [x['block_id'] for x in self.sweep_data.where('event_number == %i'%event['event_number'])]
                #involved_elements = self.sweep_data.cols.block_id[event['start_sweep_rec']:event['end_sweep_rec']]
                #event_data[col_name].append(involved_elements)
            for col_name in requested_data:
                event_data[col_name].append(event[col_name])
        #print type, start, stop
        #for event in self.event_data.where(
        return event_data
    
    @property
    def num_events(self):
        return self.event_data.nrows
    
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
