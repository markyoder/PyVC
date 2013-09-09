#!/usr/bin/env python
from . import VCSys
from . import vcexceptions

import itertools
from operator import itemgetter

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
            
    def get_event_year(self, evnum):
        return self.event_data[evnum]['event_year']
    
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
    
    def unpack_event_range(self, event_range, evid_range=False):
        #TODO: Better error checking here
        if event_range is not None:
            if evid_range and event_range['type'] != 'id':
                # convert years to event ids
                events_tmp = self.event_data.read_where('(event_year >= {}) & (event_year <= {})'.format(event_range['filter'][0], event_range['filter'][1]))
                return 'id', events_tmp[0]['event_number'], events_tmp[-1]['event_number']
            else:
                return event_range['type'], event_range['filter'][0], event_range['filter'][1]
        else:
            return 'id', 0, self.event_data.nrows-1
    
    def tmp_func(self, ev_start, ev_stop, eligible_evids_tmp):
        return sorted(list(set(itertools.ifilter(lambda x: x>=ev_start and x<=ev_stop ,list(itertools.chain.from_iterable(eligible_evids_tmp))))))
    
    def get_event_data(self, requested_data, event_range=None, magnitude_filter=None, section_filter=None):
        
        event_data = {}
        
        get_event_range_duration = False
        if 'event_range_duration' in requested_data:
            requested_data.remove('event_range_duration')
            get_event_range_duration = True
        
        for col_name in requested_data:
            event_data[col_name] = []
        
        # TODO: The section filter really slows things down. Need to find a way
        # get better performance here.
        if section_filter is not None:
            ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evid_range=True)
            
            eligible_evids_tmp = []
            for secid in section_filter['filter']:
                try:
                    eligible_evids_tmp.append(self.sim_data.file.root.events_by_section._v_children['section_{}'.format(secid)].read())
                except KeyError:
                    raise vcexceptions.BadSectionID(secid)
            #eligible_evids = sorted(list(set(itertools.ifilter(lambda x: x>=ev_start and x<=ev_stop ,list(itertools.chain.from_iterable(eligible_evids_tmp))))))
            eligible_evids = self.tmp_func(ev_start, ev_stop, eligible_evids_tmp)
            if len(eligible_evids) == 0:
                raise vcexceptions.NoEventsFound(event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
            
            if len(eligible_evids) == 1:
                event = self.event_data[eligible_evids[0]]
                if magnitude_filter is not None:
                    if eval('event_magnitude {}'.format(magnitude_filter), {'event_magnitude':event['event_magnitude']}):
                        for col_name in requested_data:
                            event_data[col_name].append(event[col_name])
                    else:
                        raise vcexceptions.NoEventsFound(event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
                else:
                    for col_name in requested_data:
                            event_data[col_name].append(event[col_name])
            else:
                ev_getter = itemgetter(*eligible_evids)
                
                #query_str = 'event_number == '
                #query_str += ' & event_number == '.join([str(x) for x in eligible_evids])
                
                #print query_str
                
                if magnitude_filter is not None:
                    exp_as_func = eval('lambda mag: ' + 'mag {}'.format(magnitude_filter))
                    for event in ev_getter(self.event_data):
                        if exp_as_func(event['event_magnitude']):
                            for col_name in requested_data:
                                event_data[col_name].append(event[col_name])
                else:
                    for event in ev_getter(self.event_data):
                        for col_name in requested_data:
                            event_data[col_name].append(event[col_name])
        else:
            ev_type, ev_start, ev_stop = self.unpack_event_range(event_range)
            
            query_str = ''
            if ev_type == 'year':
                query_str += '(event_year >= {}) & (event_year <= {})'.format(ev_start, ev_stop)
            elif ev_type == 'id':
                query_str += '(event_number >= {}) & (event_number <= {})'.format(ev_start, ev_stop)

            if magnitude_filter is not None:
                query_str += ' & (event_magnitude {})'.format(magnitude_filter)
            
            for event in self.event_data.where(query_str):
                #if get_involved_elements:
                    #involved_elements = [x['block_id'] for x in self.sweep_data.where('event_number == %i'%event['event_number'])]
                    #involved_elements = self.sweep_data.cols.block_id[event['start_sweep_rec']:event['end_sweep_rec']]
                    #event_data[col_name].append(involved_elements)
                #[x for x in self.element_data.where('event_number == {}'.format(event['event_number']))]
                for col_name in requested_data:
                    event_data[col_name].append(event[col_name])
            #print type, start, stop
            #for event in self.event_data.where(
    
        if get_event_range_duration:
            if ev_type != 'id':
                ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evid_range=True)
            event_data['event_range_duration'] = self.get_event_year(ev_stop) - self.get_event_year(ev_start)
    
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
