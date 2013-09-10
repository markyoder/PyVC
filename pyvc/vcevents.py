#!/usr/bin/env python
from . import VCSys
from . import vcexceptions

import itertools
from operator import itemgetter

#-------------------------------------------------------------------------------
# A class representing the set of events in a Virtual California simulation.
#-------------------------------------------------------------------------------
class VCEvents(VCSys):
    def __init__(self, sim_data):
        
        super(VCEvents, self).__init__(sim_data)
        
        # the tables we need from the simulation data
        self.event_data = self.sim_data.file.root.event_table
        self.sweep_data = self.sim_data.file.root.event_sweep_table
    
    @property
    def num_events(self):
        return self.event_data.nrows
    
    def get_event_year(self, evnum):
        return self.event_data[evnum]['event_year']
    
    def get_event_elements(self, evnum):
        return set([x['block_id'] for x in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]])
    
    def get_event_slip_area(self, evnum):
        areas = {}
        total_slip = 0.0
        slip_records = 0
        for sweep in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]:
            areas[sweep['block_id']] = sweep['area']
            total_slip += sweep['slip']
            slip_records += 1
        
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
            return 'id', 0, self.num_events-1

    #---------------------------------------------------------------------------
    # This method returns requested event data subject to several filters. The
    # data that will be returned is defined in requested_data which is a
    # list of strings. The allowed values for requested_data are:
    #
    #     'event_number'                 - the unique event id
    #     'event_year'                   - the decimal year in which the event
    #                                      occured
    #     'event_trigger'                - the block_in of the element that
    #                                      triggered the event
    #     'event_magnitude'              - the magnitude of the event
    #     'event_shear_init'             - the shear stress before the event
    #     'event_normal_init'            - the normal stress before the event
    #     'event_shear_final'            - the shear stress after the event
    #     'event_normal_final'           - the normal stress after the event
    #     'start_sweep_rec'              - the index of the first sweep record
    #                                      that applies to this event
    #     'end_sweep_rec'                - the index of the last sweep record
    #                                      that applies to this event
    #     'start_aftershock_rec'         - the index of the first aftershock
    #                                      record that applies to this event
    #     'end_aftershock_rec'           - the index of the last aftershock
    #                                      record that applies to this event
    #     'event_area'                   - the fault surface rupture area of the
    #                                      event in meters^2
    #     'event_average_slip'           - the average slip during the event in
    #                                      meters
    #     'event_surface_rupture_length' - the surface rupture length of the
    #                                      event in meters
    #     'event_range_duration'         - the duration in years of the
    #                                      requested event set
    #
    # The filters are as follows:
    #
    #     event_range = {'type':['year'|'id'],'filter':(start,stop)}
    #         example: {'type':'year','filter':(1000,5000)} only grabs events
    #         between year 1000 and 5000 inclusive
    #     section_filter = {'type':'section_id','filter':[list of section ids]}
    #         example: {'type':'section_id','filter':[1,2,3,10]} oly grabs
    #         events on sections 1, 2, 3, and 10
    #     magnitude_filter = '[conditional] [magnitude]'
    #         example '>= 6.0' only grabs events bigger than magnitude 6.0
    #
    # See pyvc_tests.py for examples.
    #---------------------------------------------------------------------------
    def get_event_data(self, requested_data, event_range=None, magnitude_filter=None, section_filter=None):
        
        event_data = {}
        
        get_event_range_duration = False
        if 'event_range_duration' in requested_data:
            requested_data.remove('event_range_duration')
            get_event_range_duration = True
        
        for col_name in requested_data:
            event_data[col_name] = []
        
        #-----------------------------------------------------------------------
        # Here is where the events are being grabbed. Code readability is
        # sacraficed to avoid conditionals in big loops where possible. The
        # section filter is the most expensive filter because it doesnt allow us
        # to use PyTables in kernal queries. Because of this we also need a
        # different strategy when filtering by section. The approach to applying
        # the section filter is:
        #
        # 1. Use data saved in the simulation file to compile a list of event
        #    ids that occured on the selected sections in the selected time
        #    frame. This is done first because there are always fewer sections
        #    than there are events so this is a more manageable calculation and
        #    the data in the simulation file is structured to make this easier.
        # 2. This list of event ids may not be contiguous so we cant use ranges
        #    to grab the data. Instead an itemgetter is used.
        # 3. Loop over the events applying the magnitude filter if necessary.
        #-----------------------------------------------------------------------
        # TODO: Implement a faster section filter.
        if section_filter is not None:
            ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evid_range=True)
            
            # get all of the eligible event ids from the simulation data
            eligible_evids_tmp = []
            for secid in section_filter['filter']:
                try:
                    eligible_evids_tmp.append(self.sim_data.file.root.events_by_section._v_children['section_{}'.format(secid)].read())
                except KeyError:
                    raise vcexceptions.BadSectionID(secid)
            # flatten the list, remove duplicates, get rid of events outside the
            # given range, and sort
            eligible_evids = sorted(list(set(itertools.ifilter(lambda x: x>=ev_start and x<=ev_stop ,list(itertools.chain.from_iterable(eligible_evids_tmp))))))
            # if there are no events in the remaining set throw an error
            if len(eligible_evids) == 0:
                raise vcexceptions.NoEventsFound(event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
            # if there is only one event the itemgetter doesnt work so just
            # return the requested data
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
                # an itemgetter for our (potentially) non contiguous event ids
                ev_getter = itemgetter(*eligible_evids)
                
                # only apply the magnitude filter if necessary
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
            #-------------------------------------------------------------------
            # If there is no section filter life is good. PyTables makes this
            # super fast and easy.
            #-------------------------------------------------------------------
            ev_type, ev_start, ev_stop = self.unpack_event_range(event_range)
            
            query_str = ''
            if ev_type == 'year':
                query_str += '(event_year >= {}) & (event_year <= {})'.format(ev_start, ev_stop)
            elif ev_type == 'id':
                query_str += '(event_number >= {}) & (event_number <= {})'.format(ev_start, ev_stop)

            if magnitude_filter is not None:
                query_str += ' & (event_magnitude {})'.format(magnitude_filter)
            
            for event in self.event_data.where(query_str):
                for col_name in requested_data:
                    event_data[col_name].append(event[col_name])
        
        # get the event range duration if requested
        if get_event_range_duration:
            if ev_type != 'id':
                ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evid_range=True)
            event_data['event_range_duration'] = self.get_event_year(ev_stop) - self.get_event_year(ev_start)
    
        return event_data
    
