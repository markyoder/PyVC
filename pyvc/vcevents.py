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
    
    def __getitem__(self, evnum):
        return self.event_data[evnum]
    
    @property
    def num_events(self):
        return self.event_data.nrows
    
    def get_event_year(self, evnum):
        return self.event_data[evnum]['event_year']
    
    def get_event_elements(self, evnum):
        #return set([x['block_id'] for x in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]])
        return set(self.sweep_data.read(self.event_data[evnum]['start_sweep_rec'],self.event_data[evnum]['end_sweep_rec'], field='block_id'))
    
    def get_event_element_slips(self, evnum):
        slips = {}
        for sweep in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]:
            try:
                slips[sweep['block_id']] += sweep['slip']
            except KeyError:
                slips[sweep['block_id']] = sweep['slip']
        return slips

    
    def get_event_slip_area(self, evnum):
        areas = {}
        total_slip = 0.0
        slip_records = 0
        for sweep in self.sweep_data[self.event_data[evnum]['start_sweep_rec']:self.event_data[evnum]['end_sweep_rec']]:
            areas[sweep['block_id']] = sweep['area']
            total_slip += sweep['slip']
            slip_records += 1
        
        return (total_slip, sum(areas.values()), slip_records)
    
    def unpack_event_range(self, event_range, evnum_range=False):
        #TODO: Better error checking here
        if event_range is not None:
            if evnum_range and event_range['type'] != 'id':
                # convert years to event ids
                events_tmp = self.event_data.read_where('(event_year >= {}) & (event_year <= {})'.format(event_range['filter'][0], event_range['filter'][1]))
                return 'id', events_tmp[0]['event_number'], events_tmp[-1]['event_number']
            else:
                return event_range['type'], event_range['filter'][0], event_range['filter'][1]
        else:
            return 'id', 0, self.num_events-1

    def get_event_data_from_evnums(self, evnums, requested_data, event_data=None, event_range=None, magnitude_filter=None, section_filter=None):
        
        return_event_data = False
        if event_data is None:
            event_data = {}
            return_event_data = True
            for col_name in requested_data:
                event_data[col_name] = []
        
        ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evnum_range=True)
        
        # get rid of events outside the given range, and sort
        eligible_evnums = sorted(itertools.ifilter(lambda x: x>=ev_start and x<=ev_stop ,evnums))
        # if there are no events in the remaining set throw an error
        if len(eligible_evnums) == 0:
            raise vcexceptions.NoEventsFound(event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        # if there is only one event the itemgetter doesnt work so just
        # return the requested data
        if len(eligible_evnums) == 1:
            event = self.event_data[eligible_evnums[0]]
            if magnitude_filter is not None:
                exp_as_func = eval('lambda mag: ' + 'mag {}'.format(magnitude_filter))
                self.process_events(event_data, requested_data, event, exp_as_func)
            else:
                self.process_events(event_data, requested_data, event)
        else:
            # an itemgetter for our (potentially) non contiguous event ids
            ev_getter = itemgetter(*eligible_evnums)
            
            # only apply the magnitude filter if necessary
            if magnitude_filter is not None:
                exp_as_func = eval('lambda mag: ' + 'mag {}'.format(magnitude_filter))
                self.process_events(event_data, requested_data, ev_getter(self.event_data), exp_as_func)
            else:
                self.process_events(event_data, requested_data, ev_getter(self.event_data))

        if return_event_data:
            return event_data
    
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
            ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evnum_range=True)
            
            # get all of the eligible event ids from the simulation data
            eligible_evnums_tmp = []
            for secid in section_filter['filter']:
                try:
                    eligible_evnums_tmp.append(self.sim_data.file.root.events_by_section._v_children['section_{}'.format(secid)].read())
                except KeyError:
                    raise vcexceptions.BadSectionID(secid)

            self.get_event_data_from_evnums(
                list(set(list(itertools.chain.from_iterable(eligible_evnums_tmp)))),
                requested_data,
                event_data=event_data,
                event_range=event_range,
                magnitude_filter=magnitude_filter
            )
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
                
            self.process_events(event_data, requested_data, self.event_data.where(query_str))
        
        # get the event range duration if requested
        if get_event_range_duration:
            if ev_type != 'id':
                ev_type, ev_start, ev_stop = self.unpack_event_range(event_range, evnum_range=True)
            event_data['event_range_duration'] = self.get_event_year(ev_stop) - self.get_event_year(ev_start)
    
        return event_data

    def process_events(self, event_data, requested_data, event_iterable, magnitude_filter=None):
        if magnitude_filter is not None:
            for event in event_iterable:
                if magnitude_filter(event['event_magnitude']):
                    for col_name in requested_data:
                        if col_name == 'event_elements':
                            event_data[col_name].append(self.get_event_elements(event['event_number']))
                        else:
                            event_data[col_name].append(event[col_name])
        else:
            for event in event_iterable:
                for col_name in requested_data:
                    if col_name == 'event_elements':
                        event_data[col_name].append(self.get_event_elements(event['event_number']))
                    else:
                        event_data[col_name].append(event[col_name])


    
