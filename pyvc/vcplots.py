#!/usr/bin/env python
from pyvc import *
import time

def magnitude_rupture_area(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #Instantiate VCEvents, VCGeometry classes
    with VCSimData() as sim_data:

        sim_data.open_file(sim_file)
        
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        print geometry.element_area(709)
        
        #start_time = time.time()
        #event_data = events.get_event_data(['event_magnitude', 'event_area'], event_range=event_range, magnitude_filter=magnitude_filter)
        #total_time = time.time() - start_time
        #print len(event_data['event_magnitude']), total_time
        
        #print
        
        #for i, mag in enumerate(event_data['event_magnitude']):
        #    print mag, event_data['event_area'][i]
        
        
        '''
        #time complexity test
        test_dat = []
        
        for test_num in range(100):
            ele_ids = []
            
            for id in range(1508):
                ele_ids.append(id)
                start_time = time.time()
                area = geometry.total_area([492,100])
                total_time = time.time() - start_time
                try:
                    test_dat[len(ele_ids)].append(total_time)
                except IndexError:
                    test_dat.append( [total_time] )
                #print len(ele_ids), total_time


        for n, times in enumerate(test_dat):
            print n, sum(times)/float(len(times))
            
        '''

def magnitude_average_slip(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #Instantiate VCEvents, VCGeometry classes
    with VCSimData() as sim_data:

        sim_data.open_file(sim_file)
        
        events = VCEvents(sim_data)
        #geometry = VCGeometry(sim_data)
        
        start_time = time.time()
        event_data = events.get_event_data(['event_magnitude', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter)
        total_time = time.time() - start_time
        print len(event_data['event_magnitude']), total_time
        
        print
        
        for i, mag in enumerate(event_data['event_magnitude']):
            print mag, event_data['event_average_slip'][i]