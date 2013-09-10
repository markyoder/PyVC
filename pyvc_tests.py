#!/usr/bin/env python
from pyvc import vcplots
from pyvc import *

import time

sim_data_file = 'example_simulation.h5'

kwargs1 = {
    'event_range': {'type':'year','filter':(100,900)}
}

kwargs2 = {
    'event_range': {'type':'year','filter':(100,900)},
    'magnitude_filter': '>= 6.0',
}

kwargs3 = {
    'event_range': {'type':'year','filter':(100,900)},
    'magnitude_filter': '>= 6.0',
    'section_filter': {'type':'section_id','filter':[1,2,3,4,5,6,7,8,9,10]}
}

with VCSimData() as sim_data:

    sim_data.open_file(sim_data_file)
    events = VCEvents(sim_data)
    
    print 'Getting event data \'event_year\' and \'event_magnitude\''
    
    print '    No filters'
    start_time = time.time()
    event_data = events.get_event_data(['event_year', 'event_magnitude'])
    end_time = time.time() - start_time
    print '        {} events'.format(len(event_data['event_year']))
    print '        {} seconds'.format(time.time() - start_time)
    print
    
    print '    Years 100 - 900'
    start_time = time.time()
    event_data = events.get_event_data(['event_year', 'event_magnitude'], **kwargs1)
    end_time = time.time() - start_time
    print '        {} events'.format(len(event_data['event_year']))
    print '        {} seconds'.format(time.time() - start_time)
    print
    
    print '    Years 100 - 900'
    print '    Magnitudes >= 6.0'
    start_time = time.time()
    event_data = events.get_event_data(['event_year', 'event_magnitude'], **kwargs2)
    end_time = time.time() - start_time
    print '        {} events'.format(len(event_data['event_year']))
    print '        {} seconds'.format(time.time() - start_time)
    print
    
    print '    Years 100 - 900'
    print '    Magnitudes >= 6.0'
    print '    Sections 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10'
    start_time = time.time()
    event_data = events.get_event_data(['event_year', 'event_magnitude'], **kwargs3)
    end_time = time.time() - start_time
    print '        {} events'.format(len(event_data['event_year']))
    print '        {} seconds'.format(time.time() - start_time)
    print
    print

print 'Plotting event data'

print '    No filters'
print '        Average Slip vs Surface Rupture Length'
out_file = 'test_asrl.png'
start_time = time.time()
vcplots.average_slip_surface_rupture_length(sim_data_file, out_file)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Rupture Length'
out_file = 'test_mra.png'
start_time = time.time()
vcplots.magnitude_rupture_area(sim_data_file, out_file)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Average Slip'
out_file = 'test_mas.png'
start_time = time.time()
vcplots.magnitude_average_slip(sim_data_file, out_file)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Frequency vs Magnitude'
out_file = 'test_fm.png'
start_time = time.time()
vcplots.frequency_magnitude(sim_data_file, out_file)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '    Years 100 - 900'
print '        Average Slip vs Surface Rupture Length'
out_file = 'test_asrl-year.png'
start_time = time.time()
vcplots.average_slip_surface_rupture_length(sim_data_file, out_file, **kwargs1)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Rupture Length'
out_file = 'test_mra-year.png'
start_time = time.time()
vcplots.magnitude_rupture_area(sim_data_file, out_file, **kwargs1)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Average Slip'
out_file = 'test_mas-year.png'
start_time = time.time()
vcplots.magnitude_average_slip(sim_data_file, out_file, **kwargs1)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Frequency vs Magnitude'
out_file = 'test_fm-year.png'
start_time = time.time()
vcplots.frequency_magnitude(sim_data_file, out_file, **kwargs1)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '    Years 100 - 900'
print '    Magnitudes >= 6.0'
print '        Average Slip vs Surface Rupture Length'
out_file = 'test_asrl-year_mag.png'
start_time = time.time()
vcplots.average_slip_surface_rupture_length(sim_data_file, out_file, **kwargs2)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Rupture Length'
out_file = 'test_mra-year_mag.png'
start_time = time.time()
vcplots.magnitude_rupture_area(sim_data_file, out_file, **kwargs2)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Average Slip'
out_file = 'test_mas-year_mag.png'
start_time = time.time()
vcplots.magnitude_average_slip(sim_data_file, out_file, **kwargs2)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Frequency vs Magnitude'
out_file = 'test_fm-year_mag.png'
start_time = time.time()
vcplots.frequency_magnitude(sim_data_file, out_file, **kwargs2)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '    Years 100 - 900'
print '    Magnitudes >= 6.0'
print '    Sections 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10'
print '        Average Slip vs Surface Rupture Length'
out_file = 'test_asrl-year_mag_sec.png'
start_time = time.time()
vcplots.average_slip_surface_rupture_length(sim_data_file, out_file, **kwargs3)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Rupture Length'
out_file = 'test_mra-year_mag_sec.png'
start_time = time.time()
vcplots.magnitude_rupture_area(sim_data_file, out_file, **kwargs3)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Magnitude vs Average Slip'
out_file = 'test_mas-year_mag_sec.png'
start_time = time.time()
vcplots.magnitude_average_slip(sim_data_file, out_file, **kwargs3)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

print '        Frequency vs Magnitude'
out_file = 'test_fm-year_mag_sec.png'
start_time = time.time()
vcplots.frequency_magnitude(sim_data_file, out_file, **kwargs3)
print '        plotted to {}'.format(out_file)
print '        {} seconds'.format(time.time() - start_time)
print

