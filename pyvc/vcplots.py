#!/usr/bin/env python
from pyvc import *
from pyvc import vcutils
from pyvc import vcplotutils
from pyvc import vcexceptions

import matplotlib.pyplot as mplt
import matplotlib.font_manager as mfont
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap

import numpy as np

import math
import cPickle
import time
from operator import itemgetter
import os
import sys
import gc
import itertools
import subprocess

import networkx as nx
from networkx.algorithms import bipartite

import quakelib

'''
import multiprocessing
import Queue
from PIL import Image
from mpl_toolkits.basemap import Basemap, maskoceans, interp
'''

def plot_forecast(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Get the data.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)

        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon

        fault_traces = geometry.get_fault_traces()


    intervals = np.array([   x - event_data['event_year'][n-1]
                    for n,x in enumerate(event_data['event_year'])
                    if n != 0
                ])
    
    #---------------------------------------------------------------------------
    # t vs. P(t).
    #---------------------------------------------------------------------------
    cumulative = {}
    
    cumulative['x'] = np.sort(intervals)
    cumulative['y'] = np.arange(float(intervals.size))/float(intervals.size)
    #cumulative['y'] = [float(n)/float(len(intervals)) for n,x in enumerate(cumulative['x'])]
    #mplt.plot(np.sort(intervals), [float(n)/float(len(intervals)) for n,x in enumerate(np.sort(intervals))])
    
    #---------------------------------------------------------------------------
    # t0 vs. P(t0 + dt, t0) for fixed dt.
    #---------------------------------------------------------------------------
    dt = 30
    
    conditional_dt_fixed = {'x':[],'y':[]}
    
    for t0 in np.sort(intervals):
        int_t0_dt = intervals[np.where( intervals > t0+dt)]
        int_t0 = intervals[np.where( intervals > t0)]
        
        if int_t0.size != 0:
            conditional_dt_fixed['x'].append(t0)
            conditional_dt_fixed['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))

    #---------------------------------------------------------------------------
    # t = t0 + dt vs. P(t, t0) for various t0.
    #---------------------------------------------------------------------------
    conditional = {}

    for t0 in range(0,175,25):
        conditional[t0] = {'x':[],'y':[]}
        int_t0 = intervals[np.where( intervals > t0)]
        if int_t0.size != 0:
            for dt in range(250):
                int_t0_dt = intervals[np.where( intervals > t0+dt)]
                conditional[t0]['x'].append(t0+dt)
                conditional[t0]['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))

        mplt.plot(conditional[t0]['x'],conditional[t0]['y'])


#-------------------------------------------------------------------------------
# event field animation
#-------------------------------------------------------------------------------
def event_field_animation(sim_file, output_directory, event_range,
    field_type='displacement', fringes=True, padding=0.08, cutoff=None,
    animation_target_length=60.0, animation_fps = 30.0, fade_seconds = 1.0,
    min_mag_marker = 6.5, force_plot=False):
    
    sys.stdout.write('Initializing animation :: ')
    sys.stdout.flush()
    
    # create the animation dir if needed
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    if not output_directory.endswith('/'):
        output_directory += '/'
    
    # create the animation subdirs if needed
    field_values_directory = '{}field_values/'.format(output_directory)
    frame_images_directory = '{}frame_images/'.format(output_directory)
    if not os.path.exists(field_values_directory):
        os.makedirs(field_values_directory)
    if not os.path.exists(frame_images_directory):
        os.makedirs(frame_images_directory)

    # animation properties

    #---------------------------------------------------------------------------
    # Open the data file. It needs to stay open while we do the animation.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # Get global information about the simulations geometry
        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon
        fault_traces = geometry.get_fault_traces()
        
        # Get event information
        event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], event_range=event_range)
        event_magnitudes = event_data['event_magnitude']
        event_years = event_data['event_year']
        event_numbers = event_data['event_number']
        current_year = start_year = math.floor(event_years[0])
        
        # These are the event years shifted so the begin at zero.
        _event_years = [y - start_year for y in event_years]
        
        # The large magnitudes and shifted years to be marked on the timeline.
        event_large_magnitudes = [m for m in event_magnitudes if m > min_mag_marker]
        _event_large_magnitude_years = [_event_years[i_m[0]] for i_m in enumerate(event_magnitudes) if i_m[1] >= min_mag_marker]
        event_large_magnitude_evnums = [event_numbers[i_m[0]] for i_m in enumerate(event_magnitudes) if i_m[1] > min_mag_marker]
        
        # Calculate the frames per year and the total number of frames
        total_years = math.ceil(event_years[-1]) - math.floor(event_years[0])
        fpy = math.ceil(animation_target_length*animation_fps/total_years)
        total_frames = int(fpy * total_years)
        
        # Instantiate the field and the plotter
        if field_type == 'displacement':
            EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
            EFP.calculate_look_angles(geometry[:])
        elif field_type == 'gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)

        #-----------------------------------------------------------------------
        # Find the biggest event and normalize based on these values.
        #-----------------------------------------------------------------------
        if field_type == 'displacement' and not fringes or field_type == 'gravity':
            sys.stdout.write('normalizing : ')
            sys.stdout.flush()
            max_mag_evnum = event_numbers[event_magnitudes.index(max(event_magnitudes))]
            field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, max_mag_evnum))
            if field_values_loaded:
                sys.stdout.write('event {} loaded : '.format(max_mag_evnum))
                sys.stdout.flush()
            if not field_values_loaded:
                sys.stdout.write('event {} processing : '.format(max_mag_evnum))
                sys.stdout.flush()
                event_element_slips = events.get_event_element_slips(max_mag_evnum)
                ele_getter = itemgetter(*event_element_slips.keys())
                event_element_data = ele_getter(geometry)
                if len(event_element_slips) == 1:
                    event_element_data = [event_element_data]
            
                sys.stdout.write('{} elements : '.format(len(event_element_slips)))
                sys.stdout.flush()
                
                EF.calculate_field_values(
                    event_element_data,
                    event_element_slips,
                    cutoff=cutoff,
                    save_file_prefix='{}{}_'.format(field_values_directory,max_mag_evnum)
                )
            EFP.set_field(EF)
            EFP.create_field_image()
        
        # Convert the fault traces to lat-lon
        fault_traces_latlon = {}
        for secid in fault_traces.iterkeys():
             fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(EF.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])

        # Grab all of the plot properties that we will need.
        # properties that are fringes dependent
        if fringes and field_type == 'displacement':
            cmap            = EFP.dmc['cmap_f']
            coastline_color = EFP.dmc['coastline_color_f']
            country_color   = EFP.dmc['country_color_f']
            state_color     = EFP.dmc['state_color_f']
            fault_color     = EFP.dmc['fault_color_f']
            map_tick_color  = EFP.dmc['map_tick_color_f']
            map_frame_color = EFP.dmc['map_frame_color_f']
            grid_color      = EFP.dmc['grid_color_f']
            cb_fontcolor    = EFP.dmc['cb_fontcolor_f']
        else:
            cmap            = EFP.dmc['cmap']
            coastline_color = EFP.dmc['coastline_color']
            country_color   = EFP.dmc['country_color']
            state_color     = EFP.dmc['state_color']
            fault_color     = EFP.dmc['fault_color']
            map_tick_color  = EFP.dmc['map_tick_color']
            map_frame_color = EFP.dmc['map_frame_color']
            grid_color      = EFP.dmc['grid_color']
            cb_fontcolor    = EFP.dmc['cb_fontcolor']
        
        # properties that are not fringes dependent
        boundary_width  = EFP.dmc['boundary_width']
        coastline_width = EFP.dmc['coastline_width']
        country_width   = EFP.dmc['country_width']
        state_width     = EFP.dmc['state_width']
        river_width     = EFP.dmc['river_width']
        fault_width     = EFP.dmc['fault_width']
        map_frame_width = EFP.dmc['map_frame_width']
        map_fontsize    = EFP.dmc['map_fontsize']
        arrow_inset     = EFP.dmc['arrow_inset']
        arrow_fontsize  = EFP.dmc['arrow_fontsize']
        cb_fontsize     = EFP.dmc['cb_fontsize']
        cb_height       = EFP.dmc['cb_height']
        cb_margin_t     = EFP.dmc['cb_margin_t']
        grid_width      = EFP.dmc['grid_width']
        num_grid_lines  = EFP.dmc['num_grid_lines']
        font            = EFP.dmc['font']
        font_bold       = EFP.dmc['font_bold']

        map_resolution  = EFP.dmc['map_resolution']
        map_projection  = EFP.dmc['map_projection']
        plot_resolution = EFP.dmc['plot_resolution']

        #animation specific properties
        progress_tick_color = 'k'
        progress_frame_color = 'k'
        progress_frame_width = 1
        progress_line_color = 'k'
        progress_line_width = 0
        progress_marker = '.'
        progress_marker_edge_width = 0
        progress_marker_size = 4
        progress_indicator_line_color = 'red'
        progress_indicator_linewidth = 2
        progress_indicator_fontsize = 10.0
        
        mag_color = 'k'
        current_mag_color = 'white'
        current_mag_facecolor = 'red'
        mag_facecolor = 'white'
        mag_linewidth = 0.5
        mag_fontsize = 12
        
        fm_fontsize = 9
        if field_type=='displacement':
            fm_label_color = 'white'
        else:
            fm_label_color = 'black'
            
        fm_frame_width = 1
        fm_frame_color = 'k'
        fm_line_color = '0.0'
        fm_line_width = 1
        
        mag_color_map = mcolor.LinearSegmentedColormap.from_list(
            'mag_color_map',
            [mag_color,current_mag_color],
            N=256,
            gamma=1.0
        )
        mag_line_colormap = mcolor.LinearSegmentedColormap.from_list(
            'mag_line_colormap',
            [progress_frame_color,progress_indicator_line_color],
            N=256,
            gamma=1.0
        )
        current_mag_face_colormap = mcolor.LinearSegmentedColormap.from_list(
            'current_mag_face_colormap',
            [mag_facecolor,current_mag_facecolor],
            N=256,
            gamma=1.0
        )

        # Set up the field values. In order to do smooth fading, we are gonna be
        # adding to these quantities when there are new events and subtracting a
        # bit each frame till we get back down to zero.
        EF.init_field(0.0)

        # We need to keep track of all of the magnitudes that have been plotted
        # so far for the FM plot. Also, set up some other properties of the FM
        # plot.
        cumulative_magnitudes = []
        fm_x_ticks = np.linspace(round(min(event_magnitudes)), round(max(event_magnitudes)), 5)
        fm_y_ticks = np.logspace(0, 3, 4)
        
        # Set up a cuple of dicts to maintain the state of the large magnitude
        # labels and the section lines. This will allow us to fade these items
        # out after an event.
        section_states = {}
        large_event_label_states = {}
        fm_alpha_state = 0.0
        
        sys.stdout.write('done\n')
        sys.stdout.flush()
        
        #-----------------------------------------------------------------------
        # Go through all of the frames.
        #-----------------------------------------------------------------------
        sys.stdout.write('Total frames : {}, Frames per year : {}\n'.format(total_frames, fpy))
        sys.stdout.flush()
        for the_frame in range(total_frames):
            year_frame = the_frame%fpy
            if year_frame == 0:
                current_year += 1
                events_this_year = events.get_event_data(
                    ['event_number', 'event_year', 'event_magnitude'],
                    event_range = {'type':'year', 'filter':(current_year - 1, current_year)}
                )

            progress_indicator_year = current_year - 1 - start_year + year_frame/fpy
            
            evnums_this_frame = []
            sids_this_frame = []
            for i, year in enumerate(events_this_year['event_year']):
                if math.modf(year)[0] <= float(year_frame+1)/fpy and math.modf(year)[0] > float(year_frame)/fpy:
                    evnums_this_frame.append(events_this_year['event_number'][i])
                    sids_this_frame.append(geometry.sections_with_elements(list(events.get_event_elements(events_this_year['event_number'][i]))))
                    cumulative_magnitudes.append(events_this_year['event_magnitude'][i])
            sids_this_frame = set( itertools.chain(*sids_this_frame) )

            sys.stdout.write('frame {} (year {}) of {} ({})\n'.format(the_frame, progress_indicator_year, total_frames, total_years))

            # Remove a fixed percentage from the field. This is the decay
            # that slowly fades existing field values.
            EF.shrink_field(1.0 - 1.0/(animation_fps*fade_seconds))
            
            #-------------------------------------------------------------------
            # Load or calculate all of the data for the current frame.
            #-------------------------------------------------------------------
            if len(evnums_this_frame) > 0:
                for i, evnum in enumerate(evnums_this_frame):
                    sys.stdout.write('\r Event {} :: '.format(evnum))
                    # Try and load the fields
                    field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, evnum))
                    if field_values_loaded:
                        sys.stdout.write('loaded'.format(evnum))
                    # If they havent been saved then we need to calculate them
                    elif not field_values_loaded:
                        sys.stdout.write('processing '.format(evnum))
                        sys.stdout.flush()
                        
                        event_element_slips = events.get_event_element_slips(evnum)
                        ele_getter = itemgetter(*event_element_slips.keys())
                        event_element_data = ele_getter(geometry)
                        if len(event_element_slips) == 1:
                            event_element_data = [event_element_data]
                    
                        sys.stdout.write('{} elements :: '.format(len(event_element_slips)))
                        sys.stdout.flush()
                        
                        EF.calculate_field_values(
                            event_element_data,
                            event_element_slips,
                            cutoff=cutoff,
                            save_file_prefix='{}{}_'.format(field_values_directory, evnum)
                        )
                        
                        if i < len(evnums_this_frame)-1 :
                            sys.stdout.write('\033[2K')
                        sys.stdout.flush()

            sys.stdout.write('\n')
            
            #-------------------------------------------------------------------
            # State variables that need to be maintained pre-plotting.
            #-------------------------------------------------------------------
            
            # Big magnitude event label fade state
            for elme in event_large_magnitude_evnums:
                if elme in evnums_this_frame:
                    large_event_label_states[elme] = 1.0
            
            # Active section line width
            for sid in sids_this_frame:
                section_states[sid] = 1.0
            
            #-------------------------------------------------------------------
            # If the image has not been plotted, plot it.
            #-------------------------------------------------------------------
            if not os.path.isfile('{}{}.png'.format(frame_images_directory, the_frame)) or force_plot:
                sys.stdout.write('\r Plotting :: ')
                
                # Set the plot field to be the current field
                EFP.set_field(EF)
                
                sys.stdout.write('creating the map image : ')
                sys.stdout.flush()
                # Create the field map image
                map_image = EFP.create_field_image(fringes=fringes)
                
                sys.stdout.write('finishing plot')
                sys.stdout.flush()
                #---------------------------------------------------------------
                # Plot all of the geographic data, animation timeline etc.
                #---------------------------------------------------------------
                
                # Width and height are fixed
                ph = 768.0
                pw = 1024.0
                
                # Map margin left and right
                mml = 70.0
                mmb = 70.0
                
                # Progress indicator map margin
                pimm = 70.0
                # Progress indicator margin r
                pimr = 50.0
                
                mw = EF.lons_1d.size
                mh = EF.lats_1d.size
                
                width_frac = mw/pw
                height_frac = mh/ph
                left_frac = mml/pw
                bottom_frac = mmb/ph

                pwi = pw/plot_resolution
                phi = ph/plot_resolution
                fig4 = mplt.figure(figsize=(pwi, phi), dpi=plot_resolution)

                #---------------------------------------------------------------
                # m4, fig4 is all of the boundary data.
                #---------------------------------------------------------------
                m4 = Basemap(
                    llcrnrlon=EF.min_lon,
                    llcrnrlat=EF.min_lat,
                    urcrnrlon=EF.max_lon,
                    urcrnrlat=EF.max_lat,
                    lat_0=(EF.max_lat+EF.min_lat)/2.0,
                    lon_0=(EF.max_lon+EF.min_lon)/2.0,
                    resolution=map_resolution,
                    projection=map_projection,
                    suppress_ticks=True
                )
                m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                
                # Draw a frame around the map.
                m4.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=(1,1,1,0))

                # draw coastlines, edge of map.
                m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)

                # draw countries
                m4.drawcountries(linewidth=country_width, color=country_color)

                # draw states
                m4.drawstates(linewidth=state_width, color=state_color)

                # draw parallels.
                parallels = np.linspace(EFP.lats_1d.min(), EFP.lats_1d.max(), num_grid_lines+1)
                m4_parallels = m4.drawparallels(
                    parallels,
                    labels=[1,0,0,0],
                    fontsize=map_fontsize,
                    color=grid_color,
                    fontproperties=font,
                    fmt='%.2f',
                    linewidth=grid_width,
                    dashes=[1, 10]
                )

                # draw meridians
                meridians = np.linspace(EFP.lons_1d.min(), EFP.lons_1d.max(), num_grid_lines+1)
                m4_meridians = m4.drawmeridians(
                    meridians,
                    labels=[0,0,1,0],
                    fontsize=map_fontsize,
                    color=grid_color,
                    fontproperties=font,
                    fmt='%.2f',
                    linewidth=grid_width,
                    dashes=[1, 10]
                )

                # add the displacement map image to the plot
                m4.imshow(map_image, origin='upper')
                
                #---------------------------------------------------------------
                # Plot the magnitude/progress indicator.
                #---------------------------------------------------------------
                width_frac = (pw - mw - mml - pimm - pimr) /pw
                height_frac = mh/ph
                left_frac = (mw + mml + pimm)/pw
                bottom_frac = mmb/ph
                
                mag_vs_year = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                mag_vs_year.plot(event_magnitudes, _event_years, color=progress_line_color, linewidth=progress_line_width, marker=progress_marker, mew=progress_marker_edge_width, ms=progress_marker_size, mfc=progress_line_color)
                mag_vs_year.autoscale(enable=True, axis='both', tight=True)
                mag_vs_year.set_ylim(math.floor(min(_event_years)), math.ceil(max(_event_years)))
                mag_vs_year_alt = mag_vs_year.twinx()
                mag_vs_year.set_xticks(np.linspace(min(event_magnitudes),max(event_magnitudes),3))
                
                mag_vs_year_alt.set_ylim(math.floor(min(_event_years)), math.ceil(max(_event_years)))
                
                mag_vs_year_tick_labels = ['{:0.1f}'.format(tick) for tick in np.linspace(min(event_magnitudes),max(event_magnitudes),3)]
                mag_vs_year.set_xticklabels(mag_vs_year_tick_labels)
                
                for tick in mag_vs_year.xaxis.get_major_ticks():
                    tick.label1.set_fontproperties(font)
                    tick.label1.set_fontsize(progress_indicator_fontsize)
                    tick.label1.set_color(progress_tick_color)
                
                for tl in mag_vs_year_alt.get_yticklabels():
                    tl.set_fontproperties(font)
                    tl.set_fontsize(progress_indicator_fontsize)
                    tl.set_color(progress_tick_color)
                
                for tick in mag_vs_year.yaxis.get_major_ticks():
                    tick.label1.set_alpha(0)
                
                for line in mag_vs_year.xaxis.get_ticklines() + mag_vs_year_alt.yaxis.get_ticklines() + mag_vs_year.yaxis.get_ticklines():
                    line.set_alpha(0)

                #take care of all of the frame line widths
                for spine in mag_vs_year.spines.itervalues():
                    spine.set_lw(progress_frame_width)
                    spine.set_color(progress_frame_color)
                
                mag_vs_year.set_xlabel('magnitude', fontproperties=font, size=progress_indicator_fontsize, color=progress_tick_color)
                mag_vs_year_alt.set_ylabel('year', fontproperties=font, size=progress_indicator_fontsize, color=progress_tick_color)

                #---------------------------------------------------------------
                # Add the progress indicator line
                #---------------------------------------------------------------
                mag_vs_year.axhline(y=progress_indicator_year, lw=progress_indicator_linewidth, c=progress_indicator_line_color)
            
                #---------------------------------------------------------------
                # Add the progress indicator label lines for large events.
                #---------------------------------------------------------------
                label_lines = []
                #if len(event_large_magnitudes) < 10:
                #    label_range = _event_large_magnitude_years
                #else:
                label_range = np.linspace(math.floor(min(_event_years)),math.ceil(max(_event_years)),len(event_large_magnitudes))
                for i, y in enumerate(label_range):
                    m = event_large_magnitudes[i]
                    y1 = _event_large_magnitude_years[i]
                    
                    try:
                        lels = large_event_label_states[event_large_magnitude_evnums[i]]
                    except KeyError:
                        lels = 0

                    the_color = mag_color_map(lels)
                    the_line_color = mag_line_colormap(lels)
                    the_line_width = mag_linewidth + lels * (progress_indicator_linewidth)
                    the_bbox = dict(
                        facecolor=current_mag_face_colormap(lels),
                        linewidth=the_line_width,
                        ec=the_line_color,
                        boxstyle='round4,pad=0.5'
                    )
                    mag_vs_year.text(
                        min(event_magnitudes) - 0.4, y, '{:0.1f}'.format(m),
                        fontproperties=font,
                        fontsize=mag_fontsize,
                        horizontalalignment='right',
                        verticalalignment='center',
                        color=the_color,
                        bbox=the_bbox
                    )
                    
                    label_lines.append(
                        mlines.Line2D(
                            [1.0,-0.01,-0.05,-0.085],
                            [y1/math.ceil(max(_event_years)),
                                y1/math.ceil(max(_event_years)),
                                y/math.ceil(max(_event_years)),
                                y/math.ceil(max(_event_years))
                            ],
                            linewidth=the_line_width,
                            transform=mag_vs_year.transAxes,
                            color=the_line_color,
                            solid_capstyle='round',
                            solid_joinstyle='round'
                        )
                    )
                
                mag_vs_year.lines.extend(label_lines)
                
                #---------------------------------------------------------------
                # Plot the current frequency-magnitude
                #---------------------------------------------------------------
                width_frac = 150.0/pw
                height_frac = 150.0/ph
                left_frac = 100.0/pw
                bottom_frac = 100.0/ph
                
                current_total_events = len(cumulative_magnitudes)
                if current_total_events > 1:
                    cum_freq = {}
                    fm_x = []
                    fm_y = []
                    for num, magnitude in enumerate(sorted(cumulative_magnitudes)):
                        cum_freq['{:0.10f}'.format(magnitude)] = current_total_events - (num + 1)
                    
                    for magnitude in sorted(cum_freq.iterkeys()):
                        fm_x.append(magnitude)
                        fm_y.append(float(cum_freq[magnitude]))
                    mag_vs_freq = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                    mag_vs_freq.semilogy(fm_x, fm_y, color=fm_line_color, linewidth=fm_line_width)
                    mag_vs_freq.set_ylim(bottom=min(fm_y_ticks), top=max(fm_y_ticks))
                    mag_vs_freq.set_xlim(left=round(min(event_magnitudes)), right=round(max(event_magnitudes)))
                    mag_vs_freq.set_yticks(fm_y_ticks)
                    mag_vs_freq.set_xticks(fm_x_ticks)
                    
                    mag_vs_freq_x_ticklabels = ['{:0.1f}'.format(tick) for tick in fm_x_ticks]
                    mag_vs_freq.set_xticklabels(mag_vs_freq_x_ticklabels)
                    
                    for tick in mag_vs_freq.xaxis.get_major_ticks() + mag_vs_freq.yaxis.get_major_ticks():
                        tick.label1.set_fontproperties(font)
                        tick.label1.set_fontsize(fm_fontsize)
                        tick.label1.set_color(fm_label_color)
                        tick.label1.set_alpha(fm_alpha_state)
                    
                    for line in mag_vs_freq.xaxis.get_majorticklines() + mag_vs_freq.yaxis.get_majorticklines() + mag_vs_freq.yaxis.get_minorticklines():
                        line.set_alpha(0)
                    
                    #take care of all of the frame line widths
                    for spine in mag_vs_freq.spines.itervalues():
                        spine.set_lw(fm_frame_width)
                        spine.set_color(fm_frame_color)
                        spine.set_alpha(fm_alpha_state)
                    mag_vs_freq.patch.set_alpha(fm_alpha_state)

                #---------------------------------------------------------------
                # Plot the fault traces.
                #---------------------------------------------------------------
                # Plot the sections
                for sid, sec_trace in fault_traces_latlon.iteritems():
                    trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
                    
                    try:
                        section_state = section_states[sid]
                    except KeyError:
                        section_state = 0.0

                    m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width + section_state * 2.0, solid_capstyle='round', solid_joinstyle='round')

                #---------------------------------------------------------------
                # Plot the colorbar
                #---------------------------------------------------------------
                left_frac = 70.0/pw
                bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
                width_frac = mw/pw
                height_frac = cb_height/ph
                
                cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                norm = EFP.norm
                cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
                       norm=norm,
                       orientation='horizontal')
                if field_type == 'displacement':
                    if fringes:
                        cb_title = 'Displacement [m]'
                    else:
                        cb_title = 'Total displacement [m]'

                elif field_type == 'gravity':
                    cb_title = r'Gravity changes [$\mu gal$]'

                cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

                for label in cb_ax.xaxis.get_ticklabels():
                    label.set_fontproperties(font)
                    label.set_fontsize(cb_fontsize)
                    label.set_color(cb_fontcolor)
                                      
                for line in cb_ax.xaxis.get_ticklines():
                    line.set_alpha(0)
            
                #---------------------------------------------------------------
                # If the field is a displacement field, draw the look arrows.
                #---------------------------------------------------------------
                if field_type == 'displacement':
                    # draw the azimuth look arrow
                    az_width_frac    = 50.0/pw
                    az_height_frac   = 50.0/ph
                    az_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
                    az_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac)/ph
                    az_ax = fig4.add_axes((az_left_frac,az_bottom_frac,az_width_frac,az_height_frac))

                    az_ax.set_xlim((0,1.0))
                    az_ax.set_ylim((0,1.0))
                    for item in az_ax.yaxis.get_ticklabels() + az_ax.xaxis.get_ticklabels() + az_ax.yaxis.get_ticklines() + az_ax.xaxis.get_ticklines():
                        item.set_alpha(0)

                    az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(EFP.look_azimuth)
                    az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(EFP.look_azimuth)
                    az_arrow_dx      = 0.8*math.sin(EFP.look_azimuth)
                    az_arrow_dy      = 0.8*math.cos(EFP.look_azimuth)

                    az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
                    az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                    az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - EF.convert.rad2deg(EFP.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
                    az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

                    # draw the altitude look arrow
                    al_width_frac    = 50.0/pw
                    al_height_frac   = 50.0/ph
                    al_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
                    al_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac - ph*al_height_frac)/ph
                    al_ax = fig4.add_axes((al_left_frac,al_bottom_frac,al_width_frac,al_height_frac))

                    al_ax.set_xlim((0,1.0))
                    al_ax.set_ylim((0,1.0))
                    for item in al_ax.yaxis.get_ticklabels() + al_ax.xaxis.get_ticklabels() + al_ax.yaxis.get_ticklines() + al_ax.xaxis.get_ticklines():
                        item.set_alpha(0)

                    al_arrow_start_x    = 0.1 + 0.8*math.cos(EFP.look_elevation)
                    al_arrow_start_y    = 0.1 + 0.8*math.sin(EFP.look_elevation)
                    al_arrow_dx      = -0.8*math.cos(EFP.look_elevation)
                    al_arrow_dy      = -0.8*math.sin(EFP.look_elevation)

                    al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
                    al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                    al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=EF.convert.rad2deg(EFP.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
                    al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
                #---------------------------------------------------------------
                # If the field is a gravity field change outermost tick labels
                # on colorbar.
                #---------------------------------------------------------------
                else:
                    # Want to change outermost tick labels on colorbar
                    #   from 'VALUE','-VALUE' to '>VALUE' and '<-VALUE'    
                    cb_tick_labs = [item.get_text() for item in cb_ax.get_xticklabels()]
                    cb_tick_labs[0] = '<'+cb_tick_labs[0]
                    cb_tick_labs[-1] = '>'+cb_tick_labs[-1]
                    cb_ax.set_xticklabels(cb_tick_labs)
                
                #---------------------------------------------------------------
                # Save the figure and clear out all matplotlib figures.
                #---------------------------------------------------------------
                fig4.savefig('{}{}.png'.format(frame_images_directory, the_frame), format='png', dpi=plot_resolution)
                
                fig4.clf()
                mplt.close('all')
                gc.collect()

                sys.stdout.write('\033[2K')
                sys.stdout.flush()
            
            #-------------------------------------------------------------------
            # State variables that need to be maintained post-plotting.
            #-------------------------------------------------------------------
            
            # Big magnitude event label fade state
            for k, lels in large_event_label_states.iteritems():
                if lels >= 0.04:
                    large_event_label_states[k] -= 0.04
                else:
                    large_event_label_states[k] = 0.0
            
            # Active section line width
            for sid, section_state in section_states.iteritems():
                if section_state >= 0.04:
                    section_states[sid] -= 0.04
                else:
                    section_states[sid] = 0.0

            # Frequency magnitude plot alpha
            fm_alpha_state += 0.04
            if fm_alpha_state > 1.0:
                fm_alpha_state = 1.0

            sys.stdout.write('\033[1A')
            sys.stdout.write('\033[2K')
            sys.stdout.write('\033[1A\r')
            sys.stdout.write('\033[2K')
        
        #-----------------------------------------------------------------------
        # Create the movie using ffmpeg.
        #-----------------------------------------------------------------------
        
        #ffmpeg -y -r 15 -i f_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p animation.mp4
        
        proc_args = 'ffmpeg -y -r {fps} -start_number 0 -i {dir}{inc}.png -f mp4 -vcodec h264 -pix_fmt yuv420p {out}animation.mp4'.format(
            fps=int(animation_fps),
            dir=frame_images_directory,
            inc='%d',
            out=output_directory
        )
        proc = subprocess.Popen(proc_args, shell=True)
        proc.wait()

#-------------------------------------------------------------------------------
# plots event fields
#-------------------------------------------------------------------------------
def plot_event_field(sim_file, evnum, output_file=None, field_type='displacement', fringes=True, padding=0.08, cutoff=None, save_file_prefix=None):
    
    sys.stdout.write('Initializing plot :: ')
    sys.stdout.flush()
    
    start_time = time.time()
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon

        event_data = events[evnum]
        event_element_slips = events.get_event_element_slips(evnum)
        ele_getter = itemgetter(*event_element_slips.keys())
        event_element_data = ele_getter(geometry)
        if len(event_element_slips) == 1:
            event_element_data = [event_element_data]
        fault_traces = geometry.get_fault_traces()
        event_sections = geometry.sections_with_elements(event_element_slips.keys())
    
    sys.stdout.write('{} elements in {} sections : '.format(len(event_element_slips), len(event_sections)))
    sys.stdout.flush()
    
    sys.stdout.write( '{} field : '.format(field_type))
    sys.stdout.flush()
    
    if field_type == 'displacement':
        EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
    elif field_type == 'gravity':
        EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)

    sys.stdout.write('done\n')
    sys.stdout.flush()

    field_values_loaded = EF.load_field_values(save_file_prefix)
 
    if field_values_loaded:
        sys.stdout.write('Loading event {} {} field :: '.format(evnum, field_type))
        sys.stdout.flush()
    if not field_values_loaded:
        sys.stdout.write('Processing event {} {} field :: '.format(evnum, field_type))
        sys.stdout.flush()          
        sys.stdout.write('{} elements : '.format(len(event_element_slips)))
        sys.stdout.flush()
        
        EF.calculate_field_values(
                    event_element_data,
                    event_element_slips,
                    cutoff=cutoff,
                    save_file_prefix=save_file_prefix)
                

    '''
    if field_type == 'displacement':
        np.save('local/dX.npy', EF.dX)
        np.save('local/dY.npy', EF.dY)
        np.save('local/dZ.npy', EF.dZ)
        #EF.dX = np.load('local/dX.npy')
        #EF.dY = np.load('local/dY.npy')
        #EF.dZ = np.load('local/dZ.npy')
    elif field_type == 'gravity':
        np.save('local/dG.npy', EF.dG)
        #EF.dG = np.load('local/dG.npy')
    '''
    sys.stdout.write('done\n')
    sys.stdout.flush()

    sys.stdout.write('Plotting :: initializing : ')
    sys.stdout.flush()

    if field_type == 'displacement':
        EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
    elif field_type == 'gravity':
        EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
    EFP.set_field(EF)

    if field_type == 'displacement':
        start_time = time.time()
        EFP.calculate_look_angles(event_element_data)

    sys.stdout.write('map image : ')
    sys.stdout.flush()

    map_image = EFP.create_field_image(fringes=fringes)

    sys.stdout.write('map overlay : ')
    sys.stdout.flush()
    # Convert the fault traces to lat-lon
    fault_traces_latlon = {}
    for secid in fault_traces.iterkeys():
         fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(EF.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])

    #---------------------------------------------------------------------------
    # Plot all of the geographic info on top of the displacement map image.
    #---------------------------------------------------------------------------
    
    # Grab all of the plot properties that we will need.
    # properties that are fringes dependent
    if fringes and field_type == 'displacement':
        cmap            = EFP.dmc['cmap_f']
        coastline_color = EFP.dmc['coastline_color_f']
        country_color   = EFP.dmc['country_color_f']
        state_color     = EFP.dmc['state_color_f']
        fault_color     = EFP.dmc['fault_color_f']
        map_tick_color  = EFP.dmc['map_tick_color_f']
        map_frame_color = EFP.dmc['map_frame_color_f']
        grid_color      = EFP.dmc['grid_color_f']
        cb_fontcolor    = EFP.dmc['cb_fontcolor_f']
    else:
        cmap            = EFP.dmc['cmap']
        coastline_color = EFP.dmc['coastline_color']
        country_color   = EFP.dmc['country_color']
        state_color     = EFP.dmc['state_color']
        fault_color     = EFP.dmc['fault_color']
        map_tick_color  = EFP.dmc['map_tick_color']
        map_frame_color = EFP.dmc['map_frame_color']
        grid_color      = EFP.dmc['grid_color']
        cb_fontcolor    = EFP.dmc['cb_fontcolor']
    
    # properties that are not fringes dependent
    boundary_width  = EFP.dmc['boundary_width']
    coastline_width = EFP.dmc['coastline_width']
    country_width   = EFP.dmc['country_width']
    state_width     = EFP.dmc['state_width']
    river_width     = EFP.dmc['river_width']
    fault_width     = EFP.dmc['fault_width']
    map_frame_width = EFP.dmc['map_frame_width']
    map_fontsize    = EFP.dmc['map_fontsize']
    arrow_inset     = EFP.dmc['arrow_inset']
    arrow_fontsize  = EFP.dmc['arrow_fontsize']
    cb_fontsize     = EFP.dmc['cb_fontsize']
    cb_height       = EFP.dmc['cb_height']
    cb_margin_t     = EFP.dmc['cb_margin_t']
    grid_width      = EFP.dmc['grid_width']
    num_grid_lines  = EFP.dmc['num_grid_lines']
    font            = EFP.dmc['font']
    font_bold       = EFP.dmc['font_bold']

    map_resolution  = EFP.dmc['map_resolution']
    map_projection  = EFP.dmc['map_projection']
    plot_resolution = EFP.dmc['plot_resolution']

    # The sizing for the image is tricky. The aspect ratio of the plot is fixed,
    # so we cant set all of margins to whatever we want. We will set the anchor
    # to the top, left margin position. Then scale the image based on the
    # bottom/right margin, whichever is bigger.
    
    mw = EF.lons_1d.size
    mh = EF.lats_1d.size

    if mh > mw:
        ph = 768.0
        pw = mw + 70.0 + 40.0
    else:
        pw = 790.0
        ph = mh + 70.0 + 40.0

    width_frac = mw/pw
    height_frac = mh/ph
    left_frac = 70.0/pw
    bottom_frac = 70.0/ph

    pwi = pw/plot_resolution
    phi = ph/plot_resolution

    fig4 = mplt.figure(figsize=(pwi, phi), dpi=plot_resolution)

    #---------------------------------------------------------------------------
    # m4, fig4 is all of the boundary data.
    #---------------------------------------------------------------------------
    m4 = Basemap(
        llcrnrlon=EF.min_lon,
        llcrnrlat=EF.min_lat,
        urcrnrlon=EF.max_lon,
        urcrnrlat=EF.max_lat,
        lat_0=(EF.max_lat+EF.min_lat)/2.0,
        lon_0=(EF.max_lon+EF.min_lon)/2.0,
        resolution=map_resolution,
        projection=map_projection,
        suppress_ticks=True
    )
    m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
    
    # draw coastlines, edge of map.
    m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)
    
    # draw countries
    m4.drawcountries(linewidth=country_width, color=country_color)
    
    # draw states
    m4.drawstates(linewidth=state_width, color=state_color)
    
    # draw parallels.
    parallels = np.linspace(EFP.lats_1d.min(), EFP.lats_1d.max(), num_grid_lines+1)
    m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], fontsize=map_fontsize, color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
    
    # draw meridians
    meridians = np.linspace(EFP.lons_1d.min(), EFP.lons_1d.max(), num_grid_lines+1)
    m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], fontsize=map_fontsize, color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])

    if field_type == 'displacement':
        # draw the azimuth look arrow
        az_width_frac    = 50.0/pw
        az_height_frac   = 50.0/ph
        az_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
        az_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac)/ph
        az_ax = fig4.add_axes((az_left_frac,az_bottom_frac,az_width_frac,az_height_frac))

        az_ax.set_xlim((0,1.0))
        az_ax.set_ylim((0,1.0))
        for item in az_ax.yaxis.get_ticklabels() + az_ax.xaxis.get_ticklabels() + az_ax.yaxis.get_ticklines() + az_ax.xaxis.get_ticklines():
            item.set_alpha(0)

        az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(EFP.look_azimuth)
        az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(EFP.look_azimuth)
        az_arrow_dx      = 0.8*math.sin(EFP.look_azimuth)
        az_arrow_dy      = 0.8*math.cos(EFP.look_azimuth)

        az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
        az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - EF.convert.rad2deg(EFP.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
        az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

        # draw the altitude look arrow
        al_width_frac    = 50.0/pw
        al_height_frac   = 50.0/ph
        al_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
        al_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac - ph*al_height_frac)/ph
        al_ax = fig4.add_axes((al_left_frac,al_bottom_frac,al_width_frac,al_height_frac))

        al_ax.set_xlim((0,1.0))
        al_ax.set_ylim((0,1.0))
        for item in al_ax.yaxis.get_ticklabels() + al_ax.xaxis.get_ticklabels() + al_ax.yaxis.get_ticklines() + al_ax.xaxis.get_ticklines():
            item.set_alpha(0)

        al_arrow_start_x    = 0.1 + 0.8*math.cos(EFP.look_elevation)
        al_arrow_start_y    = 0.1 + 0.8*math.sin(EFP.look_elevation)
        al_arrow_dx      = -0.8*math.cos(EFP.look_elevation)
        al_arrow_dy      = -0.8*math.sin(EFP.look_elevation)

        al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
        al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=EF.convert.rad2deg(EFP.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
        al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
        
        # draw the box with the magnitude
        mag_width_frac    = 50.0/pw
        mag_height_frac   = 10.0/ph
        mag_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
        mag_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac  - ph*az_height_frac - ph*mag_height_frac)/ph
        mag_ax = fig4.add_axes((mag_left_frac,mag_bottom_frac,mag_width_frac,mag_height_frac))

        mag_ax.set_xlim((0,1.0))
        mag_ax.set_ylim((0,1.0))
        for item in mag_ax.yaxis.get_ticklabels() + mag_ax.xaxis.get_ticklabels() + mag_ax.yaxis.get_ticklines() + mag_ax.xaxis.get_ticklines():
            item.set_alpha(0)
        
        mag_ax.text(0.5, 0.5, 'm = {:0.3f}'.format(float(event_data['event_magnitude'])), fontproperties=font_bold, size=arrow_fontsize, ha='center', va='center')

    # add the displacement map image to the plot
    m4.imshow(map_image, origin='upper')
    
    # print faults on lon-lat plot
    for sid, sec_trace in fault_traces_latlon.iteritems():
        trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
        
        if sid in event_sections:
            linewidth = fault_width + 3
        else:
            linewidth = fault_width

        m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')

    #plot the cb
    left_frac = 70.0/pw
    bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
    width_frac = mw/pw
    height_frac = cb_height/ph
    
    cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
    norm = EFP.norm
    cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
           norm=norm,
           orientation='horizontal')
    if field_type == 'displacement':
        if fringes:
            cb_title = 'Displacement [m]'
        else:
            cb_title = 'Total displacement [m]'

    elif field_type == 'gravity':
        cb_title        = r'Gravity changes [$\mu gal$]'
        # Make first and last ticks on colorbar be <MIN and >MAX
        cb_tick_labs    = [item.get_text() for item in cb_ax.get_xticklabels()]
        cb_tick_labs[0] = '<'+cb_tick_labs[0]
        cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
        cb_ax.set_xticklabels(cb_tick_labs)

    cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

    for label in cb_ax.xaxis.get_ticklabels():
        label.set_fontproperties(font)
        label.set_fontsize(cb_fontsize)
        label.set_color(cb_fontcolor)
    for line in cb_ax.xaxis.get_ticklines():
        line.set_alpha(0)

    if output_file is not None:
        # save the figure
        fig4.savefig(output_file, format='png', dpi=plot_resolution)

    sys.stdout.write('done\n')
    sys.stdout.flush()

#-------------------------------------------------------------------------------
# plots recurrence intervals
#-------------------------------------------------------------------------------
def plot_recurrence_intervals(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Plot setup
    #---------------------------------------------------------------------------
    
    num_cols = 5.0
    
    # dimensions
    simw = 270.0
    simh = 270.0
    stm = 40.0
    sbm = 40.0
    slm = 50.0
    srm = 10.0
    res = 72.0
    
    # fonts
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    titlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    subtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=8)
    
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # get the data
        section_info = geometry.get_section_info(section_filter=section_filter)
        
        #-----------------------------------------------------------------------
        # start the plot
        #-----------------------------------------------------------------------
        bins = np.linspace(0,250,50)
        # calculate the final dimensions and create the figure and axis
        spw = simw - slm - srm
        sph = simh - stm - sbm
        num_rows = math.ceil(float(len(section_info))/num_cols)
        imw = math.ceil(simw * num_cols)
        imh = math.ceil(simh * num_rows)
        imwi = imw/res
        imhi = imh/res
        fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
        
        #-----------------------------------------------------------------------
        # Calculate the recurrence intervals and plot.
        #-----------------------------------------------------------------------
        curr_row = -1.0
        for num, secid in enumerate(sorted(section_info.keys())):
            curr_col = num%num_cols
            if curr_col == 0.0:
                curr_row += 1.0
            #print curr_row, curr_col
            the_ax = fig.add_axes(((slm + curr_col * simw)/imw, (sbm + (num_rows - curr_row - 1) * simh)/imh, spw/imw, sph/imh))
            section_events = events.get_event_data_from_evnums(
                                        geometry.events_on_section(secid),
                                        ['event_magnitude', 'event_year'],
                                        event_range=event_range,
                                        magnitude_filter='>=6.5'
                                    )
            intervals = [
                x - section_events['event_year'][n-1]
                for n,x in enumerate(section_events['event_year'])
                if n != 0]
                
            intervals7 = [
                x - section_events['event_year'][n-1]
                for n,x in enumerate(section_events['event_year'])
                if n != 0 and section_events['event_magnitude'][n] >= 7.0]
            
            
            hist, bins = np.histogram(intervals, bins=bins, density=True)
            hist7, bins7 = np.histogram(intervals7, bins=bins, density=True)
            mean = np.mean(intervals)
            std = np.std(intervals)
            mean7 = np.mean(intervals7)
            std7 = np.std(intervals7)
            
            the_ax.step(bins[0:-1], hist, where='post', label='m>6.5')
            the_ax.step(bins7[0:-1], hist7, where='post', label='m>7')
            
            for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
                label.set_fontproperties(ticklabelfont)
                
            the_ax.set_ylabel('Prob. Density', fontproperties=framelabelfont)
            the_ax.set_xlabel('Recurrence Time [yr]', fontproperties=framelabelfont)
            
            the_ax.autoscale_view(tight=True)
            
            the_ax.set_title('{} {}'.format(secid,section_info[secid]['name']), position=(0.0,1.04), ha='left', fontproperties=titlefont)
            the_ax.text(0.0, 1.01, 'm>6.5: mean {mean:0.1f} std {std:0.1f}, m>7.0 mean {mean7:0.1f} std {std7:0.1f}'.format(mean=mean, std=std, mean7=mean7, std7=std7), va='bottom', ha='left', transform=the_ax.transAxes, fontproperties=subtitlefont)
    
            the_ax.legend(prop=legendfont)

    if output_file is not None:
        # Get the plot format and save the file
        plot_format = output_file.split('.')[-1]
        if plot_format != 'png' and plot_format != 'pdf':
            raise vcexceptions.PlotFormatNotSupported(plot_format)
        else:
            fig.savefig(output_file, format=plot_format, dpi=res)

#-------------------------------------------------------------------------------
# plots a matrix of a bipartite event sequence graph
#-------------------------------------------------------------------------------
def plot_bipartite_graph_matrix(graph_file, output_file=None):
    G = cPickle.load(open(graph_file, 'rb'))
    
    top_nodes = set(n for n,d in G.nodes(data=True) if d['bipartite']==0)
    bottom_nodes = [n for n in sorted(set(G) - top_nodes, key=lambda n: G.degree(n), reverse=True)]
            
    print 'calculating matrix'
    bij_matrix = bipartite.biadjacency_matrix(G, bottom_nodes, top_nodes, weight='weight')
    
    print bij_matrix.shape
    # Aspect rations are always h/w
    aspect = float(bij_matrix.shape[1])/float(bij_matrix.shape[0])
    
    print 'plotting matrix', len(top_nodes), len(bottom_nodes), bij_matrix.max()
    #print pos_sid_prob
    
    # plot parameters
    imh = 178.0        # the full image height
    imw = round(imh/aspect)    # the width is set based on the aspect ratio
    lm = 10.0
    rm = 10.0
    tm = 10.0
    bm = 20.0
    res = 72.0
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=6)
    cmap = mplt.get_cmap('hot_r')
    
    #arial14 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    #arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    #arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    #arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm # the height for both matricies
    pw = imw - lm - rm
    print imh, imw, lm/imw, bm/imh, pw/imw, ph/imh
    ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
    
    #print pos_sid_prob[::-1]
    
    ax.imshow(bij_matrix.T, interpolation='none', cmap=cmap)
    #ax.set_xticks(range(len(pos_sid_prob)))
    #ax.set_yticks(range(len(pos_sid_prob)))
    #ax.set_xticklabels(pos_sid_prob)
    #ax.set_yticklabels(pos_sid_prob)
    #for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
    #    label.set_fontproperties(ticklabelfont)

    #ax.axis('tight')
    #ax.set_ylim((15.5, -0.5))
    #ax.set_xlim((140.5, len(pos_sid_prob)-0.5))

    if output_file is not None:
        # save the figure
        fig.savefig(output_file, format='png', dpi=res)

#-------------------------------------------------------------------------------
# plots a matrix of an event graph
#-------------------------------------------------------------------------------
def plot_graph_matrix(graph_file, output_file=None):
    G = cPickle.load(open(graph_file, 'rb'))
    
    #matrix_prob, pos_sid_prob = nx.attr_matrix(G, edge_attr='weight', normalized=True)
    #matrix_mean, pos_sid_mean = nx.attr_matrix(G, edge_attr='duration_mean')
    #matrix_std, pos_sid_std = nx.attr_matrix(G, edge_attr='duration_std')
    
    for n in nx.nodes_iter(G):
        if G.node[n]['type'] == 'section':
            sequences_by_degree[n] = G.degree(n)
            
    print 'calculating matrix'
    bij_matrix = nx.attr_matrix(G, rc_order=['section', 'sequence'])
    
    print np.max(matrix_prob)
    '''
    print 'plotting matrix', matrix_prob.shape
    #print pos_sid_prob
    
    # plot parameters
    imw = 1024.0 # the full image width
    imh = 1024.0
    lm = 40.0
    rm = 50.0
    tm = 50.0
    bm = 50.0
    res = 72.0
    cbh = 20.0
    cbs = 40.0
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=6)
    
    #arial14 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    #arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    #arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    #arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm - cbh - cbs # the height for both matricies
    pw = imw - lm - rm
    ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))
    
    #print pos_sid_prob[::-1]
    
    ax.imshow(matrix_prob.T, interpolation='none')
    #ax.set_xticks(range(len(pos_sid_prob)))
    #ax.set_yticks(range(len(pos_sid_prob)))
    #ax.set_xticklabels(pos_sid_prob)
    #ax.set_yticklabels(pos_sid_prob)
    #for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
    #    label.set_fontproperties(ticklabelfont)

    #ax.axis('tight')
    #ax.set_ylim((15.5, -0.5))
    #ax.set_xlim((140.5, len(pos_sid_prob)-0.5))

    '''

#-------------------------------------------------------------------------------
# plots an event graph
#-------------------------------------------------------------------------------
def plot_graph(graph_file, output_file, degree_cut=None, label_degree_cut=0.25, self_loops=True, sequence=None):
    if type(graph_file) is str:
        G = cPickle.load(open(graph_file, 'rb'))
    else:
        G = graph_file
    
    #print(nx.clustering(nx.Graph(G), weight='weight'))
    
    # the color map for the plot
    cmap = mplt.get_cmap('GnBu_r')
    
    if degree_cut is not None:
        print 'Original Graph'
        print nx.info(G)
        degrees = G.degree(weight='weight')
        max_degree = float(max(degrees.values()))
        min_degree = float(min(degrees.values()))
        degree_cut_num = min_degree + (max_degree-min_degree)*degree_cut
        print 'max degree: {}'.format(max_degree)
        print 'min degree: {}'.format(min_degree)
        print 'degree cut: {}'.format(degree_cut_num)
        print
        print 'Cut Graph'
        sub_nodes = [n for n, d in G.degree(weight='weight').iteritems() if d > degree_cut_num]
        Gsub = G.subgraph(sub_nodes)
    else:
        Gsub = G

    if not self_loops:
        print 'Removing Self Loops'
        self_loop_edges = Gsub.selfloop_edges()
        Gsub.remove_edges_from(self_loop_edges)
    
    print nx.info(Gsub)
    
    degrees = Gsub.degree(weight='weight')
    max_degree = float(max(degrees.values()))
    min_degree = float(min(degrees.values()))
    print 'max degree: {}'.format(max_degree)
    print 'min degree: {}'.format(min_degree)
    node_min = 0.05
    node_max = 0.2
    node_line_min = 0.1
    node_line_max = 2.0
    min_label_degree = min_degree + (max_degree-min_degree)*label_degree_cut
    min_font_size = 0.5
    max_font_size = 6.0

    widths = {}
    heights = {}
    labels = {}
    styles = {}
    colors = {}
    node_line_widths = {}
    font_sizes = {}
    #print max_degree, min_degree
    for n in nx.nodes_iter(Gsub):
        try:
            sequence_node = Gsub.node[n]['bipartite']
        except KeyError:
            sequence_node = False

        if not sequence_node:
            degree = float(Gsub.degree(n, weight='weight'))
            r,g,b,a = cmap(vcutils.linear_interp(degree, min_degree, max_degree, 0.0, 1.0))
            dim = vcutils.linear_interp(degree, min_degree, max_degree, node_min, node_max)
            widths[n] = dim
            heights[n] = dim
            if degree > min_label_degree:
                labels[n] = n
                font_sizes[n] = vcutils.linear_interp(degree, min_degree, max_degree, min_font_size, max_font_size)
            else:
                labels[n] = ''
            styles[n] = 'filled'
            colors[n] = '#{r:02x}{g:02x}{b:02x}'.format(r=int(r*255.0), g=int(g*255.0), b=int(b*255.0))
            node_line_widths[n] = vcutils.linear_interp(degree, min_degree, max_degree, node_line_min, node_line_max)
        else:
            widths[n] = node_min
            heights[n] = node_min
            styles[n] = 'filled'
            colors[n] = '#FFFFF'
            node_line_widths[n] = 0.5
            labels[n] = ''

    nx.set_node_attributes(Gsub,'width',widths)
    nx.set_node_attributes(Gsub,'height',heights)
    nx.set_node_attributes(Gsub,'label',labels)
    nx.set_node_attributes(Gsub,'style',styles)
    nx.set_node_attributes(Gsub,'fillcolor',colors)
    nx.set_node_attributes(Gsub,'penwidth',node_line_widths)
    nx.set_node_attributes(Gsub,'fontsize',font_sizes)
    #print G.edges(data=True)
    
    try:
        weights = [ float(edata['weight']) for u,v,edata in Gsub.edges(data=True) ]
    except KeyError:
        weights = [ 1.0 for u,v in Gsub.edges() ]
    
    max_weight = float(max(weights))
    min_weight = float(min(weights))
    line_min = 0.1
    line_max = 5.0
    arrow_min = 0.1
    arrow_max = 0.7
    alpha_min = 10.0
    alpha_max = 255.0

    # all of the scales will be at their max
    if min_weight == max_weight:
        min_weight = 0.0
        line_max = 0.5
        arrow_max = 0.1
        alpha_max = 50.0
    
    edge_widths = {}
    arrow_sizes = {}
    edge_colors = {}

    if sequence is not None:
        sequence_parsed = []
        for n,sid in enumerate(sequence):
            if n < len(sequence) - 1:
                sequence_parsed.append('{}-{}'.format(sid, sequence[n+1]))

    for e in nx.edges_iter(Gsub):
        try:
            e_weight = float(Gsub[e[0]][e[1]]['weight'])
        except KeyError:
            e_weight = 1.0

        width = vcutils.linear_interp(e_weight, min_weight, max_weight, line_min, line_max)
        arrow_size = vcutils.linear_interp(e_weight, min_weight, max_weight, arrow_min, arrow_max)
        alpha = vcutils.linear_interp(e_weight, min_weight, max_weight, alpha_min, alpha_max)
        edge_widths[e] = width
        arrow_sizes[e] = arrow_size
        if sequence is not None:
            if '{}-{}'.format(e[0], e[1]) in sequence_parsed:
                edge_colors[e] = '#FF0000'
            else:
                edge_colors[e] = '#000000{:x}'.format(int(alpha))
        else:
            edge_colors[e] = '#000000{:x}'.format(int(alpha))
    
    nx.set_edge_attributes(Gsub, 'penwidth', edge_widths)
    nx.set_edge_attributes(Gsub, 'arrowsize', arrow_sizes)
    nx.set_edge_attributes(Gsub, 'color', edge_colors)
    #cmap = mplt.get_cmap('gray')
    
    #norm = mcolor.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
    
    A=nx.to_agraph(Gsub)        # convert to a graphviz graph
    #A.layout()            # neato layout
    A.draw(output_file, prog='sfdp', args='-Gsize="40!" -Goverlap=prism -Grepulsiveforce=1.0 -GsmoothType="graph_dist" -Goutputorder="edgesfirst" -Nfixedsize="true" -Nfontname="Helvetica"')

#-------------------------------------------------------------------------------
# space-time plot
#-------------------------------------------------------------------------------
def space_time_plot(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        # get the data
        event_data = events.get_event_data(['event_number','event_year','event_magnitude','event_elements', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        section_info = geometry.get_section_info(section_filter=section_filter)
        
        # store a sorted list of section ids
        section_ids = sorted(section_info.keys())
        
        # the section offsets determine the starting x position of each section
        section_offsets = {}
        for i, sid in enumerate(section_ids):
            section_offsets[sid] = sum([section_info[k]['blocks_along_strike'] for k in sorted(section_info.keys())[0:i]])
        
        # calculate various properties of the data set that we will need to
        # set up the plot
        min_depth = min([section_info[k]['blocks_along_dip'] for k in section_info.keys()])
        x_data_size = sum([section_info[k]['blocks_along_strike'] for k in section_info.keys()])
        max_label_len = max([len(section_info[k]['name']) for k in section_info.keys()])
        start_year = event_data['event_year'][0]
        
        # Storing all of the plot parameters here for clarity
        stp_params = {
            'output_file':output_file,
            'x_axis_data_size':x_data_size,
            'y_axis_data_size':event_data['event_range_duration'],
            'max_depth':min_depth,
            'min_mag':min(event_data['event_magnitude']),
            'max_mag':max(event_data['event_magnitude']),
            'start_year':start_year,
            'max_label_len':max_label_len,
            'geometry':geometry,
            'section_offsets':section_offsets
        }
        
        # instantiate the spacetimeplot class
        stp = vcutils.VCSpaceTimePlot(
            stp_params['output_file'],
            stp_params['x_axis_data_size'],
            stp_params['y_axis_data_size'],
            stp_params['max_depth'],
            stp_params['min_mag'],
            stp_params['max_mag'],
            stp_params['start_year'],
            stp_params['max_label_len']
        )
        
        mp = False
        #-----------------------------------------------------------------------
        # The multiprocessing stuff below is not functional. The variable "mp"
        # above should always be set to False.
        #-----------------------------------------------------------------------
        # TODO: Figure out a way to plot in parallel.
        if mp:
            '''
            num_processes = multiprocessing.cpu_count()
        
            # break the work up
            seg = int(round(float(len(event_data['event_magnitude']))/float(num_processes)))
            work_queue = multiprocessing.Queue()
            for i in range(num_processes):
                if i == num_processes - 1:
                    end_index = len(event_data['event_magnitude'])
                else:
                    end_index = seg*int(i + 1)
                work_queue.put({
                    'event_magnitude':event_data['event_magnitude'][int(i) * seg:end_index],
                    'event_elements':event_data['event_elements'][int(i) * seg:end_index],
                    'event_number':event_data['event_number'][int(i) * seg:end_index],
                    'event_year':event_data['event_year'][int(i) * seg:end_index]
                })

            # create a queue to pass to workers to store the results
            result_queue = multiprocessing.Queue()

            # spawn workers
            for i in range(num_processes):
                worker = vcutils.SpaceTimePlotter(stp_params, work_queue, result_queue)
                worker.start()
            
            # collect the results off the queue
            for i in range(num_processes):
                stp.event_lines += result_queue.get().event_lines
            '''
        else:
            # For each event in the found event set, look at the involved
            # elements, and add them to the event line array. Since the event
            # line shows only elements on the strike, elements at depths are
            # projected up to the strike: for every element along the dip the
            # strike value is incremented up to the smallest value of depth in
            # the model.
            for i, enum in enumerate(event_data['event_number']):
                event_line = np.zeros(x_data_size)
                for bid in event_data['event_elements'][i]:
                    sid = geometry[bid]['section_id']
                    try:
                        b_index = section_offsets[sid] + geometry[bid]['das_id']
                        if event_line[b_index] < min_depth:
                            event_line[b_index] += 1
                    except KeyError:
                        pass
                stp.add_event(
                    enum,
                    event_data['event_year'][i],
                    event_data['event_magnitude'][i],
                    event_line
                )

        # Add section labels
        stp.add_section_labels(section_offsets, section_info)

        # Add the title
        stp.add_title('Events from {}'.format(sim_data.filename))

        # Plot the thing
        stp.plot()
    
#-------------------------------------------------------------------------------
# magnitude rupture area plot
#-------------------------------------------------------------------------------
def magnitude_rupture_area(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_area'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_area_kmsq = [quakelib.Conversion().sqm2sqkm(x) for x in event_data['event_area']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_area_kmsq, event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    x_WC = np.linspace(2.2,5184)
    y_WC = 4.07 + 0.98 * np.log10(x_WC)
    y_error_plus_WC = 4.07+0.06 + (0.98+0.03) * np.log10(x_WC)
    y_error_minus_WC = 4.07-0.06 + (0.98-0.03) * np.log10(x_WC)
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]

    # do the standard plot
    vcutils.standard_plot(output_file, event_area_kmsq, event_data['event_magnitude'],
        axis_format='semilogx',
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'x':r'log(Rupture Area [km$^\mathsf{2}$])', 'y':'Magnitude'},
        plot_label='Magnitude-Rupture Area{}'.format(plot_label)
    )
        
#-------------------------------------------------------------------------------
# magnitude average slip plot
#-------------------------------------------------------------------------------
def magnitude_average_slip(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function

    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_data['event_average_slip'], event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)

    x_WC = np.linspace(0.05,8, num=10)
    y_WC = 6.93 + 0.82 * np.log10(x_WC)
    y_error_plus_WC = 6.93+0.05 + (0.82+0.1) * np.log10(x_WC)
    y_error_minus_WC = 6.93-0.05 + (0.82-0.1) * np.log10(x_WC)
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]
    
    # do the standard plot
    vcutils.standard_plot(output_file, event_data['event_average_slip'], event_data['event_magnitude'],
        axis_format='semilogx',
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'y':'Magnitude', 'x':'log(Average Slip [m])'},
        plot_label='Magnitude-Average Slip{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# average slip surface rupture length plot
#-------------------------------------------------------------------------------
def average_slip_surface_rupture_length(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_surface_rupture_length', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_surface_rupture_length_km = [quakelib.Conversion().m2km(x) for x in event_data['event_surface_rupture_length']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_surface_rupture_length_km, event_data['event_average_slip'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)

    x_WC = np.linspace(3.8,432, num=10)
    y_WC = 10.0**(-1.43 + 0.88 * np.log10(x_WC))
    y_error_plus_WC = 10.0**(-1.43+0.18 + (0.88+0.11) * np.log10(x_WC))
    y_error_minus_WC = 10.0**(-1.43-0.18 + (0.88-0.11) * np.log10(x_WC))
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]
    
    # do the standard plot
    vcutils.standard_plot(output_file, event_surface_rupture_length_km, event_data['event_average_slip'],
        axis_format='loglog',
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'y':'log(Average Slip [m])', 'x':'log(Surface Rupture Length [km])'},
        plot_label='Average Slip-Surface Rupture Length{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# frequency magnitude plot
#-------------------------------------------------------------------------------
def frequency_magnitude(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # initilize a dict to store the event counts and get the total number
    # of events.
    cum_freq = {}
    total_events = len(event_data['event_magnitude'])
    
    # count the number of events bigger than each magnitude
    for num, magnitude in enumerate(sorted(event_data['event_magnitude'])):
        cum_freq[magnitude] = total_events - (num + 1)
    
    # dump the counts into x and y arrays for plotting. also, divide the count
    # by the number of years so we get count per year.
    x = []
    y = []
    for magnitude in sorted(cum_freq.iterkeys()):
        x.append(magnitude)
        y.append(float(cum_freq[magnitude])/event_data['event_range_duration'])

    # create the line for b = 1
    x_b1 = np.linspace(min(x),max(x),10)
    y_b1 = 10**(math.log(y[0],10)+x[0]-x_b1)

    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # for the UCERF2 error bars
    x_UCERF = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
    y_UCERF = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
    y_error_UCERF = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
    
    # do the standard plot
    vcutils.standard_plot(output_file, x, y,
        axis_format='semilogy',
        add_lines=[{'label':'b=1', 'x':x_b1, 'y':y_b1}, {'label':'UCERF2', 'x':x_UCERF, 'y':y_UCERF, 'ls':'--', 'c':'red', 'y_error':y_error_UCERF}],
        axis_labels = {'y':'log(# events per year)', 'x':'Magnitude'},
        plot_label='Frequency-Magnitude{}'.format(plot_label),
        connect_points=True,
        legend_loc='upper right'
    )

#-------------------------------------------------------------------------------
# event return map
#-------------------------------------------------------------------------------
def event_return_map(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_trigger', 'event_year', 'event_magnitude', 'event_number'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        xs = []
        ys = []
        sections = set()
        mags = []
        max_mag_on_section = {}
        
        
        # add edges and nodes to the graph for each event
        for i, event_trigger in enumerate(event_data['event_trigger']):
            if i%round(float(len(event_data['event_year']))/100.0) == 0:
                sys.stdout.write('\r event {} of {}'.format(i, len(event_data['event_year'])))
                sys.stdout.flush()
            
            this_sid = geometry.sections_with_elements([event_trigger])[0]
            sections.add(this_sid)
            if i < len(event_data['event_trigger']) - 1:
                try:
                    if event_data['event_magnitude'][i] > max_mag_on_section[this_sid]:
                        max_mag_on_section[this_sid] = event_data['event_magnitude'][i]
                except KeyError:
                    max_mag_on_section[this_sid] = event_data['event_magnitude'][i]
                
                mags.append(max_mag_on_section[this_sid])
                
                next_sid = geometry.sections_with_elements([event_data['event_trigger'][i+1]])[0]
                xs.append(this_sid)
                ys.append(next_sid)
    
    labels = []

    #print len(xs), len(ys), len(event_data['event_magnitude'][0:-1])
    
    for index in range(max(sections)):
        if index + 1 in sections:
            labels.append('{}'.format(index+1))
        else:
            labels.append('')
    
    # plot parameters
    imw = 1024.0 # the full image width
    imh = 1024.0
    lm = 40.0
    rm = 50.0
    tm = 50.0
    bm = 50.0
    res = 72.0
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=7)
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm # the height for both matricies
    pw = imw - lm - rm
    ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))

    #ax.plot(time, output)
    #ax.set_ylim((1, max(pos_sid)))
    
    norm = mcolor.Normalize(vmin=min(mags), vmax=max(mags))
    cmap = mplt.get_cmap('YlOrRd')
    
    ax.scatter(xs,ys, c=cmap(norm(mags)), marker='s', linewidths=1, edgecolors=cmap(norm(mags)))
    ax.set_xticks(range(max(sections)))
    ax.set_yticks(range(max(sections)))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.axis('tight')

    for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)
    ax.set_ylim((0, max(sections)))
    ax.set_xlim((0, max(sections)))







