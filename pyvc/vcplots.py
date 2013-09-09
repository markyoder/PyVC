#!/usr/bin/env python
from pyvc import *
from pyvc import vcutils
from pyvc import vcexceptions
import time
import matplotlib.pyplot as plt
import matplotlib.font_manager as font
import numpy as np
import math

#-------------------------------------------------------------------------------
# standard plotting routine
#-------------------------------------------------------------------------------
def standard_plot(output_file, x, y, legend_loc='best', **kwargs):
    add_lines = kwargs.get('add_lines')
    axis_labels = kwargs.get('axis_labels')
    plot_label = kwargs.get('plot_label')
    connect_points = kwargs.get('connect_points')
    axis_format = kwargs.get('axis_format')

    plot_format = output_file.split('.')[-1]
    if plot_format != 'png' and plot_format != 'pdf' and plot_format != 'dat':
        raise vcexceptions.PlotFormatNotSupported(plot_format)
    elif plot_format == 'png' or plot_format == 'pdf':
        #-----------------------------------------------------------------------
        # plot the data using matplotlib
        #-----------------------------------------------------------------------
        
        #-----------------------------------------------------------------------
        # set up plot dimensions (all values in pixels)
        #-----------------------------------------------------------------------
        # the full image width
        imw = 501.0
        # the full image height
        imh = 501.0 
        # the left margin and bottom margin
        if axis_labels is not None:
            lm = 45.0
            bm = 45.0
        else:
            lm = 10.0
            bm = 10.0
        # the right margin
        rm = 10.0 
        # the top margin
        if plot_label is not None:
            tm = 30.0
        else:
            tm = 10.0
        # the plot resolution
        res = 72.0

        #-----------------------------------------------------------------------
        # set up the plot styles
        #-----------------------------------------------------------------------
        # the marker face color for the main plot
        mfc_main = (0,0,0,1)
        # the marker size for the main plot
        ms_main = 3.0
        # the marker type for the main plot
        mt_main = 'o'
        # the line style for the main plot
        if connect_points is not None:
            ls_main = '--'
        else:
            ls_main = 'None'
        # the line color for the main plot
        c_main = (0,0,0,1)
        # styles for additional sub plots
        if add_lines is not None:
            # the line style for the extra lines
            ls_extra = '-'
            # the line weight for the extra lines
            lw_extra = 4.0
            # the line color for the extra lines
            c_extra = (0.5,0.5,0.5,1)
    
        #-----------------------------------------------------------------------
        # set up the plot fonts
        #-----------------------------------------------------------------------
        ticklabelfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        framelabelfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        legendfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        if len(plot_label) > 100:
            title_size = 9
        elif len(plot_label) <= 100 and len(plot_label) > 89:
            title_size = 10
        elif len(plot_label) <= 89 and len(plot_label) > 50:
            title_size = 11
        else:
            title_size = 12
        titlefont = font.FontProperties(family='Arial', style='normal', variant='normal', size=title_size)
        
        #-----------------------------------------------------------------------
        # do the plot
        #-----------------------------------------------------------------------
        # calculate the final dimensions and create the figure and axis
        imwi = imw/res
        imhi = imh/res
        pw = imw - lm - rm
        ph = imh - tm - bm
        fig = plt.figure(figsize=(imwi, imhi), dpi=res)
        the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
        
        # the main plotting routines
        eval('the_ax.{}(x, y, ls=ls_main, mfc=mfc_main, ms=ms_main, marker=mt_main, c=c_main)'.format(axis_format), locals())
        # plot any additional lines
        if add_lines is not None:
            for line in add_lines:
                eval('the_ax.{}(line[\'x\'], line[\'y\'], ls=ls_extra, lw=lw_extra, c=c_extra, label=line[\'label\'])'.format(axis_format), locals())

        '''
        if axis_format is None:
            # just do linear linear plot
            pass
        elif axis_format == 'semilogx':
            # semilogx just like it says
            the_ax.semilogx(x, y, ls=ls_main, mfc=mfc_main, ms=ms_main, marker=mt_main, c=c_main)
            
            # plot any additional lines
            if add_lines is not None:
                for line in add_lines:
                    the_ax.semilogx(line['x'], line['y'], ls=ls_extra, lw=lw_extra, c=c_extra, label=line['label'])
        elif axis_format == 'loglog':
            # loglog just like it says
            the_ax.loglog(x, y, ls=ls_main, mfc=mfc_main, ms=ms_main, marker=mt_main, c=c_main)
        
            
            # plot any additional lines
            if add_lines is not None:
                for line in add_lines:
                    the_ax.loglog(line['x'], line['y'], ls=ls_extra, lw=lw_extra, c=c_extra, label=line['label'])
        elif axis_format == 'semilogy':
            # semilogx just like it says
            the_ax.semilogy(x, y, ls=ls_main, mfc=mfc_main, ms=ms_main, marker=mt_main, c=c_main)
            
            # plot any additional lines
            if add_lines is not None:
                for line in add_lines:
                    the_ax.semilogy(line['x'], line['y'], ls=ls_extra, lw=lw_extra, c=c_extra, label=line['label'])
        '''
        # set the fonts for the tick labels
        for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(ticklabelfont)
        
        # label the axes
        if axis_labels is not None:
            the_ax.set_ylabel(axis_labels['y'], fontproperties=framelabelfont)
            the_ax.set_xlabel(axis_labels['x'], fontproperties=framelabelfont)

        # label the plot
        if plot_label is not None:
            the_ax.set_title(plot_label, fontproperties=titlefont, x=0, y=1, ha='left')

        # clean up the final plot
        the_ax.autoscale_view(tight=True)
        
        # create a legend if we have extra lines
        if add_lines is not None:
            the_ax.legend(prop=legendfont, loc=legend_loc)

        # save the plot
        plt.savefig(output_file, format=plot_format, dpi=res)

#-------------------------------------------------------------------------------
# parse all of the filters and return a string for the plot title
#-------------------------------------------------------------------------------
def get_plot_label(sim_file, event_range=None, section_filter=None, magnitude_filter=None):
    # always say what file it came from
    label_str = ': from file {}'.format(sim_file)
    # do the event range
    if event_range is not None:
        label_str += ', from {} {}-{}'.format(event_range['type'], event_range['filter'][0], event_range['filter'][1])
    # do the magnitude filter
    if magnitude_filter is not None:
        label_str += ', m {}'.format(magnitude_filter)
    # do the section filter. this is tricky
    if section_filter is not None:
        label_str += ', sections '
        # section_ids are just numbers so they are small, print more of them
        if section_filter['type'] == 'section_id':
            max_sections = 5
        else:
            max_sections = 2
        # if the number of sections is less than our max defined above, print a
        # comma seperated list. if not, just print the first and last.
        if len(section_filter['filter']) < max_sections:
            label_str += ','.join([str(x) for x in section_filter['filter']])
        else:
            label_str += '{}...{}'.format(section_filter['filter'][0], section_filter['filter'][-1])

    return label_str
    
#-------------------------------------------------------------------------------
# magnitude rupture area plot
#-------------------------------------------------------------------------------
def magnitude_rupture_area(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
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
        
        start_time = time.time()
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_area'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        total_time = time.time() - start_time
        print len(event_data['event_magnitude']), total_time
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_area_kmsq = [vcutils.Converter().msq_kmsq(x) for x in event_data['event_area']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_area_kmsq, event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # do the standard plot
    standard_plot(output_file, event_area_kmsq, event_data['event_magnitude'],
        axis_format='semilogx',
        add_lines=[{'label':'binned average', 'x':x_ave, 'y':y_ave}],
        axis_labels = {'x':r'log(Rupture Area [km$^\mathsf{2}$])', 'y':'Magnitude'},
        plot_label='Magnitude-Rupture Area{}'.format(plot_label)
    )
        
#-------------------------------------------------------------------------------
# magnitude average slip plot
#-------------------------------------------------------------------------------
def magnitude_average_slip(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
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
        
        start_time = time.time()
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        total_time = time.time() - start_time
        print len(event_data['event_magnitude']), total_time
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function

    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_data['event_average_slip'], event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # do the standard plot
    standard_plot(output_file, event_data['event_average_slip'], event_data['event_magnitude'],
        axis_format='semilogx',
        add_lines=[{'label':'binned average', 'x':x_ave, 'y':y_ave}],
        axis_labels = {'y':'Magnitude', 'x':'log(Average Slip [m])'},
        plot_label='Magnitude-Average Slip{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# average slip surface rupture length plot
#-------------------------------------------------------------------------------
def average_slip_surface_rupture_length(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
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
        
        start_time = time.time()
        # get the data
        event_data = events.get_event_data(['event_surface_rupture_length', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        total_time = time.time() - start_time
        print len(event_data['event_surface_rupture_length']), total_time
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_surface_rupture_length_km = [vcutils.Converter().m_km(x) for x in event_data['event_surface_rupture_length']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_data['event_surface_rupture_length'], event_data['event_average_slip'])
    
    # get the plot label which will depend on the filters
    plot_label = get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # do the standard plot
    standard_plot(output_file, event_data['event_surface_rupture_length'], event_data['event_average_slip'],
        axis_format='loglog',
        add_lines=[{'label':'binned average', 'x':x_ave, 'y':y_ave}],
        axis_labels = {'y':'log(Average Slip [m])', 'x':'log(Surface Rupture Length [m])'},
        plot_label='Average Slip-Surface Rupture Length{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# frequency magnitude plot
#-------------------------------------------------------------------------------
def frequency_magnitude(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
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
        
        start_time = time.time()
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        total_time = time.time() - start_time
        print len(event_data['event_magnitude']), total_time
    
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
    plot_label = get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # do the standard plot
    standard_plot(output_file, x, y,
        axis_format='semilogy',
        add_lines=[{'label':'b=1', 'x':x_b1, 'y':y_b1}],
        axis_labels = {'y':'log(# events per year)', 'x':'Magnitude'},
        plot_label='Frequency-Magnitude{}'.format(plot_label),
        connect_points=True,
        legend_loc='upper right'
    )