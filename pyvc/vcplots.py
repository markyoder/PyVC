#!/usr/bin/env python
from pyvc import *
from pyvc import vcutils
from pyvc import vcexceptions
import matplotlib.pyplot as plt
import matplotlib.font_manager as font
import numpy as np
import math

#-------------------------------------------------------------------------------
# a class to manage the space-time plot
#-------------------------------------------------------------------------------
class VCSpaceTimePlot(object):
    #---------------------------------------------------------------------------
    # y_axis_data_size is the number of years in the plot, or if event_time is
    # true it is the number of events.
    #---------------------------------------------------------------------------
    def __init__(self, output_file, x_axis_data_size, y_axis_data_size, max_depth, min_mag, max_mag, event_time=False):
        # store the output file
        self.output_file = output_file
        # store the max depth. this is in number of blocks
        self.max_depth = max_depth
        # store the x and y axis data size
        self.x_axis_data_size = x_axis_data_size
        self.y_axis_data_size = y_axis_data_size
        # store the min and max magnitudes
        self.max_mag = max_mag
        self.min_mag = min_mag
        # calculate the slope of the magnitude color ramp
        self.mag_slope = 1.0/(self.max_mag - self.min_mag)
        
        #-----------------------------------------------------------------------
        # set up plot dimensions (all values in pixels)
        #-----------------------------------------------------------------------
        # the event line width and pixels per "y" are based on the y_axis_data_size
        self.elw = 1.0
        self.ppy = 1.0
        if math.ceil(self.y_axis_data_size) < 512 and math.ceil(self.y_axis_data_size) >= 256:
            if event_time:
                self.elw = 2.0
            else:
                self.ppy = 2.0
        elif math.ceil(self.y_axis_data_size) < 256:
            if event_time:
                self.elw = 4.0
            else:
                self.ppy = 4.0
        if event_time:
            self.ppy = self.elw
        # the height of the plot area depends on the pixels per "y" value
        self.ph = self.ppy * math.ceil(self.y_axis_data_size)
        
        # the event line width and pixels per "y" are based on the y_axis_data_size
        self.ppx = 1.0
        if math.ceil(self.x_axis_data_size) < 640 and math.ceil(self.x_axis_data_size) >= 320:
            self.ppx = 2.0
        elif math.ceil(self.x_axis_data_size) < 320:
            self.ppx = 4.0
        # the height of the plot area depends on the pixels per "y" value
        self.pw = self.ppx * math.ceil(self.x_axis_data_size)
        
        # the margins
        self.lm = 45.0
        self.bm = 45.0
        self.rm = 10.0
        self.tm = 150.0
        
        # the total image height
        self.imh = self.ph + self.tm + self.bm
        # the total image width
        self.imw = self.pw + self.lm + self.rm
        
        # the plot resolution
        self.res = 72.0

        #-----------------------------------------------------------------------
        # set up the plot styles
        #-----------------------------------------------------------------------
        # the color map for the plot
        self.cmap = plt.get_cmap('autumn_r',lut=1000)
        # the color of the vertical section demarcation lines
        self.sdlc = '0.5'
        # the alpha of the vertical section demarcation lines
        self.sdla = 0.5
        # the width of the vertical section demarcation lines
        self.sdlw = 0.5
        # the padding between the section names and the top of the plot frame
        self.snp = 1.05
        
        #-----------------------------------------------------------------------
        # set up the plot fonts
        #-----------------------------------------------------------------------
        self.ticklabelfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.sectionlabelfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.framelabelfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        self.legendfont = font.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        self.titlefont = font.FontProperties(family='Arial', style='normal', variant='normal', size=10)

        #-----------------------------------------------------------------------
        # start the plot
        #-----------------------------------------------------------------------
        # calculate the final dimensions and create the figure and axis
        self.imwi = self.imw/self.res
        self.imhi = self.imh/self.res
        self.fig = plt.figure(figsize=(self.imwi, self.imhi), dpi=self.res)
        self.the_ax = self.fig.add_axes((self.lm/self.imw, self.bm/self.imh, self.pw/self.imw, self.ph/self.imh))

    def add_event(self, event_number, event_year, event_magnitude, event_line):
        last_value = None
        last_slot = -1
        in_line = False
        (r,g,b,a) = self.cmap(self.mag_slope * (event_magnitude - self.min_mag))
        slots = np.nonzero(event_line)[0]
        #print slots
        #print event_line[slots]
        for k, slot in enumerate(slots):
            value = event_line[slot]
            if (slot - 1) != last_slot and k != 0:
                # we have jumped a gap
                # draw the last line
                #print 'draw', line_start, line_end + 1, last_value
                self.the_ax.axhline(
                    y=event_year,
                    xmin=float(line_start)/float(self.x_axis_data_size),
                    xmax=float(line_end + 1)/float(self.x_axis_data_size),
                    alpha=float(last_value)/float(self.max_depth),
                    linewidth=self.elw,
                    color=(r,g,b)
                )
                in_line = False
                # start a new line
                in_line = True
                line_start = slot
                line_end = slot
            else:
                # we are in a contiguous group
                if last_value != event_line[slot]:
                    if in_line:
                        # draw the last line
                        #print 'draw', line_start, line_end + 1, last_value
                        self.the_ax.axhline(
                            y=event_year,
                            xmin=float(line_start)/float(self.x_axis_data_size),
                            xmax=float(line_end + 1)/float(self.x_axis_data_size),
                            alpha=float(last_value)/float(self.max_depth),
                            linewidth=self.elw,
                            color=(r,g,b)
                        )
                        in_line = False
                    # start line
                    in_line = True
                    line_start = slot
                    line_end = slot
                else:
                    # add to current line
                    line_end = slot
            
            if k == len(slots) - 1:
                # draw the last line
                #print 'draw', line_start, line_end + 1, value
                self.the_ax.axhline(
                    y=event_year,
                    xmin=float(line_start)/float(self.x_axis_data_size),
                    xmax=float(line_end + 1)/float(self.x_axis_data_size),
                    alpha=float(value)/float(self.max_depth),
                    linewidth=self.elw,
                    color=(r,g,b)
                )
                in_line = False
            last_value = value
            last_slot = slot
        #print
        '''
        for slot, value in enumerate(event_line):
            #print slot, value
            if value != 0:
                if last_value != value:
                    if in_line:
                        # draw the last line
                        #print 'draw', line_start, line_end, last_value
                        self.the_ax.axhline(
                            y=event_year,
                            xmin=float(line_start)/float(self.x_axis_data_size),
                            xmax=float(line_end)/float(self.x_axis_data_size),
                            alpha=float(value)/float(self.max_depth),
                            linewidth=self.elw,
                            color='0.5'
                        )
                        in_line = False
                    # start line
                    in_line = True
                    line_start = slot
                    line_end = slot + 1
                else:
                    # add to current line
                    line_end = slot
            else:
                if in_line:
                    # draw the last line
                    #print 'draw', line_start, line_end, last_value
                    self.the_ax.axhline(
                        y=event_year,
                        xmin=float(line_start)/float(self.x_axis_data_size),
                        xmax=float(line_end)/float(self.x_axis_data_size),
                        alpha=float(value)/float(self.max_depth),
                        linewidth=self.elw,
                        color='0.5'
                    )
                    in_line = False
            # if we are at the end and there is an active line draw it
            if slot == len(event_line) - 1:
                if in_line:
                    # draw the last line
                    #print 'draw', line_start, line_end, last_value
                    self.the_ax.axhline(
                        y=event_year,
                        xmin=float(line_start)/float(self.x_axis_data_size),
                        xmax=float(line_end)/float(self.x_axis_data_size),
                        alpha=float(value)/float(self.max_depth),
                        linewidth=self.elw,
                        color='0.5'
                    )
                    in_line = False
    
            last_value = value
            '''
        #print
    
    def add_sections(self, section_offsets, section_info):
        #print section_info
        section_ids = sorted(section_offsets.keys())
        for snum, sid in enumerate(section_ids):
            if section_offsets[sid] != 0:
                self.the_ax.axvline(
                    x=section_offsets[sid],
                    linewidth=self.sdlw,
                    color=self.sdlc,
                    alpha=self.sdla
                )
            x_text_loc = (float(snum)+0.5)/float(len(section_ids))
            self.the_ax.text(
                x_text_loc, self.snp,
                '{} {}'.format(sid,section_info[sid]['name']),
                horizontalalignment='center',
                verticalalignment='bottom',
                rotation=90,
                fontproperties=self.sectionlabelfont,
                transform=self.the_ax.transAxes
            )
            #print section_offsets[sid], self.x_axis_data_size, section_info[sid]['name']
    
    def plot(self):
        self.the_ax.set_xlim((0, self.x_axis_data_size))
        #self.the_ax.autoscale(enable=True, axis='both', tight=True)
        
        plot_format = self.output_file.split('.')[-1]
        if plot_format != 'png' and plot_format != 'pdf':
            raise vcexceptions.PlotFormatNotSupported(plot_format)
        else:
            self.fig.savefig(self.output_file, format=plot_format, dpi=self.res)

#-------------------------------------------------------------------------------
# space-time plot
#-------------------------------------------------------------------------------
def space_time_plot(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
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
        
        section_offsets = {}
        for i, sid in enumerate(section_ids):
            section_offsets[sid] = sum([section_info[k]['blocks_along_strike'] for k in sorted(section_info.keys())[0:i]])
        
        min_depth = min([section_info[k]['blocks_along_dip'] for k in section_info.keys()])
        x_data_size = sum([section_info[k]['blocks_along_strike'] for k in section_info.keys()])
        stp = VCSpaceTimePlot(
            output_file,
            x_data_size,
            event_data['event_range_duration'],
            min_depth,
            min(event_data['event_magnitude']),
            max(event_data['event_magnitude'])
        )
        
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
        
        stp.add_sections(section_offsets, section_info)
        
        stp.plot()
    #print event_data
    
    #print math.ceil(event_data['event_range_duration']), sum([sections[k]['blocks_along_strike'] for k in sections.keys()])

    #event_elements = [events.get_event_elements(enum) for enum in event_data['event_number']]


#-------------------------------------------------------------------------------
# standard plotting routine for scaling plots
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
    #---------------------------------------------------------------------------
    # plot the data using matplotlib
    #---------------------------------------------------------------------------
        
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
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_area'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_area_kmsq = [vcutils.Converter().msq_kmsq(x) for x in event_data['event_area']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_area_kmsq, event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
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
        
        # get the data
        event_data = events.get_event_data(['event_surface_rupture_length', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_surface_rupture_length_km = [vcutils.Converter().m_km(x) for x in event_data['event_surface_rupture_length']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_data['event_surface_rupture_length'], event_data['event_average_slip'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
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
    
    # do the standard plot
    standard_plot(output_file, x, y,
        axis_format='semilogy',
        add_lines=[{'label':'b=1', 'x':x_b1, 'y':y_b1}],
        axis_labels = {'y':'log(# events per year)', 'x':'Magnitude'},
        plot_label='Frequency-Magnitude{}'.format(plot_label),
        connect_points=True,
        legend_loc='upper right'
    )