#!/usr/bin/env python
import tables
import Queue
import multiprocessing
from operator import itemgetter
import time
import itertools
import matplotlib.pyplot as mplt
import matplotlib.font_manager as mfont
import matplotlib.lines as mlines
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import numpy as np
import sys
import quakelib

import math
import calendar
from geographiclib.geodesic import Geodesic

import os
import re
import subprocess


#-------------------------------------------------------------------------------
# Comprehensive CPU availability check.
# http://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-in-python
#-------------------------------------------------------------------------------
def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # http://code.google.com/p/psutil/
    try:
        import psutil
        return psutil.NUM_CPUS
    except (ImportError, AttributeError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system')

#-------------------------------------------------------------------------------
# Given a set of maxes and mins return a linear value betweem them.
#-------------------------------------------------------------------------------
def linear_interp(x, x_min, x_max, y_min, y_max):
    return ((y_max - y_min)/(x_max - x_min) * (x - x_min)) + y_min

#-------------------------------------------------------------------------------
# A class to perform unit conversions.
#-------------------------------------------------------------------------------
'''
This is now in quakelib
class Converter:
        def __init__(self):
                self.lat0 = 31.5
                self.lon0 = -126.0
                self.earth_radius = self.km_m(6371.0)
    
        #-----------------------------------------------------------------------
        # Given two lat/lon points return the distance between them in meters.
        #-----------------------------------------------------------------------
        def surface_distance(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return conv["s12"]
    
        #-----------------------------------------------------------------------
        # Given two lat/lon point return cartesian coords with point 1 at the
        # origin and the given coords at point 2.
        #-----------------------------------------------------------------------
        def latlon_xy_2pt(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
    
        #-----------------------------------------------------------------------
        # Convert a year given in decimal format (ie. 4.234) and return the
        # year, month, day and optionaly hour, minute, second, fractional second
        #-----------------------------------------------------------------------
        def yearDecimal_YearMonthDay(self,yearDecimal,time=False):
                decimal, year = math.modf(yearDecimal)
                decimal, month = math.modf(12.0*decimal)
                decimal, day = math.modf(calendar.monthrange(1972, int(month + 1))[1]*decimal)
                decimal, hour = math.modf(24.0*decimal)
                decimal, minute = math.modf(60.0*decimal)
                decimal, second = math.modf(60.0*decimal)
                
                if time:
                    return int(year), int(month + 1), int(day + 1), int(hour), int(minute), int(second), decimal
                else:
                    return int(year), int(month + 1), int(day + 1)
    
        #-----------------------------------------------------------------------
        # Convert years to seconds.
        #-----------------------------------------------------------------------
        def year_sec(self, year):
            return 365.0 * 24.0 * 60.0 * 60.0 * year
        
        #-----------------------------------------------------------------------
        # Convert seconds to years.
        #-----------------------------------------------------------------------
        def sec_year(self, sec):
            return (1.0/365.0) * (1.0/24.0) * (1.0/60.0) * (1.0/60.0) * sec
        
        #-----------------------------------------------------------------------
        # Convert kilometers to meters.
        #-----------------------------------------------------------------------
        def km_m(self, km):
            return km * 10.0**(3.0)
        
        #-----------------------------------------------------------------------
        # Convert meters to kilometers.
        #-----------------------------------------------------------------------
        def m_km(self, m):
            return m * 10.0**(-3.0)
        
        #-----------------------------------------------------------------------
        # Convert meters squared to kilometers squared.
        #-----------------------------------------------------------------------
        def msq_kmsq(self, msq):
            return msq * 1000.0**(-2.0)
        
        #-----------------------------------------------------------------------
        # Convert degrees to radians.
        #-----------------------------------------------------------------------
        def deg_rad(self, deg):
            return deg * (math.pi/180.0)
        
        #-----------------------------------------------------------------------
        # Convert radians to degrees.
        #-----------------------------------------------------------------------
        def rad_deg(self, rad):
            return rad * (180.0/math.pi)
        
        #-----------------------------------------------------------------------
        # Set the lat and lon to serve as the zero location.
        #-----------------------------------------------------------------------
        def set_latlon0(self, lat0, lon0):
            self.lat0 = lat0
            self.lon0 = lon0
        
        #-----------------------------------------------------------------------
        # Convert a lat/lon point to cartesian coords relative to lat0/lon0.
        #-----------------------------------------------------------------------
        def latlon_xy(self, lat, lon):
            conv = Geodesic.WGS84.Inverse(self.lat0, self.lon0, lat, lon)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
        
        #-----------------------------------------------------------------------
        # Convert x/y points to lat/lon relative to lat0/lon0.
        #-----------------------------------------------------------------------
        def xy_latlon(self, x, y):
            s12 = (x**2.0 + y**2.0)**(0.5)
            azi1 = self.rad_deg(math.atan2(x,y))
            conv = Geodesic.WGS84.Direct(self.lat0, self.lon0, azi1, s12)
            
            return (conv["lat2"], conv["lon2"])
'''
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
# calculate binned averages for scatter plots with the x-axis plotted in logpace
#-------------------------------------------------------------------------------
def calculate_averages(x,y):
    num_bins = math.floor(len(x)/100)
    
    if num_bins < 20:
        num_bins = 20
    elif num_bins > 100:
        num_bins = 100
    
    x = np.array(x)
    y = np.array(y)
    
    if np.min(x) == 0:
        bin_min = 1
    else:
        bin_min = math.floor(math.log(np.min(x),10))
    bin_max = math.ceil(math.log(np.max(x),10))
    
    bins = np.logspace(bin_min,bin_max,num=num_bins)
    inds = np.digitize(x, bins)
    
    binned_data = {}

    for n, i in enumerate(inds):
        try:
            binned_data[i].append(y[n])
        except KeyError:
            binned_data[i] = [y[n]]

    x_ave = []
    y_ave = []
    for k in sorted(binned_data.keys()):
        if k != 0:
            x_ave.append(0.5*(bins[k-1]+bins[k]))
            y_ave.append(sum(binned_data[k])/float(len(binned_data[k])))

    return x_ave, y_ave

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
        ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        if len(plot_label) > 100:
            title_size = 9
        elif len(plot_label) <= 100 and len(plot_label) > 89:
            title_size = 10
        elif len(plot_label) <= 89 and len(plot_label) > 50:
            title_size = 11
        else:
            title_size = 12
        titlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=title_size)
        
        #-----------------------------------------------------------------------
        # do the plot
        #-----------------------------------------------------------------------
        # calculate the final dimensions and create the figure and axis
        imwi = imw/res
        imhi = imh/res
        pw = imw - lm - rm
        ph = imh - tm - bm
        fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
        the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
        
        # the main plotting routines
        eval('the_ax.{}(x, y, ls=ls_main, mfc=mfc_main, ms=ms_main, marker=mt_main, c=c_main)'.format(axis_format), locals())
        # plot any additional lines
        if add_lines is not None:
            for line in add_lines:
                try:
                    ls_extra = line['ls']
                except KeyError:
                    pass
                try:
                    lw_extra = line['lw']
                except KeyError:
                    pass
                try:
                    c_extra = line['c']
                except KeyError:
                    pass
                try:
                    yerr=line['y_error']
                    add_axis_format = 'errorbar'
                    y_error_key = 'yerr=line[\'y_error\'],'
                except KeyError:
                    y_error_key = ''
                    add_axis_format = axis_format
                
                eval('the_ax.{0}(line[\'x\'], line[\'y\'], ls=ls_extra, {1} lw=lw_extra, c=c_extra, label=line[\'label\'])'.format(add_axis_format, y_error_key), locals())
        
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
        
        if output_file is not None:
            # save the plot
            mplt.savefig(output_file, format=plot_format, dpi=res)

#-------------------------------------------------------------------------------
# a class to manage the space-time plot
#-------------------------------------------------------------------------------
class VCSpaceTimePlot(object):
    #---------------------------------------------------------------------------
    # y_axis_data_size is the number of years in the plot, or if event_time is
    # true it is the number of events.
    #---------------------------------------------------------------------------
    def __init__(self, output_file, x_axis_data_size, y_axis_data_size, max_depth, min_mag, max_mag, start_year, max_label_len, event_time=False):
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
        # store the max label length used to calculate the top margin
        self.max_label_len = max_label_len
        # set the start and end year
        self.start_year = start_year
        self.end_year = self.start_year + self.y_axis_data_size

        #-----------------------------------------------------------------------
        # set up the plot fonts
        #-----------------------------------------------------------------------
        self.ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.sectionlabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
        self.legendtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        self.legendticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.titlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)

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
        
        # the event line width and pixels per "y" are based on the y_axis_data_size
        self.ppx = 1.0
        if math.ceil(self.x_axis_data_size) < 640 and math.ceil(self.x_axis_data_size) >= 320:
            self.ppx = 2.0
        elif math.ceil(self.x_axis_data_size) < 320:
            self.ppx = 4.0

        # the margins
        self.lm = 50.0
        self.bm = 10.0
        self.rm = 10.0
        
        # the colorbar dimensions
        self.cbw = 200.0
        self.cbh = 10.0
        
        # the padding between the section names and the top of the plot frame
        self.snp = 40.0
        
        # the plot resolution
        self.res = 72.0
        
        #-----------------------------------------------------------------------
        # plot dimensions that are calculated based on the settings above
        #-----------------------------------------------------------------------
        # if we are using event time the pixels per "y" is just the event line
        # width
        if event_time:
            self.ppy = self.elw

        # the height of the plot area depends on the pixels per "y" value
        self.ph = self.ppy * math.ceil(self.y_axis_data_size)

        # the height of the plot area depends on the pixels per "y" value
        self.pw = self.ppx * math.ceil(self.x_axis_data_size)

        # the top margin is based on the size of the section labels and their
        # placement
        self.tm = self.snp + self.sectionlabelfont.get_size() * 0.75 * self.max_label_len + self.cbh + 20

        # the total image height
        self.imh = self.ph + self.tm + self.bm
        # the total image width
        self.imw = self.pw + self.lm + self.rm
        
        #-----------------------------------------------------------------------
        # set up the plot styles
        #-----------------------------------------------------------------------
        # the color map for the plot
        self.cmap = mplt.get_cmap('autumn_r',lut=1000)
        # the color of the vertical section demarcation lines
        self.sdlc = '0.5'
        # the alpha of the vertical section demarcation lines
        self.sdla = 0.5
        # the width of the vertical section demarcation lines
        self.sdlw = 0.5
        # the color of the section label lines
        self.sllc = '0.5'
        # the alpha of the section label lines
        self.slla = 0.5
        # the width of the section label lines
        self.sllw = 0.5

        #-----------------------------------------------------------------------
        # start the plot
        #-----------------------------------------------------------------------
        # calculate the final dimensions and create the figure and axis
        self.imwi = self.imw/self.res
        self.imhi = self.imh/self.res
        self.fig = mplt.figure(figsize=(self.imwi, self.imhi), dpi=self.res)
        self.the_ax = self.fig.add_axes((self.lm/self.imw, self.bm/self.imh, self.pw/self.imw, self.ph/self.imh))
    
    #---------------------------------------------------------------------------
    # Adds event to the plot by plotting them directly instead of storing them.
    # This behavior may have to change to storing when running in parallel.
    #---------------------------------------------------------------------------
    def add_event(self, event_number, event_year, event_magnitude, event_line):
        # initilize state variables
        last_value = None
        last_slot = -1
        in_line = False
        # the color of the event line based on the magnitude
        (r,g,b,a) = self.cmap(self.mag_slope * (event_magnitude - self.min_mag))
        # The slots are the indicies of the non-zero line elements from the
        # event line array. These correspond to x data locations where there is
        # a line.
        slots = np.nonzero(event_line)[0]
        # Go through the slots drawing lines where appropriate.
        for k, slot in enumerate(slots):
            value = event_line[slot]
            if (slot - 1) != last_slot and k != 0:
                # we have jumped a gap
                # draw the last line
                self.the_ax.axhline(
                    y=event_year,
                    xmin=float(line_start)/float(self.x_axis_data_size),
                    xmax=float(line_end + 1)/float(self.x_axis_data_size),
                    #alpha=float(last_value)/float(self.max_depth),
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
                        self.the_ax.axhline(
                            y=event_year,
                            xmin=float(line_start)/float(self.x_axis_data_size),
                            xmax=float(line_end + 1)/float(self.x_axis_data_size),
                            #alpha=float(last_value)/float(self.max_depth),
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
                self.the_ax.axhline(
                    y=event_year,
                    xmin=float(line_start)/float(self.x_axis_data_size),
                    xmax=float(line_end + 1)/float(self.x_axis_data_size),
                    #alpha=float(value)/float(self.max_depth),
                    linewidth=self.elw,
                    color=(r,g,b)
                )
                in_line = False
            last_value = value
            last_slot = slot
    
    #---------------------------------------------------------------------------
    # Add the section labels and the section demarcation lines.
    #---------------------------------------------------------------------------
    def add_section_labels(self, section_offsets, section_info):
        # an array to store all of the label lines
        label_lines = []
        
        # The sorted section ids. Sections are printed in ascending order (by
        # section id) from left to right.
        section_ids = sorted(section_offsets.keys())
        
        # Go through each section plotting the demarcation line, the section
        # label, and the section label line.
        for snum, sid in enumerate(section_ids):
            # Plot the demarcation line (except at the zero offset)
            if section_offsets[sid] != 0:
                self.the_ax.axvline(
                    x=section_offsets[sid],
                    linewidth=self.sdlw,
                    color=self.sdlc,
                    alpha=self.sdla
                )
            
            # Set the text location. This will evenly distribute the text across
            # the length of the x-axis.
            x_text_loc = (float(snum)+0.5)/float(len(section_ids))
            
            # Find the horizontal center of each section region
            if snum == len(section_ids) - 1:
                x_fault_loc = (float(section_offsets[sid] + self.x_axis_data_size)/2.0)/float(self.x_axis_data_size)
            else:
                x_fault_loc = (float(section_offsets[sid] + section_offsets[section_ids[snum + 1]])/2.0)/float(self.x_axis_data_size)
            
            # Plot the section label
            self.the_ax.text(
                x_text_loc, (self.ph + self.snp)/self.ph,
                '{} {}'.format(sid,section_info[sid]['name']),
                horizontalalignment='center',
                verticalalignment='bottom',
                rotation=90,
                fontproperties=self.sectionlabelfont,
                transform=self.the_ax.transAxes
            )

            # Create the section line
            line_xs = [x_text_loc, x_text_loc, x_fault_loc, x_fault_loc]
            line_ys = [
                (self.ph + self.snp - 5)/self.ph,
                (self.ph + self.snp * ( 3.0/4.0 ))/self.ph,
                (self.ph + self.snp * ( 1.0/4.0 ))/self.ph,
                1
            ]
            label_lines.append(mlines.Line2D(
                line_xs, line_ys,transform=self.the_ax.transAxes, solid_capstyle='round',
                    alpha=self.slla,
                    linewidth=self.sllw,
                    color=self.sllc
            ))
        
        # Add the section label lines to the plot.
        self.the_ax.lines.extend(label_lines)
    
    #---------------------------------------------------------------------------
    # Add the title.
    #---------------------------------------------------------------------------
    def add_title(self, title):
        self.the_ax.set_title(
            title,
            position=((10-self.lm)/self.pw,(self.ph + self.tm - 27)/self.ph),
            ha='left',
            fontproperties=self.titlefont,
        )
    
    #---------------------------------------------------------------------------
    # Finish off the plot and save it.
    #---------------------------------------------------------------------------
    def plot(self):
        # Set the axis limits
        self.the_ax.set_xlim((0, self.x_axis_data_size))
        self.the_ax.set_ylim((self.start_year, self.end_year))
        
        # Create the colorbar
        cb_ax    = self.fig.add_axes([(self.lm + self.pw - self.cbw)/self.imw, (self.imh - self.cbh - 10)/self.imh, self.cbw/self.imw, self.cbh/self.imh])
        norm = mcolor.Normalize(vmin=self.min_mag, vmax=self.max_mag)
        cb = mcolorbar.ColorbarBase(cb_ax, cmap=self.cmap,
                   norm=norm,
                   orientation='horizontal')
        ticks = map(float, [self.min_mag,(self.min_mag+ self.max_mag)/2.0, self.max_mag])
        cb.set_ticks([ticks[0], ticks[1], ticks[2]])
        cb.set_ticklabels(['%.2f'%ticks[0],'%.2f'%ticks[1], '%.2f'%ticks[2]])
        
        # Style and cleanup the colorbar ticks
        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(self.legendticklabelfont)
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        # Set the colorbar label
        cb.ax.set_xlabel('Magnitude',position=(1,0), ha='right', fontproperties=self.legendtitlefont)
        
        # Style and cleanup the plot ticks
        for label in self.the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(self.ticklabelfont)
        for label in self.the_ax.xaxis.get_ticklabels():
            label.set_alpha(0)
        for line in self.the_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        # Set the plot y-axis label
        self.the_ax.set_ylabel('Year',fontproperties=self.framelabelfont)
        
        if self.output_file is not None:
            # Get the plot format and save the file
            plot_format = self.output_file.split('.')[-1]
            if plot_format != 'png' and plot_format != 'pdf':
                raise vcexceptions.PlotFormatNotSupported(plot_format)
            else:
                self.fig.savefig(self.output_file, format=plot_format, dpi=self.res)

#-------------------------------------------------------------------------------
# A class for plotting spacetime events in parallel. This is currently under
# development.
#-------------------------------------------------------------------------------
# TODO: Figure out a way to plot in parallel.
class SpaceTimePlotter(multiprocessing.Process):
    def __init__(self, space_time_plot_params, work_queue, result_queue):
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        #the path to the sim file
        self.stp_params = space_time_plot_params
        
        super(SpaceTimePlotter, self).__init__()
    
    def run(self):
    
        stp = VCSpaceTimePlot(
            self.stp_params['output_file'],
            self.stp_params['x_axis_data_size'],
            self.stp_params['y_axis_data_size'],
            self.stp_params['max_depth'],
            self.stp_params['min_mag'],
            self.stp_params['max_mag'],
            self.stp_params['start_year'],
            self.stp_params['max_label_len'],
            lite_init=True
        )
        
        geometry = self.stp_params['geometry']
        x_data_size = self.stp_params['x_axis_data_size']
        min_depth = self.stp_params['max_depth']
        section_offsets = self.stp_params['section_offsets']
        
        while not self.kill_received:
            # get a task
            try:
                event_data = self.work_queue.get_nowait()
            except Queue.Empty:
                break
        
            # do the processing
            print 'processing events {} - {}'.format(event_data['event_number'][0],event_data['event_number'][-1])
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
            self.result_queue.put(stp)
        
        #self.sim_file.close()

#-------------------------------------------------------------------------------
# A class for processing event data in parallel. This is only used if a
# simulation is being analyzed for the first time. The quantities calculated
# here are saved in the simulation file.
#-------------------------------------------------------------------------------
class EventDataProcessor(multiprocessing.Process):
    def __init__(self, sim_file_path, work_queue, result_queue):
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        #the path to the sim file
        self.sim_file_path = sim_file_path
        
        super(EventDataProcessor, self).__init__()
    
    def run(self):
    
        #the simulation file
        self.sim_file = tables.open_file(self.sim_file_path, 'r')
        
        #the tables we will need to process
        sweep_table = self.sim_file.get_node('/event_sweep_table')
        event_table = self.sim_file.get_node('/event_table')
        block_info_table = self.sim_file.get_node('/block_info_table')
        
        #these getters speed up the caculation of surface rupture length
        pt1getter = itemgetter(3,4,5)
        pt4getter = itemgetter(18,19,20)
        
        while not self.kill_received:
            # get a task
            try:
                event_range = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            # do the processing
            print 'processing events {} - {}'.format(*event_range)
            results = {}
            for evnum in range(*event_range):
                areas = {}
                total_slip = 0.0
                slip_records = 0
                surface_rupture_length = 0.0
                involved_sections = []
                for sweep in sweep_table[event_table[evnum]['start_sweep_rec']:event_table[evnum]['end_sweep_rec']]:
                    eleid = sweep['block_id']
                    if block_info_table[eleid]['m_trace_flag_pt1'] > 0:
                        surface_rupture_length += (sum((x-y)**2.0 for x, y in itertools.izip(pt1getter(block_info_table[eleid]),pt4getter(block_info_table[eleid]))))**0.5
                    
                    involved_sections.append(block_info_table[eleid]['section_id'])
                    areas[eleid] = sweep['area']
                    total_slip += sweep['slip']
                    slip_records += 1
        
                results[evnum] = {'average_slip':total_slip/float(slip_records), 'area':sum(areas.values()), 'surface_rupture_length':surface_rupture_length, 'involved_sections':set(involved_sections)}
            
            self.result_queue.put(results)
        
        self.sim_file.close()

#-------------------------------------------------------------------------------
# A class that manages the connection to a simulation data file. This class
# should be used within a with statement. For example:
#
# with VCSimData() as sim_data:
#     [do stuff with sim_data]
#
# This ensures that the __exit__ method is called and the file is closed when no
# longer needed.
#-------------------------------------------------------------------------------
class VCSimData(object):
    def __init__(self, file_path=None):
        self.file = None
        self.file_path = file_path
        
        self.do_event_area = False
        self.do_event_average_slip = False
        self.do_event_surface_rupture_length = False
        self.do_events_by_section = False
        self.do_das_id = False
        self.do_depth_id = False
        
        if self.file_path is not None:
            self.open_file(self.file_path)
    
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self.file is not None:
            self.file.close()
    
    @property
    def filename(self):
        if self.file_path is not None:
            return self.file_path.split('/')[-1]
        else:
            return None

    def open_file(self, file_path):
        if self.file_path is None:
            self.file_path = file_path
        self.file = tables.open_file(self.file_path)

        #-----------------------------------------------------------------------
        # check to see if we need to calculate additional data
        #-----------------------------------------------------------------------
        if 'event_area' not in self.file.root.event_table.colnames:
            self.do_event_area = True

        if 'event_average_slip' not in self.file.root.event_table.colnames:
            self.do_event_average_slip = True
    
        if 'event_surface_rupture_length' not in self.file.root.event_table.colnames:
            self.do_event_surface_rupture_length = True
        
        if 'events_by_section' not in self.file.root._v_groups.keys():
            self.do_events_by_section = True
        
        if self.do_event_area or self.do_event_average_slip or self.do_event_surface_rupture_length or self.do_events_by_section:
            self.calculate_additional_event_data()

        if 'das_id' not in self.file.root.block_info_table.colnames:
            self.do_das_id = True
        
        if 'depth_id' not in self.file.root.block_info_table.colnames:
            self.do_depth_id = True
        
        if self.do_das_id or self.do_depth_id:
            self.calculate_additional_block_data()
        
        if 'model_extents' not in self.file.root._v_children.keys():
            self.calculate_model_extents()
    
    def calculate_model_extents(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        block_info_table = self.file.root.block_info_table
    
        #-----------------------------------------------------------------------
        # calculate the model extents
        #-----------------------------------------------------------------------
        print 'Calculating model extents'
        
        start_time = time.time()
        
        sys_max_z = -sys.float_info.max
        sys_min_z = sys.float_info.max
        sys_max_x = -sys.float_info.max
        sys_min_x = sys.float_info.max
        sys_max_y = -sys.float_info.max
        sys_min_y = sys.float_info.max
        
        for block in block_info_table:
            min_x = min((block['m_x_pt1'], block['m_x_pt2'], block['m_x_pt3'], block['m_x_pt4']))
            max_x = max((block['m_x_pt1'], block['m_x_pt2'], block['m_x_pt3'], block['m_x_pt4']))
        
            min_y = min((block['m_y_pt1'], block['m_y_pt2'], block['m_y_pt3'], block['m_y_pt4']))
            max_y = max((block['m_y_pt1'], block['m_y_pt2'], block['m_y_pt3'], block['m_y_pt4']))
            
            min_z = min((block['m_z_pt1'], block['m_z_pt2'], block['m_z_pt3'], block['m_z_pt4']))
            max_z = max((block['m_z_pt1'], block['m_z_pt2'], block['m_z_pt3'], block['m_z_pt4']))
            
            if min_x < sys_min_x:
                sys_min_x = min_x
            if max_x > sys_max_x:
                sys_max_x = max_x
                
            if min_y < sys_min_y:
                sys_min_y = min_y
            if max_y > sys_max_y:
                sys_max_y = max_y
            
            if min_z < sys_min_z:
                sys_min_z = min_z
            if max_z > sys_max_z:
                sys_max_z = max_z
        
        base_lat_lon_table = self.file.root.base_lat_lon
        conv = quakelib.Conversion(base_lat_lon_table[0], base_lat_lon_table[1])
        
        ne_corner = conv.convert2LatLon( quakelib.Vec3(sys_max_x, sys_max_y, 0.0) )
        sw_corner = conv.convert2LatLon( quakelib.Vec3(sys_min_x, sys_min_y, 0.0) )
    
        self.file.close()
    
        print 'Done! {} seconds'.format(time.time() - start_time)
        
        print 'Creating new tables'
        table_start_time = time.time()
        
        self.file = tables.open_file(self.file_path, 'a')
        
        desc = {
            'min_x':tables.Float64Col(dflt=0.0),
            'max_x':tables.Float64Col(dflt=0.0),
            'min_y':tables.Float64Col(dflt=0.0),
            'max_y':tables.Float64Col(dflt=0.0),
            'min_z':tables.Float64Col(dflt=0.0),
            'max_z':tables.Float64Col(dflt=0.0),
            'min_lat':tables.Float64Col(dflt=0.0),
            'max_lat':tables.Float64Col(dflt=0.0),
            'min_lon':tables.Float64Col(dflt=0.0),
            'max_lon':tables.Float64Col(dflt=0.0)
        }
        
        model_extents = self.file.create_table('/', 'model_extents', desc, 'Model Extents')
        
        model_extents.row.append()
        model_extents.flush()
        
        model_extents.cols.min_x[0] = sys_min_x
        model_extents.cols.max_x[0] = sys_max_x
        model_extents.cols.min_y[0] = sys_min_y
        model_extents.cols.max_y[0] = sys_max_y
        model_extents.cols.min_z[0] = sys_min_z
        model_extents.cols.max_z[0] = sys_max_z
        model_extents.cols.min_lat[0] = sw_corner.lat()
        model_extents.cols.max_lat[0] = ne_corner.lat()
        model_extents.cols.min_lon[0] = sw_corner.lon()
        model_extents.cols.max_lon[0] = ne_corner.lon()
        
        print 'Done! {} seconds'.format(time.time() - table_start_time)
        
        #-----------------------------------------------------------------------
        # close the file and reopen it with the new table
        #-----------------------------------------------------------------------
        self.file.close()
        self.file = tables.open_file(self.file_path)
        
        print 'Total time {} seconds'.format(time.time() - start_time)
        
    def calculate_additional_block_data(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        block_info_table = self.file.root.block_info_table
    
        #-----------------------------------------------------------------------
        # get the new block data
        #-----------------------------------------------------------------------
        print 'Calculating new block data'
        
        start_time = time.time()
        
        das_ids = []
        depth_ids = []
        
        curr_das = None
        curr_sec = None
        for block in block_info_table:
            sec = block['section_id']
            das = block['m_das_pt1']
            if sec != curr_sec:
                #new sec
                curr_sec = sec
                das_id = -1
                depth_id = -1
            if das != curr_das:
                #new das
                curr_das = das
                das_id += 1
                depth_id = -1
            depth_id += 1
            das_ids.append(das_id)
            depth_ids.append(depth_id)
    
        print 'Done! {} seconds'.format(time.time() - start_time)
        
        print 'Creating new tables'
        table_start_time = time.time()
        #-----------------------------------------------------------------------
        # create the new block_info_table
        #-----------------------------------------------------------------------
        #close the current file
        self.file.close()
        #reopen the file with the append flag set
        self.file = tables.open_file(self.file_path, 'a')
        block_info_table = self.file.root.block_info_table
        
        # get a description of table in dictionary format
        desc_orig = block_info_table.description._v_colObjects
        desc_new = desc_orig.copy()
        
        # add columns to description
        if self.do_das_id:
            desc_new['das_id'] = tables.UIntCol(dflt=0)
        
        if self.do_depth_id:
            desc_new['depth_id'] = tables.UIntCol(dflt=0)
        
        # create a new table with the new description
        block_info_table_new = self.file.create_table('/', 'tmp', desc_new, 'Block Info Table')
        
        # copy the user attributes
        block_info_table.attrs._f_copy(block_info_table_new)
        
        # fill the rows of new table with default values
        for i in xrange(block_info_table.nrows):
            block_info_table_new.row.append()
        
        # flush the rows to disk
        block_info_table_new.flush()
        
        # copy the columns of source table to destination
        for col in desc_orig:
            getattr(block_info_table_new.cols, col)[:] = getattr(block_info_table.cols, col)[:]

        # fill the new columns
        if self.do_das_id:
            block_info_table_new.cols.das_id[:] = das_ids
        
        if self.do_depth_id:
            block_info_table_new.cols.depth_id[:] = depth_ids
        	
        # remove the original table
        block_info_table.remove()
        
        # move table2 to table
        block_info_table_new.move('/','block_info_table')
        
        print 'Done! {} seconds'.format(time.time() - table_start_time)
        
        #-----------------------------------------------------------------------
        # close the file and reopen it with the new tables
        #-----------------------------------------------------------------------
        self.file.close()
        self.file = tables.open_file(self.file_path)
    
        print 'Total time {} seconds'.format(time.time() - start_time)
    
    def calculate_additional_event_data(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        total_events = self.file.root.event_table.nrows
        #close the current file
        self.file.close()
        
        #-----------------------------------------------------------------------
        # get the new event data
        #-----------------------------------------------------------------------
        print 'Calculating new event data'
        
        start_time = time.time()
        num_processes = multiprocessing.cpu_count()
        
        # break the work up
        seg = int(round(float(total_events)/float(num_processes)))
        work_queue = multiprocessing.Queue()
        for i in range(num_processes):
            if i == num_processes - 1:
                end_index = total_events
            else:
                end_index = seg*int(i + 1)
            work_queue.put((int(i) * seg , end_index))

        # create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()

        # spawn workers
        for i in range(num_processes):
            worker = EventDataProcessor(self.file_path, work_queue, result_queue)
            worker.start()

        # collect the results off the queue
        results = {}
        for i in range(num_processes):
            results = dict(results, **result_queue.get())
        data_process_results_sorted = [results[key] for key in sorted(results.keys())]
        
        print 'Done! {} seconds'.format(time.time() - start_time)
        
        print 'Creating new tables'
        table_start_time = time.time()
        #-----------------------------------------------------------------------
        # create the new event_table
        #-----------------------------------------------------------------------
        self.file = tables.open_file(self.file_path, 'a')
        event_table = self.file.root.event_table
        
        # get a description of table in dictionary format
        desc_orig = event_table.description._v_colObjects
        desc_new = desc_orig.copy()
        
        # add columns to description
        if self.do_event_area:
            desc_new['event_area'] = tables.Float64Col(dflt=0.0)
        
        if self.do_event_average_slip:
            desc_new['event_average_slip'] = tables.Float64Col(dflt=0.0)
        
        if self.do_event_surface_rupture_length:
            desc_new['event_surface_rupture_length'] = tables.Float64Col(dflt=0.0)
        
        # create a new table with the new description
        event_table_new = self.file.create_table('/', 'tmp', desc_new, 'Event Table')
        
        # copy the user attributes
        event_table.attrs._f_copy(event_table_new)
        
        # fill the rows of new table with default values
        for i in xrange(event_table.nrows):
            event_table_new.row.append()
        
        # flush the rows to disk
        event_table_new.flush()
        
        # copy the columns of source table to destination
        for col in desc_orig:
            getattr(event_table_new.cols, col)[:] = getattr(event_table.cols, col)[:]

        # fill the new columns
        if self.do_event_area:
            event_table_new.cols.event_area[:] = [ x['area'] for x in data_process_results_sorted ]
        
        if self.do_event_average_slip:
            event_table_new.cols.event_average_slip[:] = [ x['average_slip'] for x in data_process_results_sorted ]

        if self.do_event_surface_rupture_length:
            event_table_new.cols.event_surface_rupture_length[:] = [ x['surface_rupture_length'] for x in data_process_results_sorted ]
        	
        # remove the original table
        event_table.remove()
        
        # move table2 to table
        event_table_new.move('/','event_table')
        
        #-----------------------------------------------------------------------
        # create a new group to store the event ids for each section
        #-----------------------------------------------------------------------
        if self.do_events_by_section:
            events_by_section_group = self.file.create_group("/", 'events_by_section', 'Events on each section')
            
            # dict to store the events on each section
            events_by_section = {}
            
            # we have the sections involved in each event. we need to transpose
            # this to events on each section
            for evid, event_data in enumerate(data_process_results_sorted):
                for secid in event_data['involved_sections']:
                    try:
                        events_by_section[secid].append(evid)
                    except KeyError:
                        events_by_section[secid] = [evid]
            
            # create an array for each result
            for secid, evids in events_by_section.iteritems():
                self.file.create_array(events_by_section_group, 'section_{}'.format(secid), np.array(evids), 'Events on section {}'.format(secid))
            
        print 'Done! {} seconds'.format(time.time() - table_start_time)
        
        #-----------------------------------------------------------------------
        # close the file and reopen it with the new tables
        #-----------------------------------------------------------------------
        self.file.close()
        self.file = tables.open_file(self.file_path)
    
        print 'Total time {} seconds'.format(time.time() - start_time)
        
