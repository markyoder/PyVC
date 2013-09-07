#!/usr/bin/env python
from pyvc import *
from pyvc import vcutils
from pyvc import vcexceptions
import time
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# standard plotting routine
#-------------------------------------------------------------------------------
def standard_plot(output_file, x, y, **kwargs):
    add_lines = kwargs.get('add_lines')
    axis_labels = kwargs.get('axis_labels')
    plot_label = kwargs.get('plot_label')

    plot_format = output_file.split('.')[-1]

    if plot_format != 'png' and plot_format != 'pdf' and plot_format != 'dat':
        raise vcexceptions.PlotFormatNotSupported(plot_format)
    elif plot_format == 'png' or plot_format == 'pdf':
        #-----------------------------------------------------------------------
        # plot the data using matplotlib
        #-----------------------------------------------------------------------
        
        #-----------------------------------------------------------------------
        # set up plot dimensions
        # all values in pixels
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
            tm = 45.0
        else:
            tm = 10.0
        # the plot resolution
        res = 72.0

        # calculate the final dimensions and create the figure and axis
        imwi = imw/res
        imhi = imh/res
        pw = imw - lm - rm
        ph = imh - tm - bm
        fig = plt.figure(figsize=(imwi, imhi), dpi=res)
        the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
        
        #self.calculateAverages()
        
        x_ave = [a[0] for a in self.averages]
        y_ave = [a[1] for a in self.averages]
        
        x = map(float, self.x)
        y = map(float, self.y)

        #if connect_points:
        ls1 = '--'
        #else:
        #    ls1 = 'None'
        lw1 = 1.0
        c1 = (0.7,0.7,0.7,1)
        mfc1 = (0,0,0,1)
        ms1 = 3
        marker1 = 'o'
        
        ls_extra = '-'
        lw_extra = 4.0
        c_extra = (0.5,0.5,0.5,1)
        
        the_ax.semilogx(x, y, ls = 'None', mfc = mfc1, ms = ms1, marker = marker1)
        
        the_ax.semilogx(x_ave, y_ave, ls = ls_extra, lw = lw_extra, c = c_extra, label = 'binned average')
        
        for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(self.ticklabelfont)
            
        the_ax.set_ylabel('Magnitude', fontproperties=self.framelabelfont)
        the_ax.set_xlabel(r'log(Rupture Area [km$^\mathsf{2}$])', fontproperties=self.framelabelfont)
        
        the_ax.autoscale_view(tight=True)
        
        #the_ax.tick_params(fontproperties=arial12)
        
        #the_ax.margins(0.05,0.3)
        
        #the_ax.set_title(title, position=(0.0,1.05), ha='left', fontproperties=arial12)
        
        the_ax.legend(prop=self.legendfont, loc=2)
           
        plt.savefig(self.out_file, format=self.output_format, dpi=res)


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

        sim_data.open_file(sim_file)
        
        events = VCEvents(sim_data)
        
        start_time = time.time()
        event_data = events.get_event_data(['event_magnitude', 'event_area'], event_range=event_range, magnitude_filter=magnitude_filter)
        total_time = time.time() - start_time
        print len(event_data['event_magnitude']), total_time
        
    event_area_kmsq = [vcutils.Converter().msq_kmsq(x) for x in event_data['event_area']]
        
    x_ave, y_ave = vcutils.calculate_averages(event_area_kmsq, event_data['event_magnitude'])
    
    standard_plot(output_file, event_area_kmsq, event_data['event_magnitude'],
        add_lines=[{'label':'binned average', 'x':x_ave, 'y':y_ave}],
        axis_labels = {'x':r'log(Rupture Area [km$^\mathsf{2}$])', 'y':'Magnitude'},
        plot_label='Test Label'
    )
        
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
