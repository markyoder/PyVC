#!/usr/bin/env python
import os
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import matplotlib as mpl
mpl.use('agg')
import multiprocessing
import matplotlib.font_manager
import matplotlib.pyplot as plt
import h5py
import json
from PIL import Image, ImageEnhance
import cPickle
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import sys as py_sys
import math
import Converter
import MyVector as vec
import Queue
from matplotlib import cm
#import quakelib
import gc
'''
import matplotlib as mpl
#mpl.use('agg')
import numpy as np
import scipy.stats as stats
#import matplotlib
import matplotlib.pyplot as plt
import operator
import matplotlib.lines as mlines
import h5py
import Converter
import Utils
import pprint
import math
import sys as py_sys
import random
import cPickle
import sys
import gc

from mpl_toolkits.mplot3d import axes3d
#from matplotlib.font_manager import FontProperties
import matplotlib.font_manager
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

from PIL import Image, ImageEnhance
import PIL.ImageOps

#from scipy.interpolate import griddata as scigd

from mpl_toolkits.basemap import Basemap, maskoceans, interp

# sys.path.append('%s/classes/'%sys.path[0])
#print sys.path

#import quakelib

import scipy.sparse

#from scipy.interpolate import griddata

import MyVector as vec
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
from matplotlib.patches import Polygon
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import cm
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['mathtext.fontset'] = 'stix'
#mpl.rcParams['text.usetex']=True
#from matplotlib.backends.backend_pdf import PdfPages


import Queue
import itertools
import time
import Chaos

import subprocess

import json
import datetime

import smtplib
from email.mime.text import MIMEText

mp_status_updater = None
'''
class EventInfoPrinter(object):
    def __init__(self, sys):
        self.sys    = sys
        self.event_data = []
        
    def addEvent(self, event_data):
        self.event_data.append({'num':event_data[0],
                                'year':event_data[1],
                                'trigger':event_data[2],
                                'mag':event_data[3]
                                })
        
    def eprint(self):
        for e_dat in sorted(self.event_data,lambda a,b: cmp(b['mag'],a['mag'])):
            sid = self.sys.geometry.elements[e_dat['trigger']].sid
            s_name = self.sys.geometry.sections[sid].sname
            print  e_dat['num'], e_dat['mag'], e_dat['year'], s_name, sid, e_dat['trigger']

class InfoTheoryElementParser(multiprocessing.Process):
    def __init__(self, vc_sys, state_data, work_queue, result_queue):
        
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.vc_sys = vc_sys
        self.state_data = state_data
    
    def run(self):
        while not self.kill_received:
 
            # get a task
            try:
                element_list = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            #print element_list
            
            elements = self.vc_sys.geometry.elementsSortedByID(element_list)
            
            # the actual processing
            element_data = {}
            print '            calculating element data for %i elements '%(len(elements))
            for n, e in enumerate(elements):
                orbit_raw = self.state_data.cff_data[::,2 + e.eid]
                #orbit_shift = orbit - orbit.max()
                #orbit_normalized_raw = orbit_shift/orbit_shift.min()
                
                #remove the first 10% of values as transients
                start_t = int(float(orbit_raw.shape[0]) * 0.1)
                #orbit_normalized = orbit_normalized_raw[start_t::]
                orbit = orbit_raw[start_t::]
            
                W_L = []
                H = []
                Delta_H = []
                E = 0
                for i, L in enumerate(range(2, self.state_data.max_word_length + 1)):
                    print '%i:%i '%(e.eid, L),
                    sys.stdout.flush()
                    W_L.append(L)
                    H_L = Chaos.permutationEntropy(orbit,L)
                    H.append(H_L)
                    try:
                        Delta_H.append(H_L - H[i-1])
                    except:
                        Delta_H.append(H_L)
                
                h_mu = abs(H[-1] - H[-2])
                
                for i, H_L in enumerate(H[1::]):
                    E += H_L - h_mu
                
                element_data[e.eid] = {'W_L':W_L, 'H':H, 'Delta_H':Delta_H, 'E':E, 'h_mu':h_mu}
                
            # store the result
            self.result_queue.put(element_data)

class ElementStateData(object):
    def __init__(self, sys, cff_data, element_data, max_word_length):
        self.sys                = sys
        self.max_word_length    = max_word_length
        self.cff_data           = None
        self.element_data       = None
        parse_cff_data          = False
        parse_element_data      = False
        
        if element_data:
            try:
                print '        **** Loading element data'
                self.element_data =  cPickle.load(open('%s/element-data_WL-%i.pkl'%(self.sys.cache_dir, self.max_word_length), 'rb'))
            except:
                print '            **** Element data not found. Will calculate.'
                parse_element_data = True
                cff_data = True
                
        if cff_data:
            if self.cff_data is None:
                try:
                    print '        **** Loading CFF data'
                    self.cff_data = np.load('%s_cache/%s_CFF.npy'%(sys.name))
                except:
                    print '            **** CFF data not found. Will calculate.'
                    parse_cff_data = True
        
        if parse_cff_data:
            print '        **** Parsing CFF data'
            cff_data_raw = np.fromfile('%s_CFF.dat'%sys.name,sep=' ')
            cff_data_raw = np.reshape(cff_data_raw, (cff_data_raw.size/(len(sys.geometry.elements)+2), len(sys.geometry.elements)+2))
            
            cff_data = []
            curr_eid = 0
            in_event = False
            last_state = None
            for i, cff_state in enumerate(cff_data_raw[1::]):
                if i == 0:
                    cff_state[0] = curr_eid
                    cff_data.append(cff_state)
                else:
                    in_event = bool(cff_state[0])
                    
                    if not in_event:
                        last_state[0] = curr_eid
                        cff_data.append(last_state)
                        curr_eid += 1
                        cff_state[0] = curr_eid
                        cff_data.append(cff_state)
                    
                last_state = cff_state
            
            final_state = cff_data_raw[-1]
            
            final_state[0] = curr_eid
            cff_data.append(final_state)
                    
            self.cff_data = np.array(cff_data)
            cff_data_raw = None
            
            print '        **** Saving CFF data'
            if not os.path.exists('%s_cache'%self.sys.name):
                os.makedirs('%s_cache'%self.sys.name)
            np.save('%s_cache/%s_CFF.npy'%(sys.name,sys.name), self.cff_data)
        
        if parse_element_data:
            print '        **** Calculating element data'
            num_processes = 4
            
            all_elements = []
            for element in self.sys.geometry.elementsSortedByID():
                    all_elements.append(element.eid)
            
            seg = int(round(float(len(all_elements))/float(num_processes)))
            
            elements = []
            
            for i in range(num_processes):
                if i == num_processes - 1:
                    end_index = len(all_elements)
                else:
                    end_index = seg*int(i + 1)
                elements.append((all_elements[ int(i) * seg : end_index ]))
            
            #print elements
            
            work_queue = multiprocessing.Queue()
            for job in elements:
                work_queue.put(job)
            
            # create a queue to pass to workers to store the results
            result_queue = multiprocessing.Queue()
         
            # spawn workers
            for i in range(len(elements)):
                worker = InfoTheoryElementParser(self.sys, self, work_queue, result_queue)
                worker.start()
         
            # collect the results off the queue
            results = []
            for i in range(len(elements)):
                results.append(result_queue.get())
            
            self.element_data = {}
            
            for result in results:
                for key, value in result.iteritems():
                    self.element_data[key] = value
            
            print
            print '        **** Saving element data'
            if not os.path.exists('%s_cache'%self.sys.name):
                os.makedirs('%s_cache'%self.sys.name)
            cPickle.dump(self.element_data, open('%s/element-data_WL-%i.pkl'%(self.sys.cache_dir, self.max_word_length), 'wb'))

    
    def prepareData(self):
        
        print '        **** Removing sweeps from the raw CFF data'
        
        self.element_ids = self.cff_data_raw[0,2::]
        
        cff_data = []
        curr_eid = 0
        in_event = False
        last_state = None
        for i, cff_state in enumerate(self.cff_data_raw[1::]):
            if i == 0:
                cff_state[0] = curr_eid
                cff_data.append(cff_state)
            else:
                in_event = bool(cff_state[0])
                
                if not in_event:
                    last_state[0] = curr_eid
                    cff_data.append(last_state)
                    curr_eid += 1
                    cff_state[0] = curr_eid
                    cff_data.append(cff_state)
                
            last_state = cff_state
        
        final_state = self.cff_data_raw[-1]
        
        final_state[0] = curr_eid
        cff_data.append(final_state)
                
        self.cff_data = np.array(cff_data)
        self.cff_data_raw = None
        
        print '        **** Analyzing CFF data up to word length %i '%self.max_word_length
        
        self.element_data = {}
        
        ele_index = 0
        
        for eid, e in self.sys.geometry.elements.iteritems(): 
            print '            element %i: '%eid,
            sys.stdout.flush()
            orbit_raw = self.cff_data[::,2 + ele_index]
            #orbit_shift = orbit - orbit.max()
            #orbit_normalized_raw = orbit_shift/orbit_shift.min()
            
            #remove the first 10% of values as transients
            start_t = int(float(orbit_raw.shape[0]) * 0.1)
            #orbit_normalized = orbit_normalized_raw[start_t::]
            orbit = orbit_raw[start_t::]
        
            W_L = []
            H = []
            Delta_H = []
            E = 0
            for i, L in enumerate(range(2, self.max_word_length + 1)):
                print '%i '%L,
                sys.stdout.flush()
                W_L.append(L)
                H_L = Chaos.permutationEntropy(orbit,L)
                H.append(H_L)
                try:
                    Delta_H.append(H_L - H[i-1])
                except:
                    Delta_H.append(H_L)
            
            h_mu = abs(H[-1] - H[-2])
            
            for H_L in H:
                E += H_L - h_mu
            
            self.element_data[eid] = {'W_L':W_L, 'H':H, 'Delta_H':Delta_H, 'E':E, 'h_mu':h_mu}
            
            print
            ele_index += 1
        
    def addEvent(self, event_data):
        if self.orbits is None:
            self.orbits = []
        
        self.orbits.append(event_data)
        
        '''
        event_elements = {}
        for sweep in sweeps:
            eid = sweep[2]
            try:
                event_elements[eid]['total_slip'] += sweep[3]
                event_elements[eid]['shear_final'] = sweep[8]
                event_elements[eid]['normal_final'] = sweep[9]
                event_elements[eid]['number_of_slips'] += 1
            except KeyError:
                event_elements[eid] = {  'ele_id'                   :eid,
                                            'total_slip'            :sweep[3],
                                            'shear_init'            :sweep[6],
                                            'normal_init'           :sweep[7],
                                            'shear_final'           :sweep[8],
                                            'normal_final'          :sweep[9],
                                            'number_of_slips'       :1,
                                            'sweep_of_first_slip'   :sweep[1]
                                            }
        self.events.append([int(event_data[0]), float(event_data[1]), event_elements])
        '''
    
    def prepareOrbit(self):
        print '        **** Preparing orbit'
        self.data_changed = True
        friction_vals = {}
        orbits = []
        
        orbit_prepend = ['event_id','event_year']
        for eid in sorted(self.sys.geometry.elements.iterkeys()):
            e = self.sys.geometry.elements[eid]
            friction_vals[eid] = abs(e.dynamic_strength - e.static_strength)/155700000
            orbit_prepend.append('element_%i'%eid)
        
        if self.events[0][0] == 0:
            first_event = [-1,0]
            for eid in sorted(self.sys.geometry.elements.iterkeys()):
                e = self.sys.geometry.elements[eid]
                first_event.append(-(abs(e.dynamic_strength - e.static_strength)))
            orbits.append(first_event)
        
        for event in self.events:
            before = [event[0], event[1]]
            after = [event[0], event[1]]
            event_elements = event[2]
            for eid in sorted(self.sys.geometry.elements.iterkeys()):
                cff_before = None
                cff_after = None
                try:
                    #this element is involved in the event
                    cff_before = event_elements[eid]['shear_init'] - friction_vals[eid] * event_elements[eid]['normal_init']
                    cff_after = event_elements[eid]['shear_final'] - friction_vals[eid] * event_elements[eid]['normal_final']
                except KeyError:
                    #this element was not involved in the event get values from previous event if we can
                    try:
                        cff_before = orbits[len(orbits)-1][eid + 2]
                        cff_after = orbits[len(orbits)-1][eid + 2]
                    except IndexError:
                        e = self.sys.geometry.elements[eid]
                        cff_before = -(abs(e.dynamic_strength - e.static_strength))
                        cff_after = -(abs(e.dynamic_strength - e.static_strength))
                        
                before.append(cff_before)
                after.append(cff_after)
            orbits.append(before)
            orbits.append(after)
        
        for i in orbits:
            print i
        self.orbits = orbits
        self.events = None
    
    def analyzeOrbit(self):
        if self.orbits is None:
            self.prepareOrbit()
        
        if self.do_IT_calc or self.do_map:
            print '        **** Analyzing orbit up to word length %i '%self.max_word_length
            self.data_changed = True
            element_data = {}
            element_maps = {}
            
            np_orbits = np.array( self.orbits[1::] )
            
            ele_index = 0
            
            for eid, e in self.sys.geometry.elements.iteritems(): 
                print '            element %i: '%eid,
                sys.stdout.flush()
                orbit = np_orbits[::,2 + ele_index]
                orbit_shift = orbit - orbit.max()
                orbit_normalized_raw = orbit_shift/orbit_shift.min()
                
                #remove the first 5% of values as transients
                start_t = int(float(orbit_normalized_raw.shape[0]) * 0.05)
                orbit_normalized = orbit_normalized_raw[start_t::]
                
                if self.do_map:
                    print 'map ',
                    sys.stdout.flush()
                    X = []
                    Y = []
                    for t, val in enumerate(orbit_normalized):
                        try:
                            X.append(val)
                            Y.append(orbit_normalized[t+1])
                        except IndexError:
                            break
                    del X[-1]
                    element_maps[eid] = (X,Y)


                if self.do_IT_calc:
                    W_L = []
                    H = []
                    Delta_H = []
                    E = 0
                    for i, L in enumerate(range(2, self.max_word_length + 1)):
                        print '%i '%L,
                        sys.stdout.flush()
                        W_L.append(L)
                        H_L = Chaos.permutationEntropy(orbit_normalized,L)
                        H.append(H_L)
                        try:
                            Delta_H.append(H_L - H[i-1])
                        except:
                            Delta_H.append(H_L)
                    
                    h_mu = abs(H[-1] - H[-2])
                    
                    for H_L in H:
                        E += H_L - h_mu
                    
                    element_data[eid] = {'W_L':W_L, 'H':H, 'Delta_H':Delta_H, 'E':E, 'h_mu':h_mu}
                
                print
                ele_index += 1
            
            
            
            if self.do_IT_calc:
                self.element_data = element_data
                self.do_IT_calc = False
            if self.do_map:
                self.element_maps = element_maps
                self.do_map = False
        
class ElementStatePlotter(object):
    def __init__(self, sys):
        self.sys    = sys
        self.name = None
        self.out_file = None
        self.out_text_file = None
        self.years_per_inch = 20.0
        self.output_format = None
        self.font_family = 'Arial'
    
class ElementCFFEHmuPlotter(ElementStatePlotter):
    def __init__(self, sys):
        print '    **** Initilizing element excess entropy and entropy rate plotter'
        super(ElementCFFEHmuPlotter,self).__init__(sys)
    
    def plot(self, state_data):
        print '        **** Plotting CFF E and h_mu'
        #prepare all of the orbits
        #element_data = state_data.element_data
        #print state_data.element_data[0]['h_mu'], state_data.element_data[0]['H']
        
        #E_test = 0
        #h_mu_test = state_data.element_data[0]['h_mu']
        #print h_mu_test
        #for i, H in enumerate(state_data.element_data[0]['H']):
        #    E_test += H - h_mu_test
        #    print H - h_mu_test, H, h_mu_test
        #print E_test
        #prepare plot
        #font_family_titles = 'Arial'
        #font_family = 'Arial'
        
        fig = plt.figure(figsize=[6.5, 4],dpi=100)
        E_ax = fig.add_axes([0,0,.5,1],projection='3d')
        h_mu_ax = fig.add_axes([.5,0,.5,1],projection='3d')
        
        min_X = min_Y = py_sys.float_info.max
        max_X = max_Y = - py_sys.float_info.max
        
        ele_zs = [0.0, 0.0, 0.0, 0.0]
        
        #for f in mpl.font_manager.createFontList(mpl.font_manager.get_fontconfig_fonts()):
        #    print f
        
        #do plot
        for sid, layered_section in self.sys.geometry.layered_sections.iteritems():
            x_loc = []
            y_loc = []
            E = []
            h_mu = []
            ele_coords = []
            it = np.nditer(layered_section, flags=['multi_index'])
            #print layered_section, layered_section.shape
            while not it.finished:
                #it.multi_index[0] goes down dip is the x xais
                #it.multi_index[1] goes along strike is the y axis
                ele = self.sys.geometry.elements[int(it[0])]
                depths = map(round, ele.depths())
                das = map(round, ele.das())
                
                if min(depths) < min_X:
                    min_X = min(depths)
                if max(depths) > max_X:
                    max_X = max(depths)
                
                if min(das) < min_Y:
                    min_Y = min(das)
                if max(das) > max_Y:
                    max_Y = max(das)
                
                coords = [zip(depths,das, ele_zs)]
                #poly = Polygon(coords, zorder=0)
                ele_coords.append(coords)
                #E_ax.add_patch(poly)
                #art3d.pathpatch_2d_to_3d(poly, z=0, zdir="z")
            
                x_loc.append(0.5 * depths[0] + 0.5 * depths[1])
                y_loc.append(0.5 * das[3] + 0.5 * das[0])
                E.append(state_data.element_data[int(it[0])]['E'])
                h_mu.append(state_data.element_data[int(it[0])]['h_mu'])
                #print depths, das
                
                it.iternext()
            
            xi = np.linspace(min(x_loc),max(x_loc), num=layered_section.shape[0])
            yi = np.linspace(min(y_loc),max(y_loc), num=layered_section.shape[1])
        
            X, Y = np.meshgrid(xi, yi)
        
            Z_E = griddata(x_loc, y_loc, E, xi, yi)
            Z_h_mu = griddata(x_loc, y_loc, h_mu, xi, yi)
            
            
            #E_ax.scatter(X, Y, Z_E)
            E_ax.plot_surface(X, Y, Z_E, rstride=1, cstride=1, linewidth=0.25, edgecolor = (.4,.4,.4), zorder = 10, color=(.5,.5,.5,.5))
            h_mu_ax.plot_surface(X, Y, Z_h_mu, rstride=1, cstride=1, linewidth=0.25, edgecolor = (.4,.4,.4), zorder = 9, color=(.5,.5,.5,.5))
            
            for coords in ele_coords:
                #poly = Polygon(coords, zorder=0)
                E_ax.add_collection3d(Poly3DCollection(coords, color=(.8,.8,.8,.5), edgecolor = (.7,.7,.7),))
                h_mu_ax.add_collection3d(Poly3DCollection(coords, color=(.8,.8,.8,.5), edgecolor = (.7,.7,.7),))
                #art3d.pathpatch_2d_to_3d(poly, z=0, zdir="z")
            
            #print min_X, max_X, min_Y, max_Y
            
        E_ax.set_xlim3d((min_X, max_X))
        E_ax.set_ylim3d((min_Y, max_Y))
        E_ax.set_zlim3d((0, max(E)))
        
        E_ax.set_xlabel('depth [m]',family=self.font_family,size=8)
        E_ax.set_ylabel('distance along strike [m]',family=self.font_family,size=8)
        E_ax.set_zlabel('E',family=self.font_family,size=8)
        for tick in E_ax.xaxis.get_major_ticks():
            tick.label1.set_family(self.font_family)
            tick.label1.set_fontsize(6)
        for tick in E_ax.yaxis.get_major_ticks():
            tick.label1.set_family(self.font_family)
            tick.label1.set_fontsize(6)
        for tick in E_ax.zaxis.get_major_ticks():
            tick.label1.set_family(self.font_family)
            tick.label1.set_fontsize(6)
        
        h_mu_ax.set_xlim3d((min_X, max_X))
        h_mu_ax.set_ylim3d((min_Y, max_Y))
        h_mu_ax.set_zlim3d((0, max(h_mu)))
        
        h_mu_ax.set_xlabel('depth [m]',family=self.font_family,size=8)
        h_mu_ax.set_ylabel('distance along strike [m]',family=self.font_family,size=8)
        h_mu_ax.set_zlabel(r'h$_\mu$',family=self.font_family,size=8)
        for tick in h_mu_ax.xaxis.get_major_ticks():
            tick.label1.set_family(self.font_family)
            tick.label1.set_fontsize(6)
        for tick in h_mu_ax.yaxis.get_major_ticks():
            tick.label1.set_family(self.font_family)
            tick.label1.set_fontsize(6)
        for tick in h_mu_ax.zaxis.get_major_ticks():
            tick.label1.set_family(self.font_family)
            tick.label1.set_fontsize(6)
            
            
        fig.savefig("%s"%(self.out_file), format=self.output_format)

class ElementCFFPlotter(ElementStatePlotter):
    def __init__(self, sys):
        print '    **** Initilizing Element CFF Plotter'
        super(ElementCFFPlotter,self).__init__(sys)
        
    def plot(self, state_data, start_event, end_event):
        print '        **** Plotting CFF vs time'
        orbits = state_data.cff_data[start_event*2:(end_event+1)*2]
        
        #print state_data.cff_data.shape()
        
        #print orbits[0,0], orbits[0,1], orbits[-1,0], orbits[-1,1], state_data.cff_data.shape, orbits.shape, start_event, end_event, (orbits[-1,1] - orbits[0,1])
        
        #for i in state_data.cff_data:
        #    print int(i[0]), float(i[1])
        
        
        #fig_width = (1.0/self.years_per_inch) * (orbits[-1,1] - orbits[0,1])
        #fig_height = 8.5
        
        fig_width = 6.5
        fig_height = 4.0
        
        #print fig_width
        
        if fig_width * 100 > 32768:
            fig_width = 32768.0/100.0
        
        fig = plt.figure(figsize=[fig_width, fig_height],dpi=100)
        
        width_frac = (fig_width-1.5)/fig_width
        left_frac = 1.0/fig_width
        
        height_frac = (fig_height-1.0)/fig_height
        bottom_frac = 0.5/fig_height
        
        ax1 = fig.add_axes([left_frac, bottom_frac, width_frac, height_frac],projection='rectilinear')
        
        #np_orbits = np.array( orbits )
        X =  orbits[::,1]
        
        ax1.axhline(0,c='k',alpha=0.2,ls=':')
        
        for x in X:
            ax1.axvline(x,c='k',alpha=0.2,ls='--')
            #print x
        
        ele_index = 0
        for eid, e in self.sys.geometry.elements.iteritems():
            ax1.plot(X, orbits[::,2 + ele_index],c='k',alpha=0.3)
            ele_index += 1
        
        ax1.set_xlim((orbits[0,1],orbits[-1,1]))
        
        ax1.set_ylabel("CFF",family=self.font_family)
        ax1.set_xlabel("Year",family=self.font_family)
        
        #ax1.set_xticks( np.arange(X[0],X[-1],50))
        
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontname(self.font_family)
            tick.label1.set_fontsize(10)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontname(self.font_family)
            tick.label1.set_fontsize(10)
        
        ax1.set_title('CFF vs time for the %s system'%(self.sys.name),family=self.font_family,size=12, va='bottom', ha='left', x=0)
        
        fig.savefig("%s"%(self.out_file), format=self.output_format)

class ElementCFFMapPlotter(ElementStatePlotter):
    def __init__(self, sys):
        print '    **** Initilizing Element CFF map Plotter'
        super(ElementCFFMapPlotter,self).__init__(sys)
        
    def plot(self, state_data):
        print '        **** Plotting CFF Map'
        
        #prepare the data
        element_maps = {}
        
        ele_index = 0
        for eid, e in self.sys.geometry.elements.iteritems():
            orbit_raw = state_data.cff_data[::,2 + ele_index]
            #remove the first 30% of values as transients
            start_t = int(float(orbit_raw.shape[0]) * 0.3)
            orbit = orbit_raw[start_t::]
            
            orbit_shift = orbit - orbit.max()
            orbit_normalized = orbit_shift/orbit_shift.min()
            
            sys.stdout.flush()
            X = []
            Y = []
            for t, val in enumerate(orbit_normalized):
                try:
                    X.append(val)
                    Y.append(orbit_normalized[t+1])
                except IndexError:
                    break
            del X[-1]
            element_maps[eid] = (X,Y)
            ele_index += 1

                
        #prepare the plot
        ele_width = 2.5
        ele_height = 2.5
        
        plot_padding = 0.0
        plot_width = ele_width - 2.0*plot_padding
        plot_height = ele_height - 2.0*plot_padding
        
        total_vert_elements = 0
        max_horizontal_elements = 0
        
        for layered_section in self.sys.geometry.layered_sections.itervalues():
            total_vert_elements += layered_section.shape[0]
            if layered_section.shape[1] > max_horizontal_elements:
                max_horizontal_elements = layered_section.shape[1]
        
        fig_width = max_horizontal_elements * ele_width + 1.5
        fig_height = total_vert_elements * ele_height + (len(self.sys.geometry.layered_sections) - 1.0) * 1.0 + 1.5
        
        ele_width_frac = ele_width/fig_width
        ele_height_frac = ele_width/fig_height
        
        plot_width_frac = plot_width/fig_width
        plot_height_frac = plot_width/fig_height
                
        fig = plt.figure(figsize=[fig_width, fig_height],dpi=100)
        patch_ax = fig.add_axes([0, 0, 1, 1],projection='rectilinear')
        patch_ax.get_xaxis().set_visible(False)
        patch_ax.get_yaxis().set_visible(False)
        patch_ax.set_frame_on(False)
        
        #print patch_ax.get_xlim(), patch_ax.get_ylim()
        
        left_pad = 0.7
        top_pad = 0.9
        
        sec_index = 0
        sec_start_height = 0
        for sid, layered_section in self.sys.geometry.layered_sections.iteritems():
            sec_patches = []
            it = np.nditer(layered_section, flags=['multi_index'])
            
            while not it.finished:
                curr_ele_left_frac = (left_pad + it.multi_index[1]*ele_width)/fig_width
                curr_ele_bottom_frac = (fig_height - top_pad - ele_height*(it.multi_index[0] + 1) - sec_start_height)/fig_height
                
                curr_plot_left_frac = (left_pad + plot_padding + it.multi_index[1]*(plot_width + 2.0*plot_padding))/fig_width
                curr_plot_bottom_frac = (fig_height - top_pad + plot_padding - (plot_height + 2.0*plot_padding)*(it.multi_index[0] + 1) - sec_start_height)/fig_height
                
                #print it[0], curr_bottom_frac
                
                curr_ax = fig.add_axes([curr_plot_left_frac, curr_plot_bottom_frac, plot_width_frac, plot_height_frac],projection='rectilinear')
                
                #print it[0]
                curr_orbit = element_maps[int(it[0])]
                curr_ax.scatter(curr_orbit[0],curr_orbit[1],edgecolors='none', color='k', alpha=0.5, s=2)
                
                #print it.multi_index, layered_section.shape
                
                xpad = 0.1
                ypad = 0.1
                curr_ax.set_xlim((0 - xpad,1 + xpad))
                curr_ax.set_ylim((0 - ypad, 1 + ypad))
                
                x_on = False
                y_on = False
                
                if it.multi_index[0] == layered_section.shape[0] - 1:
                    x_on = True
                if it.multi_index[1] == 0:
                    y_on = True
                
                if x_on:
                    curr_ax.set_xlabel('CFF(t)',family=self.font_family,size=9)
                    
                    for tick in curr_ax.xaxis.get_major_ticks():
                        tick.label1.set_fontname(self.font_family)
                        tick.label1.set_fontsize(7)
                    
                    lat = self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lat
                    lon = self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lon
                    
                    patch_ax.axvline(x=curr_ele_left_frac,ymin=(fig_height - sec_start_height-top_pad+0.05)/fig_height,ymax=(fig_height - sec_start_height-top_pad+0.13)/fig_height,c='k',lw=0.5)
                    
                    patch_ax.text(curr_ele_left_frac,(fig_height - sec_start_height-top_pad+0.2)/fig_height,'%.2f N %.2f W'%(lat,lon),family=self.font_family,size=7, va='bottom', ha='center')
                    
                    #print self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lat, self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lon
                    
                else:
                    curr_ax.set_xticklabels([])
                
                if y_on:
                    curr_ax.set_ylabel('CFF(t+1)',family=self.font_family,size=9)
                    
                    for tick in curr_ax.yaxis.get_major_ticks():
                        tick.label1.set_fontname(self.font_family)
                        tick.label1.set_fontsize(7)
                    
                    depth = self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].depth
                    
                    curr_ele_top_frac = (fig_height - top_pad - ele_height*(it.multi_index[0]) - sec_start_height)/fig_height
                    
                    patch_ax.axhline(y=curr_ele_top_frac,xmin=0.95,xmax=0.955,c='k',lw=0.5)
                    
                    patch_ax.text(0.965,curr_ele_top_frac,'%.2f m'%(depth),family=self.font_family,size=7, va='center', ha='left', rotation='vertical')
                    
                    #print self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].depth
                else:
                    curr_ax.set_yticklabels([])
                
                curr_element = self.sys.geometry.elements[int(it[0])]
                
                #mu = (curr_orbit[1][-1] - curr_orbit[1][-2])/(curr_orbit[0][-1] - curr_orbit[0][-2])
                #print mu
                sigma_string = r'$\sigma_f$'
                #mu_string = r'$\mu$'
                cell_label = '%i v=%.3g %s=%.3g'%(int(it[0]), curr_element.slip, sigma_string, curr_element.static_strength)
                curr_ax.text(0,1,cell_label,family=self.font_family,size=8, va='top')
                
                sec_patches.append(mpl.patches.Rectangle((curr_ele_left_frac,curr_ele_bottom_frac), ele_width_frac, ele_height_frac))
                #print curr_left_frac, curr_bottom_frac
                it.iternext()
            
            collection = mpl.collections.PatchCollection(sec_patches,alpha=0.2)
            patch_ax.add_collection(collection)
            
            #print layered_section.shape
            sec_start_height += layered_section.shape[0]*ele_height + 1.0
            
            sec_index += 1
        
        patch_ax.text(.01,(fig_height - 0.1)/fig_height,'Map of the normalized CFF for the %s system'%(self.sys.name),family=self.font_family,size=12, va='top', ha='left')
            
        fig.savefig("%s"%(self.out_file), format=self.output_format)
        
class ElementCFFBlockEntropyPlotter(ElementStatePlotter):
    def __init__(self, sys):
        print '    **** Initilizing Element CFF block entropy Plotter'
        self.type = None
        super(ElementCFFBlockEntropyPlotter,self).__init__(sys)
        
    def plot(self, state_data):
        print '        **** Plotting CFF H/Delta H vs word length'
        #prepare all of the orbits
        #out_events = self.CFFOutEvents()
        #np_events = np.array( out_events[1::] )
        #ele_index = 0
        #element_data = state_data.element_data
        
        
        
        #prepare the plot
        if self.type == 'dbe':
            Y_label = r'$\Delta$H(L)'
        elif self.type == 'be':
            Y_label = 'H(L)'
        
        min_Y = sys.float_info.max
        max_Y = - sys.float_info.max
        
        #find the extrema Y values
        for eid, element_data in state_data.element_data.iteritems():
            if self.type == 'dbe':
                if max(element_data['Delta_H']) > max_Y:
                    max_Y = max(element_data['Delta_H'])
                if min(element_data['Delta_H']) < min_Y:
                    min_Y = min(element_data['Delta_H'])
            elif self.type == 'be':
                if max(element_data['H']) > max_Y:
                    max_Y = max(element_data['H'])
                if min(element_data['H']) < min_Y:
                    min_Y = min(element_data['H'])
            

        ele_width = 2.5
        ele_height = 2.5
        
        plot_padding = 0.0
        plot_width = ele_width - 2.0*plot_padding
        plot_height = ele_height - 2.0*plot_padding
        
        total_vert_elements = 0
        max_horizontal_elements = 0
        
        for layered_section in self.sys.geometry.layered_sections.itervalues():
            total_vert_elements += layered_section.shape[0]
            if layered_section.shape[1] > max_horizontal_elements:
                max_horizontal_elements = layered_section.shape[1]
        
        fig_width = max_horizontal_elements * ele_width + 1.5
        fig_height = total_vert_elements * ele_height + (len(self.sys.geometry.layered_sections) - 1.0) * 1.0 + 1.5
        
        ele_width_frac = ele_width/fig_width
        ele_height_frac = ele_width/fig_height
        
        plot_width_frac = plot_width/fig_width
        plot_height_frac = plot_width/fig_height
                
        fig = plt.figure(figsize=[fig_width, fig_height],dpi=100)
        patch_ax = fig.add_axes([0, 0, 1, 1],projection='rectilinear')
        patch_ax.get_xaxis().set_visible(False)
        patch_ax.get_yaxis().set_visible(False)
        patch_ax.set_frame_on(False)
        
        #print patch_ax.get_xlim(), patch_ax.get_ylim()
        
        left_pad = 0.7
        top_pad = 0.9
        
        sec_index = 0
        sec_start_height = 0
        for sid, layered_section in self.sys.geometry.layered_sections.iteritems():
            sec_patches = []
            it = np.nditer(layered_section, flags=['multi_index'])
            
            while not it.finished:
                curr_ele_left_frac = (left_pad + it.multi_index[1]*ele_width)/fig_width
                curr_ele_bottom_frac = (fig_height - top_pad - ele_height*(it.multi_index[0] + 1) - sec_start_height)/fig_height
                
                curr_plot_left_frac = (left_pad + plot_padding + it.multi_index[1]*(plot_width + 2.0*plot_padding))/fig_width
                curr_plot_bottom_frac = (fig_height - top_pad + plot_padding - (plot_height + 2.0*plot_padding)*(it.multi_index[0] + 1) - sec_start_height)/fig_height
                
                #print it[0], curr_bottom_frac
                
                curr_ax = fig.add_axes([curr_plot_left_frac, curr_plot_bottom_frac, plot_width_frac, plot_height_frac],projection='rectilinear')
                
                #print it[0]
                curr_orbit = state_data.element_data[int(it[0])]
                if self.type == 'dbe':
                    curr_ax.plot(curr_orbit['W_L'],np.around(curr_orbit['Delta_H'],3))
                elif self.type == 'be':
                    curr_ax.plot(curr_orbit['W_L'],np.around(curr_orbit['H'],3))
                
                
                #print it.multi_index, layered_section.shape
                
                xpad = (self.max_word_length - 2) * 0.05
                ypad = (max_Y - min_Y) * 0.05
                curr_ax.set_xlim((2 - xpad,self.max_word_length + xpad))
                curr_ax.set_ylim((min_Y - ypad, max_Y + ypad))
                
                x_on = False
                y_on = False
                
                if it.multi_index[0] == layered_section.shape[0] - 1:
                    x_on = True
                if it.multi_index[1] == 0:
                    y_on = True
                
                if x_on:
                    curr_ax.set_xlabel('L',family=self.font_family,size=9)
                    
                    if len(curr_orbit['W_L']) < 20:
                        stride = 1
                    else:
                        stride = int(float(len(curr_orbit['W_L']))/10.0)
                    
                    curr_ax.set_xticks(range(curr_orbit['W_L'][0], curr_orbit['W_L'][-1] + 1, stride))
                    
                    for tick in curr_ax.xaxis.get_major_ticks():
                        tick.label1.set_fontname(self.font_family)
                        tick.label1.set_fontsize(7)
                    
                    lat = self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lat
                    lon = self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lon
                    
                    patch_ax.axvline(x=curr_ele_left_frac,ymin=(fig_height - sec_start_height-top_pad+0.05)/fig_height,ymax=(fig_height - sec_start_height-top_pad+0.13)/fig_height,c='k',lw=0.5)
                    
                    patch_ax.text(curr_ele_left_frac,(fig_height - sec_start_height-top_pad+0.2)/fig_height,'%.2f N %.2f W'%(lat,lon),family=self.font_family,size=7, va='top', ha='center')
                    
                    #print self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lat, self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].lon
                    
                else:
                    curr_ax.set_xticklabels([])
                
                if y_on:
                    curr_ax.set_ylabel(Y_label,family=self.font_family,size=9)
                    
                    for tick in curr_ax.yaxis.get_major_ticks():
                        tick.label1.set_fontname(self.font_family)
                        tick.label1.set_fontsize(7)
                    
                    depth = self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].depth
                    
                    curr_ele_top_frac = (fig_height - top_pad - ele_height*(it.multi_index[0]) - sec_start_height)/fig_height
                    
                    patch_ax.axhline(y=curr_ele_top_frac,xmin=0.95,xmax=0.955,c='k',lw=0.5)
                    
                    patch_ax.text(0.965,curr_ele_top_frac,'%.2f m'%(depth),family=self.font_family,size=7, va='center', ha='left', rotation='vertical')
                    
                    #print self.sys.geometry.nodes[self.sys.geometry.elements[int(it[0])].enode_ids[0]].depth
                else:
                    curr_ax.set_yticklabels([])
                
                curr_element = self.sys.geometry.elements[int(it[0])]
                
                h_mu = curr_orbit['h_mu']
                E = curr_orbit['E']
                #print mu
                sigma_string = r'$\sigma_f$'
                h_mu_string = r'$h_{\mu}$'
                E_string = r'$\bf E$'
                cell_label = '%i v=%.3g %s=%.3g %s~%.3g'%(int(it[0]), curr_element.slip, sigma_string, curr_element.static_strength, h_mu_string, h_mu)
                cell_label_2 =  '%s~%.3g'%(E_string, E)
                curr_ax.text(2,max_Y,cell_label,family=self.font_family,size=8, va='top')
                curr_ax.text(2,max_Y-max_Y*0.05,cell_label_2,family=self.font_family,size=8, va='top')
                
                sec_patches.append(mpl.patches.Rectangle((curr_ele_left_frac,curr_ele_bottom_frac), ele_width_frac, ele_height_frac))
                #print curr_left_frac, curr_bottom_frac
                it.iternext()
            
            collection = mpl.collections.PatchCollection(sec_patches,alpha=0.2)
            patch_ax.add_collection(collection)
            
            sec_start_height += layered_section.shape[0]*ele_height + 1.0
            
            sec_index += 1
            
            
        patch_ax.text(.01,(fig_height - 0.1)/fig_height,'Permutation entropy H versus word length L for the %s'%(self.sys.name),family=self.font_family,size=12, va='top', ha='left')
            
        fig.savefig("%s"%(self.out_file), format=self.output_format)

class EventAnimationPlotter(object):
    def __init__(self, sys):
        print '    **** Initilizing Event Animator'
        self.sys    = sys
        self.name = None
        self.out_file = None
        #self.out_text_file = None
        self.events = []
        self.max_allowed_year    = 9012.0
        self.start_year          = 2012.0
        self.fade_steps = 6
        self.modified_full_year_range = None
        self.allowed_year_range = None
    
    def addEvent(self, event_data, sweeps):
        involved_elements = {}
        for sweep in sweeps:
            eid = sweep[2]
            try:
                involved_elements[eid] += sweep[3]
            except KeyError:
                involved_elements[eid] = sweep[3]
        self.events.append([event_data[0],event_data[1],event_data[3],involved_elements])
    
    def plot(self):
        print '*** Exporting KML event animation to file %s ***'%(self.out_file)
        f = open(self.out_file, 'w')
        
        self.allowed_year_range = self.max_allowed_year - self.start_year
        full_year_range = self.events[-1][1] - self.events[0][1]
        if full_year_range > self.allowed_year_range:
            self.modified_full_year_range = math.modf(full_year_range/self.allowed_year_range)[1] * self.allowed_year_range
        else:
            self.modified_full_year_range = full_year_range
        
        self.time_step = self.modified_full_year_range/500.0
        
        # write header
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        f.write('<Document>\n')
        f.write('<name>%s event animation</name>\n'%self.name)
        f.write(self.sys.geometry.KMLStyle_section())
        f.write(self.sys.geometry.KMLStyle_base())
        f.write(self.sys.geometry.KMLStyle_broke())
        f.write(self.sys.geometry.KMLStyle_elementBreak(self.fade_steps))
        f.write('<Folder>\n')
        f.write('\t<name>Fault Section Names</name>\n')
        for s in self.sys.geometry.sections.itervalues():
            s.exportKMLSectionLabel(f)
        f.write('</Folder>\n')
        f.write('<Folder>\n')
        f.write('\t<name>Fault Sections</name>\n')
        for s in self.sys.geometry.sections.itervalues():
            s.exportKMLGeometry(f)
        f.write('</Folder>\n')
        f.write('<Folder>\n')
        f.write('\t<name>Events</name>\n')
        
        for num,event in enumerate(self.events):
            if ( num%500 == 0 ):
                    print '    event:%i of %i'%(event[0], self.events[-1][0])
            self.plotEvent(f,event)
        f.write('</Folder>\n')
            
        f.write('</Document>\n')
        f.write('</kml>')
        f.close()
        print '*** Done ***'
    
    def plotEvent(self,f,event):
        curr_year = event[1] - self.events[0][1] + self.start_year
        
        if curr_year < self.start_year + self.modified_full_year_range:
            
            if curr_year > self.max_allowed_year:
                mod_curr_year = curr_year
                while mod_curr_year > self.max_allowed_year:
                    mod_curr_year -= self.allowed_year_range
            else:
                mod_curr_year = curr_year
            
            f.write('<Folder>\n')
            f.write('       <name>Event %i on Year %f (%f)</name>\n'%(event[0],curr_year,event[1]))

            f.write('       <Style>\n')
            f.write('               <ListStyle>\n')
            f.write('                       <listItemType>checkHideChildren</listItemType>\n')
            f.write('               </ListStyle>\n')
            f.write('       </Style>\n')
            
            involved_elements = event[3]
                    
            year, month, day = self.sys.geometry.converter.yearDecimalToYearMonthDay(mod_curr_year)
            year = math.modf(curr_year)[1]
            
            for eid,eslip in involved_elements.iteritems():
                    for num in range(0,self.fade_steps + 1):
                            
                        year_end, month_end, day_end = self.sys.geometry.converter.yearDecimalToYearMonthDay(mod_curr_year + self.time_step*float(num+1))
                        year_end = math.modf(curr_year + self.time_step*float(num+1))[1]
                        
                        self.sys.geometry.elements[eid].exportKML(f,'elementBreakStyle%i_%i'%(int(round(event[2])),num),(year, month, day),(year_end, month_end, day_end))
                        
                        year, month, day = year_end, month_end, day_end
            
            f.write('</Folder>\n')
        
class RuptureMapPlotter(object):
    def __init__(self, sys):
        print '    **** Initilizing Rupture Map'
        self.sys    = sys
        
        # plotting parameters for the rupture map
        imw = 1010.0 # the full image width
        imh = 1024.0 # the full image height
        lm = 70.0
        rm = 10.0
        tm = 170.0
        bm = 10.0
        res = 72.0
        
        self.ticklabelfont = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9)
        self.framelabelfont = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=10)
        self.legendfont = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9)
        self.legendticklabelfont = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9)
        self.sectionlabelfont = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9)
        
        imwi = imw/res
        imhi = imh/res
        pw = imw - lm - rm
        ph = imh - tm - bm
        self.fig    = plt.figure(figsize=(imwi, imhi), dpi=res)
        self.the_ax = self.fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
        
        self.output_format = 'pdf'
        
        
        self.sections           = None
        self.layered_sections   = None
        self.max_magnitude      = None
        self.min_magnitude      = None
        self.magnitude_slope    = None
        self.start_year         = None
        self.end_year           = None
        self.max_y_offset       = None
        self.max_x              = None
        self.y_locs             = []
        self.y_labels           = []
        self.num_of_events      = None
        self.curr_y             = None
        self.event_index        = 0
        
        #data needed to plot the events
        self.events             = None
        self.section_mag_extrema= None
        
        self.out_file           = None
        self.name               = None
        self.res                = res
        self.cmap = plt.get_cmap('autumn_r',lut=1000)
    
    def setData(self, raw_data):
        self.events = raw_data['events']
        self.section_mag_extrema = raw_data['section_mag_extrema']
        
    def loopElements(self, elements_in_event, event_magnitude):
        (r,g,b,a) = self.cmap(self.magnitude_slope * (event_magnitude - self.min_magnitude))
        section_x_offset = 0
        #the_ys = {}
        for sid in sorted(self.layered_sections.iterkeys()):
        #for sindex, s in enumerate(sections):
            it = np.nditer(self.layered_sections[sid], flags=['multi_index'])
            while not it.finished:
                #print it.multi_index, it[0]
                #try:
                    #eid = elements_in_event[int(it[0])]
                if int(it[0]) in elements_in_event:
                    the_x_start = float(it.multi_index[1] + section_x_offset)/float(self.max_x)
                    the_x_end = float(it.multi_index[1] + section_x_offset + 1)/float(self.max_x)
                    the_y = self.curr_y - it.multi_index[0]
                    #try:
                    #    the_ys[the_y].append((the_x_start, the_x_end))
                    #except KeyError:
                    #    the_ys[the_y] = [(the_x_start, the_x_end)]
                    #the_y = self.curr_y - it.multi_index[0]
                    #the_x_start = float(it.multi_index[1] + section_x_offset)/float(self.max_x)
                    #the_x_end = float(it.multi_index[1] + section_x_offset + 1)/float(self.max_x)
                    #print the_y, the_x_start, the_x_end
                    #print the_y, the_x_start, the_x_end
                    self.the_ax.axhline(y=the_y, xmin=the_x_start, xmax=the_x_end, color=(r,g,b))
                #except:
                #    pass
                
                #if e.elementInEvent(it[0]):
                #    the_y = curr_y - it.multi_index[0]
                #    the_x_start = float(it.multi_index[1] + section_x_offset)/float(max_x)
                #    the_x_end = float(it.multi_index[1] + section_x_offset + 1)/float(max_x)
                #    the_ax.axhline(y=the_y, xmin=the_x_start, xmax=the_x_end, color=(r,g,b))
                it.iternext()
            
            section_x_offset += self.layered_sections[sid].shape[1]
        
        #self.plot_data.append(the_ys)

    
    def addEvent(self, event_data, event_sweeps):
        event_number        = event_data[0]
        event_year          = event_data[1]
        #event_trigger       = event_data[2]
        event_magnitude     = event_data[3]
        #event_shear_init    = event_data[4]
        #event_normal_init   = event_data[5]
        #event_shear_final   = event_data[6]
        #event_normal_final  = event_data[7]
        #event_start_sweep   = event_data[8]
        #event_end_sweep     = event_data[9]
        #event_start_as      = event_data[10]
        #event_end_as        = event_data[11]
        
        
        
        if event_magnitude >= self.min_magnitude:
            if self.num_of_events < 200:
                self.the_ax.axhline(y=self.curr_y, alpha=0.2, linewidth=1, color='0.5')
            
            if self.num_of_events > 20:
                if (self.event_index%int(float(self.num_of_events)/10) == 0):
                    self.y_locs.append(self.curr_y)
                    self.y_labels.append('%.1f\n[%i]'%(event_year,event_number))
            else:
                self.y_locs.append(self.curr_y)
                self.y_labels.append('%.1f\n[%i]'%(event_year,event_number))
            
            elements_in_event = {}
            
            for sweep in event_sweeps:
                elements_in_event[int(sweep[2])] = int(sweep[2])
            
            self.loopElements(elements_in_event.keys(), event_magnitude)
                        
            self.curr_y += self.max_y_offset
        
        self.event_index += 1
    
    def plot(self):
        
        
        self.the_ax.autoscale(enable=True, axis='both', tight=True)
        
        label_lines = []
        section_x_start = 0
        section_x_finish = 0
        #print '    finishing'
        #for ns, s in enumerate(sections):
        for ns, sid in enumerate(sorted(self.layered_sections.iterkeys())):
            layered_section = self.layered_sections[sid]
            section_x_finish += self.layered_sections[sid].shape[1]
            
            x_text_loc = (float(ns)/float(len(self.layered_sections)))
            x_fault_loc = float(section_x_finish + section_x_start)/(2.0*float(self.max_x))
            
            section_x_start += self.layered_sections[sid].shape[1]
            
            y_text_distance = .05
            
            self.the_ax.axvline(x=section_x_finish, alpha=0.2, linewidth=1, color='0.5')
            self.the_ax.text(x_text_loc, 1+y_text_distance, self.sys.geometry.sections[sid].sname,
                horizontalalignment='center',
                verticalalignment='bottom',
                rotation=90,
                fontproperties=self.sectionlabelfont,
                transform=self.the_ax.transAxes)
            
            line_xs = [x_text_loc, x_text_loc, x_fault_loc, x_fault_loc]
            line_ys = [1+y_text_distance-0.005, 1+y_text_distance*(3.0/4.0), 1 + y_text_distance/4.0, 1]
            label_lines.append(mlines.Line2D(
                line_xs, line_ys,transform=self.the_ax.transAxes, solid_capstyle='round', alpha=0.2, linewidth=1, color='0.2'))
            
        self.the_ax.yaxis.set_ticks(self.y_locs)
        self.the_ax.yaxis.set_ticklabels(self.y_labels)
        
        x_locs = self.the_ax.yaxis.get_majorticklocs()
        
        x_tick_labels = []
        
        for x_loc in x_locs:
            x_tick_labels.append('')
        
        self.the_ax.xaxis.set_ticklabels(x_tick_labels)
        
        #mpl.artist.getp(self.the_ax.yaxis.get_clip_box())
        
        #points = self.the_ax.yaxis.get_clip_box().get_points()
        
        #points[0][0] = 60.0
        
        #self.the_ax.yaxis.set_clip_box(mpl.transforms.Bbox(points))
        
        #mpl.artist.getp(self.the_ax.yaxis.get_clip_box())
        
        #self.the_ax.yaxis.set_clip_on(False)
        
        for label in self.the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(self.ticklabelfont)
            #label.clip_on = False
            #mpl.artist.getp(label)
        for line in self.the_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        self.the_ax.autoscale(enable=True, axis='both', tight=True)
        
        self.the_ax.lines.extend(label_lines)
        
        
        self.the_ax.set_ylabel('Year [event id]',fontproperties=self.framelabelfont)
        
        #self.the_ax.tick_params(labelsize=6)
        #self.the_ax.tight_layout()
        #print self.out_file
        self.the_ax.set_title('%s'%self.sys.name, position=(-0.05,1.34), ha='left', size=12)
    
        '''
        def on_draw(event):
           bboxes = []
           for label in labels:
               bbox = label.get_window_extent()
               # the figure transform goes from relative coords->pixels and we
               # want the inverse of that
               bboxi = bbox.inverse_transformed(fig.transFigure)
               bboxes.append(bboxi)

           # this is the bbox that bounds all the bboxes, again in relative
           # figure coords
           bbox = mtransforms.Bbox.union(bboxes)
           if fig.subplotpars.left < bbox.width:
               # we need to move it over
               fig.subplots_adjust(left=1.1*bbox.width) # pad a little
               fig.canvas.draw()
           return False
   
        fig.canvas.mpl_connect('draw_event', on_draw)
        
        self.fig.subplots_adjust(left=2.0)
        '''
        
        self.fig.savefig(self.out_file, format=self.output_format, dpi=self.res)
    
    def createPlotLegend(self):
        if self.max_magnitude is None or self.min_magnitude is None:
            print '!!! Set the min and max magnitude with RuptureMapPlotter.setMinMaxMagnitude(min_magnitude, max_magnitude)!!!'
            return
        
        ax1    = self.fig.add_axes([.77, 0.98, 0.2, 0.01])
        norm = mcolor.Normalize(vmin=self.min_magnitude, vmax=self.max_magnitude)
        cb1 = mcolorbar.ColorbarBase(ax1, cmap=self.cmap,
                   norm=norm,
                   orientation='horizontal')
        ticks = map(float, [self.min_magnitude,(self.min_magnitude+ self.max_magnitude)/2.0, self.max_magnitude])
        cb1.set_ticks([ticks[0], ticks[1], ticks[2]])
        cb1.set_ticklabels(['%.2f'%ticks[0],'%.2f'%ticks[1], '%.2f'%ticks[2]])

        for label in ax1.xaxis.get_ticklabels():
            label.set_fontproperties(self.legendticklabelfont)
            
        for line in ax1.xaxis.get_ticklines():
            line.set_alpha(0)
            
        cb1.ax.set_xlabel('Magnitude',position=(1,0), ha='right', fontproperties=self.legendfont)
    
    def setLayeredSections(self, layered_sections):
        self.layered_sections = layered_sections
        
        #print layered_sections
        
        self.max_y_offset = 0
        self.max_x = 0

        for k,s in self.layered_sections.iteritems():
            #print k,s.shape
            if s.shape[0] > self.max_y_offset:
                self.max_y_offset = s.shape[0]
            self.max_x += s.shape[1]
        
        self.curr_y = self.max_y_offset
    
    def setMinMaxMagnitude(self, min_magnitude, max_magnitude):
        self.max_magnitude = max_magnitude
        self.min_magnitude = min_magnitude
        
        self.magnitude_slope = 1.0/(max_magnitude - min_magnitude)
        
        self.createPlotLegend()

class BasicPlotter(object):
    def __init__(self,sys):
        self.output_format = 'png'
        self.sys = sys
        self.event_count = 0
    
        # plotting parameters that are common to all the scaling relations
        self.imw = 501.0 # the full image width
        self.imh = 501.0 # the full image height
        self.lm = 45.0
        self.rm = 10.0
        self.tm = 10.0
        self.bm = 45.0
        self.res = 72.0
        
        self.ticklabelfont = matplotlib.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.framelabelfont = matplotlib.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        self.legendfont = matplotlib.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        self.titlefont = matplotlib.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
        #arial8_light = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=8, weight='light')
    
    def addEvent(self):
        self.event_count += 1
    
    def calculateAverages(self):
        num_bins = math.floor(len(self.x)/100)
        
        if num_bins < 20:
            num_bins = 20
        elif num_bins > 100:
            num_bins = 100
        
        x = np.array(self.x)
        y = np.array(self.y)
        
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
        
        for k in sorted(binned_data.keys()):
            self.averages.append(( 0.5*(bins[k-1]+bins[k]),sum(binned_data[k])/float(len(binned_data[k])) ))
        
        
    
    def basicPlot(self, x_label = None, y_label = None, title = None, out_file = None, axis_scale = None, extra_x = None, extra_y = None, extra_title = None, connect_points = False, extra_x_error = None, extra_y_error = None):
        
    # set up the plot
        imw = 503.0 # the full image width
        imh = 503.0 # the full image height
        lm = 40.0
        rm = 10.0
        tm = 10.0
        bm = 40.0
        res = 72.0
        
        con_line = '--'
        
        arial12 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=12)
        arial10 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=10)
        arial8_light = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=8, weight='light')
        
    # do the plot
        imwi = imw/res
        imhi = imh/res
        pw = imw - lm - rm
        ph = imh - tm - bm
        fig = plt.figure(figsize=(imwi, imhi), dpi=res)
        the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
        
        x = self.x
        y = self.y
        
        if connect_points:
            ls1 = '--'
        else:
            ls1 = 'None'
        lw1 = 0.5
        c1 = (0.5,0.5,0.5,0.5)
        mfc1 = (0,0,0,1)
        ms1 = 3
        marker1 = 'o'
        
        ls_extra = '-'
        lw_extra = 0.5
        c_extra = (0.0,0.0,1,1)
        
        if axis_scale is None:
            the_ax.plot(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                the_ax.plot(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                
        elif axis_scale is 'semilogy':
            the_ax.semilogy(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                #print extra_x_error, extra_y_error
                if extra_x_error is not None or extra_y_error is not None:
                    the_ax.errorbar(extra_x, extra_y, yerr = extra_y_error, xerr = extra_x_error, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                else:
                    the_ax.semilogy(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                
        elif axis_scale is 'semilogx':
            the_ax.semilogx(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                the_ax.semilogx(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
        
        elif axis_scale is 'loglog':
            the_ax.loglog(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                the_ax.loglog(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
        
        for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(arial10)
            
        the_ax.set_ylabel(y_label, fontproperties=arial12)
        the_ax.set_xlabel(x_label, fontproperties=arial12)
        
        #the_ax.tick_params(fontproperties=arial12)
        
        the_ax.margins(0.05,0.3)
        
        #the_ax.set_title(title, position=(0.0,1.05), ha='left', fontproperties=arial12)
        
        the_ax.legend()
           
        if out_file is not None:
            plt.savefig(out_file, format=self.output_format)

class  FrequencyMagnitudePlotter(BasicPlotter):
    def __init__(self, sys):
        print '    **** Initilizing Frequency Magnitude'
        self.number_of_years = None
        self.name = None
        self.out_file = None
        self.out_text_file = None
        self.magnitudes = []
        super(FrequencyMagnitudePlotter,self).__init__(sys)
        
    def addEvent(self, magnitude):
        super(FrequencyMagnitudePlotter,self).addEvent()
        self.magnitudes.append(magnitude)
    
    def plot(self):
        self.x = []
        self.y = []
        
        # for the UCERF2 error bars
        extra_x = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        extra_y = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]    
        extra_y_error = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        
        cum_freq = {}
        total_events = len(self.magnitudes)
        
        if self.out_text_file is not None:
            f = open(self.out_text_file,'w')
        
        #print self.number_of_years
        
        for num, magnitude in enumerate(sorted(self.magnitudes)):
            cum_freq['%.10f'%magnitude] = total_events - (num + 1)
        
        for magnitude in sorted(cum_freq.iterkeys()):
            if self.out_text_file is not None:
                f.write('%f %f\n'%(float(magnitude), float(cum_freq[magnitude])/self.number_of_years))
            #print magnitude, float(cum_freq[magnitude])/self.number_of_years
            if self.out_file is not None:
                self.x.append(magnitude)
                self.y.append(float(cum_freq[magnitude])/self.number_of_years)
        
        if self.out_text_file is not None:
            f.close()
        
        # do the plot
        if self.out_file is not None:
            title = '%s Frequency-Magnitude'%self.name
            
            # plotting parameters that are common to all the scaling relations
            imw = self.imw
            imh = self.imh
            lm = self.lm
            rm = self.rm
            tm = self.tm
            bm = self.bm
            res = self.res

        
            imwi = imw/res
            imhi = imh/res
            pw = imw - lm - rm
            ph = imh - tm - bm
            fig = plt.figure(figsize=(imwi, imhi), dpi=res)
            the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
            
            x = map(float, self.x)
            y = map(float, self.y)
            
            # create the line for b = 1
            x_b1 = np.linspace(np.min(x),np.max(x),10)
            y_b1 = 10**(math.log(y[0],10)+x[0]-x_b1)
            
            #y_b1 *= 10**np.mean(y)
            
            #if connect_points:
            ls1 = '--'
            #else:
            #    ls1 = 'None'
            lw1 = 1.0
            c1 = (0.7,0.7,0.7,1)
            mfc1 = (0,0,0,1)
            ms1 = 3
            marker1 = 'o'
            
            ls_extra = ':'
            lw_extra = 1.0
            c_extra = (0.5,0.5,0.5,1)
            
            the_ax.semilogy(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            the_ax.semilogy(x_b1, y_b1, ls = ls_extra, lw = lw_extra, c = c_extra, label = 'b = 1')
            
            for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
                label.set_fontproperties(self.ticklabelfont)
                
            the_ax.set_ylabel('Log(Cumulative Number of Events)', fontproperties=self.framelabelfont)
            the_ax.set_xlabel('Magnitude', fontproperties=self.framelabelfont)
            
            the_ax.autoscale_view(tight=True)
            
            #the_ax.tick_params(fontproperties=arial12)
            
            #the_ax.margins(0.05,0.3)
            
            #the_ax.set_title(title, position=(0.0,1.05), ha='left', fontproperties=arial12)
            
            the_ax.legend(prop=self.legendfont)
               
            plt.savefig(self.out_file, format=self.output_format, dpi=res)

            #self.basicPlot(x_label='Magnitude', y_label='Log(Cumulative Number of Events)', title=title, out_file = self.out_file, axis_scale = 'semilogy', connect_points = True, extra_x = extra_x, extra_y = extra_y, extra_y_error = extra_y_error, extra_title = 'UCERF2 Observed 95% Confidence')

class MagnitudeRuptureAreaPlotter(BasicPlotter):
    def __init__(self, sys):
        print '    **** Initilizing Magnitude Rupture Area'
        self.name = None
        self.out_file = None
        self.out_text_file = None
        
        self.x = []
        self.y = []
        self.averages = []
        super(MagnitudeRuptureAreaPlotter,self).__init__(sys)
    
    def addEvent(self, magnitude, rupture_area):
        '''
        self.y.append(event_data[3])
        
        ruptured_elements = {}
        
        rupture_area = 0.0
        
        for sweep in event_sweeps:
            eid = sweep[2]
            try:
                tmp = ruptured_elements[eid]
            except KeyError:
                rupture_area += VCGeometry.elements[eid].sideLength() ** 2.0
                ruptured_elements[eid] = True
        
        self.x.append(rupture_area  * 1.0e-6)
        '''
        super(MagnitudeRuptureAreaPlotter,self).addEvent()
        self.y.append(magnitude)
        self.x.append(rupture_area)
        
    def plot(self):
        
        self.calculateAverages()
                
        if self.out_text_file is not None:
            f = open(self.out_text_file,'w')
            for index, x in enumerate(self.x):
                if index < len(self.averages):
                    f.write('%f %f %f %f\n'%(x, self.y[index], self.averages[index][0], self.averages[index][1]))
                else:
                    f.write('%f %f\n'%(x, self.y[index]))
            f.close()
        
        if self.out_file is not None:
            #wc = 4.07 + 0.98 * np.log10(self.x)
            
            title = '%s Magnitude-Rupture_Area'%self.name
            #self.basicPlot(x_label=r'log(Rupture Area [km$^\mathsf{2}$])', y_label='Magnitude', title=title, out_file = self.out_file, axis_scale = 'semilogx', extra_x = self.x, extra_y = wc, extra_title = 'Wells-Coppersmith')
            # plotting parameters that are common to all the scaling relations
            imw = self.imw
            imh = self.imh
            lm = self.lm
            rm = self.rm
            tm = self.tm
            bm = self.bm
            res = self.res

        
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

class MagnitudeAverageSlipPlotter(BasicPlotter):
    def __init__(self, sys):
        print '    **** Initilizing Magnitude Average Slip'
        self.name = None
        self.out_file = None
        self.out_text_file = None
        
        self.x = []
        self.y = []
        self.averages = []
        super(MagnitudeAverageSlipPlotter,self).__init__(sys)
    
    def addEvent(self, magnitude, average_slip):
        '''
        self.y.append(event_data[3])
        
        ruptured_elements = {}
        
        total_slip = 0.0
        
        for sweep in event_sweeps:
            total_slip += sweep[3]
            eid = sweep[2]
            try:
                tmp = ruptured_elements[eid]
            except KeyError:
                ruptured_elements[eid] = True
        
        self.x.append(total_slip/float(len(ruptured_elements)))
        '''
        super(MagnitudeAverageSlipPlotter,self).addEvent()
        self.y.append(magnitude)
        self.x.append(average_slip)
    
    def plot(self):
        
        self.calculateAverages()

        #output the data
        if self.out_text_file is not None:
            f = open(self.out_text_file,'w')
            for index, x in enumerate(self.x):
                if index < len(self.averages):
                    f.write('%f %f %f %f\n'%(x, self.y[index], self.averages[index][0], self.averages[index][1]))
                else:
                    f.write('%f %f\n'%(x, self.y[index]))
            f.close()
        
        if self.out_file is not None:
            #wc = 6.93 + 0.82 * np.log10(self.x)
            
            title = '%s Magnitude-Average_Slip'%self.name
            #self.basicPlot(x_label='log(Average Slip [m])', y_label='Magnitude', title=title, out_file = self.out_file, axis_scale = 'semilogx', extra_x = self.x, extra_y = wc, extra_title = 'Wells-Coppersmith')
            # plotting parameters that are common to all the scaling relations
            imw = self.imw
            imh = self.imh
            lm = self.lm
            rm = self.rm
            tm = self.tm
            bm = self.bm
            res = self.res

        
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
            the_ax.set_xlabel('log(Average Slip [m])', fontproperties=self.framelabelfont)
            
            the_ax.autoscale_view(tight=True)
            
            #the_ax.tick_params(fontproperties=arial12)
            
            #the_ax.margins(0.05,0.3)
            
            #the_ax.set_title(title, position=(0.0,1.05), ha='left', fontproperties=arial12)
            
            the_ax.legend(prop=self.legendfont, loc=2)

            plt.savefig(self.out_file, format=self.output_format, dpi=res)

class AverageSlipRuptureLengthPlotter(BasicPlotter):
    def __init__(self, sys):
        print '    **** Initilizing Average Slip Rupture Length'
        self.name = None
        self.out_file = None
        self.out_text_file = None
        
        self.x = []
        self.y = []
        self.averages = []
        super(AverageSlipRuptureLengthPlotter,self).__init__(sys)
    
    def addEvent(self, average_slip, rupture_length):
        '''
        ruptured_elements = {}
        
        rupture_length = 0.0
        total_slip = 0.0
        
        for sweep in event_sweeps:
            total_slip += sweep[3]
            eid = sweep[2]
            try:
                tmp = ruptured_elements[eid]
            except KeyError:
                if VCGeometry.elements[eid].onTrace():
                    rupture_length += VCGeometry.elements[eid].sideLength()
                ruptured_elements[eid] = True
        
        if rupture_length != 0.0:
            self.x.append(rupture_length * 1.0e-3)
            self.y.append(total_slip/float(len(ruptured_elements)))
        '''
        super(AverageSlipRuptureLengthPlotter,self).addEvent()
        self.y.append(average_slip)
        self.x.append(rupture_length)
    
    def plot(self):
        
        self.calculateAverages()

        #output the data
        if self.out_text_file is not None:
            f = open(self.out_text_file,'w')
            for index, x in enumerate(self.x):
                if index < len(self.averages):
                    f.write('%f %f %f %f\n'%(x, self.y[index], self.averages[index][0], self.averages[index][1]))
                else:
                    f.write('%f %f\n'%(x, self.y[index]))
            f.close()
        
        if self.out_file is not None:
            #wc = 10 ** (-1.43 + 0.88 * np.log10(self.x))
        
            title = '%s Average_Slip-Rupture_Length'%self.name
            #self.basicPlot(x_label='log(Rupture Length [km])', y_label='log(Average Slip [m])', title=title, out_file = self.out_file, axis_scale = 'loglog', extra_x = self.x, extra_y = wc, extra_title = 'Wells-Coppersmith')
            # plotting parameters that are common to all the scaling relations
            imw = self.imw
            imh = self.imh
            lm = self.lm
            rm = self.rm
            tm = self.tm
            bm = self.bm
            res = self.res

        
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
            
            the_ax.loglog(x, y, ls = 'None', mfc = mfc1, ms = ms1, marker = marker1)
            
            the_ax.loglog(x_ave, y_ave, ls = ls_extra, lw = lw_extra, c = c_extra, label = 'binned average')
            
            for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
                label.set_fontproperties(self.ticklabelfont)
                
            the_ax.set_ylabel('log(Average Slip [m])', fontproperties=self.framelabelfont)
            the_ax.set_xlabel('log(Rupture Length [km])', fontproperties=self.framelabelfont)
            
            the_ax.autoscale_view(tight=True)
            
            #the_ax.tick_params(fontproperties=arial12)
            
            #the_ax.margins(0.05,0.3)
            
            #the_ax.set_title(title, position=(0.0,1.05), ha='left', fontproperties=arial12)
            
            the_ax.legend(prop=self.legendfont, loc=2)

            plt.savefig(self.out_file, format=self.output_format, dpi=res)

class SlipInSlipOutPlotter(BasicPlotter):
    def __init__(self, sys):
        print '    **** Initilizing Slip In Slip Out'
        self.name = None
        self.out_file = None
        self.out_text_file = None
        
        self.slip_in = {}
        self.slip_out = {}
        self.duration = 0.0
        #self.averages = []
        super(SlipInSlipOutPlotter,self).__init__(sys)

        self.calulateInputSlip()

    def calulateInputSlip(self):
    # get the slip rate from the elements. this is the input slip.
        
        print '        *** Getting input slip rate for all elements. ***'
        slip_in = {}
        
        for eid, ele in self.sys.geometry.elements.iteritems():
            if self.sys.section_filter is None or ele.sid in self.sys.section_filter:
                try:
                    slip = self.slip_in[ele.eid]
                    print '!!! Duplicate element id %i !!!'%eid
                    break
                except KeyError:
                    self.slip_in[ele.eid] = ele.slip

    def addSlip(self, eid, slip):
    # get the total slip from the elements. divide this by the length of the simulations and this is the output slip.
        
        if self.sys.section_filter is None or self.sys.geometry.elements[eid].sid in self.sys.section_filter:
            try:
                self.slip_out[eid] += slip
            except KeyError:
                self.slip_out[eid] = slip

    def plot(self):
        
        # sort by the input slip
        sorted_eids = sorted(self.slip_in, key=lambda k: self.slip_in[k])
        
        f = open(self.out_text_file,'w')
        f.write('element slip_in slip_out\n')
        
        for eid in sorted_eids:
            try:
                i = float(self.slip_in[eid])
                o = float(self.slip_out[eid])
            except KeyError:
                i = float(self.slip_in[eid])
                o = 0.0
            f.write('%i %g %g\n'%(eid, i, o/self.duration))
        f.close()

class  RecurrenceIntervalPlotter(BasicPlotter):
    def __init__(self, sys):
        print '    **** Initilizing Recurrence Intervals'
        self.name = None
        self.out_file_prepend = None
        self.out_text_file_prepend = None
        self.sections = {}
        super(RecurrenceIntervalPlotter,self).__init__(sys)
        
    def addEvent(self, magnitude, year, involved_sections):
        super(RecurrenceIntervalPlotter,self).addEvent()
        for sid in involved_sections:
            try:
                self.sections[sid].append((magnitude,year))
            except KeyError:
                self.sections[sid] = [(magnitude,year)]
    
    def plot(self):
        
        bins = np.linspace(0,250,50)
        for sid in sorted(self.sections.keys()):
            years = [x[1] for x in self.sections[sid]]
            years7 = [x[1] for x in self.sections[sid] if x[0] >= 7.0]
            intervals = []
            intervals7 = []
            
            for n, year in enumerate(years):
                if n != 0:
                    intervals.append(year - years[n-1])
            for n, year in enumerate(years7):
                if n != 0:
                    intervals7.append(year - years7[n-1])
            
            # do the plot
            if self.out_file_prepend is not None:
                title = '%s %s Recurrence Intervals'%(self.name, self.sys.geometry.sections[sid].sname)
                
                # plotting parameters that are common to all the scaling relations
                imw = self.imw
                imh = self.imh + 30.0
                lm = self.lm
                rm = self.rm
                tm = self.tm + 30.0
                bm = self.bm
                res = self.res
            
                imwi = imw/res
                imhi = imh/res
                pw = imw - lm - rm
                ph = imh - tm - bm
                fig = plt.figure(figsize=(imwi, imhi), dpi=res)
                the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
                
                # some data crunching
                hist7, bins7 = np.histogram(intervals7, bins=bins, density=True)
                hist, bins = np.histogram(intervals, bins=bins, density=True)
                mean = np.mean(intervals)
                std = np.std(intervals)
                mean7 = np.mean(intervals7)
                std7 = np.std(intervals7)
                    
                #print hist
                # plotting commands
                the_ax.step(bins[0:-1], hist, where='post', label='m>6.5')
                the_ax.step(bins7[0:-1], hist7, where='post', label='m>7')
                
                for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
                    label.set_fontproperties(self.ticklabelfont)
                    
                the_ax.set_ylabel('Prob. Density', fontproperties=self.framelabelfont)
                the_ax.set_xlabel('Recurrence Time [yr]', fontproperties=self.framelabelfont)
                
                the_ax.autoscale_view(tight=True)
                
                the_ax.set_title(title, position=(0.0,1.04), ha='left', fontproperties=self.titlefont)
                the_ax.text(0.0, 1.01, 'm>6.5: mean %.1f std %.1f, m>7.0 mean %.1f std %.1f'%(mean, std, mean7, std7), va='bottom', ha='left', transform=the_ax.transAxes, fontproperties=self.framelabelfont)
        
                the_ax.legend(prop=self.legendfont)
                   
                plt.savefig('%ss-%i-%s_recurrence-interval.%s'%(self.out_file_prepend, sid, self.sys.geometry.sections[sid].sname, self.output_format), format=self.output_format, dpi=res)
            # output the text file
            if self.out_text_file_prepend is not None:
                
                f = open('%ss-%i-%s_recurrence-interval.dat'%(self.out_text_file_prepend, sid, self.sys.geometry.sections[sid].sname),'w')
                
                f.write('m>6.5 m>7.0\n')
                
                if len(intervals7) > len(intervals):
                    for i, interval in enumerate(intervals7):
                        try:
                            f.write('%f %f\r'%(intervals[i], interval))
                        except IndexError:
                            f.write('0.0 %f\r'%(interval))
                else:
                    for i, interval in enumerate(intervals):
                        try:
                            f.write('%f %f\r'%(interval, intervals7[i]))
                        except IndexError:
                            f.write('%f 0.0\r'%(interval))
                
                f.close()

            

            
        
        '''
        self.x = []
        self.y = []
        
        # for the UCERF2 error bars
        extra_x = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        extra_y = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]    
        extra_y_error = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        
        cum_freq = {}
        total_events = len(self.magnitudes)
        
        if self.out_text_file is not None:
            f = open(self.out_text_file,'w')
        
        #print self.number_of_years
        
        for num, magnitude in enumerate(sorted(self.magnitudes)):
            cum_freq['%.10f'%magnitude] = total_events - (num + 1)
        
        for magnitude in sorted(cum_freq.iterkeys()):
            if self.out_text_file is not None:
                f.write('%f %f\n'%(float(magnitude), float(cum_freq[magnitude])/self.number_of_years))
            #print magnitude, float(cum_freq[magnitude])/self.number_of_years
            if self.out_file is not None:
                self.x.append(magnitude)
                self.y.append(float(cum_freq[magnitude])/self.number_of_years)
        
        if self.out_text_file is not None:
            f.close()
        
        # do the plot
        if self.out_file is not None:
            title = '%s Frequency-Magnitude'%self.name
            
            # plotting parameters that are common to all the scaling relations
            imw = self.imw
            imh = self.imh
            lm = self.lm
            rm = self.rm
            tm = self.tm
            bm = self.bm
            res = self.res

        
            imwi = imw/res
            imhi = imh/res
            pw = imw - lm - rm
            ph = imh - tm - bm
            fig = plt.figure(figsize=(imwi, imhi), dpi=res)
            the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
            
            x = map(float, self.x)
            y = map(float, self.y)
            
            # create the line for b = 1
            x_b1 = np.linspace(np.min(x),np.max(x),10)
            y_b1 = 10**(math.log(y[0],10)+x[0]-x_b1)
            
            #y_b1 *= 10**np.mean(y)
            
            #if connect_points:
            ls1 = '--'
            #else:
            #    ls1 = 'None'
            lw1 = 1.0
            c1 = (0.7,0.7,0.7,1)
            mfc1 = (0,0,0,1)
            ms1 = 3
            marker1 = 'o'
            
            ls_extra = ':'
            lw_extra = 1.0
            c_extra = (0.5,0.5,0.5,1)
            
            the_ax.semilogy(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            the_ax.semilogy(x_b1, y_b1, ls = ls_extra, lw = lw_extra, c = c_extra, label = 'b = 1')
            
            for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
                label.set_fontproperties(self.ticklabelfont)
                
            the_ax.set_ylabel('Log(Cumulative Number of Events)', fontproperties=self.framelabelfont)
            the_ax.set_xlabel('Magnitude', fontproperties=self.framelabelfont)
            
            the_ax.autoscale_view(tight=True)
            
            #the_ax.tick_params(fontproperties=arial12)
            
            #the_ax.margins(0.05,0.3)
            
            #the_ax.set_title(title, position=(0.0,1.05), ha='left', fontproperties=arial12)
            
            the_ax.legend(prop=self.legendfont)
               
            plt.savefig(self.out_file, format=self.output_format, dpi=res)

            #self.basicPlot(x_label='Magnitude', y_label='Log(Cumulative Number of Events)', title=title, out_file = self.out_file, axis_scale = 'semilogy', connect_points = True, extra_x = extra_x, extra_y = extra_y, extra_y_error = extra_y_error, extra_title = 'UCERF2 Observed 95% Confidence')
            '''
        
        

        

class EqSimEventExporter(object):
    def __init__(self,sys):
        print '    **** Initilizing EqSim Event Export'
        self.eqsim_events = quakelib.EQSimEventWriter()
        self.eqsim_event_id = int(1)
        self.sys = sys
        self.out_file = None
        self.min_slip_map_mag = 5.0
        self.eid_remap = None
    
    def setEidRemap(self, eid_remap_file):
        self.eid_remap = {}
    
        for ln, line in enumerate(eid_remap_file):
            if ln != 0:
                orig_id, vc_id = map(int, line.split())
                self.eid_remap[vc_id] = orig_id
    
    def addEvent(self, event_data, sweeps):
        event_sections = {}
            
        event_elements = {}
        #print event_data[0], event_data[1]
        for sweep in sweeps:
            # ele_info is:
            # {ele_id, total_slip, shear_init, normal_init, shear_final, normal_final, number_of_slips}
            #total_slip += VCEvents.sweeps[sweep_id].slip
            eid = int(sweep[2])
            #print '   ',eid
            try:
                event_elements[eid]['total_slip'] += float(sweep[3])
                event_elements[eid]['shear_final'] = float(sweep[8])
                event_elements[eid]['normal_final'] = float(sweep[9])
                event_elements[eid]['area'] = float(sweep[4])
                event_elements[eid]['number_of_slips'] += 1
            except KeyError:
                event_elements[eid] = {  'ele_id'                   :eid,
                                            'total_slip'            :float(sweep[3]),
                                            'shear_init'            :float(sweep[6]),
                                            'normal_init'           :float(sweep[7]),
                                            'shear_final'           :float(sweep[8]),
                                            'normal_final'          :float(sweep[9]),
                                            'number_of_slips'       :1,
                                            'sweep_of_first_slip'   :int(sweep[1]),
                                            'area'                  :float(sweep[4])
                                            }

        for ele_info in [e for e in sorted(event_elements.values(),lambda a,b: cmp(a['sweep_of_first_slip'],b['sweep_of_first_slip']))]:
            #print ele_info['ele_id']
            sid = self.sys.geometry.elements[ele_info['ele_id']].sid
            try:
                event_sections[sid].append(ele_info)
            except:
                event_sections[sid] = [ele_info]
        
        #print eid, event.trigger_id
        for sid, event_section in event_sections.iteritems():
            #print sid
            eqsim_event = quakelib.EQSimEventSummary()
            
            eqsim_event.set_event_id(self.eqsim_event_id)
            eqsim_event.set_magnitude(event_data[3])
            eqsim_event.set_time(self.sys.geometry.converter.year_sec(event_data[1]))
            eqsim_event.set_duration(1.0)
            eqsim_event.set_sid(int(sid))
            eqsim_event.set_shear_before(event_data[4])
            eqsim_event.set_shear_after(event_data[6])
            eqsim_event.set_normal_before(event_data[5])
            eqsim_event.set_normal_after(event_data[7])
            
            depth_lo    = 0
            depth_hi    = -py_sys.float_info.max
            das_lo      = py_sys.float_info.max
            das_hi      = 0
            area        = 0
            all_slips  = []
            
            for ele_num, ele_info in enumerate(event_section):
                curr_element = self.sys.geometry.elements[ele_info['ele_id']]
                node0 = self.sys.geometry.nodes[curr_element.enode_ids[0]]
                node1 = self.sys.geometry.nodes[curr_element.enode_ids[1]]
                node2 = self.sys.geometry.nodes[curr_element.enode_ids[2]]
                node3 = self.sys.geometry.nodes[curr_element.enode_ids[3]]
                
                if ele_num == 0:
                    lame_mu = curr_element.lame_mu
                    eqsim_event.set_hypo_depth( (node0.depth + node1.depth)/2.0 )
                    eqsim_event.set_hypo_das( (max((node0.das, node3.das)) + min((node0.das, node3.das)))/2.0 )
                
                ele_area = ele_info['area']
                
                area += ele_area
                all_slips.append(ele_info['total_slip'])
                
                if max((node0.das, node3.das)) > das_hi:
                    das_hi = max((node0.das, node3.das))
                
                if min((node0.das, node3.das)) < das_lo:
                     das_lo = min((node0.das, node3.das))
                
                if node0.depth > depth_hi:
                    depth_hi = node0.depth
                    
                if node1.depth < depth_lo:
                    depth_lo = node1.depth
                
                #print curr_element.lame_mu
                
                if event_data[3] >= self.min_slip_map_mag:
                    slip_map = quakelib.EQSimEventSlipMap()
                    slip_map.set_depth_lo(node1.depth)
                    slip_map.set_depth_hi(node0.depth)
                    slip_map.set_das_lo(min((node0.das, node3.das)))
                    slip_map.set_das_hi(max((node0.das, node3.das)))
                    slip_map.set_area(ele_area)
                    slip_map.set_mean_slip(ele_info['total_slip'])
                    slip_map.set_moment(ele_info['total_slip']*ele_area*curr_element.lame_mu)
                    slip_map.set_shear_before(ele_info['shear_init']) 
                    slip_map.set_shear_after(ele_info['shear_final'])
                    slip_map.set_normal_before(ele_info['normal_init']) 
                    slip_map.set_normal_after(ele_info['normal_final'])
                    if self.eid_remap is None:
                        slip_map.add_slip_entry(quakelib.EQSimEventSlipElement(int(curr_element.eid)))
                    else:
                        slip_map.add_slip_entry(quakelib.EQSimEventSlipElement(self.eid_remap[int(curr_element.eid)]))
                    
                    eqsim_event.add_slip_map(slip_map)
    
                    #print self.eid_remap[int(curr_element.eid)], int(curr_element.eid)
                    
                
                
            
            eqsim_event.set_depth_lo(depth_lo)
            eqsim_event.set_depth_hi(depth_hi)
            eqsim_event.set_das_lo(das_lo)
            eqsim_event.set_das_hi(das_hi)
            eqsim_event.set_area(area)
            eqsim_event.set_mean_slip(np.mean(all_slips))
            eqsim_event.set_moment(np.mean(all_slips)*area*lame_mu)
            
            self.eqsim_events.add_event_summary(eqsim_event)
        
        self.eqsim_event_id += 1
    
    def write(self):
        self.eqsim_events.open(self.out_file)
        self.eqsim_events.write()
        self.eqsim_events.flush()
        self.eqsim_events.close()
        
        #fix the first line
        os.rename(self.out_file,'tmp.dat')
        in_file = open('tmp.dat','r')
        out_file = open(self.out_file,'w')
        
        lines = in_file.readlines()
        lines[0] = '101 EQSim_Output_Event_2 2\n'
        out_file.writelines(lines)
        in_file.close()
        out_file.close()
        os.remove('tmp.dat')
        
class SectionLayerer(multiprocessing.Process):
    def __init__(self, vc_sys, work_queue, result_queue):
        
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.vc_sys = vc_sys
    
    def run(self):
        while not self.kill_received:
 
            # get a task
            try:
                section_list = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            sections = self.vc_sys.geometry.sectionsSortedByID(section_list)
            
            # the actual processing
            layered_sections = {}
            print '    layering %i sections '%(len(sections))
            for n, s in enumerate(sections):
                sid = s.sid
                #print sid
                curr_section = []
                for eid in s.selement_ids:
                    if self.vc_sys.geometry.elements[eid].onTrace():
                        #print s.elementsAlongDip(eid)
                        curr_section.append( s.elementsAlongDip(eid) )
                
                curr_section_np = np.array(curr_section)
                
                layered_sections[s.sid] = curr_section_np.T

            # store the result
            self.result_queue.put(layered_sections)
'''            
class DisplacementGridProcessor(multiprocessing.Process):
    def __init__(self, vc_sys, work_queue, result_queue, Xs, Ys):
        
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.vc_sys = vc_sys
        self.Xs = Xs
        self.Ys = Ys
    
    def run(self):
        while not self.kill_received:
 
            # get a task
            try:
                elements = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            sections = self.vc_sys.geometry.sectionsSortedByID(section_list)
            
            # the actual processing
            layered_sections = {}
            print '    layering %i sections '%(len(sections))
            for n, s in enumerate(sections):
                sid = s.sid
                #print sid
                curr_section = []
                for eid in s.selement_ids:
                    if self.vc_sys.geometry.elements[eid].onTrace():
                        #print s.elementsAlongDip(eid)
                        curr_section.append( s.elementsAlongDip(eid) )
                
                curr_section_np = np.array(curr_section)
                
                layered_sections[s.sid] = curr_section_np.T

            # store the result
            self.result_queue.put(layered_sections)
'''
class DisplacementGridProcessor(multiprocessing.Process):
    def __init__(self, vc_sys, work_queue, result_queue, Xs, Ys, counter, total_tasks):#, min_lat, min_lon, max_lat, max_lon):
        
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.vc_sys = vc_sys
        self.Xs = Xs
        self.Ys = Ys
        self.counter = counter
        self.total_tasks = total_tasks
        #self.m = Basemap(   llcrnrlon=min_lon, llcrnrlat=min_lat, urcrnrlon=max_lon, urcrnrlat=max_lat,\
        #                lat_0=(max_lat+min_lat)/2.0, lon_0=(max_lon+min_lon)/2.0,\
        #                resolution='i', projection='lcc', suppress_ticks=True
        #                )
    
    def run(self):
        while not self.kill_received:
 
            # get a task
            try:
                element_list = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            #sections = self.vc_sys.geometry.sectionsSortedByID(section_list)
            
            # the actual processing
            processed_displacements = {}
            #print '        processing displacements for %i elements '%(len(element_list))
            
            dX = np.empty(self.Xs.shape)
            dY = np.empty(self.Xs.shape)
            dZ = np.empty(self.Xs.shape)
        
            okada = quakelib.Okada()
            
            for eindex, element_data in enumerate(element_list):
                eid = element_data['eid']
                slip = element_data['slip']
                
                element = self.vc_sys.geometry.elements[eid]
                c = self.vc_sys.geometry.nodes[element.enode_ids[1]].depth
                US = slip * math.cos( element.rake )
                UD = slip * math.sin( element.rake )
                UT = 0
                dip = element.dip()
                vector0 = self.vc_sys.geometry.nodes[element.enode_ids[0]].cartVector()
                vector3 = self.vc_sys.geometry.nodes[element.enode_ids[3]].cartVector()
                vector1 = self.vc_sys.geometry.nodes[element.enode_ids[1]].cartVector()
                L = (vector0 - vector3).length()
                W = (vector0 - vector1).length()
                lame_lambda = element.lame_lambda
                lame_mu = element.lame_mu
                
                strike = element.strike()
                
                lat0 = self.vc_sys.geometry.nodes[element.enode_ids[0]].lat
                lon0 = self.vc_sys.geometry.nodes[element.enode_ids[0]].lon
                xp0, yp0 = self.vc_sys.geometry.converter.latlon_xy(lat0,lon0)
                #xp0, yp0 = self.m(lon0, lat0)
                
                #lat3 = self.vc_sys.geometry.nodes[element.enode_ids[3]].lat
                #lon3 = self.vc_sys.geometry.nodes[element.enode_ids[3]].lon
                #xp3, yp3 = self.vc_sys.geometry.converter.latlon_xy(lat3,lon3)
                
                # rotation matrix
                #
                # | sin -cos |
                # | cos  sin |
                
                it = np.nditer(self.Xs, flags=['multi_index'])
                while not it.finished:
                    #x, y = self.geometry.converter.latlon_xy(lats[it.multi_index],lons[it.multi_index])
                    
                    x, y = self.Xs[it.multi_index], self.Ys[it.multi_index]
                    
                    xp = (x-xp0) * round(math.sin( strike ),3) - (y-yp0) * round(math.cos( strike ),3)
                    yp = (x-xp0) * round(math.cos( strike ),3) + (y-yp0) * round(math.sin( strike ),3)
                    
                    loc = quakelib.Vec3(xp,yp,0)
                    dxp, dyp, dzp = okada.calc_displacement_vector(loc, c, dip, L, W, US, UD, UT, lame_lambda, lame_mu)
                    
                    dx =  dxp * round(math.sin( strike ),3) + dyp * round(math.cos( strike ),3)
                    dy = -dxp * round(math.cos( strike ),3) + dyp * round(math.sin( strike ),3)
                    dz =  dzp
                    
                    if eindex == 0:
                        dX[it.multi_index] = dx
                        dY[it.multi_index] = dy
                        dZ[it.multi_index] = dz
                    else:
                        dX[it.multi_index] += dx
                        dY[it.multi_index] += dy
                        dZ[it.multi_index] += dz
                    
                    it.iternext()
                
                self.counter.value += 1
                py_sys.stdout.write('\r')
                py_sys.stdout.flush()
                py_sys.stdout.write('        %i of %i completed'%(self.counter.value, self.total_tasks))
                py_sys.stdout.flush()
            
            processed_displacements['dX'] = dX
            processed_displacements['dY'] = dY
            processed_displacements['dZ'] = dZ
                
            '''
            for n, s in enumerate(sections):
                curr_section = []
                for eid in s.selement_ids:
                    if self.vc_sys.geometry.elements[eid].onTrace():
                        curr_section.append( s.elementsAlongDip(eid) )
                
                curr_section_np = np.array(curr_section)
                
                layered_sections[s.sid] = curr_section_np.T
            '''
            # store the result
            self.result_queue.put(processed_displacements)

class multiThreadStatusUpdater(object):
    def __init__(self, num_to_process):
        self.num_to_process = num_to_process
        self.num_processed = 0
        
    def update(self):
        self.num_processed += 1
        #sys.stdout.write('\r')
        #sys.stdout.flush()
        py_sys.stdout.write('    %i of %i completed'%(self.num_processed, self.num_to_process))
        py_sys.stdout.flush()

class VCSys(object):
    def __init__(self, sys_name, data_file):
        self._geometry = None
        self._events = VCEvents(self)
        self._layered_sections = None
        #self.greens = None
        self.name = sys_name
        self.data_file = data_file
        self.root_dir = os.path.dirname(data_file)
        self.cache_dir = '%s/%s_cache'%(self.root_dir, self.name)
        self.image_dir = '%s/%s_images'%(self.root_dir, self.name)
        self.dat_dir = '%s/%s_dat'%(self.root_dir, self.name)
        self.eqsim_dir = '%s/%s_eqsim'%(self.root_dir, self.name)
        
        self.output_format = 'pdf'
        
        #self.data_file = None
        #self.plot_faults = False
        #self.start_event = 0
        #self.end_event = None
        self.single_event_animation = False
        self.event_animation = False
        self.rupture_map = False
        self.section_filter = None
        self.rupture_map_min_mag = None
        #self.rupture_map_outfile = None
        self.frequency_magnitude = False
        self.magnitude_rupture_area = False
        self.magnitude_average_slip = False
        self.average_slip_rupture_length = False
        self.area_average_slip = False
        self.shear_stress_drop_magnitude = False
        self.export_eqsim_events = False
        self.export_eqsim_events_min_slip_map_mag = 5.0
        self.slip_in_slip_out = False
        self.recurrence_intervals = False
        #self.export_eqsim_geometry = False
        #self.trace_file = None
        self.export_eqsim_events = False
        self.print_event_info = False
        
        self.element_CFF                        = False
        self.element_CFF_map                    = False
        self.element_CFF_block_entropy          = None
        self.element_CFF_delta_block_entropy    = None
        self.element_CFF_E_Hmu                  = None
        
        # values for the fringes map are denoted by a {value}_f
        self.displacementMapConfig = {
            'font':               matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', weight='normal'),
            'font_bold':          matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', weight='bold'),
            'cmap':               plt.get_cmap('YlOrRd'),
            'cmap_f':             plt.get_cmap('jet'),
        #water
            'water_color':          '#4eacf4',
            'water_color_f':        '#4eacf4',
        #map boundaries
            'boundary_color':       '#000000',
            'boundary_color_f':     '#ffffff',
            'boundary_width':       1.0,
            'coastline_color':      '#000000',
            'coastline_color_f':    '#ffffff',
            'coastline_width':      1.0,
            'country_color':        '#000000',
            'country_color_f':      '#ffffff',
            'country_width':        1.0,
            'state_color':          '#000000',
            'state_color_f':        '#ffffff',
            'state_width':          1.0,
        #rivers
            'river_width':          0.25,
        #faults
            'fault_color':          '#000000',
            'fault_color_f':        '#ffffff',
            'event_fault_color':    '#ff0000',
            'event_fault_color_f':  '#ffffff',
            'fault_width':          0.5,
        #lat lon grid
            'grid_color':           '#000000',
            'grid_color_f':         '#ffffff',
            'grid_width':           0.0,
            'num_grid_lines':       5,
        #map props
            'map_resolution':       'i',
            'plot_resolution':      72.0,
            'map_tick_color':       '#000000',
            'map_tick_color_f':     '#000000',
            'map_frame_color':      '#000000',
            'map_frame_color_f':    '#000000',
            'map_frame_width':      1,
            'map_fontsize':         12,
            'arrow_inset':          10.0,
            'arrow_fontsize':       9.0,
            'cb_fontsize':          10.0,
            'cb_fontcolor':         '#000000',
            'cb_fontcolor_f':       '#000000',
            'cb_height':            20.0,
            'cb_margin_t':          10.0
        }
        #test
        self.checkCaches()
    
    @property
    def geometry(self):
        if self._geometry is None:
            self._geometry = cPickle.load(open('%s/geometry.pkl'%(self.cache_dir), 'rb'))
        return self._geometry
    
    @property
    def sections(self):
        if self._geometry is None:
            self._geometry = cPickle.load(open('%s/geometry.pkl'%(self.cache_dir), 'rb'))
        return self._geometry.sections

    @property
    def elements(self):
        if self._geometry is None:
            self._geometry = cPickle.load(open('%s/geometry.pkl'%(self.cache_dir), 'rb'))
        return self._geometry.elements
    
    @property
    def layered_sections(self):
        if self._layered_sections is None:
            self._layered_sections = cPickle.load(open('%s/layered-sections.pkl'%(self.cache_dir), 'rb'))
        return self._layered_sections

    @property
    def events(self):
        return self._events
    
    def webgui_checkCaches(self):
        print '*** Checking the webgui caches. ***'
        # checks for webgui specific caches and creates them if nessecary.
        
        webgui_cacheDir = '%s/webgui'%self.cache_dir
        
        # create the webgui cache dir if needed
        if not os.path.exists(webgui_cacheDir):
            os.makedirs(webgui_cacheDir)
    
        # create the full sys displacement grid
        print '    *** Displacement map grid. ***'
        self.loadAndSaveDisplacementGrid('full', True)
        
        # create the fault map data
        print '    *** Fault map data. ***'
        try:
            fault_map_data = cPickle.load(open('%s/fault-map-data.pkl'%(webgui_cacheDir),'rb'))
        except IOError:
            fault_map_data = self.webgui_getFaultTraces()
            cPickle.dump(fault_map_data, open('%s/fault-map-data.pkl'%(webgui_cacheDir), 'wb'))
        
        # create the system info
        print '    *** System info. ***'
        try:
            system_info = cPickle.load(open('%s/system-info.pkl'%(webgui_cacheDir),'rb'))
        except IOError:
            f = h5py.File(self.data_file, 'r')
            events = f['event_table']
            system_info =   {   'sections': len(self.geometry.sections),
                                'elements': len(self.geometry.elements),
                                'events':len(events),
                                'years':events[-1][1]
                            }
            cPickle.dump(system_info, open('%s/system-info.pkl'%(webgui_cacheDir), 'wb'))
            f.close()

        # create the event_ids, years and mags arrays.
        #print '    *** Event arrays. ***'
        #try:
        #    years, event_ids = cPickle.load(open('%s/years-event_ids.pkl'%(self.cache_dir), 'rb'))
        #    years, mags = cPickle.load(open('%s/years-mags.pkl'%(self.cache_dir), 'rb'))
        #except IOError:
        #    self.loadHdf5Events(self.data_file)
            #self.eventYearMappingHdf5(self.data_file)

        # create the sorted event lists.
        print '    *** Event lists. ***'
        dtypes =    [
                        ('evid', 'i4'),
                        ('evyear', 'f4'),
                        ('evmag', 'f4'),
                        ('evtriggersec', 'a30'),
                        ('evtriggerele', 'i4'),
                        ('evinvolvedsec', 'i4'),
                        ('evinvolvedele', 'i4'),
                        ('evaveslip', 'f4')
                    ]
        event_lists_cache_dir = '%s/event-lists'%webgui_cacheDir
        if not os.path.exists(event_lists_cache_dir):
            os.makedirs(event_lists_cache_dir)

        try:
            for col_id, col in enumerate(dtypes):
                f = h5py.File('%s/event-list_sort-%i-ascending.h5'%(event_lists_cache_dir, col_id), 'r')
                f.close()
                f = h5py.File('%s/event-list_sort-%i-descending.h5'%(event_lists_cache_dir, col_id), 'r')
                f.close()
        except IOError:
            print '        *** Reading events. ***'
            f = h5py.File(self.data_file, 'r')
                
            events = f['event_table']
            sweeps = f['event_sweep_table']

            #create the datatype for the fields
            vc_dtype = np.dtype(dtypes)

            number_of_events = len(events)

            event_list = np.empty((number_of_events,), dtype=vc_dtype)
            
            for ind, raw_event in enumerate(events):
                
                #event number
                event_list['evid'][ind] = raw_event[0]
                #event year
                event_list['evyear'][ind] = raw_event[1]
                #event magnitude
                event_list['evmag'][ind] = raw_event[3]
                #trigger section
                event_list['evtriggersec'][ind] = self.geometry.sections[self.geometry.elements[raw_event[2]].sid].sname
                #trigger element
                event_list['evtriggerele'][ind] = raw_event[2]
                #involved sections and elements and total slip
                involved_sections = {}
                involved_elements = {}
                total_slip = 0.0
                
                for sweep in sweeps[raw_event[8]:raw_event[9]]:
                    eid = sweep[2]
                    sid = self.geometry.elements[eid].sid
                    total_slip += sweep[3]
                    involved_sections[sid] = sid
                    involved_elements[eid] = eid
            
                event_list['evinvolvedsec'][ind] = int(len(involved_sections))
                event_list['evinvolvedele'][ind] = int(len(involved_elements))
                event_list['evaveslip'][ind] = total_slip/float(len(involved_elements))
                
                if ( ind%10000 == 0 ):
                    print '            event:%i of %i'%(ind, number_of_events)

            f.close()
            
            print '        *** Saving event cache files. ***'
            for col_id, col in enumerate(dtypes):
                data = np.sort(event_list, order=col[0])
                f = h5py.File('%s/event-list_sort-%i-ascending.h5'%(event_lists_cache_dir, col_id), 'w')
                f.create_dataset('event_table', data=data)
                f.close()
                
                f = h5py.File('%s/event-list_sort-%i-descending.h5'%(event_lists_cache_dir, col_id), 'w')
                f.create_dataset('event_table', data=data[::-1])
                f.close()

    def checkCaches(self):
        # a function to check if the caches exist.
        # if they do they are loaded here.
        # if they dont they are created.
        
        # create the cache dir if needed
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        
        # create the image dir if needed. this is where all plots will go.
        if not os.path.exists(self.image_dir):
            os.makedirs(self.image_dir)
    
        # create the dat dir if needed. this is where all dat format exports will go.
        if not os.path.exists(self.dat_dir):
            os.makedirs(self.dat_dir)
            
        # create the eqsim dir if needed. this is where all eqsim format exports will go.
        if not os.path.exists(self.eqsim_dir):
            os.makedirs(self.eqsim_dir)
        
        # create the geometry cache. This contains the model and the layered sections.
        if not os.path.isfile('%s/geometry.pkl'%(self.cache_dir)):
            self._geometry = VCGeometry(self)
            self.loadHdf5Geometry(self.data_file)
            cPickle.dump(self.geometry, open('%s/geometry.pkl'%(self.cache_dir), 'wb'))
        
        if not os.path.isfile('%s/layered-sections.pkl'%(self.cache_dir)):
            self.layerSections()
            cPickle.dump(self._layered_sections, open('%s/layered-sections.pkl'%(self.cache_dir), 'wb'))

        #if not os.path.isfile('%s/layered-sections.pkl'%(self.cache_dir)):


        '''
        try:
            self.geometry = cPickle.load(open('%s/geometry.pkl'%(self.cache_dir), 'rb'))
        except IOError:
            self.geometry = VCGeometry(self)
            self.loadHdf5Geometry(self.data_file)
            self.layerSections()
            cPickle.dump(self.geometry, open('%s/geometry.pkl'%(self.cache_dir), 'wb'))
        '''
        #self.events = VCEvents(self)
        
        # create the events cache. this contains the event year mapping.
        '''
        try:
            self.events = cPickle.load(open('%s/events.pkl'%(self.cache_dir), 'rb'))
        except IOError:
            self.events = VCEvents(self)
            self.loadHdf5Events(self.data_file)
            #self.eventYearMappingHdf5(self.data_file)
            cPickle.dump(self.events, open('%s/events.pkl'%(self.cache_dir), 'wb'))
        '''
        
    
    def fileNamePrepend(self, type='image', event_range=None, event=None, supress_sections=False):
        
        if type == 'image':
            directory = self.image_dir
        elif type == 'cache':
            directory = self.cache_dir
        elif type == 'dat':
            directory = self.dat_dir
        elif type == 'eqsim':
            directory = self.eqsim_dir
        else:
            directory = self.root_dir
        
        if self.section_filter is None or supress_sections:
            file_name_prepend = ''
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                """ Return list of consecutive lists of sids.
                    Cribbed from here: http://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy
                """
                step = 1
                run = []
                clustered_sids = [run]
                expect = None
                for sid in sorted(self.section_filter):
                    if (sid == expect) or (expect is None):
                        run.append(sid)
                    else:
                        run = [sid]
                        clustered_sids.append(run)
                    expect = sid + step
                
                for cluster in clustered_sids:
                    if len(cluster) > 1:
                        section_str += '%i--%i-'%(cluster[0],cluster[-1])
                    else:
                        section_str += '%i-'%(cluster[0])
                section_str = section_str.rpartition('-')[0]
            file_name_prepend = '%s_'%(section_str)
        
        if event_range is not None:
            file_name_prepend = '%s%i_%i_'%(file_name_prepend, event_range[0], event_range[1])
        elif event is not None:
            file_name_prepend = '%s%i_'%(file_name_prepend, event)
        return '%s/%s'%(directory, file_name_prepend)

    def webgui_getActionPlots(self, start_evid=0, end_evid=None, fm=False, mra=False, mas=False, asrl=False, st=False):
        event_file = self.data_file
        
        # returns a string with html for the view scaling plot page
        # also creates the images if ness
        
        #os.chdir(os.path.dirname(event_file))
        
        fmap_name = None
        ts_name = None
        fm_name = None
        mra_name = None
        mas_name = None
        asrl_name = None
        st_name = None
        
        fm_plot = False
        mra_plot = False
        mas_plot = False
        asrl_plot = False
        st_plot = False
        
        # if there is a section filter we need a fault map plot.
        if self.section_filter is not None:
            section_filter = {}
            for sid in self.section_filter:
                section_filter[sid] = self.geometry.sections[sid].sname
            fmap_name = '%sfault-map.png'%(self.fileNamePrepend())
            try:
                fmap_w, fmap_h = Image.open(fmap_name).size
            except IOError:
                fmap_w, fmap_h = self.plotFaultMap()
        
        f = h5py.File(event_file, 'r')
        events = f['event_table']
        if end_evid is None:
            end_evid = events[-1][0]
            end_year = events[-1][1]
            include_last_event = True
        else:
            end_year = events[end_evid][1]
            include_last_event = False
        
        start_year = events[start_evid][1]
        
        f.close()
        
        self.output_format = 'png'
        
        file_name_prepend_events = self.fileNamePrepend(event_range=[start_evid, end_evid])

        #### Frequency Magnitude
        if fm:
            fm_name = '%sfrequency-magnitude.png'%(file_name_prepend_events)
            try:
                fm_w, fm_h = Image.open(fm_name).size
            except IOError:
                fm_plot = True
                self.frequency_magnitude = 'plot'
            
        #### Magnitude Rupture Area
        if mra:
            mra_name = '%smagnitude-rupture-area.png'%(file_name_prepend_events)
            try:
                mra_w, mra_h = Image.open(mra_name).size
            except IOError:
                mra_plot = True
                self.magnitude_rupture_area = 'plot'

        #### Magnitude Average Slip
        if mas:
            mas_name = '%smagnitude-average-slip.png'%(file_name_prepend_events)
            try:
                mas_w, mas_h = Image.open(mas_name).size
            except IOError:
                mas_plot = True
                self.magnitude_average_slip = 'plot'

        #### Average Slip Rupture Length
        if asrl:
            asrl_name = '%saverage-slip-rupture-length.png'%(file_name_prepend_events)
            try:
                asrl_w, asrl_h = Image.open(asrl_name).size
            except IOError:
                asrl_plot = True
                self.average_slip_rupture_length = 'plot'
        
        #### Space Time Plot
        if st:
            st_name = '%sspace-time.png'%(file_name_prepend_events)
            try:
                st_w, st_h = Image.open(st_name).size
            except IOError:
                st_plot = True
                self.rupture_map = True
    
        try:
            event_count = self.eventActions(start_evid, end_evid)
        except ValueError:
            return {'error':'No events found on the selected sections in the selected event range'}

        if fm_plot:
            fm_w, fm_h = Image.open(fm_name).size
        if mra_plot:
            mra_w, mra_h = Image.open(mra_name).size
        if mas_plot:
            mas_w, mas_h = Image.open(mas_name).size
        if asrl_plot:
            asrl_w, asrl_h = Image.open(asrl_name).size
        if st_plot:
            st_w, st_h = Image.open(st_name).size

        #section_filter_tmp = self.section_filter
        #self.section_filter = None
        # create event time series if it doesnt exist
        if start_evid == 0 and include_last_event:
            # its all events, get the complete timeseries
            ts_name = '%sevent-time-series.png'%(self.fileNamePrepend(supress_sections=True))
            try:
                ts_w, ts_h = Image.open(ts_name).size
            except IOError:
                ts_w, ts_h = self.eventTimeSeries()
        else:
            # its a subset of events
            ts_name = '%sevent-time-series.png'%(self.fileNamePrepend(supress_sections=True, event_range=[start_evid,end_evid]))
            try:
                ts_w, ts_h = Image.open(ts_name).size
            except IOError:
                ts_w, ts_h = self.eventTimeSeries(start_event=start_evid, end_event=end_evid)

        ret = {}
        
        ret['start_year'] = start_year
        ret['end_year'] = end_year
        
        ret['start_evid'] = start_evid
        ret['end_evid'] = end_evid
        
        try:
            ret['section_filter'] = section_filter
        except UnboundLocalError:
            ret['section_filter'] = None
        
        try:
            ret['fmap_image'] = {'name':fmap_name.split('vc_data/')[-1], 'width':fmap_w, 'height':fmap_h}
        except (UnboundLocalError, AttributeError):
            ret['fmap_image'] = None
    
        try:
            ret['ts_image'] = {'name':ts_name.split('vc_data/')[-1], 'width':ts_w, 'height':ts_h}
        except (UnboundLocalError, AttributeError):
            ret['ts_image'] = None
    
        try:
            ret['fm_image'] = {'name':fm_name.split('vc_data/')[-1], 'width':fm_w, 'height':fm_h}
        except (UnboundLocalError, AttributeError):
            ret['fm_image'] = None
    
        try:
            ret['mas_image'] = {'name':mas_name.split('vc_data/')[-1], 'width':mas_w, 'height':mas_h}
        except (UnboundLocalError, AttributeError):
            ret['mas_image'] = None
    
        try:
            ret['mra_image'] = {'name':mra_name.split('vc_data/')[-1], 'width':mra_w, 'height':mra_h}
        except (UnboundLocalError, AttributeError):
            ret['mra_image'] = None
    
        try:
            ret['asrl_image'] = {'name':asrl_name.split('vc_data/')[-1], 'width':asrl_w, 'height':asrl_h}
        except (UnboundLocalError, AttributeError):
            ret['asrl_image'] = None

        try:
            ret['st_image'] = {'name':st_name.split('vc_data/')[-1], 'width':st_w, 'height':st_h}
        except (UnboundLocalError, AttributeError):
            ret['st_image'] = None

        return ret
    
    def webgui_getEventDetail(self, evid):
        event_file = self.data_file
        # returns a string with html for the event detail page
        # also creates the image for the event detail page if ness
        
        #os.chdir(os.path.dirname(event_file))
        
        file_name_prepend = self.fileNamePrepend(event=evid, supress_sections=True)
        fmap_name = '%sevent-fault-map.png'%(file_name_prepend)
        rm_name = '%sevent-rupture-map.png'%(file_name_prepend)
        ts_name = '%sevent-time-series.png'%(file_name_prepend)
        
        # create fault map image if it doesnt exist
        try:
            fmap_w, fmap_h = Image.open(fmap_name).size
        except IOError:
            fmap_w, fmap_h = self.plotFaultMap(evid=evid, event_file=self.data_file)
        
        # create the event rupture map image if it doesnt exist
        # old way #if not os.path.isfile('%s_Event-Rupture-Map_%i.png'%(self.name, evid)):
        try:
            rm_w, rm_h = Image.open(rm_name).size
        except IOError:
            rm_w, rm_h = self.eventRuptureMap(evid, title=False)
            
        # create the event time series image if it doesnt exist
        try:
            ts_w, ts_h = Image.open(ts_name).size
        except IOError:
            ts_w, ts_h = self.eventTimeSeries(marked_event=evid)
        
        # find the info for the event
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        event = events[evid]
        event_sweeps = sweeps[event[8]:event[9]]
        
        rupture_area = 0.0
        total_slip = 0.0
        rupture_length = 0.0
        
        involved_sections = {}
        involved_elements = {}
        
        for sweep in event_sweeps:
            eid = sweep[2]
            sid = self.geometry.elements[eid].sid
            total_slip += sweep[3]
            
            try:
                tmp = involved_elements[eid]
            except KeyError:
                side_length = self.geometry.elements[eid].sideLength()
                rupture_area += side_length ** 2.0
                if self.geometry.elements[eid].onTrace():
                    rupture_length += side_length
                involved_elements[eid] = True
                
            involved_sections[sid] = self.geometry.sections[sid].sname
        
        trigger_eid = event[2]
        trigger_lat, trigger_lon = self.geometry.elements[trigger_eid].pointLatLon()
        trigger_section = self.geometry.elements[trigger_eid].sid
        trigger_section_name = self.geometry.sections[trigger_section].sname
        
        f.close()
        
        year, month, day, hour, minute, second, decimal_second = self.geometry.converter.yearDecimalToYearMonthDay(event[1], time=True)
        
        average_slip = total_slip/float(len(involved_elements))
        
        ret = {}
        
        ret['number'] = evid
        ret['year'] = year
        ret['month'] = month
        ret['day'] = day
        ret['hour'] = hour
        ret['minute'] = minute
        ret['second'] = second
        ret['year_decimal'] = event[1]
        ret['magnitude'] = event[3]
        ret['rupture_area'] = rupture_area * 1.0e-6
        ret['rupture_length'] = rupture_length * 1.0e-3
        ret['average_slip'] = average_slip
        ret['trigger_eid'] = trigger_eid
        ret['trigger_lat'] = trigger_lat
        ret['trigger_lon'] = trigger_lon
        ret['trigger_section_sid'] = trigger_section
        ret['trigger_section_name'] = trigger_section_name
        ret['involved_sections'] = involved_sections
        ret['ts_image'] = {'name':ts_name.split('vc_data/')[-1], 'width':ts_w, 'height':ts_h }
        ret['fmap_image'] = {'name':fmap_name.split('vc_data/')[-1], 'width':fmap_w, 'height':fmap_h }
        ret['rm_image'] = {'name':rm_name.split('vc_data/')[-1], 'width':rm_w, 'height':rm_h }

        return ret
        
        '''
        
        involved_sections_string = ''
        
        for n, sid in enumerate(sorted(involved_sections.keys())):
            involved_sections_string += '%i %s'%(sid, self.geometry.sections[sid].sname)
            if n < len(involved_sections) - 1:
                involved_sections_string += ', '
        
        ret = []
        
        ret.append('<div id="event_info">\n')
        ret.append('    <div id="event_number">Number: <b>%i</b></div>\n'%(evid))
        ret.append('    <div id="event_year">Time: <b>%i-%i-%i %i:%i:%i (%.5f)</b></div>\n'%(year, month, day, hour, minute, second, event[1]))
        ret.append('    <div id="event_magnitude">Magnitude: <b>%.5f</b></div>\n'%(event[3]))
        ret.append('    <div id="event_rupture_area">Rupture Area: <b>%.5f [km<sup>2</sup>]</b></div>\n'%(rupture_area * 1.0e-6))
        ret.append('    <div id="event_surface_ruptue_length">Surface Rupture Length: <b>%.5f [km]</b></div>\n'%(rupture_length * 1.0e-3))
        ret.append('    <div id="event_average_slip">Average Slip: <b>%.5f [m]</b></div>\n'%(average_slip))
        ret.append('    <div id="event_trigger_location">Trigger Element: <b>%i</b> Lat: <b>%.5f</b> Lon: <b>%.5f</b></div>\n'%(trigger_eid, trigger_lat, trigger_lon))
        ret.append('    <div id="event_trigger_section">Trigger Section: <b>%i %s</b></div>\n'%(trigger_section, trigger_section_name))
        ret.append('    <div id="event_involved_sections">Involved Sections: <b>%s</b></div>\n'%(involved_sections_string))
        ret.append('</div>\n')

        ret.append('<div id="event_time_series">\n')
        ret.append('    <img src="images/%s_Event-Time-Series_%i.png" width="%i" height="%i" />\n'%(self.name, evid, ts_w, ts_h))
        ret.append('</div>\n')

        ret.append('<div id="event_fault_map">\n')
        ret.append('    <img src="images/%s_Fault-Map_%i.png" width="%i" height="%i" />\n'%(self.name, evid, fmap_w, fmap_h))
        ret.append('</div>\n')

        ret.append('<div id="event_rupture_map">\n')
        ret.append('    <img src="images/%s_Event-Rupture-Map_%i.png" width="%i" height="%i" />\n'%(self.name, evid, rm_w, rm_h))
        ret.append('</div>\n')

        return ''.join(ret)
        '''
    
    
    def webgui_getFaultTraces(self):
        ret_obj = {}
        for sid, sec in self.geometry.sections.iteritems():
            trace = []
            lats, lons = sec.traceLocs()
            for lat_index, lat in enumerate(lats):
                trace.append([lat, lons[lat_index]])
            ret_obj['%i'%sid] = {'name':'%s'%sec.sname, 'trace':trace}
        
        min_lat, max_lat, min_lon, max_lon =  self.geometry.minMaxLatLon()
        #ret_string += '%f, %f'%(min_lat + (max_lat - min_lat)/2.0, min_lon + (max_lon - min_lon)/2.0)
        
        ret_obj['system_center'] = [ min_lat + (max_lat - min_lat)/2.0, min_lon + (max_lon - min_lon)/2.0 ]
        '''
        # returns a string with all of the fault traces for display in the webgui
        ret_string = '{'
        sec_index = 0
        for sid, sec in self.geometry.sections.iteritems():
            ret_string += '%i:{name:"%s", trace:['%(sid,sec.sname)
            lats, lons = sec.traceLocs()
            for lat_index, lat in enumerate(lats):
                ret_string += '[%f, %f]'%(lat, lons[lat_index])
                if lat_index < len(lats) - 1:
                    ret_string +=', '
            ret_string += ']}'
            if sec_index < len(self.geometry.sections) - 1:
                ret_string += ', '
            sec_index += 1
        
        ret_string += '}'
        '''
        return json.dumps(ret_obj)
    
    def webgui_getEventList(self, page_number, number_per_page, sort_field=0, sort_direction='ascending'):
        print '*** Creating event list. Sort field:%i, Number per page:%i, Page number:%i, Sort direction:%s ***'%(sort_field, number_per_page, page_number, sort_direction)
        # returns a json string
        # sort_field is:
        # 0 : event number
        # 1 : event year
        # 2 : magnitude
        # 3 : trigger section
        # 4 : trigger element
        # 5 : involved sections
        # 6 : involved elements
        # 7 : average slip
        # page number should be 0 based
        
        cache_dir = '%s/webgui/event-lists/'%self.cache_dir
        
        # open the presorted list.
        f = h5py.File('%sevent-list_sort-%i-%s.h5'%(cache_dir, sort_field, sort_direction), 'r')
        events = f['event_table']
        
        ret_obj = {'Result':'OK', 'Records':[], 'TotalRecordCount':len(events)}
        
        for n, item in enumerate(events[page_number*(number_per_page):page_number*(number_per_page) + number_per_page]):
            #year, month, day, hour, minute, second, decimal_second = self.geometry.converter.yearDecimalToYearMonthDay(item[1], time=True)
            #print year, month, day, hour, minute, second, decimal_second
            ret_obj['Records'].append( {
                'evid':             int(item[0]),
                'evyear':           float(item[1]),
                'evmag':            float(item[2]),
                'evtriggersec':     str(item[3]),
                'evtriggerele':     int(item[4]),
                'evinvolvedsec':    int(item[5]),
                'evinvolvedele':    int(item[6]),
                'evaveslip':        float(item[7])
                
                } )
        f.close()
        print '*** Done ***'
        return json.dumps(ret_obj)
        
        
    '''
    def webgui_getSystemCenter(self):
        # returns a string with the lat lon of the system center for use in the webgui
        ret_string = '['
        
        min_lat, max_lat, min_lon, max_lon =  self.geometry.minMaxLatLon()
        
        ret_string += '%f, %f'%(min_lat + (max_lat - min_lat)/2.0, min_lon + (max_lon - min_lon)/2.0)
        
        ret_string += ']'
        
        return ret_string
    '''

    def magsBySection(self, start_event, end_event, num_bins = 100, landscape=False):
        event_file = self.data_file
        print '*** Plotting the number of events at specific magnitudes by section ***'
        print '    Event file %s'%event_file
        
        mags = self.events.mags
        
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        if end_event is not None:
            if end_event < start_event:
                print '!!! The end event must be bigger than the start event !!!'
                return
        else:
            end_event = events[-1][0]
        
        '''
        if self.section_filter is None:
            file_name_prepend = '%s_%i-%i_'%(self.name,start_event,end_event)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_%i-%i_'%(self.name,section_str,start_event,end_event)
        '''
        
        file_name_prepend = self.fileNamePrepend()
        
    # prepare the plot data
        #mags = []
        #for event in events[start_event:end_event]:
        #    mags.append(event[3])
        
        hist, mag_bins = np.histogram(mags, bins = num_bins)

        if self.section_filter is None:
            section_ids = np.array(sorted(self.geometry.sections.keys()))
        else:
            section_ids = np.array(sorted(self.section_filter))
        
        section_vs_mag = np.zeros( (mag_bins.size, len(section_ids)) )
        section_id_index_mapping = {}
        section_names = []

        for i, sid in enumerate(section_ids):
            section_id_index_mapping[sid] = i
            section_names.append(self.geometry.sections[sid].sname)
    
        print '    start event: %i, end event: %i'%(start_event, end_event)
        print '    *** Analyzing events ***'
        for event_num, event in enumerate(events[start_event:end_event]):
            if ( event_num%5000 == 0 ):
                py_sys.stdout.write('\r')
                py_sys.stdout.flush()
                py_sys.stdout.write('        event:%i of %i'%(event_num, end_event - start_event))
                py_sys.stdout.flush()
            
            involved_sections = {}
            
            event_sweeps = sweeps[event[8]:event[9]]
            
            for sweep in event_sweeps:
                eid = sweep[2]
                sid = self.geometry.elements[eid].sid
                involved_sections[sid] = sid

            event_mag = event[3]
            involved_section_ids = involved_sections.keys()
            
            mag_bin_num = np.digitize(np.array([event_mag]), mag_bins) - 1

            for sid in involved_section_ids:
                if self.section_filter is None or sid in self.section_filter:
                    section_vs_mag[mag_bin_num,section_id_index_mapping[sid]] += 1
        
        #print section_vs_mag
        #if not landscape:
        #    np.swapaxes(section_vs_mag,0,1)
        #print section_vs_mag
        print
        print '    *** Plotting ***'
        #for i in section_vs_mag.T:
        #    print np.array_str(i, max_line_width=600)
    # plot parameters
        if landscape:
            imw = 2642.0 # the full image width
        else:
            imh = 2642.0 # the full image height
        lm = 150.0
        rm = 20.0
        tm = 20.0
        bm = 100.0
        res = 300.0
        cbh = 25.0
        cbs = 150.0

        arial7 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=4)
        arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=6)

        cmap = cm.YlOrRd
    
    # do the plot
        if landscape:
            imh = imw * (float(len(mag_bins))/float(len(section_ids))) + cbh + cbs
        else:
            imw = imh * (float(len(mag_bins))/float(len(section_ids))) + cbh + cbs
        imwi = imw/res
        imhi = imh/res
        fig = plt.figure(figsize=(imwi, imhi), dpi=res)
        ph = imh - tm - bm - cbh - cbs # the height for the matrix
        pw = imw - lm - rm
        fig_ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))

        if section_vs_mag.max() - section_vs_mag.min() < 31:
            norm = mpl.colors.Normalize(vmin=section_vs_mag.min(), vmax=section_vs_mag.max())
        else:
            norm = mpl.colors.LogNorm(vmin=1, vmax=section_vs_mag.max())

        if landscape:
            fig_ax.pcolor(section_vs_mag, cmap=cmap, norm=norm )
        else:
            fig_ax.pcolor(section_vs_mag.T, cmap=cmap, norm=norm )

        fig_ax.axis('tight')

        #print len(section_names), section_vs_mag.T.shape
        if landscape:
            fig_ax.set_xticks(range(0, len(section_names)))
            for i in range(0, len(section_names)):
                fig_ax.axvline(x=float(i), lw=0.2, c='0.7')
            fig_ax.set_xticklabels(section_names)

            fig_ax.set_yticks(range(0, len(mag_bins)))
            mag_bin_labels = []
            for i in range(0, len(mag_bins)):
                fig_ax.axhline(y=float(i)-0.1, lw=0.2, c='0.7')
                mag_bin_labels.append('%f'%(mag_bins[i]))
            fig_ax.set_yticklabels(mag_bin_labels)
        else:
            fig_ax.set_yticks(range(0, len(section_names)))
            for i in range(0, len(section_names)):
                fig_ax.axhline(y=float(i) - 0.05, lw=0.2, c='0.7')
            fig_ax.set_yticklabels(section_names)

            fig_ax.set_xticks(range(0, len(mag_bins)))
            mag_bin_labels = []
            for i in range(0, len(mag_bins)):
                fig_ax.axvline(x=float(i), lw=0.2, c='0.7')
                mag_bin_labels.append('%f'%(mag_bins[i]))
            fig_ax.set_xticklabels(mag_bin_labels)

        for line in fig_ax.xaxis.get_ticklines() + fig_ax.yaxis.get_ticklines():
            line.set_alpha(0)

        for label in fig_ax.xaxis.get_ticklabels():
            label.set_fontproperties(arial7)
            label.set_rotation(90)
            label.set_ha('left')
        for label in fig_ax.yaxis.get_ticklabels():
            label.set_fontproperties(arial7)
            label.set_va('bottom')

        if landscape:
            fig_ax.set_ylabel('Magnitude', fontproperties=arial10)
        else:
            fig_ax.set_xlabel('Magnitude', fontproperties=arial10)

        for spine in fig_ax.spines.itervalues():
            spine.set_lw(0.5)

    # plot the colorbar
        cb_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, cbh/imh))
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='horizontal')

        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(arial7)
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)

        cb_ax.set_xlabel('Number of events', fontproperties=arial10)

        cb.outline.set_lw(0.5)

        #for spine in cb_ax.spines.itervalues():
        #    spine.set_lw(0.5)

        fig.savefig('%s%i_%i_mag_by_section.png'%(file_name_prepend, start_event, end_event), format='png', dpi=res)

        f.close()
    
    def plotGreensInteractions(self, output_type):
        data_file = self.data_file
        print '*** Plotting greens interactions. Output type: %s ***'%output_type
    
    # make sure the greens file exists
        try:
            f = h5py.File('%s/%s_Greens.h5'%(os.path.dirname(data_file), self.name), 'r')
        except:
            print '!!! Cannot find the greens file. It should be in the same directory as the data dile and named [sys-name]_Greens.h5 !!!'

    # plot parameters
        imw = 1024.0 # the full image width
        lm = 40.0
        rm = 50.0
        tm = 50.0
        bm = 50.0
        res = 72.0
        cbh = 20.0
        cbs = 40.0
        vss = 50.0 # vertical section spacing
        
        arial14 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=14)
        arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
        arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
        
        greens_normal = f['greens_normal']
        greens_shear = f['greens_shear']

        greens_normal_np = np.array(greens_normal)
        greens_shear_np = np.array(greens_shear)

    # create colormaps so that zero values map to white
        shear_white_pt = -greens_shear_np.min()/(greens_shear_np.max() - greens_shear_np.min())
        shear_color_dict = {'red':      (   (0.0, 0.0, 1.0),
                                            (shear_white_pt, 1.0, 1.0),
                                            (1.0, 0.0, 0.0)),
            
                            'green':    (   (0.0, 0.0, 0.0),
                                            (shear_white_pt, 1.0, 1.0),
                                            (1.0, 0.0, 0.0)),
            
                            'blue':     (   (0.0, 0.0, 0.0),
                                            (shear_white_pt, 1.0, 1.0),
                                            (1.0, 1.0, 0.0))
                            }

        normal_white_pt = -greens_normal_np.min()/(greens_normal_np.max() - greens_normal_np.min())
        normal_color_dict = {'red':      (   (0.0, 0.0, 1.0),
                                            (normal_white_pt, 1.0, 1.0),
                                            (1.0, 0.0, 0.0)),
            
                            'green':    (   (0.0, 0.0, 0.0),
                                            (normal_white_pt, 1.0, 1.0),
                                            (1.0, 0.0, 0.0)),
            
                            'blue':     (   (0.0, 0.0, 0.0),
                                            (normal_white_pt, 1.0, 1.0),
                                            (1.0, 1.0, 0.0))
                            }
            

        shear_cmap = LinearSegmentedColormap('shear_cmap', shear_color_dict, N=256, gamma=1.0)
        normal_cmap = LinearSegmentedColormap('normal_cmap', normal_color_dict, N=256, gamma=1.0)
    
    # do the plot
        imh = imw/2.0 + cbh + cbs
        imwi = imw/res
        imhi = imh/res
        fig = plt.figure(figsize=(imwi, imhi), dpi=res)
        ph = imh - tm - bm - cbh - cbs # the height for both matricies
        pw = ph
        shear_ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))
        normal_ax = fig.add_axes(((imw - pw - rm)/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))
        
        shear_ax.pcolor(greens_shear_np, cmap=shear_cmap )
        normal_ax.pcolor(greens_normal_np, cmap=normal_cmap  )

        shear_ax.invert_yaxis()
        normal_ax.invert_yaxis()

        shear_ax.axis('tight')
        normal_ax.axis('tight')

        for tick in shear_ax.xaxis.get_major_ticks() + normal_ax.xaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
        
        for label in shear_ax.xaxis.get_ticklabels() + shear_ax.yaxis.get_ticklabels() + normal_ax.xaxis.get_ticklabels() + normal_ax.yaxis.get_ticklabels():
            label.set_fontproperties(arial12)
        
        shear_ax.set_title('Shear interactions', fontproperties=arial14, color='k', va='bottom', ha='left', position=(0,1.05))
        normal_ax.set_title('Normal interactions', fontproperties=arial14, color='k', va='bottom', ha='left', position=(0,1.05))
    
    # create the color bars
        shear_cb_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, cbh/imh))
        normal_cb_ax = fig.add_axes(((imw - pw - rm)/imw, bm/imh, pw/imw, cbh/imh))

        norm = mpl.colors.Normalize(vmin=greens_shear_np.min(), vmax=greens_shear_np.max())
        cb = mpl.colorbar.ColorbarBase(shear_cb_ax, cmap=shear_cmap, norm=norm, orientation='horizontal')

        norm = mpl.colors.Normalize(vmin=greens_normal_np.min(), vmax=greens_normal_np.max())
        cb = mpl.colorbar.ColorbarBase(normal_cb_ax, cmap=normal_cmap, norm=norm, orientation='horizontal')

        for label in shear_cb_ax.xaxis.get_ticklabels() + normal_cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(arial10)
        for line in shear_cb_ax.xaxis.get_ticklines() + normal_cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        shear_cb_ax.set_title('shear stress per unit slip [Pa/m]', fontproperties=arial10, color='k', va='top', ha='left', position=(0,-1.1))
        normal_cb_ax.set_title('normal stress per unit slip [Pa/m]', fontproperties=arial10, color='k', va='top', ha='left', position=(0,-1.1))

    # mark what rows correspond to which sections
        for sid, section in self.geometry.sections.iteritems():
            min_id = min(section.selement_ids)
            max_id = max(section.selement_ids)
            
            shear_ax.add_line(mlines.Line2D((0,greens_shear_np.shape[0] + 1), (max_id + 1,max_id + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))
            shear_ax.add_line(mlines.Line2D((max_id + 1,max_id + 1), (0,greens_shear_np.shape[0] + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))
            
            shear_ax.text((max_id + min_id)/2, max_id, '%i %s'%(sid, section.sname), ha='center', va='bottom', fontproperties=arial7_light)

            normal_ax.add_line(mlines.Line2D((0,greens_shear_np.shape[0] + 1), (max_id + 1,max_id + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))
            normal_ax.add_line(mlines.Line2D((max_id + 1,max_id + 1), (0,greens_shear_np.shape[0] + 1), lw=0.5, ls=':', c='k', dashes=(2.0,1.0)))

            normal_ax.text((max_id + min_id)/2, max_id, '%i %s'%(sid, section.sname), ha='center', va='bottom', fontproperties=arial7_light)
        
        

    # output the results
        file_name_prepend = self.fileNamePrepend(supress_sections=True)
    
        if output_type == 'plot' or output_type == 'both':
            fig.savefig('%sgreens.png'%(file_name_prepend), format='png')

        if output_type == 'text' or output_type == 'both':
            tfn = open('%sgreens-normal.dat'%(file_name_prepend), 'w')
            tfs = open('%sgreens-shear.dat'%(file_name_prepend), 'w')

            for (x,y), val in np.ndenumerate(greens_normal_np):
                tfn.write('%i %i %f\n'%(x,y,val))
            for (x,y), val in np.ndenumerate(greens_shear_np):
                tfs.write('%i %i %f\n'%(x,y,val))
            
            tfn.close()
            tfs.close()

        f.close()
        print '*** Done ***'
    # done
    '''
    def slipInSlipOut(self, start_event, end_event):
        data_file = self.data_file
        print '*** Comparing input slip vs output slip ***'
        print '    Data file %s'%data_file
        
        f = h5py.File(data_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        if end_event is not None:
            if end_event < start_event:
                print '!!! The end event must be bigger than the start event !!!'
                return
        else:
            end_event = events[-1][0]
        
        if self.section_filter is None:
            file_name_prepend = '%s_%i-%i_'%(self.name,start_event,end_event)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_%i-%i_'%(self.name,section_str,start_event,end_event)

    # get the slip rate from the elements. this is the input slip.
        
        print '    *** Getting input slip rate for all elements. ***'
        slip_in = {}
        
        for eid, ele in self.geometry.elements.iteritems():
            if self.section_filter is None or ele.sid in self.section_filter:
                try:
                    slip = slip_in[ele.eid]
                    print '!!! Duplicate element id %i !!!'%eid
                    break
                except KeyError:
                    slip_in[ele.eid] = ele.slip


    # get the total slip from the elements. divide this by the length of the simulations and this is the output slip.

        print '    *** Getting output slip rate for all elements. ***'
        slip_out = {}
    
        if start_event == 0:
            duration = self.geometry.converter.year_sec(events[end_event][1])
        else:
            duration = self.geometry.converter.year_sec(events[end_event][1] - events[start_event][1])
    
        print '    start event: %i, end event: %i. %.3f years'%(start_event, end_event, self.geometry.converter.sec_year(duration))
        print '    **** Scanning events'
        for event_num, event_data in enumerate(events[start_event:end_event]):
            for sweep in sweeps[event_data[8]:event_data[9]]:
                eid = sweep[2]
                slip = sweep[3]
                if self.section_filter is None or self.geometry.elements[eid].sid in self.section_filter:
                    try:
                        slip_out[eid] += slip
                    except KeyError:
                        slip_out[eid] = slip

        # sort by the input slip
        sorted_eids = sorted(slip_in, key=lambda k: slip_in[k])
            
        
        for eid in sorted_eids:
            try:
                i = float(slip_in[eid])
                o = float(slip_out[eid])
            except KeyError:
                i = float(slip_in[eid])
                o = 0.0
            print eid, i, o/duration

        f.close()
        print '*** Done ***'
    '''
    def plotFaultMap(self, event_file=None, evid=None):
        print '*** Plotting fault map ***'
    
        # <Font 'Arial' (arial.ttf) normal normal 400 normal>
        ticklabelfont = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9)
    
    # set up plot
        pw = 501.0
        # ph set by the map aspect ratio
        lm = 40.0
        rm = 10.0
        tm = 10.0
        bm = 0.0
        res = 72.0
        
        water_color = '#4eacf4'
        cmap = plt.get_cmap('YlOrRd')
        land_color = cmap(0)
        
        boundary_color = 'k'
        boundary_width = '1'
        
        grid_color = '#5f84a0'
        grid_width = '0.1'
        num_grid_lines = 5
        
        coastline_color = 'k'
        coastline_width = '1'
        
        country_color = 'k'
        country_width = '1'
        
        state_color = 'k'
        state_width = '1'
        
        river_width = '0.25'
        
        fault_color = 'k'
        selected_fault_color = 'white'
        trigger_fault_color = 'white'
        fault_width = 0.5
        selected_fault_width = 6.0
        trigger_fault_width = 10.0
        
        map_resolution = 'i'
        
        map_tick_color = 'k'
        map_frame_color = 'k'
        map_frame_width = 1
        map_fontsize = 12
        
        '''
        if self.section_filter is None:
            file_name_prepend = '%s_'%(self.name)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_'%(self.name,section_str)
        '''
    
    # if the evid and the event_file are not none, figure out what sections are involved in the event
        if evid is not None and event_file is not None:
            f = h5py.File(event_file, 'r')
            
            events = f['event_table']
            sweeps = f['event_sweep_table']
            
            event = events[evid]
            
            event_sweeps = sweeps[event[8]:event[9]]
            
            involved_sections = {}
            
            for sweep in event_sweeps:
                eid = sweep[2]
                sid = self.geometry.elements[eid].sid
                involved_sections[sid] = sid
            
            involved_sections = involved_sections.keys()
        
            trigger_eid = event[2]
            trigger_section = self.geometry.elements[trigger_eid].sid
        
            f.close()
    
    # get the plot area and the trace vectors
        trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = self.loadAndSaveDisplacementGrid('full', True)

    # do plot
        m = Basemap(    llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                                lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                                resolution=map_resolution, projection='cyl'
                                )
        
        ph = pw*m.aspect
        
        pwi = pw/res
        phi = ph/res
        fig = plt.figure(figsize=(pwi, phi), dpi=res)
        m.ax = fig.add_axes((lm/pw, bm/ph, (pw-lm-rm)/pw, (ph-tm-bm)/ph))
        
        # draw the map boundary
        m.drawmapboundary(color=boundary_color, linewidth=1, fill_color=water_color)
        # fill the continents
        m.fillcontinents(color=land_color, lake_color=water_color)
        # draw coastlines, edge of map.
        m.drawcoastlines(color=coastline_color, linewidth=coastline_width)
        # draw countries
        m.drawcountries(linewidth=country_width, color=country_color)
        # draw states
        m.drawstates(linewidth=state_width, color=state_color)
        # draw parallels.
        parallels = np.linspace(lats.min(), lats.max(), num_grid_lines+1)
        m_parallels = m.drawparallels(parallels, labels=[1,0,0,0], fontsize=map_fontsize, color=grid_color, fontproperties=ticklabelfont, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
        # draw meridians
        meridians = np.linspace(lons.min(), lons.max(), num_grid_lines+1)
        m_meridians = m.drawmeridians(meridians, labels=[0,0,1,0], fontsize=map_fontsize, color=grid_color, fontproperties=ticklabelfont, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
        
        # print faults on lon-lat plot
        for sid, sec_trace in trace_vectors.iteritems():
            trace_lons = []
            trace_lats = []
            for trace_ele in sec_trace:
                trace_lons.append(trace_ele.lon)
                trace_lats.append(trace_ele.lat)
            trace_Xs, trace_Ys =m(trace_lons, trace_lats)
        
            #try:
            #    section_state = section_states[sid]
            #except KeyError:
            #    section_state = 0
            if self.section_filter is None:
                if evid is not None and event_file is not None:
                    if sid in involved_sections:
                        if sid == trigger_section:
                            the_color = trigger_fault_color
                            the_linewidth = trigger_fault_width
                        else:
                            the_color = selected_fault_color
                            the_linewidth = selected_fault_width
                        m.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=the_linewidth, solid_capstyle='round', solid_joinstyle='round')
                        m.plot(trace_Xs, trace_Ys, color=the_color, linewidth=the_linewidth - 3.0, solid_capstyle='round', solid_joinstyle='round')
                    else:
                        m.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width, solid_capstyle='round', solid_joinstyle='round')
                else:
                    m.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width, solid_capstyle='round', solid_joinstyle='round')
            else:
                if sid in self.section_filter:
                    m.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=selected_fault_width, solid_capstyle='round', solid_joinstyle='round')
                    m.plot(trace_Xs, trace_Ys, color=selected_fault_color, linewidth=selected_fault_width - 3.0, solid_capstyle='round', solid_joinstyle='round')
                else:
                    m.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width, solid_capstyle='round', solid_joinstyle='round')
                    
        
        #save the image
        if evid is not None and event_file is not None:
            file_name_prepend = self.fileNamePrepend(event=evid)
            fig.savefig('%sevent-fault-map.png'%(file_name_prepend), format='png', dpi=res)
        else:
            file_name_prepend = self.fileNamePrepend()
            fig.savefig('%sfault-map.png'%(file_name_prepend), format='png', dpi=res)
        
        print '*** Done ***'
        
        return pw, ph
        
    def eventTimeSeries(self, marked_event=None, start_event=0, end_event=None):
        event_file = self.data_file
        print '*** Plotting event time series ***'
        print '    Event file %s'%event_file
        
        if start_event == 0 and end_event is None:
            mark_event_range = False
        else:
            mark_event_range = True
        
        #print self.events.simple_event_list[0]
        
        #years = [e[1] for e in self.events.simple_event_list]
        #mags = [e[2] for e in self.events.simple_event_list]
        #print
        ##years = years_mags[0]
        #mags = years_mags[1]
        
        
        #try:
        #    years, mags = cPickle.load(open('%s/years-mags.pkl'%(self.cache_dir), 'rb'))
        #    print '    *** Data loaded from cache. ***'
        #except:
        #    self.loadHdf5Events(self.data_file)
        years = self.events.years
        mags = self.events.mags
        
    # set up the plot
        pw = 501.0
        ph = 150.0
        lm = 35.0
        rm = 20.0
        tm = 10.0
        bm = 35.0
        res = 72.0
        
        #low_mag_color = '#92C8FC'
        #high_mag_color = 'k'
        
        arial7 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9)
        #event_color_map = LinearSegmentedColormap.from_list('event_color_map', [low_mag_color,high_mag_color], N=256, gamma=1.0)
        
        event_color = '0.5'
        current_event_facecolor = 'white'
        current_event_edgecolor = 'k'
        event_range_color = 'k'
        
    # plot the time series
        pwi = pw/res
        phi = ph/res
        ts = plt.figure(figsize=(pwi, phi), dpi=res)
        ts_ax = ts.add_axes((lm/pw, bm/ph, (pw-lm-rm)/pw, (ph-tm-bm)/ph))
        
        
        
        #ts_ax.plot(years, mags, ls='None', marker=',', mfc='k')

        if mark_event_range:
            ts_ax.scatter(years[0:start_event-1], mags[0:start_event-1], marker=',', c=event_color, s=0.2, edgecolors='none')
            ts_ax.scatter(years[start_event:end_event], mags[start_event:end_event], marker=',', c=event_range_color, s=0.2, edgecolors='none')
            ts_ax.scatter(years[end_event+1:-1], mags[end_event+1:-1], marker=',', c=event_color, s=0.2, edgecolors='none')
        else:
            ts_ax.scatter(years, mags, marker=',', c=event_color, s=0.2, edgecolors='none')
        
        if marked_event is not None:
            ts_ax.scatter(years[marked_event], mags[marked_event], marker='*', c=current_event_facecolor, s=120.0, edgecolors=current_event_edgecolor, lw=1)
        
        ts_ax.set_xlim(0,max(years))
        ts_ax.set_ylim(min(mags),max(mags))
        
        mag_ticks = np.linspace(min(mags),max(mags),7)
        year_ticks = np.linspace(0,max(years),11)
        
        mag_ticklabels = []
        for i in mag_ticks:
            mag_ticklabels.append('%.1f'%i)
        year_ticklabels = []
        for i in year_ticks:
            year_ticklabels.append('%.1f'%i)
        
        ts_ax.set_xticks(year_ticks)
        ts_ax.set_yticks(mag_ticks)
        ts_ax.set_xticklabels(year_ticklabels)
        ts_ax.set_yticklabels(mag_ticklabels)
        
        
        
        for label in ts_ax.xaxis.get_ticklabels() + ts_ax.yaxis.get_ticklabels():
            label.set_fontproperties(arial7)
        
        for line in ts_ax.xaxis.get_ticklines() + ts_ax.yaxis.get_ticklines():
            line.set_alpha(0)
        
        ts_ax.set_xlabel('year', fontproperties=arial7)
        ts_ax.set_ylabel('magnitude', fontproperties=arial7)

        if marked_event is not None:
            file_name_prepend = self.fileNamePrepend(event=marked_event, supress_sections=True)
        elif mark_event_range:
            file_name_prepend = self.fileNamePrepend(event_range=[start_event,end_event], supress_sections=True)
        else:
            file_name_prepend = self.fileNamePrepend(supress_sections=True)
        
        
        
        ts.savefig('%sevent-time-series.png'%(file_name_prepend), format='png', dpi=res)
        
        print '*** Done ***'
        
        return pw, ph
    
    def slipVsTime(self, start_event, end_event):
        event_file = self.data_file
        print '*** Plotting slip vs time ***'
        print '    Event file %s'%event_file
        
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        if end_event is not None:
            if end_event < start_event:
                print '!!! The end event must be bigger than the start event !!!'
                return
        else:
            end_event = events[-1][0]
        
        if self.section_filter is None:
            file_name_prepend = '%s_%i-%i_'%(self.name,start_event,end_event)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_%i-%i_'%(self.name,section_str,start_event,end_event)
        
        element_ids = []

        for eid in sorted(self.geometry.elements.iterkeys()):
            ele = self.geometry.elements[eid]
            if self.section_filter is None or ele.sid in self.section_filter:
                element_ids.append(eid)

        
        #section_elements = sorted(self.geometry.sections[sid].selement_ids)
        
    # initilize the slips at zero and count the slip
        element_slip_history = {}
        event_years = [0.0,0.0]
        for eid in element_ids:
            element_slip_history[eid] = [0.0,0.0]
        
        if start_event == 0:
            duration = self.geometry.converter.year_sec(events[end_event][1])
        else:
            duration = self.geometry.converter.year_sec(events[end_event][1] - events[start_event][1])
        
        print '    start event: %i, end event: %i. %.3f years'%(start_event, end_event, self.geometry.converter.sec_year(duration))
        print '    **** Scanning events'
        
        for index, event_data in enumerate(events):
            curr_year = event_data[1]
            
            curr_slip_history_index = 2*(index + 1)
            curr_year = event_data[1]
            prev_year = event_years[curr_slip_history_index - 1]
            
            backslip_duration = self.geometry.converter.year_sec(curr_year - prev_year)
            
            slip_before_event = {}
            slip_after_event = {}
            
            
            # figure out the slip before the event. this will be the previous events slip minus backslip.
            for eid in element_ids:
                curr_slip = element_slip_history[eid][curr_slip_history_index - 1]
                slip_before_event[eid] = curr_slip - backslip_duration * self.geometry.elements[eid].slip
                slip_after_event[eid] = slip_before_event[eid]
            
            for sweep in sweeps[event_data[8]:event_data[9]]:
                eid = sweep[2]
                try:
                    slip_after_event[eid] += sweep[3]
                except KeyError:
                    pass

            for eid in element_ids:
                element_slip_history[eid].append(slip_before_event[eid])
                element_slip_history[eid].append(slip_after_event[eid])

            # append twice for before and after values
            event_years.append(curr_year)
            event_years.append(curr_year)
    
    # set up the plot
        pw = 1024
        ph = 768
        lm = 40.0
        rm = 20.0
        tm = 20.0
        bm = 20.0
        res = 72.0
        
        arial10 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=10)
        
        base_color = (0.5,0.5,0.5,0.5)
        base_width = 0.5
        
        hl_color = 'k'
        hl_width = 2
        
    # plot the time series
        pwi = pw/res
        phi = ph/res
        ts = plt.figure(figsize=(pwi, phi), dpi=res)
        ts_ax = ts.add_axes((lm/pw, bm/ph, (pw-lm-rm)/pw, (ph-tm-bm)/ph))
        
        for eid in element_ids:
            if eid in [37,149]:
            #if eid in [160,161,162,163,168,169,170,171,176,177,178,179,184,185,186,187,192,193,194,195]:
                ts_ax.plot(event_years[start_event:end_event], element_slip_history[eid][start_event:end_event], c=hl_color, lw=hl_width)
            else:
                ts_ax.plot(event_years[start_event:end_event], element_slip_history[eid][start_event:end_event], c=base_color, lw=base_width)
        
        for label in ts_ax.xaxis.get_ticklabels() + ts_ax.yaxis.get_ticklabels():
            label.set_fontproperties(arial10)
        
        ts_ax.axis('tight')

        ts.savefig('%sSlip-vs-Time.png'%(file_name_prepend), format='png')
        '''
        header_str = 'year '
        for eid in element_ids:
            header_str += 'e%i '%eid
        print header_str

        for index, year in enumerate(event_years):
            out_str = '%f '%year
            for eid in element_ids:
                out_str += '%f '%element_slip_history[eid][index]
            print out_str
        '''

    def sectionSlipMap(self, sid, start_event, end_event):
        event_file = self.data_file
        print '*** Creating slip map for section %i %s ***'%(sid, self.geometry.sections[sid].sname)
        print '    Event file %s'%event_file
        
    # figure out what elements are involved and find the total slip on each
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        if end_event is not None:
            if end_event < start_event:
                print '!!! The end event must be bigger than the start event !!!'
                return
        else:
            end_event = events[-1][0]
        
        involved_sections = [sid]
        
        layered_sections = self.geometry.getLayeredSections(involved_sections)
        
        total_vertical_elements = 0.0
        max_horizontal_elements = 0.0
        
        for k in sorted(layered_sections.keys()):
            total_vertical_elements += layered_sections[k].shape[0]
            if layered_sections[k].shape[1] > max_horizontal_elements:
                max_horizontal_elements = layered_sections[k].shape[1]
    
    # get the slip rate from the elements. this is the input slip.
        
        print '    *** Getting input slip rate for all elements. ***'
        slip_in = {}
        
        for eid in self.geometry.sections[sid].selement_ids:
            try:
                slip = slip_in[eid]
                print '!!! Duplicate element id %i !!!'%eid
                break
            except KeyError:
                slip_in[eid] = self.geometry.elements[eid].slip

    # get the total slip from the elements. divide this by the length of the simulations and this is the output slip.

        print '    *** Getting output slip rate for all elements. ***'
        slip_out = {}
    
        if start_event == 0:
            duration = self.geometry.converter.year_sec(events[end_event][1])
        else:
            duration = self.geometry.converter.year_sec(events[end_event][1] - events[start_event][1])
    
        print '    start event: %i, end event: %i. %.3f years'%(start_event, end_event, self.geometry.converter.sec_year(duration))
        print '    **** Scanning events'
        for event_num, event_data in enumerate(events[start_event:end_event]):
            for sweep in sweeps[event_data[8]:event_data[9]]:
                eid = sweep[2]
                if eid in self.geometry.sections[sid].selement_ids:
                    slip = sweep[3]
                    try:
                        slip_out[eid] += slip/duration
                    except KeyError:
                        slip_out[eid] = slip/duration

        #print slips

    # set up the plot
        imw = 1024.0 # the full image width
        vss = 50.0 # vertical section spacing
        cbh = 50.0 # colorbar height
        lm = 20.0
        rm = 20.0
        tm = 40.0
        bm = 40.0
        res = 72.0
        
        arial12 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=12)
        arial8_light = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=8, weight='light')
        
        not_involved_color = 'white'
        min_slip_color = 'yellow'
        max_slip_color = 'red'
        slip_color_map = LinearSegmentedColormap.from_list('slip_color_map', [min_slip_color,max_slip_color], N=256, gamma=1.0)
        
    # plot the rupture map
        pw = imw - lm - rm + 1.0 # the plot width for all axes
        ppe = (pw-1.0)/float(max_horizontal_elements) # pixels per element
        ph = 2.0*total_vertical_elements*ppe + vss + 1.0 # the height of the rupture map axis
        imh = ph + cbh + tm + bm # the full image height
        
        imwi = imw/res
        imhi = imh/res
        ssmr = plt.figure(figsize=(imwi, imhi), dpi=res)
        ssmr_ax = ssmr.add_axes((lm/imw, (bm+cbh+1.0)/imh, pw/imw, ph/imh))
        ssmr_ax.get_xaxis().set_visible(False)
        ssmr_ax.get_yaxis().set_visible(False)
        for spine in ssmr_ax.spines.itervalues():
            spine.set_alpha(0)
        
        vsl = 1.0 - ppe/ph # vertical start location
        min_slip = min((min(slip_out.values()), min(slip_in.values())))
        max_slip = max((max(slip_out.values()), max(slip_in.values())))
        if min_slip == max_slip:
            min_slip -= min_slip*0.1
            max_slip += max_slip*0.1
        
        do_trigger_loc = False
        
    # plot the output slip
        #print k, len(self.geometry.sections[k].selement_ids), layered_sections[k].shape[1]*layered_sections[k].shape[0]
        ssmr_patches = []
        
        element_index = 0
        for hi in range(layered_sections[sid].shape[1]):
            for vi in range(layered_sections[sid].shape[0]):
                eid = self.geometry.sections[sid].selement_ids[element_index]
                #print eid, event[2]
                try:
                    current_slip = slip_out[eid]
                    #print current_slip, min_slip, max_slip, (current_slip-min_slip)/(max_slip-min_slip)
                    element_color = slip_color_map((current_slip-min_slip)/(max_slip-min_slip))
                except KeyError:
                    #print element_index, self.geometry.sections[k].selement_ids[element_index], 'no slip'
                    element_color = not_involved_color
                ssmr_patches.append(mpatches.Rectangle(((hi*ppe)/pw, vsl - (vi*ppe)/ph), ppe/pw, ppe/ph, fc=element_color))
                #if eid == event[2]:
                #    trigger_loc = ((hi*ppe)/pw, vsl - (vi*ppe)/ph)
                #    do_trigger_loc = True
                element_index += 1
    
        #if do_trigger_loc:
            #ppe_frac = ppe*0.05
            #erm_patches.append(mpatches.Rectangle((trigger_loc[0]+ppe_frac/pw,trigger_loc[1]+ppe_frac/ph), (ppe-2.0*ppe_frac)/pw, (ppe-2.0*ppe_frac)/ph, fc='none', lw=ppe_frac))
            #erm_ax.text(trigger_loc[0]+0.5*ppe/pw, vsl - ((layered_sections[k].shape[0]-1.0)*ppe)/ph - 6.0/ph, 'Trigger', ha='left', va='top', fontproperties=arial8_light, transform=erm_ax.transAxes)
            #erm_ax.lines.extend([mlines.Line2D([trigger_loc[0]+0.5*ppe/pw, trigger_loc[0]+0.5*ppe/pw],[trigger_loc[1], vsl - ((layered_sections[k].shape[0]-1.0)*ppe)/ph - 4.0/ph],
            #                        linewidth=1.0, transform=erm_ax.transAxes, color='k', solid_capstyle='round', solid_joinstyle='round')])
            #do_trigger_loc = False
            #print vsl
    
        collection = mpl.collections.PatchCollection(ssmr_patches,match_original=True)
        ssmr_ax.add_collection(collection)
    
        ssmr_ax.text(0.0, vsl + ppe/ph + 10.0/ph, '%s output'%self.geometry.sections[sid].sname, ha='left', va='bottom', fontproperties=arial12, transform=ssmr_ax.transAxes)
    
        vsl -= (layered_sections[sid].shape[0]*ppe)/ph + vss/ph

    # plot the input slip
        #print k, len(self.geometry.sections[k].selement_ids), layered_sections[k].shape[1]*layered_sections[k].shape[0]
        ssmr_patches = []
        
        element_index = 0
        for hi in range(layered_sections[sid].shape[1]):
            for vi in range(layered_sections[sid].shape[0]):
                eid = self.geometry.sections[sid].selement_ids[element_index]
                #print eid, event[2]
                try:
                    current_slip = slip_in[eid]
                    #print current_slip, min_slip, max_slip, (current_slip-min_slip)/(max_slip-min_slip)
                    element_color = slip_color_map((current_slip-min_slip)/(max_slip-min_slip))
                except KeyError:
                    #print element_index, self.geometry.sections[k].selement_ids[element_index], 'no slip'
                    element_color = not_involved_color
                ssmr_patches.append(mpatches.Rectangle(((hi*ppe)/pw, vsl - (vi*ppe)/ph), ppe/pw, ppe/ph, fc=element_color))
                #if eid == event[2]:
                #    trigger_loc = ((hi*ppe)/pw, vsl - (vi*ppe)/ph)
                #    do_trigger_loc = True
                element_index += 1
    
        #if do_trigger_loc:
            #ppe_frac = ppe*0.05
            #erm_patches.append(mpatches.Rectangle((trigger_loc[0]+ppe_frac/pw,trigger_loc[1]+ppe_frac/ph), (ppe-2.0*ppe_frac)/pw, (ppe-2.0*ppe_frac)/ph, fc='none', lw=ppe_frac))
            #erm_ax.text(trigger_loc[0]+0.5*ppe/pw, vsl - ((layered_sections[k].shape[0]-1.0)*ppe)/ph - 6.0/ph, 'Trigger', ha='left', va='top', fontproperties=arial8_light, transform=erm_ax.transAxes)
            #erm_ax.lines.extend([mlines.Line2D([trigger_loc[0]+0.5*ppe/pw, trigger_loc[0]+0.5*ppe/pw],[trigger_loc[1], vsl - ((layered_sections[k].shape[0]-1.0)*ppe)/ph - 4.0/ph],
            #                        linewidth=1.0, transform=erm_ax.transAxes, color='k', solid_capstyle='round', solid_joinstyle='round')])
            #do_trigger_loc = False
            #print vsl
    
        collection = mpl.collections.PatchCollection(ssmr_patches,match_original=True)
        ssmr_ax.add_collection(collection)
    
        ssmr_ax.text(0.0, vsl + ppe/ph + 10.0/ph, '%s input'%self.geometry.sections[sid].sname, ha='left', va='bottom', fontproperties=arial12, transform=ssmr_ax.transAxes)
    
        vsl -= (layered_sections[sid].shape[0]*ppe)/ph + vss/ph

    
        #print round(min_slip,1), round(max_slip,1)
    #plot the cb
        print min_slip, max_slip
        cb_ax = ssmr.add_axes((lm/imw, bm/imh, pw/imw, (cbh-30.0)/imh))
        norm = mpl.colors.Normalize(vmin=min_slip, vmax=max_slip, clip=False)
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=slip_color_map, norm=norm, orientation='horizontal')
        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(arial12)
            label.set_color('k')
            #label.set_va('bottom')
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)
            
        cb_ax.set_title('Total output slip rate [m/s]', fontproperties=arial12, color='k', va='top', ha='left', position=(0,-1.2))

        ssmr.savefig('%s_%i_Section-Slip-Rate-Map.png'%(self.name, sid), format='png')
        
        f.close()
        
    def eventRuptureMap(self, evid, title=True):
        event_file = self.data_file
        print '*** Creating rupture map for event %i ***'%evid
        print '    Event file %s'%event_file
        
    # figure out what elements are involved and find the total slip on each
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        event = events[evid]
        
        event_sweeps = sweeps[event[8]:event[9]]
        
        involved_elements = {}
        involved_sections = {}
        
        for sweep in event_sweeps:
            eid = sweep[2]
            sid = self.geometry.elements[eid].sid
            try:
                involved_elements[eid] += sweep[3]
            except KeyError:
                involved_elements[eid] = sweep[3]
            
            involved_sections[sid] = sid
        
        involved_sections = involved_sections.keys()
        slips = involved_elements.values()
        
        layered_sections = self.geometry.getLayeredSections(involved_sections)
        
        total_vertical_elements = 0.0
        max_horizontal_elements = 0.0
        
        for k in sorted(layered_sections.keys()):
            total_vertical_elements += layered_sections[k].shape[0]
            if layered_sections[k].shape[1] > max_horizontal_elements:
                max_horizontal_elements = layered_sections[k].shape[1]
        
    # set up the plot
        imw = 1010.0 # the full image width
        vss = 60.0 # vertical section spacing
        cbh = 50.0 # colorbar height
        lm = 20.0
        rm = 20.0
        if title:
            tm = 80.0
        else:
            tm = 40.0
        bm = 40.0
        res = 72.0
        
        arial14 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=14)
        arial12 = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=12)
        arial8_light = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9, weight='light')
        arial7_light = matplotlib.font_manager.FontProperties(family='FreeSans', style='normal', variant='normal', size=9, weight='light')
        
        not_involved_color = 'white'
        min_slip_color = 'yellow'
        max_slip_color = 'red'
        slip_color_map = LinearSegmentedColormap.from_list('slip_color_map', [min_slip_color,max_slip_color], N=256, gamma=1.0)
        
    # plot the rupture map
        pw = imw - lm - rm + 1.0 # the plot width for all axes
        ppe = (pw-1.0)/float(max_horizontal_elements) # pixels per element
        ph = total_vertical_elements*ppe + (len(involved_sections)-1)*vss + 1.0 # the height of the rupture map axis
        imh = ph + cbh + tm + bm # the full image height
        
        imwi = imw/res
        imhi = imh/res
        erm = plt.figure(figsize=(imwi, imhi), dpi=res)
        erm_ax = erm.add_axes((lm/imw, (bm+cbh+1.0)/imh, pw/imw, ph/imh))
        erm_ax.get_xaxis().set_visible(False)
        erm_ax.get_yaxis().set_visible(False)
        for spine in erm_ax.spines.itervalues():
            spine.set_alpha(0)
        
        if title:
            erm_ax.set_title('Event %i. Magnitude %.3f'%(evid, event[3]), fontproperties=arial14, color='k', va='top', ha='left', position=(0,(ph + 55)/ph))
        
        vsl = 1.0 - ppe/ph # vertical start location
        min_slip = min(slips)
        max_slip = max(slips)
        if round(min_slip,1) == round(max_slip,1):
            min_slip -= min_slip*0.1
            max_slip += max_slip*0.1
        
        do_trigger_loc = False
        
        for k in sorted(layered_sections.keys()):
            #print k, len(self.geometry.sections[k].selement_ids), layered_sections[k].shape[1]*layered_sections[k].shape[0]
            erm_patches = []
            
            trigger_ele = False
            
            element_index = 0
            for hi in range(layered_sections[k].shape[1]):
                for vi in range(layered_sections[k].shape[0]):
                    eid = self.geometry.sections[k].selement_ids[element_index]
                    
                    if eid == event[2]:
                        trigger_loc = ((hi*ppe)/pw, vsl - (vi*ppe)/ph)
                        do_trigger_loc = True
                        trigger_ele = True
                    
                    try:
                        current_slip = involved_elements[eid]
                        #print current_slip, min_slip, max_slip, (current_slip-min_slip)/(max_slip-min_slip)
                        element_color = slip_color_map((current_slip-min_slip)/(max_slip-min_slip))
                    except KeyError:
                        #print element_index, self.geometry.sections[k].selement_ids[element_index], 'no slip'
                        element_color = not_involved_color
                    erm_patches.append(mpatches.Rectangle(((hi*ppe)/pw, vsl - (vi*ppe)/ph), ppe/pw, ppe/ph, fc=element_color))
                    if trigger_ele:
                        erm_ax.text((hi*ppe + 8)/pw, vsl - (vi*ppe - 8)/ph, '%i'%eid, fontproperties=arial7_light)
                        trigger_ele = False
                    else:
                        erm_ax.text((hi*ppe + 2)/pw, vsl - (vi*ppe - 2)/ph, '%i'%eid, fontproperties=arial7_light)
                
                    #print self.geometry.elements[eid].onTrace(), element_index, layered_sections[k].size - layered_sections[k].shape[0],
                    if self.geometry.elements[eid].onTrace() and element_index == 0:
                        erm_ax.text((hi*ppe)/pw, vsl - (vi*ppe - ppe - 5)/ph, '(%.3f, %.3f)'%self.geometry.elements[eid].pointLatLon(), fontproperties=arial7_light)
                    if self.geometry.elements[eid].onTrace() and element_index == layered_sections[k].size - layered_sections[k].shape[0]:
                        erm_ax.text((hi*ppe + ppe)/pw, vsl - (vi*ppe - ppe - 5)/ph, '(%.3f, %.3f)'%self.geometry.elements[eid].pointLatLon(), fontproperties=arial7_light, ha='right')
                    
                    element_index += 1
        
            if do_trigger_loc:
                ppe_frac = ppe*0.05
                erm_patches.append(mpatches.Rectangle((trigger_loc[0]+ppe_frac/pw,trigger_loc[1]+ppe_frac/ph), (ppe-2.0*ppe_frac)/pw, (ppe-2.0*ppe_frac)/ph, fc='none', lw=ppe_frac))
                erm_ax.text(trigger_loc[0]+0.5*ppe/pw, vsl - ((layered_sections[k].shape[0]-1.0)*ppe)/ph - 6.0/ph, 'Trigger', ha='left', va='top', fontproperties=arial8_light, transform=erm_ax.transAxes)
                erm_ax.lines.extend([mlines.Line2D([trigger_loc[0]+0.5*ppe/pw, trigger_loc[0]+0.5*ppe/pw],[trigger_loc[1], vsl - ((layered_sections[k].shape[0]-1.0)*ppe)/ph - 4.0/ph],
                                        linewidth=1.0, transform=erm_ax.transAxes, color='k', solid_capstyle='round', solid_joinstyle='round')])
                do_trigger_loc = False
                #print vsl
        
            collection = mpl.collections.PatchCollection(erm_patches,match_original=True)
            erm_ax.add_collection(collection)
        
            erm_ax.text(0.0, vsl + ppe/ph + 20.0/ph, '%i %s'%(k, self.geometry.sections[k].sname), ha='left', va='bottom', fontproperties=arial12, transform=erm_ax.transAxes)
        
            vsl -= (layered_sections[k].shape[0]*ppe)/ph + vss/ph
    
        #print round(min_slip,1), round(max_slip,1)
    #plot the cb
        #print round(min_slip,2), round(max_slip,2)
        cb_ax = erm.add_axes((lm/imw, bm/imh, pw/imw, (cbh-30.0)/imh))
        norm = mpl.colors.Normalize(vmin=round(min_slip,1), vmax=round(max_slip,1), clip=False)
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=slip_color_map, norm=norm, orientation='horizontal')
        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(arial12)
            label.set_color('k')
            #label.set_va('bottom')
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)
            
        cb_ax.set_title('Total slip [m]', fontproperties=arial12, color='k', va='top', ha='left', position=(0,-1.2))
        
        file_name_prepend = self.fileNamePrepend(event=evid, supress_sections=True)
        
        erm.savefig('%sevent-rupture-map.png'%(file_name_prepend), format='png', dpi=res)
        
        f.close()

        print '*** Done ***'

        return imw, imh
        
    def calculateEventDisplacements(self, Xs, Ys, involved_elements, evid):#,min_lat, min_lon, max_lat, max_lon):
        # involved elements is a list [{'eid':eid,'slip':slip}, {'eid':eid,'slip':slip} ...]
        print '    *** Calculating event %i displacements. %i involved elements. %i sample points. ***'%(evid, len(involved_elements), Xs.shape[0]*Xs.shape[1])
        
        #global mp_status_updater
        #mp_status_updater = multiThreadStatusUpdater(len(involved_elements))
        
        counter = multiprocessing.Value('i',0)
        
        num_processes = multiprocessing.cpu_count()
            
        seg = int(round(float(len(involved_elements))/float(num_processes)))
        
        if seg < 1:
            seg = 1
        
        segmented_elements_init = []
            
        for i in range(num_processes):
            if i == num_processes - 1:
                end_index = len(involved_elements)
            else:
                end_index = seg*int(i + 1)
            segmented_elements_init.append((involved_elements[ int(i) * seg : end_index ]))
        
        segmented_elements = []
        for i, n in enumerate(segmented_elements_init):
            if len(n) != 0:
                segmented_elements.append(n)
        
        #print len(segmented_elements)
        
        work_queue = multiprocessing.Queue()
        for job in segmented_elements:
            #print len(job)
            work_queue.put(job)
        
        # create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()
         
        py_sys.stdout.write('        %i of %i completed'%(0, len(involved_elements)))
        py_sys.stdout.flush()
        # spawn workers
        for i in range(len(segmented_elements)):
            worker = DisplacementGridProcessor(self, work_queue, result_queue, Xs, Ys, counter, len(involved_elements))#, min_lat, min_lon, max_lat, max_lon)
            worker.start()
        
        # collect the results off the queue
        results = []
        for i in range(len(segmented_elements)):
            results.append(result_queue.get())
        
        print
        
        dX = None
        dY = None
        dZ = None
        
        for result_num, result in enumerate(results):
            if dX is None:
                dX = result['dX']
            else:
                dX += result['dX']
            
            #print result['dX']
                
            if dY is None:
                dY = result['dY']
            else:
                dY += result['dY']
                
            if dZ is None:
                dZ = result['dZ']
            else:
                dZ += result['dZ']
        
        return dX, dY, dZ
    
    def calculateDisplacementMapGrid(self, area_pad = 0.01, map_resolution = 'i', map_area = None):
        print '    *** Calculating displacement map grid ***'
        
        # these are constrained this way so we can plot on 1024x780 for the animations
        max_plot_width = 690.0
        max_plot_height = 658.0
        
        if map_area is None:
            min_lat, max_lat, min_lon, max_lon = self.geometry.minMaxLatLon()
        else:
            min_lat, max_lat, min_lon, max_lon = map_area[1].lat, map_area[0].lat, map_area[1].lon, map_area[0].lon
        
        #trace_lats = []
        #trace_lons = []
        
        trace_vectors = {}
        
        #for sid in involved_sections.itervalues():
        for sid, section in self.geometry.sections.iteritems():
            lats, lons =  self.geometry.sections[sid].traceLocs()
            sec_trace_vecs = []
            for i in range(len(lats)):
                x1, y1 = self.geometry.converter.latlon_xy(lats[i],lons[i])
                try:
                    x2, y2 = self.geometry.converter.latlon_xy(lats[i+1],lons[i+1])
                    vec1 = vec.MyVector(x1, y1, 0)
                    vec2 = vec.MyVector(x2, y2, 0)
                    l = (vec2 - vec1).length()
                    if l > 0.0:
                        #sec_trace_vecs.append({ 'v1':vec1,
                        #                        'v2':vec2,
                        #                        'v2-v1':vec2-vec1,
                        #                        'length':l,
                        #                        '#x2-x1':vec2.x-vec1.x,
                        #                        'y2-y1':vec2.y-vec1.y,
                        #                        'p1':vec.LatLonPoint(lats[i],lons[i]),
                        #                        'p2':vec.LatLonPoint(lats[i+1],lons[i+1])
                        #                       })
                        sec_trace_vecs.append(vec.LatLonPoint(lats[i],lons[i]))
                        sec_trace_vecs.append(vec.LatLonPoint(lats[i+1],lons[i+1]))
                except IndexError:
                    break
                    
            #trace_lats.append(lats)
            #trace_lons.append(lons)
            trace_vectors[sid] = sec_trace_vecs
        
        lon_range = max_lon - min_lon
        lat_range = max_lat - min_lat
        max_range = max((lon_range, lat_range))
        #print lon_range, lat_range
        padded_min_lon = min_lon - lon_range*area_pad
        padded_min_lat = min_lat - lat_range*area_pad
        padded_max_lon = max_lon + lon_range*area_pad
        padded_max_lat = max_lat + lat_range*area_pad
        
        #print min_lat, max_lat, min_lon, max_lon
        #print padded_min_lat, padded_max_lat, padded_min_lon, padded_max_lon
        #m1, fig1 is the oceans and the continents. This will lie behind the masked data image
        #    initilizing this here so we can use it for various measurements
        m1 = Basemap(   llcrnrlon=padded_min_lon, llcrnrlat=padded_min_lat, urcrnrlon=padded_max_lon, urcrnrlat=padded_max_lat,\
                        lat_0=(padded_max_lat+padded_min_lat)/2.0, lon_0=(padded_max_lon+padded_min_lon)/2.0,\
                        resolution=map_resolution, projection='cyl', suppress_ticks=True
                        )
        # aspect is height/width
        if m1.aspect > 1.0:
            plot_height = max_plot_height
            plot_width = max_plot_height/m1.aspect
        else:
            plot_width = max_plot_width
            plot_height = max_plot_width*m1.aspect
        
        #print plot_width, plot_height, m1.aspect
        
        lon_divisions = int(plot_width)
        lat_divisions = int(plot_height)
        
        lons_1d = np.linspace(padded_min_lon,padded_max_lon,lon_divisions)
        lats_1d = np.linspace(padded_min_lat,padded_max_lat,lat_divisions)
        
        lons,lats = np.meshgrid(lons_1d,lats_1d)
        
        #print lons.min(), lons.max(), lats.min(),lats.max()
        #print padded_min_lon, padded_max_lon, padded_min_lat, padded_max_lat
        
        #M = np.zeros(lons.shape, dtype='bool')
        
        Xs = np.empty(lons.shape)
        Ys = np.empty(lons.shape)
        
        #Xs, Ys = m1(lons, lats)
        
        #found_point = False
        it = np.nditer(lons, flags=['multi_index'])
        while not it.finished:
            Xs[it.multi_index], Ys[it.multi_index] = self.geometry.converter.latlon_xy(lats[it.multi_index],lons[it.multi_index])
            it.iternext()
            
        return trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys
    
    def loadAndSaveDisplacementGrid(self, plot_area, cache_results, involved_sections = None, map_resolution = 'i', evid = None):
        # load the lat lon/x y points and the fault traces
        save_displacement_grid = False
        try:
            if plot_area == 'full':
                trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = cPickle.load(open('%s/displacements_full-sys/displacement_map_grid.pkl'%(self.cache_dir), 'rb'))
            elif plot_area == 'event':
                trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = cPickle.load(open('%s/displacements_events/displacement_map_grid-event_%i.pkl'%(self.cache_dir, evid), 'rb'))
            else:
                 trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = cPickle.load(open('%s/displacements_%.3f-%.3f_%.3f-%.3f/displacement_map_grid.pkl'%(self.cache_dir,plot_area[0].lat,plot_area[0].lon,plot_area[1].lat,plot_area[1].lon), 'rb'))
            print '    *** Displacement map grid loaded from cache. ***'
        except:
            if plot_area == 'full':
                map_area = None
                area_pad = 0.05
            elif plot_area == 'event':
                #calc map area for the specific event
                minLat, maxLat, minLon, maxLon = self.geometry.areaOfSections(involved_sections)
                NE_pt = vec.LatLonPoint( maxLat,maxLon)
                SW_pt = vec.LatLonPoint( minLat,minLon)
                map_area = [NE_pt, SW_pt]
                area_pad = 0.1
            else:
                map_area = plot_area
                area_pad = 0.1
                
            trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = self.calculateDisplacementMapGrid(map_resolution=map_resolution, map_area=map_area, area_pad=area_pad)
            
            save_displacement_grid = True
        
        #print minLat, maxLat, minLon, maxLon
        #print lats.min(), lats.max(), lons.min(), lons.max()
        
        #save the displacement grid
        
        if save_displacement_grid and cache_results:
            if plot_area == 'full':
                # save the full sys displacement grid
                if not os.path.exists('%s/displacements_full-sys'%(self.cache_dir)):
                    os.makedirs('%s/displacements_full-sys'%(self.cache_dir))
                cPickle.dump((trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys), open('%s/displacements_full-sys/displacement_map_grid.pkl'%(self.cache_dir), 'wb'))
            elif plot_area == 'event':
                # save displacement grid for this specific event
                if not os.path.exists('%s/displacements_events'%(self.cache_dir)):
                    os.makedirs('%s/displacements_events'%(self.cache_dir))
                cPickle.dump((trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys), open('%s/displacements_events/displacement_map_grid-event_%i.pkl'%(self.cache_dir, evid), 'wb'))
            else:
                # save displacement grid for the specified area
                if not os.path.exists('%s/displacements_%.3f-%.3f_%.3f-%.3f'%(self.cache_dir,map_area[0].lat,map_area[0].lon,map_area[1].lat,map_area[1].lon)):
                    os.makedirs('%s/displacements_%.3f-%.3f_%.3f-%.3f'%(self.cache_dir,map_area[0].lat,map_area[0].lon,map_area[1].lat,map_area[1].lon))
                cPickle.dump((trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys), open('%s/displacements_%.3f-%.3f_%.3f-%.3f/displacement_map_grid.pkl'%(self.cache_dir,plot_area[0].lat,plot_area[0].lon,plot_area[1].lat,plot_area[1].lon), 'wb'))
        
        return trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys
    
    def loadAndSaveDisplacements(self, plot_area, cache_results, Xs, Ys, involved_elements, evid):
        # load the calculated displacements
        save_displacements = False
        try:
            if plot_area == 'full':
                dX, dY, dZ = cPickle.load(open('%s/displacements_full-sys/displacement_map-event_%i.pkl'%(self.cache_dir,evid), 'rb'))
            elif plot_area == 'event':
                dX, dY, dZ = cPickle.load(open('%s/displacements_events/displacement_map-event_%i.pkl'%(self.cache_dir,evid), 'rb'))
            else:
                dX, dY, dZ = cPickle.load(open('%s/displacements_%.3f-%.3f_%.3f-%.3f/displacement_map-event_%i.pkl'%(self.cache_dir,map_area[0].lat,map_area[0].lon,map_area[1].lat,map_area[1].lon,evid), 'rb'))
            print '    *** Event %i displacements loaded from cache. ***'%evid
        except:
            dX, dY, dZ = self.calculateEventDisplacements(Xs,Ys,involved_elements,evid)#, lats.min(), lons.min(), lats.max(), lons.max())
            save_displacements = True

        # save the displacements
        if save_displacements and cache_results:
            if plot_area == 'full':
                cPickle.dump((dX, dY, dZ), open('%s/displacements_full-sys/displacement_map-event_%i.pkl'%(self.cache_dir,evid), 'wb'))
            elif plot_area == 'event':
                cPickle.dump((dX, dY, dZ), open('%s/displacements_events/displacement_map-event_%i.pkl'%(self.cache_dir,evid), 'wb'))
            else:
                cPickle.dump((dX, dY, dZ), open('%s/displacements_%.3f-%.3f_%.3f-%.3f/displacement_map-event_%i.pkl'%(self.cache_dir,map_area[0].lat,map_area[0].lon,map_area[1].lat,map_area[1].lon,evid), 'wb'))
        
        return dX, dY, dZ
    
    def eventDisplacementMapAnimation(self, start_evid, end_evid, movie_only = False, fringes = False, look_angles = None):
        event_file = self.data_file
        
        print '*** Creating displacement map animation. Start event: %i, end event: %i***'%(start_evid, end_evid)
        print '    Event file %s'%event_file
    
        arial = self.displacementMapConfig['font']
        arial_bold = self.displacementMapConfig['font_bold']
    
    # properties that are fringes dependent
        if fringes:
            cmap                = self.displacementMapConfig['cmap_f']
            water_color         = self.displacementMapConfig['water_color_f']
            boundary_color      = self.displacementMapConfig['boundary_color_f']
            coastline_color     = self.displacementMapConfig['coastline_color_f']
            country_color       = self.displacementMapConfig['country_color_f']
            state_color         = self.displacementMapConfig['state_color_f']
            fault_color         = self.displacementMapConfig['fault_color_f']
            event_fault_color   = self.displacementMapConfig['event_fault_color_f']
            map_tick_color      = self.displacementMapConfig['map_tick_color_f']
            map_frame_color     = self.displacementMapConfig['map_frame_color_f']
            grid_color          = self.displacementMapConfig['grid_color_f']
            cb_fontcolor        = self.displacementMapConfig['cb_fontcolor_f']
        else:
            cmap                = self.displacementMapConfig['cmap']
            water_color         = self.displacementMapConfig['water_color']
            boundary_color      = self.displacementMapConfig['boundary_color']
            coastline_color     = self.displacementMapConfig['coastline_color']
            country_color       = self.displacementMapConfig['country_color']
            state_color         = self.displacementMapConfig['state_color']
            fault_color         = self.displacementMapConfig['fault_color']
            event_fault_color   = self.displacementMapConfig['event_fault_color']
            map_tick_color      = self.displacementMapConfig['map_tick_color']
            map_frame_color     = self.displacementMapConfig['map_frame_color']
            grid_color          = self.displacementMapConfig['grid_color']
            cb_fontcolor        = self.displacementMapConfig['cb_fontcolor']
        
    # properties that are not fringes dependent
        land_color      = cmap(0)
        boundary_width  = self.displacementMapConfig['boundary_width']
        coastline_width = self.displacementMapConfig['coastline_width']
        country_width   = self.displacementMapConfig['country_width']
        state_width     = self.displacementMapConfig['state_width']
        river_width     = self.displacementMapConfig['river_width']
        fault_width     = self.displacementMapConfig['fault_width']
        map_resolution  = self.displacementMapConfig['map_resolution']
        plot_resolution = self.displacementMapConfig['plot_resolution']
        map_frame_width = self.displacementMapConfig['map_frame_width']
        map_fontsize    = self.displacementMapConfig['map_fontsize']
        arrow_inset     = self.displacementMapConfig['arrow_inset']
        arrow_fontsize  = self.displacementMapConfig['arrow_fontsize']
        cb_fontsize     = self.displacementMapConfig['cb_fontsize']
        cb_height       = self.displacementMapConfig['cb_height']
        cb_margin_t     = self.displacementMapConfig['cb_margin_t']
        grid_width      = self.displacementMapConfig['grid_width']
        num_grid_lines  = self.displacementMapConfig['num_grid_lines']        
        fault_color_map = LinearSegmentedColormap.from_list('fault_color_map', [fault_color,event_fault_color], N=256, gamma=1.0)
        
    #animation specific properties
        progress_tick_color = 'k'
        progress_frame_color = 'k'
        progress_frame_width = 1
        progress_line_color = 'k'
        progress_line_width = 1
        progress_indicator_line_color = 'red'
        progress_indicator_linewidth = 2
        progress_indicator_fontsize = 10.0
        
        mag_color = 'k'
        current_mag_color = 'white'
        current_mag_facecolor = 'red'
        mag_facecolor = 'white'
        mag_linewidth = 0.5
        mag_fontsize = 12
        
        mag_color_map = LinearSegmentedColormap.from_list('mag_color_map', [mag_color,current_mag_color], N=256, gamma=1.0)
        mag_line_colormap = LinearSegmentedColormap.from_list('mag_line_colormap', [progress_frame_color,progress_indicator_line_color], N=256, gamma=1.0)
        current_mag_face_colormap = LinearSegmentedColormap.from_list('current_mag_face_colormap', [mag_facecolor,current_mag_facecolor], N=256, gamma=1.0)
        
        look_arrow_fill_color = '#d0e2f7'
        look_arrow_line_color = land_color
        look_arrow_line_width = 1
        
        fm_fontsize = 9
        fm_label_color = 'white'
        fm_frame_width = 1
        fm_frame_color = 'k'
        fm_line_color = '0.0'
        fm_line_width = 1
        
        anim_target_length = 60.0
        anim_fps = 30.0
        fade_seconds = 1.0
        
        if not movie_only:
        # load the lat lon/x y points and the fault traces
            trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = self.loadAndSaveDisplacementGrid('full', True)
            
        # get the mags and the years
            f = h5py.File(event_file, 'r')
            
            events = f['event_table']
            sweeps = f['event_sweep_table']
            
            mags = [e[3] for e in events[start_evid:end_evid]]
            years_un = [e[1] for e in events[start_evid:end_evid]]
            years = [e[1] - math.floor(min(years_un)) for e in events[start_evid:end_evid]]
            big_mags = [(e[3], e[1] - math.floor(min(years_un)), e[0]) for e in events[start_evid:end_evid] if e[3] >= 6.5]
            
        # set the look angles
            if look_angles is None:
                look_azimuth = self.geometry.converter.deg_rad(135)
                look_elevation = self.geometry.converter.deg_rad(20)
            else:
                look_azimuth = self.geometry.converter.deg_rad(look_angles['azimuth'])
                look_elevation = self.geometry.converter.deg_rad(look_angles['elevation'])
            
            total_years = math.ceil(years_un[-1]) - math.floor(years_un[0])
            fpy = math.ceil(anim_target_length*anim_fps/total_years)
            current_year = math.floor(years_un[0])
            total_frames = int(fpy * total_years)
            
            section_states = {}
            mag_states = {}
            mag_images = {}
            
        # each item in displacement_map_states is {state: , mag: , image:}
            displacement_map_states = []
            
            current_mags = []
            total_events = len(mags)
            
            fm_x_ticks = np.linspace(round(min(mags)), round(max(mags)), 5)
            fm_y_ticks = np.logspace(0, 3, 4)
            fm_alpha = 0.0
            
        #create the directory to store the animation
            if not os.path.exists('%s_Displacement-Animation_%i-%i'%(self.name,start_evid,end_evid)):
                os.makedirs('%s_Displacement-Animation_%i-%i'%(self.name,start_evid,end_evid))
            
        #create the directory to store the animation images
            if fringes:
                if not os.path.exists('%s_Displacement-Animation_%i-%i/frames_f'%(self.name,start_evid,end_evid)):
                    os.makedirs('%s_Displacement-Animation_%i-%i/frames_f'%(self.name,start_evid,end_evid))
            else:
                if not os.path.exists('%s_Displacement-Animation_%i-%i/frames'%(self.name,start_evid,end_evid)):
                    os.makedirs('%s_Displacement-Animation_%i-%i/frames'%(self.name,start_evid,end_evid))
            
        # this is the background with no events. we need this for blending
            m1 = Basemap(   llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                            lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                            resolution=map_resolution, projection='cyl', suppress_ticks=True
                            )

        # set up all the plot dimensions in inches
            mw = lons_1d.shape[0]
            mh = lats_1d.shape[0]
            mwi = mw/plot_resolution
            mhi = mh/plot_resolution
            
            fig1 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
            m1.ax = fig1.add_axes((0,0,1,1))
            m1.drawmapboundary(color=boundary_color, linewidth=0, fill_color=water_color)
            m1.fillcontinents(color=land_color, lake_color=water_color)
            
            # FIGURE 1 draw the renderer
            fig1.canvas.draw()
            
            # FIGURE 1 Get the RGBA buffer from the figure
            w,h = fig1.canvas.get_width_height()
            buf = np.fromstring ( fig1.canvas.tostring_argb(), dtype=np.uint8 )
            buf.shape = ( w, h,4 )
         
            # FIGURE 1 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
            buf = np.roll ( buf, 3, axis = 2 )
            generic_bg = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
            
        #Set up the deformations. In order to do smooth fading, we are gonna be adding to these quantities
        #when there are new events and subtracting a bit each frame till we get back down to zero.
            dX = 1e-4*np.ones(Xs.shape)
            dY = 1e-4*np.ones(Xs.shape)
            dZ = 1e-4*np.ones(Xs.shape)
            
        # create the frames
            for the_frame in range(total_frames):
                year_frame = the_frame%fpy
                if year_frame == 0:
                    current_year += 1
                    events_this_year = []
                    for event in events[start_evid:end_evid]:
                        if event[1] > current_year - 1 and event[1] <= current_year:
                            events_this_year.append(event)

                progress_indicator_year = current_year - 1 - math.floor(years_un[0]) + year_frame/fpy
                print '      - frame %i of %i'%(the_frame + 1, total_frames)
                
                events_this_frame = []
                for event in events_this_year:
                    if math.modf(event[1])[0] <= float(year_frame+1)/fpy and math.modf(event[1])[0] > float(year_frame)/fpy:
                        events_this_frame.append(event)

            # remove a fixed amount from the displacements. this is the decay that happens for deformations
            # that have already happened

                dX *= 1.0 - 1.0/(anim_fps*fade_seconds)
                dY *= 1.0 - 1.0/(anim_fps*fade_seconds)
                dZ *= 1.0 - 1.0/(anim_fps*fade_seconds)
                
                if len(events_this_frame) > 0:
                    #there are events in this frame. get their net displacements.
                    for event in events_this_frame:
                        # load the calculated displacements
                        event_sweeps = sweeps[event[8]:event[9]]
                        
                        involved_elements_tmp = {}
                        involved_sections = {}
                        involved_elements = []
                        
                        for sweep in event_sweeps:
                            eid = sweep[2]
                            sid = self.geometry.elements[eid].sid
                            try:
                                involved_elements_tmp[eid] += sweep[3]
                            except KeyError:
                                involved_elements_tmp[eid] = sweep[3]
                            
                            involved_sections[sid] = sid
                            section_states[sid] = 1.0
                        
                        involved_sections = involved_sections.keys()
                        
                        for eid, slip in involved_elements_tmp.iteritems():
                            involved_elements.append({'eid':eid,'slip':slip})
                        
                        dX_tmp, dY_tmp, dZ_tmp = self.loadAndSaveDisplacements('full', True, Xs, Ys, involved_elements, event[0])
                        
                        dX += dX_tmp
                        dY += dY_tmp
                        dZ += dZ_tmp

                        current_mags.append(event[3])
                    
            #now plot the net displacements
                dMags = -dX * math.sin(look_azimuth) * math.cos(look_elevation) - dY * math.cos(look_azimuth) * math.cos(look_elevation) + dZ * math.sin(look_elevation)
            
            #m2, fig2 is the plotted data
                m2 = Basemap(   llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                                lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                                resolution=map_resolution, projection='cyl', suppress_ticks=True
                                )
                fig2 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
                m2.ax = fig2.add_axes((0,0,1,1))
                
            #prepare the colors for the plot
                wavelength = 0.03
                dMags_transformed = m2.transform_scalar(dMags, lons_1d, lats_1d, Xs.shape[1], Xs.shape[0])
                dMags_colors = np.empty((dMags_transformed.shape[0],dMags_transformed.shape[1],4))
                
                if fringes:
                    it = np.nditer(dMags_transformed, flags=['multi_index'])
                    while not it.finished:
                        r,g,b,a = cmap(math.modf(abs(dMags_transformed[it.multi_index])/wavelength)[0])
                        dMags_colors[it.multi_index[0], it.multi_index[1], 0] = r
                        dMags_colors[it.multi_index[0], it.multi_index[1], 1] = g
                        dMags_colors[it.multi_index[0], it.multi_index[1], 2] = b
                        dMags_colors[it.multi_index[0], it.multi_index[1], 3] = a
                        it.iternext()
                    im = m2.imshow(dMags_colors, interpolation='spline36')
                else:
                    dMags_colors = np.fabs(dMags_transformed)
                    mod_vmax = 1000
                    im = m2.imshow(dMags_colors, cmap=cmap, norm=mpl.colors.LogNorm(vmin=1e-4, vmax=mod_vmax, clip=True))

            #m3, fig3 is the ocean land mask
                m3 = Basemap(   llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                                lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                                resolution=map_resolution, projection='cyl', suppress_ticks=True
                                )
                fig3 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
                m3.ax = fig3.add_axes((0,0,1,1))
                m3.fillcontinents(color='#000000', lake_color='#ffffff')
            
            # composite image 2 and 3 together with the generic bg
                # FIGURE 2 draw the renderer
                fig2.canvas.draw()
                
                # FIGURE 2 Get the RGBA buffer from the figure
                w,h = fig2.canvas.get_width_height()
                buf = np.fromstring ( fig2.canvas.tostring_argb(), dtype=np.uint8 )
                buf.shape = ( w, h,4 )
             
                # FIGURE 2 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
                buf = np.roll ( buf, 3, axis = 2 )
                im2 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
                
                # FIGURE 3 draw the renderer
                fig3.canvas.draw()
                
                # FIGURE 3 Get the RGBA buffer from the figure
                w,h = fig3.canvas.get_width_height()
                buf = np.fromstring ( fig3.canvas.tostring_argb(), dtype=np.uint8 )
                buf.shape = ( w, h,4 )
             
                # FIGURE 3 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
                buf = np.roll ( buf, 3, axis = 2 )
                im3 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
                
                mask = im3.convert('L')
                
                map_image = Image.composite(generic_bg, im2, mask)
                fig2.clf()
                fig3.clf()

            #do other plotting tasks. blending, fading etc.
                
            #now plot all of the geo data to go over the composited rupture
                #this is the map to overlay over the map image
                m4 = Basemap(    llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                                lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                                resolution=map_resolution, projection='cyl'
                                )
                
                ph = 768.0
                pw = 1024.0
                
                full_pwi = pw/plot_resolution
                full_phi = ph/plot_resolution
                
                fig4 = plt.figure(figsize=(full_pwi, full_phi), dpi=plot_resolution)

                width_frac = mw/pw
                height_frac = mh/ph
                left_frac = 70.0/pw
                bottom_frac = 70.0/ph
                
                m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                
                m4.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=(1,1,1,0))
                # draw coastlines, edge of map.
                m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)
                # draw countries
                m4.drawcountries(linewidth=country_width, color=country_color)
                # draw states
                m4.drawstates(linewidth=state_width, color=state_color)
                # draw parallels.
                parallels = np.linspace(lats.min(), lats.max(), num_grid_lines+1)
                m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], fontsize=map_fontsize, color=grid_color, fontproperties=arial, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
                # draw meridians
                meridians = np.linspace(lons.min(), lons.max(), num_grid_lines+1)
                m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], fontsize=map_fontsize, color=grid_color, fontproperties=arial, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
                
                m4.imshow(map_image, origin='upper')
                
            # plot the magnitude/progress indicator
                width_frac = 184.0/pw
                height_frac = mh/ph
                left_frac = 790.0/pw
                bottom_frac = 70.0/ph
                
                mag_vs_year = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                mag_vs_year.plot(mags, years, color=progress_line_color, linewidth=0, marker='.', mew=0, ms=4, mfc=progress_line_color)
                mag_vs_year.autoscale(enable=True, axis='both', tight=True)
                mag_vs_year.set_ylim(math.floor(min(years)), math.ceil(max(years)))
                mag_vs_year_alt = mag_vs_year.twinx()
                mag_vs_year.set_xticks(np.linspace(min(mags),max(mags),3))
                
                mag_vs_year_alt.set_ylim(math.floor(min(years)), math.ceil(max(years)))
                
                mag_vs_year_tick_labels = []
                for tick in np.linspace(min(mags),max(mags),3):
                    mag_vs_year_tick_labels.append('%.1f'%tick)
                mag_vs_year.set_xticklabels(mag_vs_year_tick_labels)
                
                for tick in mag_vs_year.xaxis.get_major_ticks():
                    tick.label1.set_fontproperties(arial)
                    tick.label1.set_fontsize(progress_indicator_fontsize)
                    tick.label1.set_color(progress_tick_color)
                
                for tl in mag_vs_year_alt.get_yticklabels():
                    tl.set_fontproperties(arial)
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
                
                mag_vs_year.set_xlabel('magnitude', fontproperties=arial, size=progress_indicator_fontsize, color=progress_tick_color)
                mag_vs_year_alt.set_ylabel('year', fontproperties=arial, size=progress_indicator_fontsize, color=progress_tick_color)
            
            # plot the current frequency-magnitude
                width_frac = 150.0/pw
                height_frac = 150.0/ph
                left_frac = 100.0/pw
                bottom_frac = 100.0/ph
                
                current_total_events = len(current_mags)
                if current_total_events > 1:
                    cum_freq = {}
                    fm_x = []
                    fm_y = []
                    for num, magnitude in enumerate(sorted(current_mags)):
                        cum_freq['%.10f'%magnitude] = current_total_events - (num + 1)
                    #print cum_freq
                    for magnitude in sorted(cum_freq.iterkeys()):
                        fm_x.append(magnitude)
                        fm_y.append(float(cum_freq[magnitude]))
                    mag_vs_freq = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                    mag_vs_freq.semilogy(fm_x, fm_y, color=fm_line_color, linewidth=fm_line_width)
                    mag_vs_freq.set_ylim(bottom=min(fm_y_ticks), top=max(fm_y_ticks))
                    mag_vs_freq.set_xlim(left=round(min(mags)), right=round(max(mags)))
                    mag_vs_freq.set_yticks(fm_y_ticks)
                    mag_vs_freq.set_xticks(fm_x_ticks)
                    
                    mag_vs_freq_x_ticklabels = []
                    for tick in fm_x_ticks:
                        mag_vs_freq_x_ticklabels.append('%.1f'%tick)
                    mag_vs_freq.set_xticklabels(mag_vs_freq_x_ticklabels)
                    
                    #mag_vs_freq_y_ticklabels = []
                    #for tick in fm_y_ticks:
                    #    mag_vs_freq_y_ticklabels.append('%.1g'%tick)
                    #mag_vs_freq.set_yticklabels(mag_vs_freq_y_ticklabels)
                    
                    for tick in mag_vs_freq.xaxis.get_major_ticks() + mag_vs_freq.yaxis.get_major_ticks():
                        tick.label1.set_fontproperties(arial)
                        tick.label1.set_fontsize(fm_fontsize)
                        tick.label1.set_color(fm_label_color)
                        tick.label1.set_alpha(fm_alpha)
                    
                    for line in mag_vs_freq.xaxis.get_majorticklines() + mag_vs_freq.yaxis.get_majorticklines() + mag_vs_freq.yaxis.get_minorticklines():
                        line.set_alpha(0)
                    
                    #take care of all of the frame line widths
                    for spine in mag_vs_freq.spines.itervalues():
                        spine.set_lw(fm_frame_width)
                        spine.set_color(fm_frame_color)
                        spine.set_alpha(fm_alpha)
                    mag_vs_freq.patch.set_alpha(fm_alpha)
                    if fm_alpha < 1:
                        fm_alpha += 0.04
                    else:
                        fm_alpha = 1

            # progress indicator 
                mag_vs_year.axhline(y=progress_indicator_year, lw=progress_indicator_linewidth, c=progress_indicator_line_color)
            
            # print label lines for the progress meter
                label_lines = []
                if len(big_mags) < 10:
                    label_range = [bm[1] for bm in big_mags]
                else:
                    label_range = np.linspace(math.floor(min(years)),math.ceil(max(years)),len(big_mags))
                for i, x in enumerate(label_range):
                    m = big_mags[i][0]
                    y = big_mags[i][1]
                    if big_mags[i][2] in [e[0] for e in events_this_frame]:
                        mag_states[big_mags[i][2]] = 1.0
                        mag_images[big_mags[i][2]] = map_image
                    
                    try:
                        mag_state = mag_states[int(big_mags[i][2])]
                    except KeyError:
                        mag_state = 0
                    
                    
                    the_color = mag_color_map(mag_state)
                    the_line_color = mag_line_colormap(mag_state)
                    the_line_width = mag_linewidth + mag_state * (progress_indicator_linewidth)
                    the_bbox = dict(facecolor=current_mag_face_colormap(mag_state), linewidth=the_line_width, ec=the_line_color, boxstyle='round4,pad=0.5')
                    mag_vs_year.text(min(mags) - 0.4, x, '%.1f'%m, fontproperties=arial,
                                            fontsize=mag_fontsize, horizontalalignment='right', 
                                            verticalalignment='center', color=the_color, bbox=the_bbox)
                    
                    label_lines.append(mlines.Line2D([1.0,-0.01,-0.05,-0.085],[y/math.ceil(max(years)),y/math.ceil(max(years)),x/math.ceil(max(years)),x/math.ceil(max(years))],
                                        linewidth=the_line_width, transform=mag_vs_year.transAxes, color=the_line_color, solid_capstyle='round', solid_joinstyle='round'))
                    
                    if mag_state >= 0.04:
                        mag_states[int(big_mags[i][2])] -= 0.04
                    else:
                        mag_states[int(big_mags[i][2])] = 0.0
                
                mag_vs_year.lines.extend(label_lines)
                
            # print faults on lon-lat plot
                for sid, sec_trace in trace_vectors.iteritems():
                    trace_lons = []
                    trace_lats = [] 
                    for trace_ele in sec_trace:
                        trace_lons.append(trace_ele.lon)
                        trace_lats.append(trace_ele.lat)
                    trace_Xs, trace_Ys =m4(trace_lons, trace_lats)
                    try:
                        section_state = section_states[sid]
                    except KeyError:
                        section_state = 0
                    #if sid in involved_sections:
                    m4.plot(trace_Xs, trace_Ys, color=fault_color_map(section_state), linewidth=fault_width + section_state * 2.0, solid_capstyle='round', solid_joinstyle='round')
                    if section_state >= 0.04:
                        section_states[sid] -= 0.04
                    else:
                        section_states[sid] = 0.0
            
            #plot the cb
                left_frac = 70.0/pw
                bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
                width_frac = mw/pw
                height_frac = cb_height/ph
                
                
                cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                if fringes:
                    norm = mpl.colors.Normalize(vmin=0, vmax=wavelength)
                    cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap,
                       norm=norm,
                       orientation='horizontal')
                    cb_ax.set_title('Displacement [m]', fontproperties=arial, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )
                else:
                    norm = mpl.colors.LogNorm(vmin=1e-4, vmax=mod_vmax, clip=True)
                    cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap,
                               norm=norm,
                               orientation='horizontal')
                    cb_ax.set_title('Total displacement [m]', fontproperties=arial, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )
                
                for label in cb_ax.xaxis.get_ticklabels():
                    label.set_fontproperties(arial)
                    label.set_fontsize(cb_fontsize)
                    label.set_color(cb_fontcolor)
                for line in cb_ax.xaxis.get_ticklines():
                    line.set_alpha(0)
            
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

                az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(look_azimuth)
                az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(look_azimuth)
                az_arrow_dx      = 0.8*math.sin(look_azimuth)
                az_arrow_dy      = 0.8*math.cos(look_azimuth)

                az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
                az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - self.geometry.converter.rad_deg(look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
                az_ax.text(1.0, 1.0, 'az = %.1f%s'%(self.geometry.converter.rad_deg(look_azimuth),r'$^{\circ}$'), fontproperties=arial_bold, size=arrow_fontsize, ha='right', va='top')

        
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

                al_arrow_start_x    = 0.1 + 0.8*math.cos(look_elevation)
                al_arrow_start_y    = 0.1 + 0.8*math.sin(look_elevation)
                al_arrow_dx      = -0.8*math.cos(look_elevation)
                al_arrow_dy      = -0.8*math.sin(look_elevation)

                al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
                al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=self.geometry.converter.rad_deg(look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
                al_ax.text(1.0, 1.0, 'al = %.1f%s'%(self.geometry.converter.rad_deg(look_elevation),r'$^{\circ}$'), fontproperties=arial_bold, size=arrow_fontsize, ha='right', va='top')

                if fringes:
                    fig4.savefig('%s_Displacement-Animation_%i-%i/frames_f/%i.png'%(self.name,start_evid,end_evid,the_frame), format='png', dpi=plot_resolution)
                else:
                    fig4.savefig('%s_Displacement-Animation_%i-%i/frames/%i.png'%(self.name,start_evid,end_evid,the_frame), format='png', dpi=plot_resolution)
                
                fig4.clf()
                plt.close()
                
                gc.collect()
            
            f.close()
        print '    *** Creating movie ***'
        if fringes:
            proc_args = "ffmpeg -y -r %i -f image2 -i %s_Displacement-Animation_%i-%i/frames_f/%s.png %s_Displacement-Animation_%i-%i/%s_Fringes-Displacement-Animation_%i-%i.mp4"%(anim_fps,self.name,start_evid,end_evid,'%d',self.name,start_evid,end_evid,self.name,start_evid,end_evid)
        else:
            proc_args = "ffmpeg -y -r %i -f image2 -i %s_Displacement-Animation_%i-%i/frames/%s.png %s_Displacement-Animation_%i-%i/%s_Displacement-Animation_%i-%i.mp4"%(anim_fps,self.name,start_evid,end_evid,'%d',self.name,start_evid,end_evid,self.name,start_evid,end_evid)
        proc = subprocess.Popen(proc_args, shell=True)
        proc.wait()
        print '*** Done ***'
    
    def eventDisplacementMapAll(self):
        event_file = self.data_file
        # figure out what elements are involved and find the total slip on each
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        
        num_events = len(events)
        
        f.close()
    
                
        startid = 1911
        numevents = 1
        for evid in range(num_events):
            if evid >= startid and evid <= startid+numevents:
                self.eventDisplacementMap(evid, plot_area = 'full', cache_results=False)
                #self.eventDisplacementMap(event_file, evid, plot_area = 'event', cache_results=False)
            elif evid > startid+numevents:
                break
        
        msg = MIMEText('events %i through %i complete'%(startid, startid+numevents))

        msg['Subject'] = 'events %i through %i complete'%(startid, startid+numevents)
        msg['From'] = 'msachs@me.com'
        msg['To'] = 'mksachs@ucdavis.edu'
    
        s = smtplib.SMTP('smtp.me.com',25)
        s.ehlo()
        s.starttls()
        s.ehlo()
        s.login('msachs','password')
        s.sendmail('msachs@me.com', 'mksachs@ucdavis.edu', msg.as_string())
        s.quit()

        
    def eventDisplacementMap(self, evid, plot_area = 'event', fringes=False, cache_results=True, look_angles = None):
        event_file = self.data_file
        print '*** Plotting event %i displacement map ***'%evid
        if plot_area == 'event':
            print '    Map area defined by event size'
        elif plot_area == 'full':
            print '    Map area for the full fault system'
        else:
            print '    Map area: NE corner (%.3f, %.3f) to SW corner (%.3f, %.3f)'%(plot_area[0].lat, plot_area[0].lon, plot_area[1].lat, plot_area[1].lon)
        print '    Event file %s'%event_file
    
        arial = self.displacementMapConfig['font']
        arial_bold = self.displacementMapConfig['font_bold']
    
    # properties that are fringes dependent
        if fringes:
            cmap            = self.displacementMapConfig['cmap_f']
            water_color     = self.displacementMapConfig['water_color_f']
            boundary_color  = self.displacementMapConfig['boundary_color_f']
            coastline_color = self.displacementMapConfig['coastline_color_f']
            country_color   = self.displacementMapConfig['country_color_f']
            state_color     = self.displacementMapConfig['state_color_f']
            fault_color     = self.displacementMapConfig['fault_color_f']
            map_tick_color  = self.displacementMapConfig['map_tick_color_f']
            map_frame_color = self.displacementMapConfig['map_frame_color_f']
            grid_color      = self.displacementMapConfig['grid_color_f']
            cb_fontcolor    = self.displacementMapConfig['cb_fontcolor_f']
        else:
            cmap            = self.displacementMapConfig['cmap']
            water_color     = self.displacementMapConfig['water_color']
            boundary_color  = self.displacementMapConfig['boundary_color']
            coastline_color = self.displacementMapConfig['coastline_color']
            country_color   = self.displacementMapConfig['country_color']
            state_color     = self.displacementMapConfig['state_color']
            fault_color     = self.displacementMapConfig['fault_color']
            map_tick_color  = self.displacementMapConfig['map_tick_color']
            map_frame_color = self.displacementMapConfig['map_frame_color']
            grid_color      = self.displacementMapConfig['grid_color']
            cb_fontcolor    = self.displacementMapConfig['cb_fontcolor']
        
    # properties that are not fringes dependent
        land_color      = cmap(0)
        boundary_width  = self.displacementMapConfig['boundary_width']
        coastline_width = self.displacementMapConfig['coastline_width']
        country_width   = self.displacementMapConfig['country_width']
        state_width     = self.displacementMapConfig['state_width']
        river_width     = self.displacementMapConfig['river_width']
        fault_width     = self.displacementMapConfig['fault_width']
        map_resolution  = self.displacementMapConfig['map_resolution']
        plot_resolution = self.displacementMapConfig['plot_resolution']
        map_frame_width = self.displacementMapConfig['map_frame_width']
        map_fontsize    = self.displacementMapConfig['map_fontsize']
        arrow_inset     = self.displacementMapConfig['arrow_inset']
        arrow_fontsize  = self.displacementMapConfig['arrow_fontsize']
        cb_fontsize     = self.displacementMapConfig['cb_fontsize']
        cb_height       = self.displacementMapConfig['cb_height']
        cb_margin_t     = self.displacementMapConfig['cb_margin_t']
        grid_width      = self.displacementMapConfig['grid_width']
        num_grid_lines  = self.displacementMapConfig['num_grid_lines']
        
    # create filename prepend
        '''
        if self.section_filter is None:
            file_name_prepend = '%s_%i_'%(self.name,evid)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_%i_'%(self.name,section_str,evid)
        '''
        
        file_name_prepend = self.fileNamePrepend()
        
    # figure out what elements are involved and find the total slip on each
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        event = events[evid]
        
        event_sweeps = sweeps[event[8]:event[9]]
        
        involved_elements_tmp = {}
        involved_sections = {}
        involved_elements = []
        
        for sweep in event_sweeps:
            eid = sweep[2]
            sid = self.geometry.elements[eid].sid
            try:
                involved_elements_tmp[eid] += sweep[3]
            except KeyError:
                involved_elements_tmp[eid] = sweep[3]
            
            involved_sections[sid] = sid
        
        involved_sections = involved_sections.keys()
    
    # set the look angles
        strikes = []
        rakes = []
        
        for eid, slip in involved_elements_tmp.iteritems():
            strikes.append(self.geometry.elements[eid].strike())
            rakes.append(self.geometry.elements[eid].get_rake())
            involved_elements.append({'eid':eid,'slip':slip})
        
        if look_angles is None:
            look_azimuth = -sum(strikes)/len(strikes)
            average_rake = sum(rakes)/len(rakes)
            
            if average_rake >= math.pi/2.0:
                average_rake = math.pi - average_rake
            look_elevation = abs(average_rake)
        else:
            look_azimuth = self.geometry.converter.deg_rad(look_angles['azimuth'])
            look_elevation = self.geometry.converter.deg_rad(look_angles['elevation'])
        
    # load the lat lon/x y points and the fault traces
        trace_vectors, lons, lats, lons_1d, lats_1d, Xs, Ys = self.loadAndSaveDisplacementGrid(plot_area, cache_results, involved_sections, map_resolution, evid=evid)

    # load the calculated displacements
        dX, dY, dZ = self.loadAndSaveDisplacements(plot_area, cache_results, Xs, Ys, involved_elements, evid)
        
        print '    *** Preparing plot ***'

        dMags = -dX * math.sin(look_azimuth) * math.cos(look_elevation) - dY * math.cos(look_azimuth) * math.cos(look_elevation) + dZ * math.sin(look_elevation)
        #dMags = np.empty(dX.shape)
        
        #it = np.nditer(dX, flags=['multi_index'])
        #while not it.finished:
        #    dMags[it.multi_index] = -dX[it.multi_index] * math.sin(look_azimuth) * math.cos(look_elevation) - dY[it.multi_index] * math.cos(look_azimuth) * math.cos(look_elevation) + dZ[it.multi_index] * math.sin(look_elevation)
        #    it.iternext()
        
        print '    *** Plotting ***'

    #m1, fig1 is the oceans and the continents. This will lie behind the masked data image
    #    initilizing this here so we can use it for various measurements
        m1 = Basemap(   llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                        lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                        resolution=map_resolution, projection='cyl', suppress_ticks=True
                        )

        # set up all the plot dimensions in inches
        mw = lons_1d.shape[0]
        mh = lats_1d.shape[0]
        mwi = mw/plot_resolution
        mhi = mh/plot_resolution
        
        fig1 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        m1.ax = fig1.add_axes((0,0,1,1))
        m1.drawmapboundary(color=boundary_color, linewidth=0, fill_color=water_color)
        m1.fillcontinents(color=land_color, lake_color=water_color)
        
    #m2, fig2 is the plotted deformation data
        m2 = Basemap(   llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                        lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                        resolution=map_resolution, projection='cyl', suppress_ticks=True
                        )
        fig2 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        m2.ax = fig2.add_axes((0,0,1,1))
        
    #prepare the colors for the plot
        wavelength = 0.03
        dMags_transformed = m2.transform_scalar(dMags, lons_1d, lats_1d, Xs.shape[1], Xs.shape[0])
        dMags_colors = np.empty((dMags_transformed.shape[0],dMags_transformed.shape[1],4))
        
        if fringes:
            it = np.nditer(dMags_transformed, flags=['multi_index'])
            while not it.finished:
                r,g,b,a = cmap(math.modf(abs(dMags_transformed[it.multi_index])/wavelength)[0])
                dMags_colors[it.multi_index[0], it.multi_index[1], 0] = r
                dMags_colors[it.multi_index[0], it.multi_index[1], 1] = g
                dMags_colors[it.multi_index[0], it.multi_index[1], 2] = b
                dMags_colors[it.multi_index[0], it.multi_index[1], 3] = a
                it.iternext()
            im = m2.imshow(dMags_colors, interpolation='spline36')
        else:
            dMags_colors = np.fabs(dMags_transformed)
            vmax = np.amax(dMags_colors)
            if vmax <= 1:
                mod_vmax = 1
            elif vmax > 1 and vmax <= 10:
                mod_vmax = 10
            elif vmax > 10 and vmax <= 100:
                mod_vmax = 100
            elif vmax > 100 and vmax <= 1000:
                mod_vmax = 1000
            elif vmax > 1000:
                mod_vmax = 1000
            im = m2.imshow(dMags_colors, cmap=cmap, norm=mpl.colors.LogNorm(vmin=1e-4, vmax=mod_vmax, clip=True))

    #m3, fig3 is the ocean land mask
        m3 = Basemap(   llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                        lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                        resolution=map_resolution, projection='cyl', suppress_ticks=True
                        )
        fig3 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        m3.ax = fig3.add_axes((0,0,1,1))
        m3.fillcontinents(color='#000000', lake_color='#ffffff')
    
    # composite image 1 - 3 together
        # FIGURE 1 draw the renderer
        fig1.canvas.draw()
        
        # FIGURE 1 Get the RGBA buffer from the figure
        w,h = fig1.canvas.get_width_height()
        buf = np.fromstring ( fig1.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 1 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im1 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        # FIGURE 2 draw the renderer
        fig2.canvas.draw()
        
        # FIGURE 2 Get the RGBA buffer from the figure
        w,h = fig2.canvas.get_width_height()
        buf = np.fromstring ( fig2.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 2 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im2 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        # FIGURE 3 draw the renderer
        fig3.canvas.draw()
        
        # FIGURE 3 Get the RGBA buffer from the figure
        w,h = fig3.canvas.get_width_height()
        buf = np.fromstring ( fig3.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 3 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im3 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        mask = im3.convert('L')
        
        map_image = Image.composite(im1, im2, mask)
        
    # now plot all of the geo data over the composited rupture
        #this is the map to overlay over the map image
        m4 = Basemap(    llcrnrlon=lons.min(), llcrnrlat=lats.min(), urcrnrlon=lons.max(), urcrnrlat=lats.max(),\
                        lat_0=(lats.max()+lats.min())/2.0, lon_0=(lons.max()+lons.min())/2.0,\
                        resolution=map_resolution, projection='cyl'
                        )

        # the sizing for the image is tricky. the aspect ratio of the plot is fixed, so we cant set all of margins to whatever we want.
        # we will set the anchor to the top, left margin position. then scale the image based on the bottom/right margin, whichever is bigger.
        
        if mh > mw:
            ph = 768.0
            pw = mw + 70.0 + 40.0
        else:
            pw = 790.0
            ph = mh + 70.0 + 40.0
            
        
        full_pwi = pw/plot_resolution
        full_phi = ph/plot_resolution
        
        fig4 = plt.figure(figsize=(full_pwi, full_phi), dpi=plot_resolution)

        width_frac = mw/pw
        height_frac = mh/ph
        left_frac = 70.0/pw
        bottom_frac = 70.0/ph
        
        m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
        
        involved_sections = {}
        
        # draw coastlines, edge of map.
        m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)
        # draw countries
        m4.drawcountries(linewidth=country_width, color=country_color)
        # draw states
        m4.drawstates(linewidth=state_width, color=state_color)
        # draw parallels.
        parallels = np.linspace(lats.min(), lats.max(), num_grid_lines+1)
        m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], fontsize=map_fontsize, color=grid_color, fontproperties=arial, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
        # draw meridians
        meridians = np.linspace(lons.min(), lons.max(), num_grid_lines+1)
        m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], fontsize=map_fontsize, color=grid_color, fontproperties=arial, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
        
        #for m in m4_meridians.itervalues():
        #    if len(m[1]) != 0:
        #        m[1][0].set_color(map_tick_color)
        
        #for p in m4_parallels.itervalues():
        #    if len(p[1]) != 0:
        #        p[1][0].set_color(map_tick_color)
        
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

        az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(look_azimuth)
        az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(look_azimuth)
        az_arrow_dx      = 0.8*math.sin(look_azimuth)
        az_arrow_dy      = 0.8*math.cos(look_azimuth)

        az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
        az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - self.geometry.converter.rad_deg(look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
        az_ax.text(1.0, 1.0, 'az = %.1f%s'%(self.geometry.converter.rad_deg(look_azimuth),r'$^{\circ}$'), fontproperties=arial_bold, size=arrow_fontsize, ha='right', va='top')
        
        
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

        al_arrow_start_x    = 0.1 + 0.8*math.cos(look_elevation)
        al_arrow_start_y    = 0.1 + 0.8*math.sin(look_elevation)
        al_arrow_dx      = -0.8*math.cos(look_elevation)
        al_arrow_dy      = -0.8*math.sin(look_elevation)

        al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
        al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=self.geometry.converter.rad_deg(look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
        al_ax.text(1.0, 1.0, 'al = %.1f%s'%(self.geometry.converter.rad_deg(look_elevation),r'$^{\circ}$'), fontproperties=arial_bold, size=arrow_fontsize, ha='right', va='top')

        
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
        
        mag_ax.text(0.5, 0.5, 'm = %.3f'%(float(event[3])), fontproperties=arial_bold, size=arrow_fontsize, ha='center', va='center')

        m4.imshow(map_image, origin='upper')
        
    # print faults on lon-lat plot
        section_states = {}
        for sweep in event_sweeps:
            eid = sweep[2]
            sid = self.geometry.elements[eid].sid
            section_states[sid] = 1.0
        for sid, sec_trace in trace_vectors.iteritems():
            trace_lons = []
            trace_lats = [] 
            for trace_ele in sec_trace:
                trace_lons.append(trace_ele.lon)
                trace_lats.append(trace_ele.lat)
            trace_Xs, trace_Ys =m4(trace_lons, trace_lats)
            
            try:
                section_state = section_states[sid]
            except KeyError:
                section_state = 0
            
            m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width + section_state * 0.0, solid_capstyle='round', solid_joinstyle='round')
    
    #plot the cb
    
        cb_left_frac = 70.0/pw
        cb_bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
        cb_width_frac = width_frac
        cb_height_frac = cb_height/ph
        
        
        cb_ax = fig4.add_axes([cb_left_frac, cb_bottom_frac, cb_width_frac, cb_height_frac])
        if fringes:
            norm = mpl.colors.Normalize(vmin=0, vmax=wavelength)
            cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap,
                       norm=norm,
                       orientation='horizontal')
            cb_ax.set_title('Displacement [m]', fontproperties=arial, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )
        else:
            norm = mpl.colors.LogNorm(vmin=1e-4, vmax=mod_vmax, clip=True)
            cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap,
                       norm=norm,
                       orientation='horizontal')
            cb_ax.set_title('Total displacement [m]', fontproperties=arial, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )
        

        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(arial)
            label.set_fontsize(cb_fontsize)
            label.set_color(cb_fontcolor)
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        file_name_prepend += '%i_'%evid

        if fringes:
            file_name_prepend += 'fringes-'

        #save the image
        if plot_area == 'event':
            fig4.savefig('%sdisplacement-map.png'%file_name_prepend, format='png', dpi=plot_resolution)
        elif plot_area == 'full':
            fig4.savefig('%sfull-sys-displacement-map.png'%file_name_prepend, format='png', dpi=plot_resolution)
        else:
            fig4.savefig('%s%.2f-%.2f_%.2f-%.2f_displacement-map.png'%(file_name_prepend,plot_area[0].lat,plot_area[0].lon,plot_area[1].lat,plot_area[1].lon), format='png', dpi=plot_resolution)

        f.close()

        fig1.clf()
        fig2.clf()
        fig3.clf()
        fig4.clf()
        plt.close()
            
        gc.collect()
        
        #it = np.nditer(dMags, flags=['multi_index'])
        #while not it.finished:
        #    print dMags[it.multi_index]
        #    it.iternext()
        
        
        print '*** Done ***'
        
    def infoTheoryActions(self, start_event, end_event):
        event_file = self.data_file
        print '*** Performing information theory actions ***'
        print '    Event file %s'%event_file
        
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        if end_event is not None:
            if end_event < start_event:
                print '!!! The end event must be bigger than the start event !!!'
                return
        else:
            end_event = events[-1][0]
        
        if self.section_filter is None:
            file_name_prepend = '%s_%i-%i_'%(self.name,start_event,end_event)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_%i-%i_'%(self.name,section_str,start_event,end_event)
        
        #aftershocks = f['aftershock_table']
        #bgevents = f['background_event_table']
        
        print '    start event: %i, end event: %i'%(start_event, end_event)
        
        cff_data = False
        element_data = False
        
        if self.element_CFF is not False:
            cff_data = True
            element_CFF_plotter                      = ElementCFFPlotter(self)
            element_CFF_plotter.name                 = self.name
            if self.element_CFF == 'both' or self.element_CFF == 'plot':
                element_CFF_plotter.output_format    = self.output_format
                element_CFF_plotter.out_file         = '%sElement-CFF.%s'%(file_name_prepend,self.output_format)
            if self.element_CFF == 'both' or self.element_CFF == 'text':
                element_CFF_plotter.out_text_file    = '%sElement-CFF.dat'%(file_name_prepend)
        
        if self.element_CFF_map:
            cff_data = True
            element_CFF_map_plotter                 = ElementCFFMapPlotter(self)
            element_CFF_map_plotter.output_format   = self.output_format
            element_CFF_map_plotter.out_file        = '%sElement-CFF-map.%s'%(file_name_prepend,self.output_format)
        
        if self.element_CFF_block_entropy is not None:
            element_data = True
            element_CFF_block_entropy_plotter                   = ElementCFFBlockEntropyPlotter(self)
            element_CFF_block_entropy_plotter.output_format     = self.output_format
            element_CFF_block_entropy_plotter.max_word_length   = self.element_CFF_block_entropy
            element_CFF_block_entropy_plotter.type              = 'be'
            element_CFF_block_entropy_plotter.out_file          = '%sElement-CFF-block-entropy.%s'%(file_name_prepend,self.output_format)
        
        if self.element_CFF_delta_block_entropy is not None:
            element_data = True
            element_CFF_delta_block_entropy_plotter                 = ElementCFFBlockEntropyPlotter(self)
            element_CFF_delta_block_entropy_plotter.output_format   = self.output_format
            element_CFF_delta_block_entropy_plotter.max_word_length = self.element_CFF_delta_block_entropy
            element_CFF_delta_block_entropy_plotter.type            = 'dbe'
            element_CFF_delta_block_entropy_plotter.out_file        = '%sElement-CFF-delta-block-entropy.%s'%(file_name_prepend,self.output_format)
        
        if self.element_CFF_E_Hmu is not None:
            element_data = True
            element_CFF_E_Hmu_plotter                 = ElementCFFEHmuPlotter(self)
            element_CFF_E_Hmu_plotter.output_format   = self.output_format
            element_CFF_E_Hmu_plotter.max_word_length = self.element_CFF_E_Hmu
            #element_CFF_E_Hmu_plotter.type            = 'dbe'
            element_CFF_E_Hmu_plotter.out_file        = '%sElement-E-Hmu.%s'%(file_name_prepend,self.output_format)
        
        state_data = ElementStateData(self, cff_data, element_data, max((self.element_CFF_E_Hmu, self.element_CFF_delta_block_entropy, self.element_CFF_block_entropy, 15)))
        
        if self.element_CFF is not False:
            element_CFF_plotter.plot(state_data, start_event, end_event)
        
        if self.element_CFF_map:
            element_CFF_map_plotter.plot(state_data)
        
        if self.element_CFF_block_entropy is not None:
            element_CFF_block_entropy_plotter.plot(state_data)
        
        if self.element_CFF_delta_block_entropy is not None:
            element_CFF_delta_block_entropy_plotter.plot(state_data)
        
        if self.element_CFF_E_Hmu is not None:
            element_CFF_E_Hmu_plotter.plot(state_data)
        
        f.close()
        print '*** Done ***'
        
        return file_name_prepend
    
    def eventActions(self, start_event, end_event):
        event_file = self.data_file
        print '*** Performing event actions ***'
        print '    Event file %s'%event_file
        
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        if end_event is not None:
            if end_event < start_event:
                print '!!! The end event must be bigger than the start event !!!'
                return
        else:
            end_event = events[-1][0]
        
        file_name_prepend = self.fileNamePrepend(event_range=[start_event, end_event])
        file_name_prepend_dat = self.fileNamePrepend(event_range=[start_event, end_event], type='dat')
        file_name_prepend_eqsim = self.fileNamePrepend(event_range=[start_event, end_event], type='eqsim')
        
        '''
        if self.section_filter is None:
            file_name_prepend = '%s_%i-%i_'%(self.name,start_event,end_event)
        else:
            section_str = ''
            if len(self.section_filter) < 5:
                for sec_id in self.section_filter:
                    section_str += '%i-'%sec_id
                section_str = section_str.rpartition('-')[0]
            else:
                section_str = '%i--%i'%(self.section_filter[0],self.section_filter[-1])
            file_name_prepend = '%s_%s_%i-%i_'%(self.name,section_str,start_event,end_event)
        '''
        #aftershocks = f['aftershock_table']
        #bgevents = f['background_event_table']
        
        print '    start event: %i, end event: %i'%(start_event, end_event)
        
        scan_events = True
        cff_data = False
        
        if self.element_CFF_E_Hmu                       is not None or \
                self.element_CFF_delta_block_entropy    is not None or \
                self.element_CFF_block_entropy          is not None or \
                self.element_CFF_map                                or \
                self.element_CFF is not False:
            cff_data = True
            try:
                print '    **** Loading state data'
                state_data = cPickle.load(open('%s/state-data_WL-%i.pkl'%(self.cache_dir, max((self.element_CFF_E_Hmu, self.element_CFF_delta_block_entropy, self.element_CFF_block_entropy,15))), 'rb'))
                scan_events = False
            except:
                print '        **** No state data saved. Will calculate.'
                state_data = ElementStateData(self)
                state_data.max_word_length = max((self.element_CFF_E_Hmu, self.element_CFF_delta_block_entropy, self.element_CFF_block_entropy, 15))
                scan_events = True
        else:
            state_data = None
        
        if self.rupture_map:
            scan_events = True
            min_magnitude = py_sys.float_info.max
            max_magnitude = 0
            
            for event_data in events[start_event:end_event]:
                mag = event_data[3]
                if self.section_filter is None:
                    if mag < min_magnitude:
                        min_magnitude = mag
                    if mag > max_magnitude:
                        max_magnitude = mag
                else:
                    section_in_event = False
                    event_sweeps = sweeps[event_data[8]:event_data[9]]
                    for event_sweep in event_sweeps:
                        if self.geometry.elements[int(event_sweep[2])].sid in self.section_filter:
                            section_in_event = True
                            break
                    if section_in_event:
                        if mag < min_magnitude:
                            min_magnitude = mag
                        if mag > max_magnitude:
                            max_magnitude = mag
            
            rupture_map_plotter = RuptureMapPlotter(self)
            rupture_map_plotter.setMinMaxMagnitude(min_magnitude, max_magnitude)
            
            if self.geometry.layered_sections is None:
                self.layerSections()
            
            if self.section_filter is None:
                layered_sections = self.geometry.layered_sections
            else:
                layered_sections = {}
                for sid in self.section_filter:
                    layered_sections[sid] = self.geometry.layered_sections[sid]
            
            rupture_map_plotter.setLayeredSections(layered_sections)
            
            rupture_map_plotter.start_year      = events[start_event][1]
            rupture_map_plotter.end_year        = events[end_event][1]
            rupture_map_plotter.sections        = self.geometry.sectionsSortedByID(self.section_filter)
            rupture_map_plotter.num_of_events   = len(events[start_event:end_event])
            rupture_map_plotter.name            = self.name
            rupture_map_plotter.output_format   = self.output_format
            rupture_map_plotter.out_file        = '%sspace-time.%s'%(file_name_prepend,self.output_format)
            
        if self.event_animation:
            scan_events = True
            event_animation_plotter = EventAnimationPlotter(self)
            event_animation_plotter.number_of_years = events[end_event][1] - events[start_event][1]
            event_animation_plotter.name            = self.name
            event_animation_plotter.out_file        = '%sevents-animation.kml'%(file_name_prepend_dat)
            
        if self.frequency_magnitude is not False:
            scan_events = True
            frequency_magnitude_plotter = FrequencyMagnitudePlotter(self)
            frequency_magnitude_plotter.number_of_years = events[end_event][1] - events[start_event][1]
            frequency_magnitude_plotter.name            = self.name
            if self.frequency_magnitude == 'both' or self.frequency_magnitude == 'plot':
                frequency_magnitude_plotter.output_format   = self.output_format
                frequency_magnitude_plotter.out_file        = '%sfrequency-magnitude.%s'%(file_name_prepend,self.output_format)
            if self.frequency_magnitude == 'both' or self.frequency_magnitude == 'text':
                frequency_magnitude_plotter.out_text_file = '%sfrequency-magnitude.dat'%(file_name_prepend_dat)
                    
        if self.magnitude_rupture_area is not False:
            scan_events = True
            magnitude_rupture_area_plotter = MagnitudeRuptureAreaPlotter(self)
            magnitude_rupture_area_plotter.name            = self.name
            if self.magnitude_rupture_area == 'both' or self.magnitude_rupture_area == 'plot':
                magnitude_rupture_area_plotter.output_format   = self.output_format
                magnitude_rupture_area_plotter.out_file        = '%smagnitude-rupture-area.%s'%(file_name_prepend,self.output_format)
            if self.magnitude_rupture_area == 'both' or self.magnitude_rupture_area == 'text':
                magnitude_rupture_area_plotter.out_text_file = '%smagnitude-rupture-area.dat'%(file_name_prepend_dat)
        
        if self.magnitude_average_slip is not False:
            scan_events = True
            magnitude_average_slip_plotter = MagnitudeAverageSlipPlotter(self)
            magnitude_average_slip_plotter.name            = self.name
            if self.magnitude_average_slip == 'both' or self.magnitude_average_slip == 'plot':
                magnitude_average_slip_plotter.output_format   = self.output_format
                magnitude_average_slip_plotter.out_file        = '%smagnitude-average-slip.%s'%(file_name_prepend,self.output_format)
            if self.magnitude_average_slip == 'both' or self.magnitude_average_slip == 'text':
                magnitude_average_slip_plotter.out_text_file = '%smagnitude-average-slip.dat'%(file_name_prepend_dat)

        if self.average_slip_rupture_length is not False:
            scan_events = True
            average_slip_rupture_length_plotter = AverageSlipRuptureLengthPlotter(self)
            average_slip_rupture_length_plotter.name            = self.name
            if self.average_slip_rupture_length == 'both' or self.average_slip_rupture_length == 'plot':
                average_slip_rupture_length_plotter.output_format   = self.output_format
                average_slip_rupture_length_plotter.out_file        = '%saverage-slip-rupture-length.%s'%(file_name_prepend,self.output_format)
            if self.average_slip_rupture_length == 'both' or self.average_slip_rupture_length == 'text':
                average_slip_rupture_length_plotter.out_text_file = '%saverage-slip-rupture-length.dat'%(file_name_prepend_dat)
    
        if self.slip_in_slip_out is not False:
            scan_events = True
            slip_in_slip_out_plotter = SlipInSlipOutPlotter(self)
            slip_in_slip_out_plotter.name            = self.name
            if start_event == 0:
                slip_in_slip_out_plotter.duration = self.geometry.converter.year_sec(events[end_event][1])
            else:
                slip_in_slip_out_plotter.duration = self.geometry.converter.year_sec(events[end_event][1] - events[start_event][1])
            slip_in_slip_out_plotter.out_text_file = '%sslip_in_slip_out.dat'%(file_name_prepend_dat)

        if self.recurrence_intervals is not False:
            scan_events = True
            recurrence_interval_plotter         = RecurrenceIntervalPlotter(self)
            recurrence_interval_plotter.name    = self.name
            if self.recurrence_intervals == 'both' or self.recurrence_intervals == 'plot':
                recurrence_interval_plotter.output_format       = self.output_format
                recurrence_interval_plotter.out_file_prepend    = file_name_prepend
            if self.recurrence_intervals == 'both' or self.recurrence_intervals == 'text':
                recurrence_interval_plotter.out_text_file_prepend = file_name_prepend_dat
        
        if self.export_eqsim_events:
            scan_events = True
            eqsim_event_exporter                    = EqSimEventExporter(self)
            eqsim_event_exporter.min_slip_map_mag   = self.export_eqsim_events_min_slip_map_mag
            eqsim_event_exporter.out_file           = '%sevents_slip-map-%.1f.dat'%(file_name_prepend_eqsim,self.export_eqsim_events_min_slip_map_mag)
            try:
                eqsim_event_exporter.setEidRemap(open('%s%s_element-remap.dat'%(self.root_dir, self.name),'r'))
            except:
                eqsim_event_exporter.eid_remap = None
        
        if self.element_CFF is not False:
            element_CFF_plotter                      = ElementCFFPlotter(self)
            element_CFF_plotter.name                 = self.name
            if self.element_CFF == 'both' or self.element_CFF == 'plot':
                element_CFF_plotter.output_format    = self.output_format
                element_CFF_plotter.out_file         = '%selement-CFF.%s'%(file_name_prepend,self.output_format)
            if self.element_CFF == 'both' or self.element_CFF == 'text':
                element_CFF_plotter.out_text_file    = '%selement-CFF.dat'%(file_name_prepend_dat)
        
        if self.element_CFF_map:
            element_CFF_map_plotter                 = ElementCFFMapPlotter(self)
            element_CFF_map_plotter.output_format   = self.output_format
            element_CFF_map_plotter.out_file        = '%selement-CFF-map.%s'%(file_name_prepend,self.output_format)
        
        if self.element_CFF_block_entropy is not None:
            element_CFF_block_entropy_plotter                   = ElementCFFBlockEntropyPlotter(self)
            element_CFF_block_entropy_plotter.output_format     = self.output_format
            element_CFF_block_entropy_plotter.max_word_length   = self.element_CFF_block_entropy
            element_CFF_block_entropy_plotter.type              = 'be'
            element_CFF_block_entropy_plotter.out_file          = '%selement-CFF-block-entropy.%s'%(file_name_prepend,self.output_format)
        
        if self.element_CFF_delta_block_entropy is not None:
            element_CFF_delta_block_entropy_plotter                 = ElementCFFBlockEntropyPlotter(self)
            element_CFF_delta_block_entropy_plotter.output_format   = self.output_format
            element_CFF_delta_block_entropy_plotter.max_word_length = self.element_CFF_delta_block_entropy
            element_CFF_delta_block_entropy_plotter.type            = 'dbe'
            element_CFF_delta_block_entropy_plotter.out_file        = '%selement-CFF-delta-block-entropy.%s'%(file_name_prepend,self.output_format)
        
        if self.element_CFF_E_Hmu is not None:
            element_CFF_E_Hmu_plotter                 = ElementCFFEHmuPlotter(self)
            element_CFF_E_Hmu_plotter.output_format   = self.output_format
            element_CFF_E_Hmu_plotter.max_word_length = self.element_CFF_E_Hmu
            #element_CFF_E_Hmu_plotter.type            = 'dbe'
            element_CFF_E_Hmu_plotter.out_file        = '%selement-E-Hmu.%s'%(file_name_prepend,self.output_format)
        
        if self.print_event_info:
            event_info_printer                 = EventInfoPrinter(self)
        
        if scan_events and cff_data:
            print '    **** Scanning events from cff data file %s_CFF.dat'%self.name
            
            cff_data_raw = np.fromfile('%s_CFF.dat'%self.name,sep=' ')
            
            #print cff_data.shape, (cff_data.size/(len(self.geometry.elements)+1), len(self.geometry.elements)+1)
            
            cff_data_raw = np.reshape(cff_data_raw, (cff_data_raw.size/(len(self.geometry.elements)+2), len(self.geometry.elements)+2))
            
            state_data.cff_data_raw = cff_data_raw
            state_data.prepareData()
            
            print '        **** Saving state data'
            cPickle.dump(state_data, open('%s/state-data_WL-%i.pkl'%(self.cache_dir, max((self.element_CFF_E_Hmu, self.element_CFF_delta_block_entropy, self.element_CFF_block_entropy, 15))), 'wb'))
        elif scan_events:
            print '    **** Scanning events'
            for event_num, event_data in enumerate(events[start_event:end_event]):
                if ( event_num%5000 == 0 ):
                    py_sys.stdout.write('\r')
                    py_sys.stdout.flush()
                    py_sys.stdout.write('        event:%i of %i'%(event_num, end_event - start_event))
                    py_sys.stdout.flush()
                if state_data is not None:
                    state_data.addEvent(event_data, sweeps[event_data[8]:event_data[9]])
                    
                if self.export_eqsim_events:
                    eqsim_event_exporter.addEvent(event_data, sweeps[event_data[8]:event_data[9]])
                if self.rupture_map:
                    rupture_map_plotter.addEvent(event_data, sweeps[event_data[8]:event_data[9]])
                if self.event_animation:
                    event_animation_plotter.addEvent(event_data, sweeps[event_data[8]:event_data[9]])
                if self.frequency_magnitude is not False:
                    if self.section_filter is not None:
                        if self.geometry.elements[int(event_data[2])].sid in self.section_filter:
                            frequency_magnitude_plotter.addEvent(event_data[3])
                    else:
                        frequency_magnitude_plotter.addEvent(event_data[3])
                if self.print_event_info:
                    if self.section_filter is not None:
                        if self.geometry.elements[int(event_data[2])].sid in self.section_filter:
                            event_info_printer.addEvent(event_data)
                    else:
                        event_info_printer.addEvent(event_data)
                    
                analyze_sweeps = False
                
                if self.magnitude_rupture_area is not False or self.magnitude_average_slip is not False or self.average_slip_rupture_length is not False or self.slip_in_slip_out is not False or self.recurrence_intervals is not False:
                    if self.section_filter is not None:
                        if self.geometry.elements[int(event_data[2])].sid in self.section_filter:
                            analyze_sweeps = True
                    else:
                        analyze_sweeps = True
                
                if analyze_sweeps:
                    ruptured_elements = {}
                    
                    if self.magnitude_rupture_area is not False:
                        rupture_area = 0.0
                    
                    if self.magnitude_average_slip is not False or self.average_slip_rupture_length is not False:
                        total_slip = 0.0
                    
                    if self.average_slip_rupture_length is not False:
                        rupture_length = 0.0
                    
                    for sweep in sweeps[event_data[8]:event_data[9]]:
                    #for sweepid in range(event_data[8],event_data[9]):
                        #sweep = sweeps[sweepid]
                        #print sweep[1]
                        eid = sweep[2]
                        if self.magnitude_average_slip is not False or self.average_slip_rupture_length is not False:
                            total_slip += sweep[3]
                        if self.slip_in_slip_out is not False:
                            slip_in_slip_out_plotter.addSlip(eid, sweep[3])
                        try:
                            tmp = ruptured_elements[eid]
                        except KeyError:
                            if self.magnitude_rupture_area is not False:
                                rupture_area += self.geometry.elements[eid].sideLength() ** 2.0
                            if self.average_slip_rupture_length is not False:
                                if self.geometry.elements[eid].onTrace():
                                    rupture_length += self.geometry.elements[eid].sideLength()
                            ruptured_elements[eid] = True
                    
                    if self.magnitude_average_slip is not False or self.average_slip_rupture_length is not False:
                        average_slip = total_slip/float(len(ruptured_elements))
                        
                    if self.magnitude_rupture_area is not False:
                        magnitude_rupture_area_plotter.addEvent(event_data[3], rupture_area * 1.0e-6)
                    if self.magnitude_average_slip is not False:
                        magnitude_average_slip_plotter.addEvent(event_data[3], average_slip)
                    if self.average_slip_rupture_length is not False and rupture_length != 0.0:
                        average_slip_rupture_length_plotter.addEvent(average_slip, rupture_length * 1.0e-3)
                    if self.recurrence_intervals is not False and event_data[3] >= 6.5:
                        involved_sections = {}
                        for eid in ruptured_elements.keys():
                            sid = self.geometry.elements[eid].sid
                            involved_sections[sid] = sid
                        recurrence_interval_plotter.addEvent(event_data[3], event_data[1], involved_sections.keys())
                        
    
        event_count = {}
        print
        print '    **** Begining plot'
        if self.rupture_map:
            rupture_map_plotter.plot()
            event_count['st'] = rupture_map_plotter.num_of_events
        
        if self.event_animation:
            event_animation_plotter.plot()
        
        if self.frequency_magnitude is not False:
            frequency_magnitude_plotter.plot()
            event_count['fm'] = frequency_magnitude_plotter.event_count
        
        if self.magnitude_rupture_area is not False:
            magnitude_rupture_area_plotter.plot()
            event_count['mra'] = magnitude_rupture_area_plotter.event_count
        
        if self.magnitude_average_slip is not False:
            magnitude_average_slip_plotter.plot()
            event_count['mas'] = magnitude_average_slip_plotter.event_count
        
        if self.average_slip_rupture_length is not False:
            average_slip_rupture_length_plotter.plot()
            event_count['asrl'] = average_slip_rupture_length_plotter.event_count

        if self.slip_in_slip_out is not False:
            slip_in_slip_out_plotter.plot()
        
        if self.recurrence_intervals is not False:
            recurrence_interval_plotter.plot()
    
        if self.export_eqsim_events:
            eqsim_event_exporter.write()
        '''
        if state_data is not None:
            if self.element_CFF_E_Hmu                       is not None or \
                    self.element_CFF_delta_block_entropy    is not None or \
                    self.element_CFF_block_entropy          is not None:
                if state_data.element_data is None:
                    state_data.do_IT_calc = True
            
            if self.element_CFF_map:
                if state_data.element_maps is None:
                    state_data.do_map = True
            
            state_data.analyzeOrbit()
        
        if state_data is not None:
            if state_data.data_changed:
                state_data.data_changed = False
                print '        **** Saving state data'
                cPickle.dump(state_data, open('%sState-Data_WL-%i.pkl'%(file_name_prepend, max((self.element_CFF_E_Hmu, self.element_CFF_delta_block_entropy, self.element_CFF_block_entropy, 10))), 'wb'))
        '''
        if self.element_CFF is not False:
            element_CFF_plotter.plot(state_data, start_event, end_event)
        
        if self.element_CFF_map:
            element_CFF_map_plotter.plot(state_data)
        
        if self.element_CFF_block_entropy is not None:
            element_CFF_block_entropy_plotter.plot(state_data)
        
        if self.element_CFF_delta_block_entropy is not None:
            element_CFF_delta_block_entropy_plotter.plot(state_data)
        
        if self.element_CFF_E_Hmu is not None:
            element_CFF_E_Hmu_plotter.plot(state_data)
            
        if self.print_event_info:
            event_info_printer.eprint()
        
        f.close()
        print '*** Done ***'
        
        return event_count
        
    def layerSections(self, section_list = None):
        if section_list is None:
            all_sections = []
            for section in self.geometry.sectionsSortedByID():
                all_sections.append(section.sid)
        else:
            all_sections = section_list
        
        #print multiprocessing.cpu_count()
        
        num_processes = multiprocessing.cpu_count()
        
        print '*** Layering sections. Threads: %s ***'%num_processes
        
        seg = int(round(float(len(all_sections))/float(num_processes)))
        
        sections = []
        
        for i in range(num_processes):
            if i == num_processes - 1:
                end_index = len(all_sections)
            else:
                end_index = seg*int(i + 1)
            sections.append((all_sections[ int(i) * seg : end_index ]))
            
        work_queue = multiprocessing.Queue()
        for job in sections:
            work_queue.put(job)
        
        # create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()
     
        # spawn workers
        for i in range(len(sections)):
            worker = SectionLayerer(self, work_queue, result_queue)
            worker.start()
     
        # collect the results off the queue
        results = []
        for i in range(len(sections)):
            results.append(result_queue.get())
        
        self._layered_sections = {}
        
        for result in results:
            for key, value in result.iteritems():
                self._layered_sections[key] = value
        

        print '*** Done ***'
    
    '''
    def saveRuptureMapData(self, event_file):
        print '*** Saving Rupture Map Data ***'
        print '    Event file %s'%event_file
        
        #print self.geometry.layered_sections[1].shape
        
        rupture_map_events = []
        section_mag_extrema = {}
                
        f = h5py.File(event_file, 'r')
        
        events = f['event_table']
        sweeps = f['event_sweep_table']
        
        #x = scipy.sparse.lil_matrix( (10,N) )
        section_sparse = None
        
        end_event = len(events)-1
        
        for event_num, event_data in enumerate(events):
            if ( event_num%5000 == 0 ):
                print '        event:%i of %i'%(event_num, end_event)
            event_magnitude = event_data[3]
            event_year = event_data[1]
            
            #print event_magnitude
            
            elements_in_event_dict = {}
            
            
            for sweep in sweeps[event_data[8]:event_data[9]]:
                elements_in_event_dict[int(sweep[2])] = int(sweep[2])
            
            elements_in_event = elements_in_event_dict.keys()
            
            event_sections = {}
            #event_data = {'event_number':event_num, 'event_magnitude':event_magnitude}
            #(r,g,b,a) = self.cmap(self.magnitude_slope * (event_magnitude - self.min_magnitude))
            #section_x_offset = 0
            #the_ys = {}
            for sid in sorted(self.geometry.layered_sections.iterkeys()):
            #for sindex, s in enumerate(sections):
                #print self.geometry.layered_sections[sid].shape
                it = np.nditer(self.geometry.layered_sections[sid], flags=['multi_index'])
                while not it.finished:
                    #print it.multi_index, it[0]
                    #try:
                        #eid = elements_in_event[int(it[0])]
                    if int(it[0]) in elements_in_event:
                        if section_sparse is not None:
                            section_sparse[it.multi_index[0], it.multi_index[1]] = True
                        else:
                            section_sparse = scipy.sparse.lil_matrix(self.geometry.layered_sections[sid].shape, dtype=np.dtype(bool))
                            section_sparse[it.multi_index[0], it.multi_index[1]] = True
                        
                        #the_x_start = float(it.multi_index[1] + section_x_offset)/float(self.max_x)
                        #the_x_end = float(it.multi_index[1] + section_x_offset + 1)/float(self.max_x)
                        #the_y = self.curr_y - it.multi_index[0]
                    it.iternext()
                if section_sparse is not None:
                    event_sections[sid] = section_sparse
                    section_sparse = None
                    try:
                        extrema = section_mag_extrema[sid]
                    except KeyError:
                        extrema = [py_sys.float_info.max, 0]
                    
                    if event_magnitude < extrema[0]:
                        extrema[0] = event_magnitude
                    if event_magnitude > extrema[1]:
                        extrema[1] = event_magnitude
                        
                    section_mag_extrema[sid] = extrema
            
            rupture_map_events.append({'sections':event_sections, 'event_number':event_num, 'event_magnitude':event_magnitude, 'event_year':event_year})
            
            #for item in event_dict.iteritems():
            #    print item[0], item[1]
            #print
        cPickle.dump({'section_mag_extrema':section_mag_extrema,'events':rupture_map_events}, open('%s_Rupture-Map-Data.pkl'%self.name, 'wb'))
        
        f.close()
        print '*** Done ***'
    '''
    
    def sectionWithElement(self, element_id):
        return self.geometry.sectionWithElement(element_id)
        
    def geometryFromTraceFile(self, trace_file):
        print '*** Adding geometry from trace file %s ***'%(trace_file)
                    
        f = open(trace_file, 'r')
        self.geometry.geometryFromTraceFile(f)
        f.close()
        print '*** Done ***'
    
    def proximityTable(self):
        print '*** Building proximity table ***'
        self.geometry.proximityTable()
        print '*** Done ***'
        
    def totalEvents(self):
        events = self.events.eventsSortedByYear()
        return events[-1].eid
    
    def basicPlot(self, x, y, x_label = None, y_label = None, title = None, out_file = None, axis_scale = None, extra_x = None, extra_y = None, extra_title = None, connect_points = False, extra_x_error = None, extra_y_error = None):
        fig = plt.figure(figsize=[11, 7.33723],dpi=100)
        the_ax = fig.add_axes([0.06, 0.08, 0.9, 0.8])
        
        if connect_points:
            ls1 = '--'
        else:
            ls1 = 'None'
        lw1 = 0.5
        c1 = (0.5,0.5,0.5,0.5)
        mfc1 = (0,0,0,1)
        ms1 = 3
        marker1 = 'o'
        
        ls_extra = '-'
        lw_extra = 0.5
        c_extra = (0.0,0.0,1,1)
        
        if axis_scale is None:
            the_ax.plot(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                the_ax.plot(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                
        elif axis_scale is 'semilogy':
            the_ax.semilogy(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                #print extra_x_error, extra_y_error
                if extra_x_error is not None or extra_y_error is not None:
                    the_ax.errorbar(extra_x, extra_y, yerr = extra_y_error, xerr = extra_x_error, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                else:
                    the_ax.semilogy(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                
        elif axis_scale is 'semilogx':
            the_ax.semilogx(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                the_ax.semilogx(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
        
        elif axis_scale is 'loglog':
            the_ax.loglog(x, y, ls = ls1, lw = lw1, c = c1, mfc = mfc1, ms = ms1, marker = marker1)
            
            if extra_x is not None and extra_y is not None:
                the_ax.loglog(extra_x, extra_y, ls = ls_extra, lw = lw_extra, c = c_extra, label = extra_title)
                        
        the_ax.set_ylabel(y_label,size=9)
        the_ax.set_xlabel(x_label,size=9)
        
        the_ax.tick_params(labelsize=9)
        
        the_ax.margins(0.05,0.3)
        
        the_ax.set_title(title, position=(0.0,1.05), ha='left', size=12)
        
        the_ax.legend()
           
        if out_file is not None:
            plt.savefig(out_file, format='pdf')

    def loadHdf5Events(self, event_file):
        print '*** Creating event cache from file %s ***'%(event_file)
        #self.events = VCEvents()
        
        self.events.loadHdf5Events(event_file)

        print '*** Done ***'
    
    def sectionsSortedByID(self):
        return self.geometry.sectionsSortedByID()
    
    def sectionsInArea(self, NE_pt, SW_pt):
        #NE_pt is the LatLonPoint of the north east corner of the area
        #SW_pt is the LatLonPoint of the south west corner of the area
        return self.geometry.sectionsInArea(NE_pt, SW_pt)
        
    def setLayeredSections(self, layered_sections):
        self.geometry.layered_sections = layered_sections
    
    def saveLayeredSections(self, out_file):
        cPickle.dump(self.geometry.layered_sections, open(out_file, 'wb'))
            
    def loadEqsimEvents(self, event_file, start_event = 0, end_event = None, event_unit = 'event_num'):
        pass
    
    def loadHdf5Geometry(self, geometry_file, section_list = None):
        print '*** Creating geometry cache from file %s ***'%(geometry_file)
        #if section_list is not None:
        #    print '    sections: ',
        #    print section_list
            
        #self.geometry = VCGeometry()
        f = h5py.File(self.data_file, 'r')
        self.geometry.loadHdf5Geometry(f, section_list)
        f.close()
        print '*** Done ***'
    
    '''
    def eventYearMappingHdf5(self, event_file):
        print '*** Creating event-year map from hdf5 file %s ***'%(event_file)
        
        f = h5py.File(event_file, 'r')
        events = f['event_table']
        number_of_events = len(events)
        self.events.years = np.zeros(number_of_events)
        self.events.event_ids = np.zeros(number_of_events,dtype=int)
        
        for ind, event in enumerate(events):
            self.events.years[ind] = event[1]
            self.events.event_ids[ind] = event[0]
            
            if ( ind%10000 == 0 ):
                print '    event:%i of %i'%(ind, number_of_events)
        
        f.close()
        print '*** Done ***'
    '''
    
    def eventForYear(self, year):
        return self.events.eventForYear(year)
    
    def loadEqsimGeometry(self, section_list = None):
        print '*** Loading Eqsim geometry from file %s ***'%(geometry_file)
        if section_list is not None:
            print '    sections: ',
            print section_list
            
        #self.geometry = VCGeometry()
        f = open(geometry_file, 'r')
        self.geometry.loadEqsimGeometry(f, section_list)
        f.close()
        print '*** Done ***'
    
    def loadGreens(self, greens_file, section_list = None):
        pass
    
    def exportHdf5Geometry(self, out_geometry_file, out_friction_file, section_list = None):
        pass
        
    def exportEqsimGeometry(self, geometry_file, friction_file, section_list = None):
        print '*** Exporting EqSim geometry to file %s and %s ***'%(geometry_file, friction_file)
        self.geometry.exportEqsimGeometry(geometry_file, friction_file)
        print '*** Done ***'
    
    def exportKMLGeometry(self, out_file, section_list = None):
        print '*** Exporting KML geometry to file %s ***'%(out_file)
        min_lat, max_lat, min_lon, max_lon = self.geometry.minMaxLatLon()
        #print 'test'
        #print min_lat, max_lat, min_lon, max_lon
        
        dh = self.geometry.converter.distanceOnEarthsSurface(max_lat - (max_lat - min_lat)/2, min_lon, max_lat - (max_lat - min_lat)/2, max_lon)
        
        dv = self.geometry.converter.distanceOnEarthsSurface(min_lat, max_lon - (max_lon - min_lon)/2, max_lat, max_lon - (max_lon - min_lon)/2)
        
        #print dv, dh
        
        range = (max((dv,dh))/2.0) * (1/math.tan(self.geometry.converter.deg_rad(30.0)))
        
        #print range
        f = open(out_file, 'w')
        # write header
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        f.write('<Document>\n')
        f.write('<LookAt>\n')
        f.write('   <longitude>%f</longitude>\n'%(max_lon - (max_lon - min_lon)/2))
        f.write('   <latitude>%f</latitude>\n'%(max_lat - (max_lat - min_lat)/2))
        f.write('   <altitude>0</altitude>\n')
        f.write('   <range>%f</range>\n'%(range + range/10.0))
        f.write('   <tilt>0</tilt>\n')   
        f.write('   <heading>0</heading>\n')   
        f.write('   <altitudeMode>absolute</altitudeMode>\n') 
        f.write('</LookAt>\n')
        f.write('<name>%s</name>\n'%self.name)
        f.write(self.geometry.KMLStyle_section())
        f.write(self.geometry.KMLStyle_base())
        f.write('<Folder id="fault_names">\n')
        f.write('   <name>Fault Section Names</name>\n')
        for s in self.geometry.sections.itervalues():
            s.exportKMLSectionLabel(f)
        f.write('</Folder>\n')
        
        f.write('<Folder id="faults">\n')
        f.write('   <name>Fault Sections</name>\n')
        for s in self.geometry.sections.itervalues():
            s.exportKMLGeometry(f)
        f.write('</Folder>\n')
        f.write('</Document>\n')
        f.write('</kml>')
        f.close()
        print '*** Done ***'
    
    def exportHdf5Events(self, out_file, start_event = 0, end_event = None, event_unit = 'event_num'):
        pass
    
    def exportEqsimEvents(self, out_file, min_slip_map_mag = 5.0):
        print '*** Exporting EqSim events with minimum slip map magnitude %.2f to file %s ***'%(min_slip_map_mag, out_file)
        self.events.exportEqsimEvents(out_file, min_slip_map_mag)
        print '*** Done ***'

    '''
    def ruptureMap(self, out_file = None, section_list = None, min_mag = None):
        
        print '*** Printing rupture map ***'
        
        print '    preparing'
        
        try:
            events_by_year = self.events.eventsSortedByYear()
            events_by_magnitude = self.events.eventsSortedByMagnitude()
        except:
            print '!!! Events have not been loaded !!!'
            return
        
        
        if min_mag is None:
            min_magnitude = round(events_by_magnitude[0].magnitude)
        else:
            min_magnitude = min_mag
        max_magnitude = round(events_by_magnitude[-1].magnitude)
        
        magnitude_slope = 1.0/(max_magnitude - min_magnitude)
        
        cmap = plt.get_cmap('autumn_r')
        
        end_year = events_by_year[-1].year
        start_year = events_by_year[0].year
        
        print '    %i events from year %f to year %f'%(len(events_by_year), start_year, end_year)
        
        
        
        if self.geometry.layered_sections is None:
            try:
                self.setLayeredSections(cPickle.load(open('%s_layered-sections.pkl'%self.name, 'rb')))
            except:
                self.layerSections()
                self.saveLayeredSections('%s_layered-sections.pkl'%vc_sys.name)
            
            #self.layerSections(section_list)
        
        layered_sections = self.geometry.getLayeredSections(section_list)
        
        max_y_offset = 0
        max_x = 0
        
        for s in layered_sections.itervalues():
            if s.shape[0] > max_y_offset:
                max_y_offset = s.shape[0]
            max_x += s.shape[1]
        
        fig = plt.figure(figsize=[11, 7.33723],dpi=100)
        the_ax = fig.add_axes([0.08, 0.03, 0.89, 0.7])
        
        ax1 = fig.add_axes([.77, 0.98, 0.2, 0.01])
        
        norm = mcolor.Normalize(vmin=min_magnitude, vmax=max_magnitude)
        
        cb1 = mcolorbar.ColorbarBase(ax1, cmap=cmap,
                   norm=norm,
                   orientation='horizontal')
        
        cb1.set_ticks([min_magnitude,(min_magnitude+ max_magnitude)/2.0, max_magnitude])

        for label in ax1.xaxis.get_ticklabels():
            label.set_fontsize(6)
            
        for line in ax1.xaxis.get_ticklines():
            line.set_alpha(0)
            
        cb1.ax.set_xlabel('Magnitude',size=6,position=(1,0), ha='right')
        
        years = []
        curr_y = max_y_offset
        for n, e in enumerate(events_by_year):
            if e.magnitude >= min_magnitude:
                if len(events_by_year) < 200:
                    the_ax.axhline(y=curr_y, alpha=0.2, linewidth=1, color='0.5')
                years.append((e.year, e.eid, curr_y))
                
                (r,g,b,a) = cmap(magnitude_slope * (e.magnitude - min_magnitude))
                
                if ( n%500 == 0 ):
                    print '    event:%i of %i'%(n, len(events_by_year))
                
                section_x_offset = 0
                for sid in sorted(layered_sections.iterkeys()):
                #for sindex, s in enumerate(sections):
                    it = np.nditer(layered_sections[sid], flags=['multi_index'])
                    while not it.finished:
                        if e.elementInEvent(it[0]):
                            the_y = curr_y - it.multi_index[0]
                            the_x_start = float(it.multi_index[1] + section_x_offset)/float(max_x)
                            the_x_end = float(it.multi_index[1] + section_x_offset + 1)/float(max_x)
                            the_ax.axhline(y=the_y, xmin=the_x_start, xmax=the_x_end, color=(r,g,b))
                        it.iternext()
                    
                    section_x_offset += layered_sections[sid].shape[1]
                
                
                curr_y += max_y_offset
        
        the_ax.autoscale(enable=True, axis='both', tight=True)
        
        label_lines = []
        section_x_start = 0
        section_x_finish = 0
        print '    finishing'
        #for ns, s in enumerate(sections):
        for ns, sid in enumerate(sorted(layered_sections.iterkeys())):
            layered_section = layered_sections[sid]
            section_x_finish += layered_sections[sid].shape[1]
            
            x_text_loc = (float(ns)/float(len(layered_sections)))
            x_fault_loc = float(section_x_finish + section_x_start)/(2.0*float(max_x))
            
            section_x_start += layered_sections[sid].shape[1]
            
            y_text_distance = .08
            
            the_ax.axvline(x=section_x_finish, alpha=0.2, linewidth=1, color='0.5')
            the_ax.text(x_text_loc, 1+y_text_distance, self.geometry.sections[sid].sname,
                horizontalalignment='center',
                verticalalignment='bottom',
                rotation=90,
                size=6,
                transform=the_ax.transAxes)
            
            line_xs = [x_text_loc, x_text_loc, x_fault_loc, x_fault_loc]
            line_ys = [1+y_text_distance-0.005, 1+y_text_distance*(3.0/4.0), 1 + y_text_distance/4.0, 1]
            label_lines.append(mlines.Line2D(
                line_xs, line_ys,transform=the_ax.transAxes, solid_capstyle='round', alpha=0.2, linewidth=1, color='0.2'))
        
        y_locs = []
        y_labels = []
        
        for n, year in enumerate(years):
            if len(years) > 20:
                if (n%int(float(len(years))/10) == 0):
                    y_locs.append(year[2])
                    y_labels.append('%.1f [%i]'%(year[0],year[1]))
            else:
                y_locs.append(year[2])
                y_labels.append('%.1f [%i]'%(year[0],year[1]))
            
        the_ax.yaxis.set_ticks(y_locs)
        the_ax.yaxis.set_ticklabels(y_labels)
        
        x_locs = the_ax.yaxis.get_majorticklocs()
        
        x_tick_labels = []
        
        for x_loc in x_locs:
            x_tick_labels.append('')
        
        the_ax.xaxis.set_ticklabels(x_tick_labels)
        
        the_ax.autoscale(enable=True, axis='both', tight=True)
        
        the_ax.lines.extend(label_lines)
        
        
        the_ax.set_ylabel('Year [event id]',size=9)
        
        the_ax.tick_params(labelsize=6)
        
        if out_file is not None:
            the_ax.set_title('%s'%self.name, position=(-0.05,1.34), ha='left', size=12)
            plt.savefig(out_file, format='pdf')
        
        print '*** Done ***'
    '''
    '''
    def frequencyMagnitude(self, out_file = None, section_list = None):
        print '*** Printing Frequency-Magnitude plot to file %s ***'%(out_file)
        
        try:
            events_by_magnitude = self.events.eventsSortedByMagnitude(order = 'r', section_list = section_list)
            events_by_year = self.events.eventsSortedByYear()
            end_year = events_by_year[-1].year
            start_year = events_by_year[0].year
        except:
            print '!!! Events have not been loaded !!!'
            return
        
        magnitudes = np.asarray([e.magnitude for e in events_by_magnitude])
        
        cum_freq = {}
        total_events = len(magnitudes)
        
        for num, mag in enumerate(reversed(magnitudes)):
            cum_freq['%.10f'%mag] = total_events - (num + 1)
            #print '%.10f'%mag, cum_freq['%.10f'%mag]
        
        #print
        
        x = []
        y = []
        
        extra_x = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        extra_y = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
        #extra_y_error_plus = [1.50, 0.43, 0.28, 0.11, 0.06, 0.035]       
        extra_y_error = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        
        number_of_years = end_year - start_year
        
        # for key in sorted(mydict.iterkeys()):
        for mag in sorted(cum_freq.iterkeys()):
            #print mag, float(cum_freq[mag])/number_of_years
            x.append(mag)
            y.append(float(cum_freq[mag])/number_of_years)
        
        #print end_year, start_year, end_year - start_year
        #return [e for e in sorted(self.events.values(),lambda a,b: cmp(a.magnitude,b.magnitude))]
        
        ''
        magnitudes_mod = magnitudes[0] - magnitudes
        
        cumfreqs, lowlim, binsize, extrapoints = stats.cumfreq(magnitudes_mod, numbins=100)
        
        x = []
        y = []
        
        extra_x = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        extra_y = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
        #extra_y_error_plus = [1.50, 0.43, 0.28, 0.11, 0.06, 0.035]       
        extra_y_error = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        
        for bin,freq in enumerate(cumfreqs):
            x.append(max(magnitudes) + abs(lowlim) - binsize*(bin+1))
            y.append(freq/number_of_years)
            #print max(magnitudes) + abs(lowlim) - binsize*(bin+1), freq/number_of_years
        ''
        if section_list is None:
            title = '%s Frequency-Magnitude'%self.name
        else:
            title = '%s %s Frequency-Magnitude'%(self.name, section_list)
            
        self.basicPlot(x, y, x_label='Magnitude', y_label='Log(Cumulative Number of Events)', title=title, out_file = out_file, axis_scale = 'semilogy', connect_points = True, extra_x = extra_x, extra_y = extra_y, extra_y_error = extra_y_error, extra_title = 'UCERF2 Observed 95% Confidence')
        
        print '*** Done ***'
    '''
    '''
    def magnitudeRuptureArea(self, out_file = None, section_list = None):
        print '*** Printing Magnitude-Rupture Area plot to file %s ***'%(out_file)
        
        try:
            events_by_magnitude = self.events.eventsSortedByMagnitude(order = 'r', section_list = section_list)
        except:
            print '!!! Events have not been loaded !!!'
            return
        
        x = []
        y = []
        
        for event in events_by_magnitude:
            x.append(event.ruptureArea() * 1.0e-6 )
            y.append(event.magnitude)
            #print event.ruptureArea() * 1.0e-6,event.magnitude
        
        wc = 4.07 + 0.98 * np.log10(x)
        
        if section_list is None:
            title = '%s Magnitude-Rupture_Area'%self.name
        else:
            title = '%s %s Magnitude-Rupture_Area'%(self.name, section_list)
            
        self.basicPlot(x, y, x_label=r'log(Rupture Area [km$^\mathsf{2}$])', y_label='Magnitude', title=title, out_file = out_file, axis_scale = 'semilogx', extra_x = x, extra_y = wc, extra_title = 'Wells-Coppersmith')
        
        print '*** Done ***'
    '''
    '''
    def magnitudeAverageSlip(self, out_file = None, section_list = None):
        print '*** Printing Magnitude-Average Slip plot to file %s ***'%(out_file)
        
        try:
            events_by_magnitude = self.events.eventsSortedByMagnitude(order = 'r', section_list = section_list)
        except:
            print '!!! Events have not been loaded !!!'
            return
        
        x = []
        y = []
        
        for event in events_by_magnitude:
            x.append(event.averageSlip() )
            y.append(event.magnitude)
            #print event.averageSlip(),event.magnitude
        
        wc = 6.93 + 0.82 * np.log10(x)
        
        if section_list is None:
            title = '%s Magnitude-Average_Slip'%self.name
        else:
            title = '%s %s Magnitude-Average_Slip'%(self.name, section_list)
            
        self.basicPlot(x, y, x_label='log(Average Slip [m])', y_label='Magnitude', title=title, out_file = out_file, axis_scale = 'semilogx', extra_x = x, extra_y = wc, extra_title = 'Wells-Coppersmith')
        
        print '*** Done ***'
    '''
    '''
    def averageSlipRuptureLength(self, out_file = None, section_list = None, start_event = 0, end_event = None, event_unit = 'event_num'):
        print '*** Printing Average Slip-Rupture Length plot to file %s ***'%(out_file)
        
        try:
            events_by_magnitude = self.events.eventsSortedByMagnitude(order = 'r', section_list = section_list)
        except:
            print '!!! Events have not been loaded !!!'
            return
        
        x = []
        y = []
        
        for event in events_by_magnitude:
            rupture_length = event.ruptureLength() * 1.0e-3
            if rupture_length != 0.0:
                x.append(rupture_length)
                y.append(event.averageSlip())
                #print rupture_length,event.averageSlip()
            
        wc = 10 ** (-1.43 + 0.88 * np.log10(x))
        
        if section_list is None:
            title = '%s Average_Slip-Rupture_Length'%self.name
        else:
            title = '%s %s Average_Slip-Rupture_Length'%(self.name, section_list)
            
        self.basicPlot(x, y, x_label='log(Rupture Length [km])', y_label='log(Average Slip [m])', title=title, out_file = out_file, axis_scale = 'loglog', extra_x = x, extra_y = wc, extra_title = 'Wells-Coppersmith')
        
        print '*** Done ***'
    '''

    def eventAnimation(self, out_file, section_list = None, fade_steps = 6, event_number = None):
        print '*** Exporting KML event animation to file %s ***'%(out_file)
        f = open(out_file, 'w')
        
        # write header
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        f.write('<Document>\n')
        if event_number is None:
            f.write('<name>%s event animation</name>\n'%self.name)
        else:
            f.write('<name>%s event animation of event %i</name>\n'%(self.name, event_number))
        f.write(self.geometry.KMLStyle_section())
        f.write(self.geometry.KMLStyle_base())
        f.write(self.geometry.KMLStyle_broke())
        f.write(self.geometry.KMLStyle_elementBreak(fade_steps))
        f.write('<Folder>\n')
        f.write('\t<name>Fault Section Names</name>\n')
        for s in self.geometry.sections.itervalues():
            s.exportKMLSectionLabel(f)
        f.write('</Folder>\n')
        f.write('<Folder>\n')
        f.write('\t<name>Fault Sections</name>\n')
        for s in self.geometry.sections.itervalues():
            s.exportKMLGeometry(f)
        f.write('</Folder>\n')
        f.write('<Folder>\n')
        f.write('\t<name>Events</name>\n')
        if event_number is None:
            events = self.events.eventsSortedByYear()
            e_end = events[-1].eid
            for event in events:
                    if ( event.eid%500 == 0 ):
                            print '    event:%i of %i'%(event.eid, e_end)
                    event.exportKML(f,self,events[0].year, events[-1].year,fade_steps)
        else:
            self.events.events[event_number].exportSweepKML(f, fade_steps)
        f.write('</Folder>\n')
        f.write('</Document>\n')
        f.write('</kml>')
        f.close()
        print '*** Done ***'
    
    def deleteSection(self, sid):
        print '*** Deleting section %s ***'%(sid)
    
        target_eids = self.geometry.sections[sid].selement_ids
        target_nids = self.geometry.sections[sid].snode_ids
    
        for eid in target_eids:
            del self.geometry.elements[eid]
    
        for nid in target_nids:
            del self.geometry.nodes[nid]

        del self.geometry.sections[sid]
    
        print '*** Done ***'
        
    def __str__(self):
        str = 'VCSys: %s'%self.name
        return str

class VCGeometry(object):
    def __init__(self,sys):
        self.base_lat_lon = None
        self.sections = None
        self.elements = None
        self.nodes = None
        self.converter = Converter.Converter()
        self.element_section_map = None
        self.proximity_table = None
        self.distance_table = None
        
        self.max_lat = -90.0
        self.min_lat = 90.0
        self.max_lon = -180.0
        self.min_lon = 180.0
        self.max_depth = -py_sys.float_info.max
        self.min_depth = 0.0
        
        self.sys = sys
    
    @property
    def layered_sections(self):
        return self.sys.layered_sections
    
    def getLayeredSections(self, section_list):
        layered_sections = {}
        for sid in section_list:
            layered_sections[sid] = self.layered_sections[sid]
        return layered_sections
    
    def sectionsInArea(self, NE_pt, SW_pt):
        #NE_pt is the LatLonPoint of the north east corner of the area
        #SW_pt is the LatLonPoint of the south west corner of the area
        
        max_lat = NE_pt.lat
        max_lon = NE_pt.lon
        min_lat = SW_pt.lat
        min_lon = SW_pt.lon
        
        sections = []
        
        for sec in self.sections.itervalues():
            #curr_min_lat = sec.minLat()
            #curr_max_lat = sec.maxLat()
            #curr_min_lon = sec.minLon()
            #curr_max_lon = sec.maxLon()
            
            self.sys.geometry.nodes[sec.traceStartNode()]
            self.sys.geometry.nodes[sec.traceEndNode()]
            
            max_in = False
            min_in = False
            
            if (    curr_max_lat <= max_lat and
                    curr_max_lat >= min_lat and
                    curr_max_lon <= max_lon and
                    curr_max_lon >= min_lon
                ):
                max_in = True
            
            if (    curr_min_lat <= max_lat and
                    curr_min_lat >= min_lat and
                    curr_min_lon <= max_lon and
                    curr_min_lon >= min_lon
                ):
                min_in = True
            
            if max_in or min_in:
                print max_in, min_in, sec.sname, curr_max_lat, curr_max_lon, curr_min_lat, curr_min_lon
    
    def areaOfSections(self, sids):
        minLat = 90.0
        maxLat = -90.0
        minLon = 180.0
        maxLon = -180.0
        for sid in sids:
            sMinLat = self.sections[sid].minLat()
            sMaxLat = self.sections[sid].maxLat()
            sMinLon = self.sections[sid].minLon()
            sMaxLon = self.sections[sid].maxLon()
            
            if sMinLat < minLat:
                minLat = sMinLat
            if sMaxLat > maxLat:
                maxLat = sMaxLat
            if sMinLon < minLon:
                minLon = sMinLon
            if sMaxLon > maxLon:
                maxLon = sMaxLon
            
        return minLat, maxLat, minLon, maxLon
        
    def proximityTable(self):
        proximity_table = []
        #distance_table = []
        for eid, e in self.elements.iteritems():
            #print 'stress drop:', e.dynamic_strength - e.static_strength
            my_center = e.centerVector()
            my_distances = []
            for teid, te in self.elements.iteritems():
                target_center = te.centerVector()
                my_distances.append((teid, (my_center - target_center).length()))
            
            proximity_table.append( [e for e in sorted(my_distances,lambda a,b: cmp(a[1],b[1]))] )
        
        return np.array(proximity_table)
                
    
    
        #if section_list is None:
        #    return self.layered_sections
        #else:
        #    layered_sections = {}
        #    for sid in section_list:
        #        layered_sections[sid] = self.layered_sections[sid]
        #    return layered_sections
    '''
    def layerSections(self, section_list = None):
        try:
            sections = self.sectionsSortedByID(section_list)
        except:
            print '!!! Sections have not been loaded !!!'
            return
        
        #max_y_offset = 0
        #max_x = 0
        
        self.layered_sections = {}
        
        print '    layering %i sections '%(len(sections)),
        
        for n, s in enumerate(sections):
            curr_section = []
            for eid in s.selement_ids:
                if self.sys.geometry.elements[eid].onTrace():
                    curr_section.append( s.elementsAlongDip(eid) )
            
            curr_section_np = np.array(curr_section)
            
            self.layered_sections[s.sid] = curr_section_np.T
            
            #if layered_sections[s.sid].shape[0] > max_y_offset:
                #max_y_offset = layered_sections[s.sid].shape[0]
            
            #print layered_sections[s.sid].shape, s.sid, s.sname
            
            #max_x += layered_sections[s.sid].shape[1]
            
            if ( n%10 == 0 ):
                print '%.1f%s'%((float(n)/float(len(sections)))* 100.0,'%' ),
            else:
                print '.',
            py_sys.stdout.flush()
        
        print
    '''
    def exportEqsimGeometry(self, geometry_file, friction_file, section_list = None):
        geom = quakelib.EQSimGeometryWriter()
        fric = quakelib.EQSimFrictionWriter()
        conv = quakelib.Conversion(quakelib.LatLonDepth(self.base_lat_lon.lat, self.base_lat_lon.lon, 0))
        
        #print self.converter.lat0, conv.get_base_lat_lon_depth().lat(), self.converter.lon0, conv.get_base_lat_lon_depth().lon()
        
        #if self.lame_lambda is None or self.lame_mu is None:
        #    fric.set_lame_lambda_mu(3.200000e010, 3.200000e010)
        #else:
        #    fric.set_lame_lambda_mu(self.lame_lambda, self.lame_mu)
        
        #print quakelib.Vec3(0,0,10)
        
        for sid, section in self.sections.iteritems():
            eqsim_section = geom.new_section()
            eqsim_section.set_sid(int(sid))
            eqsim_section.set_name(section.sname)
            #eqsim_element_list = []
            for eid in section.selement_ids:
                element = self.elements[eid]
                
                #eqsim_element = quakelib.ElementRect()
                eqsim_element = eqsim_section.new_rectangle()
                eqsim_element.set_index(int(element.eid)+1)
                eqsim_element.set_rake(element.rake)
                eqsim_element.set_slip_rate(element.slip)
                eqsim_element.set_perfect_flag(1)
                eqsim_element.set_dip(self.converter.rad_deg(element.dip()))
                #print '%f pi, %f deg'%(element.strike()/math.pi, self.converter.rad_deg(element.strike()))
                eqsim_element.set_strike(self.converter.rad_deg(element.strike()))
                
                for index, nid in enumerate(element.enode_ids):
                    node = self.nodes[nid]
                    
                    #x, y = self.converter.latlon_xy(node.lat, node.lon)
                    #z = node.depth
                    
                    #eqsim_xyz = conv.convert2xyz(quakelib.LatLonDepth(node.lat, node.lon, node.depth))
                    
                    #print abs(x - eqsim_xyz[0]), abs(y - eqsim_xyz[1]), abs(z - eqsim_xyz[2])
                    
                    #eqsim_element.set_vert(index, quakelib.Vec3(x,y,z))
                    
                    eqsim_vertex = eqsim_section.new_vertex()
                    eqsim_vertex.set_index(node.nid)
                    eqsim_vertex.set_loc(quakelib.LatLonDepth(node.lat, node.lon, node.depth))
                    eqsim_vertex.set_das(node.das)
                    #print x,eqsim_element.vert(index)[0], y, eqsim_element.vert(index)[1], z, eqsim_element.vert(index)[2]
                    
                    #my_conv = self.converter.xy_latlon(x,y)
                    #their_conv =  conv.convert2LatLon(quakelib.Vec3(x,y,z))
                    
                    #print abs(my_conv[0] - their_conv.lat()) , abs(my_conv[1] - their_conv.lon())
                    
                    if node.trace_flag == 1:
                        eqsim_vertex.set_trace_flag(quakelib.MIDDLE_TRACE)
                        #eqsim_element.set_trace_flag(index, quakelib.MIDDLE_TRACE)
                    elif node.trace_flag == 2:
                        eqsim_vertex.set_trace_flag(quakelib.BEGINNING_TRACE)
                        #eqsim_element.set_trace_flag(index, quakelib.BEGINNING_TRACE)
                    elif node.trace_flag == 3:
                        eqsim_vertex.set_trace_flag(quakelib.END_TRACE)
                        #eqsim_element.set_trace_flag(index, quakelib.END_TRACE)
                    else:
                        eqsim_vertex.set_trace_flag(quakelib.NOT_ON_TRACE)
                        #eqsim_element.set_trace_flag(index, quakelib.NOT_ON_TRACE)
                    
                    eqsim_element.set_vertex(index, eqsim_vertex.index())
                    #eqsim_element.set_das(index, node.das)
                
                #print eqsim_element.min_depth(), eqsim_element.max_depth()
                #ind = eqsim_section.add_element(eqsim_element, conv)
                #ind = eqsim_element.index()
                fric.set_strengths(eqsim_element.index(), element.static_strength, element.dynamic_strength)
                fric.set_lame_lambda_mu(element.lame_lambda, element.lame_mu)
                #eqsim_element_list.append(eqsim_element)
        
        geom.open(geometry_file)
        geom.write()
        geom.close()
        fric.open(friction_file)
        fric.write()
        fric.close()
        
    def KMLStyle_base(self):
        # base style for KML export
        basePolyColor = (int(255* 0.5),0,126,255)
        baseLineColor = (int(255* 1.0),0,200,255)
        baseLineWidth = 1
    
        return '<Style id="baseStyle">\n\t<LineStyle>\n\t\t<width>%i</width>\n\t\t<color>%02x%02x%02x%02x</color>\n\t</LineStyle>\n\t<PolyStyle>\n\t\t<color>%02x%02x%02x%02x</color>\n\t</PolyStyle>\n</Style>\n'%(baseLineWidth,baseLineColor[0],baseLineColor[3],baseLineColor[2],baseLineColor[1],basePolyColor[0],basePolyColor[3],basePolyColor[2],basePolyColor[1])
    
    def KMLStyle_section(self):
        return '<Style id="sectionLabel">\n\t<IconStyle>\n\t\t<Icon>\n\t\t\t<href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>\n\t\t</Icon>\n\t</IconStyle>\n</Style>\n'
    
    def KMLStyle_broke(self):
        brokePolyColor = (int(255* 0.5),255,0,0)
        brokeLineColor = (int(255* 1.0),255,0,0)
        brokeLineWidth = 1
        
        return '<Style id="brokeStyle">\n\t<LineStyle>\n\t\t<width>%i</width>\n\t\t<color>%02x%02x%02x%02x</color>\n\t</LineStyle>\n\t<PolyStyle>\n\t\t<color>%02x%02x%02x%02x</color>\n\t</PolyStyle>\n</Style>\n'%(brokeLineWidth, brokeLineColor[0],brokeLineColor[3],brokeLineColor[2],brokeLineColor[1], brokePolyColor[0],brokePolyColor[3],brokePolyColor[2],brokePolyColor[1])
    
    def KMLStyle_elementBreak(self, fade_steps):
        m3PolyColor = (int(255* 1),255,255,0)
        m3LineColor = (int(255* 1),255,255,0)
        m3LineWidth = 3
        
        m4PolyColor = (int(255* 1),255,212,0)
        m4LineColor = (int(255* 1),255,212,0)
        m4LineWidth = 4
        
        m5PolyColor = (int(255* 1),255,170,0)
        m5LineColor = (int(255* 1),255,170,0)
        m5LineWidth = 5
        
        m6PolyColor = (int(255* 1),255,127,0)
        m6LineColor = (int(255* 1),255,127,0)
        m6LineWidth = 6
        
        m7PolyColor = (int(255* 1),255,85,0)
        m7LineColor = (int(255* 1),255,85,0)
        m7LineWidth = 7
        
        m8PolyColor = (int(255* 1),255,42,0)
        m8LineColor = (int(255* 1),255,42,0)
        m8LineWidth = 8
        
        m9PolyColor = (int(255* 1),255,0,0)
        m9LineColor = (int(255* 1),255,0,0)
        m9LineWidth = 9
        
        str = ''
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle3_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m3LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m3LineColor[0])*frac),m3LineColor[3],m3LineColor[2],m3LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m3PolyColor[0])*frac),m3PolyColor[3],m3PolyColor[2],m3PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle4_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m4LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m4LineColor[0])*frac),m4LineColor[3],m4LineColor[2],m4LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m4PolyColor[0])*frac),m4PolyColor[3],m4PolyColor[2],m4PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle5_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m5LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m5LineColor[0])*frac),m5LineColor[3],m5LineColor[2],m5LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m5PolyColor[0])*frac),m5PolyColor[3],m5PolyColor[2],m5PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle6_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m6LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m6LineColor[0])*frac),m6LineColor[3],m6LineColor[2],m6LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m6PolyColor[0])*frac),m6PolyColor[3],m6PolyColor[2],m6PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle7_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m7LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m7LineColor[0])*frac),m7LineColor[3],m7LineColor[2],m7LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m7PolyColor[0])*frac),m7PolyColor[3],m7PolyColor[2],m7PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle8_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m8LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m8LineColor[0])*frac),m8LineColor[3],m8LineColor[2],m8LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m8PolyColor[0])*frac),m8PolyColor[3],m8PolyColor[2],m8PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        for num in range(0,fade_steps + 1):
            frac = 1.0 - float(num)/float(fade_steps)
            str += '<Style id="elementBreakStyle9_%i">\n'%num
            str += '    <LineStyle>\n'
            str += '        <width>%f</width>\n'%(float(m9LineWidth)*frac)
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m9LineColor[0])*frac),m9LineColor[3],m9LineColor[2],m9LineColor[1])
            str += '    </LineStyle>\n'
            str += '    <PolyStyle>\n'
            str += '        <color>%02x%02x%02x%02x</color>\n'%(int(float(m9PolyColor[0])*frac),m9PolyColor[3],m9PolyColor[2],m9PolyColor[1])
            str += '    </PolyStyle>\n'
            str += '</Style>\n'
        
        return str
    
    def minMaxLatLon(self):
        min_lat = 90.0
        max_lat = -90.0
        min_lon = 180.0
        max_lon = -180.0
        for sec in self.sections.itervalues():
            curr_min_lat = sec.minLat()
            curr_max_lat = sec.maxLat()
            curr_min_lon = sec.minLon()
            curr_max_lon = sec.maxLon()
            #print sec.minLat(), sec.maxLat(), sec.minLon(), sec.maxLon()
            
            if curr_min_lat < min_lat:
                min_lat = curr_min_lat
            
            if curr_max_lat > max_lat:
                max_lat = curr_max_lat
            
            if curr_min_lon < min_lon:
                min_lon = curr_min_lon
            
            if curr_max_lon > max_lon:
                max_lon = curr_max_lon
        
        return min_lat, max_lat, min_lon, max_lon
     
    def sectionsSortedBySectionSize(self):
        return [s for s in sorted(self.sections.values(),lambda a,b: cmp(len(a.selement_ids),len(b.selement_ids)))]
           
    def sectionsSortedByID(self, section_list=None):
        if section_list is not None:
            return [self.sections[k] for k in sorted(self.sections.keys()) if k in section_list]
        else:
            return [self.sections[k] for k in sorted(self.sections.keys())]
    
    def elementsSortedByID(self, element_list=None):
        if element_list is not None:
            return [self.elements[k] for k in sorted(self.elements.keys()) if k in element_list]
        else:
            return [self.elements[k] for k in sorted(self.elements.keys())]
    
    def nodesSortedByID(self):
        return [self.nodes[k] for k in sorted(self.nodes.keys())]
    
    def setBaseLatLon(self, base_lat_lon):
        self.converter.setLatLon0(base_lat_lon.lat, base_lat_lon.lon)
        self.base_lat_lon = base_lat_lon
    
    def addSection(self, sid):
        if self.sections is None:
            self.sections = {}
        section = VCSection(self.sys)
        section.sid = sid
        self.sections[sid] = section
        del section
    
    def addQuad(self, eid):
        if self.elements is None:
            self.elements = {}
        element = VCQuad(self.sys)
        element.eid = eid
        self.elements[eid] = element
        del element
    
    def addNode(self, nid):
        if self.nodes is None:
            self.nodes = {}
        node = VCNode(self.sys)
        node.nid = nid
        self.nodes[nid] = node
        del node
    
    def geometryFromTraceFile(self, trace_file):
        f = trace_file
        
        if self.base_lat_lon is None:
            self.setBaseLatLon( vec.LatLonPoint(36.0, -120.5) )
        
        if self.sections is None:
            sid = 0
        else:
            sid = self.sectionsSortedByID()[-1].sid
        
        current_section_info = None
        current_trace_points = None
        
        for line in f:
            if not line.startswith("#"):
                split_line = line.split()
                if len(split_line) != 0:
                    # This is the start of a section
                    if split_line[0] == "section":
                        if current_trace_points is not None:
                            # add the geometry to the section using the trace points
                            self.addGeometryToSectionFromTracePoints(sid, current_trace_points, current_section_info[0], current_section_info[1], current_section_info[2], current_section_info[3], current_section_info[4], current_section_info[5], current_section_info[6], current_section_info[7], current_section_info[8], current_section_info[9], current_section_info[10], current_section_info[11])
                            current_trace_points = None
                        sid += 1
                        self.addSection(sid)
                        self.sections[sid].sname  = split_line[1]
                        self.sections[sid].fid    = int(split_line[2])
                        current_section_info = map( float, split_line[3:])
                    else:
                        # this is a lat lon point defining a fault trace
                        if current_trace_points is None:
                            current_trace_points = []
                        current_trace_points.append( vec.LatLonPoint(float(split_line[0]), float(split_line[1])) )
        
        #print sid, current_trace_points, current_section_info[0], 0, current_section_info[1], current_section_info[2], current_section_info[3], current_section_info[4], current_section_info[5], current_section_info[6], current_section_info[7], current_section_info[8], current_section_info[9]
        # add the geometry to the section using the trace points. this is for the last section only
        self.addGeometryToSectionFromTracePoints(sid, current_trace_points, current_section_info[0], current_section_info[1], current_section_info[2], current_section_info[3], current_section_info[4], current_section_info[5], current_section_info[6], current_section_info[7], current_section_info[8], current_section_info[9], current_section_info[10], current_section_info[11])                
                
    def addGeometryToSectionFromTracePoints(self, sid, trace_points, target_element_size, top_depth, depth_along_dip, dip_angle, rake_angle, slip_rate, slip_rate_noise_percentage, static_strength, dynamic_strength, aseismicity_factor, lame_lambda, lame_mu):
        
        if len(trace_points) == 1:
            print "at least two trace points are needed to create a section"
            return
            
        # so we dont get any duplicate nodes
        tmp_nodes = {}
        
        try:
            nid = self.nodesSortedByID()[-1].nid + 1
        except AttributeError:
            nid = 0
        try:
            eid = self.elementsSortedByID()[-1].eid + 1
        except AttributeError:
            eid = 0
        
        for point_num, trace_point in enumerate(trace_points):
            this_x, this_y = self.converter.latlon_xy(trace_point.lat, trace_point.lon)
            this_z = -top_depth
            if point_num + 1 == len(trace_points):
                # this is the last point
                #print "last point"
                pass
            else:
                next_x, next_y = self.converter.latlon_xy(trace_points[point_num+1].lat, trace_points[point_num+1].lon)
                next_z = -top_depth
                horiz_distance = ((this_x-next_x)**2.0+(this_y-next_y)**2.0)**(0.5)
                horiz_mesh_elements = int(round(horiz_distance/target_element_size))
                element_size = horiz_distance/float(horiz_mesh_elements)
                vert_mesh_elements = int(round(depth_along_dip/element_size))
                #print element_size, horiz_mesh_elements, vert_mesh_elements
                #element_size = horiz_distance/25.0
                #horiz_mesh_elements = int(round(horiz_distance/element_size))
                #vert_distance = abs(depth_top - depth_bottom)
                #vert_mesh_elements = int(round(vert_distance/element_size))
                
                this_vec = vec.MyVector(this_x,this_y,this_z)
                next_vec = vec.MyVector(next_x,next_y,next_z)
                
                horiz_sep_vec = (next_vec - this_vec).norm()
                
                y_axis = vec.MyVector(0,1,0)
                x_axis = vec.MyVector(1,0,0)
                
                strike_angle = self.converter.rad_deg(y_axis < horiz_sep_vec)
                strike_angle_compliment = self.converter.rad_deg(x_axis < horiz_sep_vec)
                
                #print self.sections[sid].sname
                
                #print strike_angle, strike_angle_compliment
                if strike_angle_compliment > 90.0:
                    strike_angle = 360.0 - strike_angle
                #print strike_angle
                
                d = vert_mesh_elements * element_size
                dip = self.converter.deg_rad(dip_angle)
                strike = self.converter.deg_rad(strike_angle)
                
                #print horiz_sep_vec
                vert_sep_vec = (vec.MyVector(this_vec.x - d * math.cos(dip) * math.cos(strike), this_vec.y + d * math.cos(dip) * math.sin(strike), this_vec.z - d * math.sin(dip) ) - this_vec).norm()
                #print vert_sep_vec
                
                for horiz_n in range(horiz_mesh_elements):
                    for vert_n in range(vert_mesh_elements):
                        
                        # Set up the element
                        self.sections[sid].addQuad(eid)
                        
                        self.elements[eid].sid               = sid
                        self.elements[eid].fid               = self.sections[sid].fid
                        self.elements[eid].slip              = slip_rate + random.uniform(-slip_rate_noise_percentage,slip_rate_noise_percentage) * slip_rate
                        self.elements[eid].aseis             = aseismicity_factor
                        self.elements[eid].rake              = rake_angle
                        self.elements[eid].dynamic_strength  = dynamic_strength
                        self.elements[eid].static_strength   = static_strength
                        self.elements[eid].lame_lambda       = lame_lambda
                        self.elements[eid].lame_mu           = lame_mu
                        
                        # Set up the nodes
                        root_trace_flag     = 0
                        two_trace_flag      = 0
                        three_trace_flag    = 0
                        four_trace_flag     = 0
                        if horiz_n == 0 and vert_n == 0:
                            root_trace_flag = 2
                            four_trace_flag = 1
                        elif horiz_n == horiz_mesh_elements - 1 and vert_n == 0:
                            root_trace_flag = 1
                            four_trace_flag = 3
                        elif vert_n == 0:
                            root_trace_flag = 1
                            four_trace_flag = 1
                        
                        root_vec    = this_vec + horiz_sep_vec * element_size * horiz_n + vert_sep_vec * element_size * vert_n
                        two_vec     = root_vec + vert_sep_vec * element_size
                        three_vec   = root_vec + vert_sep_vec * element_size + horiz_sep_vec * element_size
                        four_vec    = root_vec + horiz_sep_vec * element_size
                        
                        root_das    = element_size * horiz_n
                        two_das     = root_das
                        three_das   = element_size * (horiz_n + 1.0)
                        four_das    = three_das
                        
                        self.elements[eid].addNode(nid)
                        self.sections[sid].addNode(nid)
                        self.nodes[nid].eid         = eid
                        lat, lon = self.converter.xy_latlon(root_vec.x, root_vec.y)
                        self.nodes[nid].lat         = lat
                        self.nodes[nid].lon         = lon
                        self.nodes[nid].depth       = root_vec.z
                        self.nodes[nid].das         = root_das
                        self.nodes[nid].trace_flag  = root_trace_flag
                        nid += 1
                        
                        self.elements[eid].addNode(nid)
                        self.sections[sid].addNode(nid)
                        self.nodes[nid].eid         = eid
                        lat, lon = self.converter.xy_latlon(two_vec.x, two_vec.y)
                        self.nodes[nid].lat         = lat
                        self.nodes[nid].lon         = lon
                        self.nodes[nid].depth       = two_vec.z
                        self.nodes[nid].das         = two_das
                        self.nodes[nid].trace_flag  = two_trace_flag
                        nid += 1
                        
                        self.elements[eid].addNode(nid)
                        self.sections[sid].addNode(nid)
                        self.nodes[nid].eid         = eid
                        lat, lon = self.converter.xy_latlon(three_vec.x, three_vec.y)
                        self.nodes[nid].lat         = lat
                        self.nodes[nid].lon         = lon
                        self.nodes[nid].depth       = three_vec.z
                        self.nodes[nid].das         = three_das
                        self.nodes[nid].trace_flag  = three_trace_flag
                        nid += 1
                        
                        self.elements[eid].addNode(nid)
                        self.sections[sid].addNode(nid)
                        self.nodes[nid].eid         = eid
                        lat, lon = self.converter.xy_latlon(four_vec.x, four_vec.y)
                        self.nodes[nid].lat         = lat
                        self.nodes[nid].lon         = lon
                        self.nodes[nid].depth       = four_vec.z
                        self.nodes[nid].das         = four_das
                        self.nodes[nid].trace_flag  = four_trace_flag
                        nid += 1
                        
                        eid += 1
                        
    def loadHdf5Geometry(self, geometry_file, section_list = None):
        f = geometry_file
        
        latlon = f['base_lat_lon']
        
        self.setBaseLatLon( vec.LatLonPoint(float(latlon[0]), float(latlon[1])) )
        
        #c = quakelib.Conversion(quakelib.LatLonDepth(float(latlon[0]), float(latlon[1]), 0))
        
        vc_blocks = f['block_info_table']
        
        nid = 1
        block_num = 0
        total_blocks = len(vc_blocks)
        
        for block in vc_blocks:
            block_num += 1
            sid = block[2]

            if section_list is None or sid in section_list:
                eid = block[0]
                
                try:
                    self.sections[sid].addQuad(eid)
                except (KeyError,TypeError):
                    self.addSection(sid)
                    self.sections[sid].sname    = block[31]
                    self.sections[sid].fid      = block[1]
                    self.sections[sid].addQuad(eid)
                
                self.elements[eid].sid                = sid
                self.elements[eid].fid                = block[1]
                self.elements[eid].slip               = block[23]
                self.elements[eid].aseis              = block[24]
                self.elements[eid].rake               = self.converter.rad_deg(block[25])
                self.elements[eid].dynamic_strength   = block[27]
                self.elements[eid].static_strength    = block[28]
                self.elements[eid].lame_mu            = block[29]
                self.elements[eid].lame_lambda        = block[30]
                
                for i in range(4):
                    self.elements[eid].addNode(nid)
                    #self.sections[sid].addNode(nid)
                    self.nodes[nid].eid         = eid
                    #lat_lon_depth = c.convert2LatLon(quakelib.Vec3(block[3 + i*5],block[4 + i*5],block[5 + i*5]))
                    lat, lon = self.converter.xy_latlon(block[3 + i*5],block[4 + i*5])
                    #self.nodes[nid].lat         = lat_lon_depth.lat()
                    #self.nodes[nid].lon         = lat_lon_depth.lon()
                    self.nodes[nid].lat         = lat
                    self.nodes[nid].lon         = lon
                    self.nodes[nid].depth       = block[5 + i*5]
                    self.nodes[nid].das         = block[6 + i*5]
                    self.nodes[nid].trace_flag  = block[7 + i*5]
                    #print eid, nid, self.nodes[nid].trace_flag
                    nid += 1
                
                self.elements[eid].setSideLength()
                
            #print eid, self.elements[eid].dynamic_strength, self.elements[eid].static_strength
            
            if ( block_num%1000 == 0 ):
                    print '    block:%i of %i'%(block_num, total_blocks)
        
        self.finishGeometryLoad()
    
    def finishGeometryLoad(self):
        for node in self.nodes.itervalues():
            lat = node.lat
            lon = node.lon
            depth = node.depth
            
            if lat > self.max_lat: self.max_lat = lat
            if lat < self.min_lat: self.min_lat = lat
            
            if lon > self.max_lon: self.max_lon = lon
            if lon < self.min_lon: self.min_lon = lon
            
            if depth > self.max_depth: self.max_depth = depth
            if depth < self.min_depth: self.min_depth = depth
    
    def loadEqsimGeometry(self, geometry_file, friction_file, section_list = None):
        pass
    
    def exportHdf5Geometry(self, out_file, section_list = None):
        pass
        
    '''
    def __str__(self):
        str = ''
        
        if self.base_lat_lon is not None:
            str += 'base lat: ' + self.base_lat_lon.lat + ', base lon: ' + self.base_lat_lon.lon + '\n'
        
        if self.sections is not None:
            str += len(self.sections) + ' sections\n'
        
        if self.elements is not None:
            str += len(self.elements) + ' elements\n'
        
        if self.nodes is not None:
            str += len(self.nodes) + ' nodes\n'
        
        if self.sections is not None:
            for section in [ self.sections[k] for k in sorted(self.sections.keys())]:
                str += '    ' + section.__str__()
        
        return str
    '''

class VCSection(object):
    def __init__(self,sys):
        self.sid = None
        self.sname = None
        self.fid = None
        self.selement_ids = None
        self.snode_ids = None
        self.sys = sys
        #super(VCSection,self).__init__()
    
    @property
    def surface_area(self):
        surface_area = 0.0
        for eid in self.selement_ids:
            surface_area += self.sys.elements[eid].area()
        return surface_area

    def exportKMLSectionLabel(self, out_file):
        trace_start_node_id = self.traceStartNode()
        if trace_start_node_id is None:
            trace_start_node_id = self.firstOnTraceNode()
        if trace_start_node_id is None:
            #we are screwed
            print '!! Could not find a node to attach the flag to !!'
            return
        f = out_file
        f.write('   <Placemark id="section_%i_label">\n'%(self.sid))
        f.write('       <name>%i %s</name>\n'%(self.sid, self.sname))
        #f.write('       <description>Tethered to the ground by a customizable &quot;tail&quot;</description>\n')
        f.write('       <styleUrl>#sectionLabel</styleUrl>\n')
        f.write('       <Point>\n')
        f.write('           <extrude>1</extrude>\n')
        f.write('           <altitudeMode>relativeToGround</altitudeMode>\n')
        f.write('           <coordinates>%f,%f,%f</coordinates>\n'%(self.sys.geometry.nodes[trace_start_node_id].lon,self.sys.geometry.nodes[trace_start_node_id].lat,self.sys.geometry.nodes[trace_start_node_id].depth+abs(self.sys.geometry.min_depth)))
        f.write('       </Point>\n')
        f.write('   </Placemark>\n')
    
    def exportKMLGeometry(self, out_file):
        f = out_file
        
        f.write('   <Folder id="section_%i">\n'%(self.sid))
        f.write('       <name>%i %s</name>\n'%(self.sid, self.sname))
        
        for eid in self.selement_ids:
            self.sys.geometry.elements[eid].exportKML(f)
        
        f.write('   </Folder>\n')
        
    def traceStartNode(self):
        for nid in self.snode_ids:
            #print nid, self.sys.geometry.nodes[nid].eid, self.sys.geometry.nodes[nid].trace_flag
            if self.sys.geometry.nodes[nid].trace_flag == 2:
                #print nid
                return nid
        return None
    
    def traceEndNode(self):
        for nid in self.snode_ids:
            #print nid, self.sys.geometry.nodes[nid].eid, self.sys.geometry.nodes[nid].trace_flag
            if self.sys.geometry.nodes[nid].trace_flag == 3:
                #print nid
                return nid
        return None
    
    def firstOnTraceNode(self):
        for nid in self.snode_ids:
            #print nid, self.sys.geometry.nodes[nid].eid, self.sys.geometry.nodes[nid].trace_flag
            if self.sys.geometry.nodes[nid].trace_flag == 1:
                #print nid
                return nid
        return None
    
    def traceLocs(self):
        lats = []
        lons = []
        #print self.snode_ids
        for nid in self.snode_ids:
            #print self.sys.geometry.nodes[nid].trace_flag
            if self.sys.geometry.nodes[nid].trace_flag == 1 or self.sys.geometry.nodes[nid].trace_flag == 2 or self.sys.geometry.nodes[nid].trace_flag == 3:
                lat = self.sys.geometry.nodes[nid].lat
                lon = self.sys.geometry.nodes[nid].lon
                #if lat not in lats and lon not in lons:
                lats.append(lat)
                lons.append(lon)
        return lats, lons
                
    #returns (element id, distance in meters)
    def nearestElements(self, source_eid):
        element_distances = {}
        
        source_vector = self.sys.geometry.nodes[self.sys.geometry.elements[source_eid].enode_ids[0]].cartVector()
            
        for eid in self.selement_ids:
            test_vector = self.sys.geometry.nodes[self.sys.geometry.elements[eid].enode_ids[0]].cartVector()
            element_distances[eid] = (test_vector - source_vector).length()
        
        return_elements = []
        
        source_side_length = self.sys.geometry.elements[source_eid].sideLength()
        
        for eid,distance in element_distances.iteritems():
            target_side_length = self.sys.geometry.elements[eid].sideLength()
            if distance < 1.5 * (source_side_length + target_side_length) * 0.5 and eid != source_eid:
                return_elements.append((eid,distance))
        
        #print return_elements
        return return_elements
    
    #returns (next element id, next element distance, next element node[0] depth)
    def nextElementDownDip(self, source_eid):
        neighbor_elements_and_distances = self.nearestElements(source_eid)
        neighbor_elements_and_distances_and_depths = [(eid,distance,self.sys.geometry.nodes[self.sys.geometry.elements[eid].enode_ids[0]].depth) for eid,distance in neighbor_elements_and_distances]
        
        #print source_eid
        #print neighbor_elements_and_distances
        #print sorted(neighbor_elements_and_distances_and_depths,lambda a,b: cmp((a[2],a[1]),(b[2],b[1])))
        
        #curr_depth = VCGeometry.nodes[VCGeometry.elements[source_eid].enode_ids[0]].depth
        #ret_ele = []
        #min_distance = sys.float_info.max
        
        #for ele in neighbor_elements_and_distances_and_depths:
        #        if ele[2] < curr_depth:
        #            min_distance = ele[1]
                    
        
        #try:
        if len(neighbor_elements_and_distances_and_depths) > 1:
            #print sorted(neighbor_elements_and_distances_and_depths,lambda a,b: cmp((a[2],a[1]),(b[2],b[1])))
            return sorted(neighbor_elements_and_distances_and_depths,lambda a,b: cmp((a[2],a[1]),(b[2],b[1])))[0]
        else:
            return None
        #except:
            
            
        
    def elementsAlongDip(self, start_eid = None):
        if start_eid is None:
            start_element = self.sys.geometry.nodes[self.traceStartNode()].eid
        else:
            start_element = start_eid
        elements = [start_element]
        previous_depth = self.sys.geometry.nodes[self.sys.geometry.elements[start_element].enode_ids[0]].depth
        previous_element = start_element
        for i in range(50):
            current_element = self.nextElementDownDip(previous_element)
            if current_element is not None:
                current_depth = current_element[2]
                
                #print previous_depth, current_depth
                if current_depth >= previous_depth:
                    break
                else:
                    elements.append(current_element[0])
                    previous_depth = current_depth
                    previous_element = current_element[0]
            else:
                break
        
        #print elements
        #print 
        return elements
         
    def addQuad(self, eid):
        if self.selement_ids is None:
            self.selement_ids = []
        self.selement_ids.append(eid)
        self.sys.geometry.addQuad(eid)
    
    def addNode(self, nid):
        if self.snode_ids is None:
            self.snode_ids = []
        self.snode_ids.append(nid)
    
    def layerCount(self):
        return len(self.elementsAlongDip())
    
    def nodeCount(self):
        pass
    
    def elementCount(self):
        pass
    
    def minLat(self):
        lats = sorted([self.sys.geometry.nodes[nid].lat for nid in self.snode_ids])
        return lats[0]
    
    def maxLat(self):
        lats = sorted([self.sys.geometry.nodes[nid].lat for nid in self.snode_ids])
        return lats[-1]
    
    def minLon(self):
        lons = sorted([self.sys.geometry.nodes[nid].lon for nid in self.snode_ids])
        return lons[0]
    
    def maxLon(self):
        lons = sorted([self.sys.geometry.nodes[nid].lon for nid in self.snode_ids])
        return lons[-1]
    
    def minDepth(self):
        depths = sorted([self.sys.geometry.nodes[nid].depth for nid in self.snode_ids])
        return depths[0]
    
    def maxDepth(self):
        depths = sorted([self.sys.geometry.nodes[nid].depth for nid in self.snode_ids])
        return depths[-1]
        
    
    def __str__(self):
        str = ''
        
        str += '%i %i %s'%(self.sid, self.fid, self.sname)
        
        return str
    
class VCElement(object):
    def __init__(self,sys):
        self.eid                = None
        self.sid                = None
        self.fid                = None
        self.slip               = None
        self.aseis              = None
        self.rake               = None
        self.enode_ids          = None
        self.dynamic_strength   = None
        self.static_strength    = None
        self.lame_lambda        = None
        self.lame_mu            = None
        self.side_length        = None
        self.sys = sys
        #super(VCElement,self).__init__()
    
    def exportKML(self, out_file, kml_style = 'baseStyle', time_span_start = None, time_span_end = None):
        f = out_file
        
        offset = abs(self.sys.geometry.min_depth)
        
        lat1    = self.sys.geometry.nodes[self.enode_ids[0]].lat
        lon1    = self.sys.geometry.nodes[self.enode_ids[0]].lon
        z1      = self.sys.geometry.nodes[self.enode_ids[0]].depth + offset
        lat2    = self.sys.geometry.nodes[self.enode_ids[1]].lat
        lon2    = self.sys.geometry.nodes[self.enode_ids[1]].lon
        z2      = self.sys.geometry.nodes[self.enode_ids[1]].depth + offset
        lat3    = self.sys.geometry.nodes[self.enode_ids[2]].lat
        lon3    = self.sys.geometry.nodes[self.enode_ids[2]].lon
        z3      = self.sys.geometry.nodes[self.enode_ids[2]].depth + offset
        lat4    = self.sys.geometry.nodes[self.enode_ids[3]].lat
        lon4    = self.sys.geometry.nodes[self.enode_ids[3]].lon
        z4      = self.sys.geometry.nodes[self.enode_ids[3]].depth + offset
        
        f.write('       <Placemark>\n')
        if time_span_start is not None and time_span_end is not None:
                f.write('           <TimeSpan>\n')
                f.write('               <begin>%i-%02i-%02i</begin>\n'%(time_span_start[0], time_span_start[1], time_span_start[2]))
                f.write('               <end>%i-%02i-%02i</end>\n'%(time_span_end[0], time_span_end[1], time_span_end[2]))
                f.write('           </TimeSpan>\n')
                
        f.write('           <styleUrl>#%s</styleUrl>\n'%kml_style)
        f.write('           <Polygon>\n')
        f.write('               <extrude>0</extrude>\n')
        f.write('               <altitudeMode>relativeToGround</altitudeMode>\n')
        f.write('               <outerBoundaryIs>\n')
        f.write('                   <LinearRing>\n')
        f.write('                       <coordinates>\n')
        f.write('                           %f,%f,%f\n'%(lon2,lat2,z2))
        f.write('                           %f,%f,%f\n'%(lon1,lat1,z1))
        f.write('                           %f,%f,%f\n'%(lon4,lat4,z4))
        f.write('                           %f,%f,%f\n'%(lon3,lat3,z3))
        f.write('                           %f,%f,%f\n'%(lon2,lat2,z2))
        f.write('                       </coordinates>\n')
        f.write('                   </LinearRing>\n')
        f.write('               </outerBoundaryIs>\n')
        f.write('           </Polygon>\n')
        f.write('       </Placemark>\n')
    
    def setSideLength(self):
        vector0 = self.sys.geometry.nodes[self.enode_ids[0]].cartVector()
        vector1 = self.sys.geometry.nodes[self.enode_ids[1]].cartVector()
        self.side_length = (vector0 - vector1).length()
        #if self.eid == 2550:
        #    print [[a.depth,a.nid] for a in sorted([self.sys.geometry.nodes[nid] for nid in self.enode_ids],lambda a,b: cmp(a.depth,b.depth))]
        #depths = sorted([self.sys.geometry.nodes[nid].depth for nid in self.enode_ids])
        #top = self.topNorthNode()
        #bottom = self.bottomNode()
        #print top, self.sys.geometry.nodes[top].depth, bottom, self.sys.geometry.nodes[bottom].depth
        #return abs(depths[0] - depths[2])
    
    def sideLength(self):
        return self.side_length
    
    def area(self):
        return self.side_length**2.0
    
    def onTrace(self):
        if self.sys.geometry.nodes[self.enode_ids[0]].trace_flag != 0:
            return True
        else:
            return False

    def pointLatLon(self, index=0):
        return (self.sys.geometry.nodes[self.enode_ids[index]].lat, self.sys.geometry.nodes[self.enode_ids[index]].lon)
    
    def centerVector(self):
        vec0 = self.sys.geometry.nodes[self.enode_ids[0]].cartVector()
        vec1 = self.sys.geometry.nodes[self.enode_ids[1]].cartVector()
        vec3 = self.sys.geometry.nodes[self.enode_ids[3]].cartVector()
        a = vec1 - vec0
        b = vec3 - vec0
        
        #print vec0, vec1, vec3, a, b
        #print (a+b) * 0.5, 0.5 * (a+b)
        return vec0 + 0.5 * (a+b)
    
    def normalVector(self):
        # this is a left handed normal so the conventions are consistant with Okada
        vector0 = self.sys.geometry.nodes[self.enode_ids[0]].cartVector()
        vector1 = self.sys.geometry.nodes[self.enode_ids[1]].cartVector()
        #vector2 = self.sys.geometry.nodes[self.enode_ids[2]].cartVector()
        vector3 = self.sys.geometry.nodes[self.enode_ids[3]].cartVector()
        
        vector01 = vector1 - vector0
        vector03 = vector3 - vector0
        
        return (vector01 * vector03).norm()
        
    def dip(self):
        return self.normalVector() < vec.MyVector(0,0,-1)
    
    def get_rake(self):
        return self.sys.geometry.converter.deg_rad(self.rake)
    
    def strike(self):
        #strike is positive when rotating west, negative east
        norm = self.normalVector()
        strike = (vec.MyVector(norm.x,norm.y,0) < vec.MyVector(0,1,0)) - math.pi/2.0
        return strike
        
    def bottomNode(self):
        pass
        
    def topNode(self):
        pass
        
    def topNorthNode(self):
        max_depth = -sys.float_info.max
        max_lat = -90.0
        max_nid = None
        for nid in self.enode_ids:
            if self.sys.geometry.nodes[nid].lat > max_lat and self.sys.geometry.nodes[nid].depth > max_depth:
                max_depth = self.sys.geometry.nodes[nid].depth
                max_lat = self.sys.geometry.nodes[nid].lat
                max_nid = nid
        return max_nid
    
    def das(self):
        das0 = self.sys.geometry.nodes[self.enode_ids[0]].das
        das1 = self.sys.geometry.nodes[self.enode_ids[1]].das
        das2 = self.sys.geometry.nodes[self.enode_ids[2]].das
        das3 = self.sys.geometry.nodes[self.enode_ids[3]].das
        return (das0, das1, das2, das3)
    
    def depths(self):
        depth0 = self.sys.geometry.nodes[self.enode_ids[0]].depth
        depth1 = self.sys.geometry.nodes[self.enode_ids[1]].depth
        depth2 = self.sys.geometry.nodes[self.enode_ids[2]].depth
        depth3 = self.sys.geometry.nodes[self.enode_ids[3]].depth
        return (depth0, depth1, depth2, depth3)
        
    def addNode(self, nid):
        if self.enode_ids is None:
            self.enode_ids = []
        self.enode_ids.append(nid)
        self.sys.geometry.addNode(nid)
        try:
            self.sys.geometry.sections[self.sid].addNode(nid)
        except:
            print '!!! Sections have been loaded improperly !!!'
            
            
    '''     
    def __str__(self):
        pass
    '''

class VCQuad(VCElement):
    def __init__(self,sys):
        super(VCQuad,self).__init__(sys)
    '''
    def __str__(self):
        pass
    '''

class VCTri(VCElement):
    def __init__(self,sys):
        super(VCTri,self).__init__(sys)
    '''
    def __str__(self):
        pass
    '''

class VCNode(object):
    def __init__(self,sys):
        self.nid         = None
        self.eid         = None
        self.lat         = None
        self.lon         = None
        self.depth       = None
        self.das         = None
        self.trace_flag  = None
        self.sys = sys
        #super(VCNode,self).__init__()
    
    def cartVector(self):
        x,y = self.sys.geometry.converter.latlon_xy(self.lat, self.lon)
        return vec.MyVector(x,y,self.depth)
    '''
    def __str__(self):
        pass      
    '''

class VCEvents(object):
    def __init__(self,sys):
        self.sys = sys
        self._years = None
        self._event_ids = None
        self._mags = None
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
    
    @property
    def years(self):
        f = h5py.File(self.sys.data_file, 'r')
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
    
    '''
    def loadHdf5Events(self, event_file):
                
        f = h5py.File(event_file, 'r')
            
        events = f['event_table']
        #sweeps = f['event_sweep_table']
        
        #number_of_events = len(events)
        #self.years = np.zeros(number_of_events)
        #self.event_ids = np.zeros(number_of_events,dtype=int)
        #self.mags = np.zeros(number_of_events)
        
        self.years = events['event_year']
        self.event_ids = events['event_number']
        self.mags = events['event_magnitude']
        
        #self.simple_event_list = []
        
        #for ind, raw_event in enumerate(events):
            #for the event year map
            #self.years[ind] = raw_event[1]
            #self.event_ids[ind] = raw_event[0]
            #self.mags[ind] = raw_event[3]
            #event = []
            #event number
            #event.append(raw_event[0])
            #event year
            #event.append(raw_event[1])
            #event magnitude
            #event.append(raw_event[3])
            #trigger section
            #event.append(self.sys.geometry.sections[self.sys.geometry.elements[raw_event[2]].sid].sname)
            #trigger element
            #event.append(raw_event[2])
            #involved sections and elements and total slip
            #involved_sections = {}
            #involved_elements = {}
            #total_slip = 0.0
            
            #for sweep in sweeps[raw_event[8]:raw_event[9]]:
            #    eid = sweep[2]
            #    sid = self.sys.geometry.elements[eid].sid
            #    total_slip += sweep[3]
            #    involved_sections[sid] = sid
            #    involved_elements[eid] = eid
        
            #event.append(len(involved_sections))
            #event.append(len(involved_elements))
            #event.append(total_slip/float(len(involved_elements)))
            
            #self.simple_event_list.append(event)
    
            #if ( ind%10000 == 0 ):
            #    print '    event:%i of %i'%(ind, number_of_events)
    
        cPickle.dump([self.years, self.event_ids], open('%s/years-event_ids.pkl'%(self.sys.cache_dir), 'wb'))
        cPickle.dump([self.years, self.mags], open('%s/years-mags.pkl'%(self.sys.cache_dir), 'wb'))
        '''
    def eventForYear(self, year):
        '''if self.years is None:
            try:
                self.years, self.event_ids = cPickle.load(open('%s/years-event_ids.pkl'%(self.sys.cache_dir), 'rb'))
            except IOError:
                self.sys.loadHdf5Events(self.sys.data_file)
                #self.eventYearMappingHdf5(self.data_file)
        '''       
        it = np.nditer(self.years, [], ['readonly'])
        event = None
        while not it.finished:
            if it[0] >= year:
                event = self.event_ids[it.iterindex] - 1
                break
            it.iternext()
        return event


'''
class VCEvents(object):
    def __init__(self,sys):
        self.sys = sys
        self.events = None
        self.sweeps = None
        self.bg_events = None
        self.as_events = None
        self.years = None
        self.event_ids = None
        
        self.max_allowed_year    = 9012.0
        self.start_year          = 2012.0
    
    def eventForYear(self, year):
        it = np.nditer(self.years, [], ['readonly'])
        event = None
        while not it.finished:
            if it[0] >= year:
                event = self.event_ids[it.iterindex] - 1
                break
            it.iternext()
        return event
        
    def minMaxMagnitude(self):
        events = [e for e in sorted(self.events.values(),lambda a,b: cmp(a.magnitude,b.magnitude))]
        
        return (events[0].magnitude, events[-1].magnitude)
    
    def eventsSortedByYear(self):
        return [e for e in sorted(self.events.values(),lambda a,b: cmp(a.year,b.year))]
    
    def eventsSortedByMagnitude(self, order = None, section_list = None):
        if section_list is None:
            if order is None:
                return [e for e in sorted(self.events.values(),lambda a,b: cmp(a.magnitude,b.magnitude))]
            else:
                return [e for e in sorted(self.events.values(),lambda a,b: cmp(b.magnitude,a.magnitude))]
        else:
            return_events = []
            
            for e in self.events.values():
                for sid in section_list:
                    for eid in self.sys.geometry.sections[sid].selement_ids:
                        if e.elementInEvent(eid):
                            return_events.append(e)
            
            if order is None:
                return [e for e in sorted(return_events,lambda a,b: cmp(a.magnitude,b.magnitude))]
            else:
                return [e for e in sorted(return_events,lambda a,b: cmp(b.magnitude,a.magnitude))]
    
    def eventsSortedByID(self, start_event = 0, end_event = None):
        return [self.events[k] for k in sorted(self.events.keys())]
    
    def addEvent(self, eid):
        if self.events is None:
            self.events = {}
        event = VCEvent(self.sys)
        event.eid = eid
        self.events[eid] = event
        del event
    
    def addSweep(self, sid):
        if self.sweeps is None:
            self.sweeps = {}
        sweep = VCSweep(self.sys)
        sweep.sid = sid
        self.sweeps[sid] = sweep
        del sweep
    
    def loadHdf5Events(self, event_file, start_event = 0, end_event = None, event_unit = 'event_num'):
        f = event_file
            
        events = f['event_table']
        sweeps = f['event_sweep_table']
        aftershocks = f['aftershock_table']
        bgevents = f['background_event_table']
        
        if event_unit == 'event_num':
            if end_event == None:
                end_event = len(events) - 1
            elif end_event > len(events) - 1:
                end_event = len(events) - 1
                
        curr_event = start_event
        
        print '    start event: %i, end event: %i'%(start_event, end_event)
        
        sid = 0
        
        while curr_event <= end_event:
            event = events[curr_event]
            self.addEvent(curr_event)
            
            self.events[curr_event].year            = event[1]
            self.events[curr_event].trigger_id      = event[2]
            self.events[curr_event].magnitude       = event[3]
            self.events[curr_event].shear_init      = event[4]
            self.events[curr_event].normal_init     = event[5]
            self.events[curr_event].shear_final     = event[6]
            self.events[curr_event].normal_final    = event[7]
            
            for sweep in sweeps[event[8]:event[9]]:
                self.events[curr_event].addSweep(sid)
                
                self.sweeps[sid].snum           = sweep[1]
                self.sweeps[sid].seid           = sweep[2]
                self.sweeps[sid].slip           = sweep[3]
                self.sweeps[sid].shear_init     = sweep[6]
                self.sweeps[sid].normal_init    = sweep[7]
                self.sweeps[sid].shear_final    = sweep[8]
                self.sweeps[sid].normal_final   = sweep[9]
                
                sid += 1
            
            # do aftershocks
            
            if ( curr_event%1000 == 0 ):
                    print '    event:%i of %i'%(curr_event, end_event)
            
            curr_event += 1
            
            del event

class VCEvent(object):
    def __init__(self,sys):
        self.eid            = None
        self.year           = None
        self.trigger_id     = None
        self.magnitude      = None
        self.shear_init     = None
        self.normal_init    = None
        self.shear_final    = None
        self.normal_final   = None
        self.esweep_ids     = None
        self.eas_ids        = None
        self.sys            = sys
        #super(VCEvent,self).__init__()
    
    def ruptureLength(self, surface_only = True):
        ruptured_elements = {}
        
        rupture_length = 0.0
        
        for sweep_id in self.esweep_ids:
            eid = self.sys.events.sweeps[sweep_id].seid
            try:
                tmp = ruptured_elements[eid]
            except KeyError:
                if surface_only and self.sys.geometry.elements[eid].onTrace():
                    rupture_length += self.sys.geometry.elements[eid].sideLength()
                ruptured_elements[eid] = True
        
        return rupture_length
        
    def ruptureArea(self):
        ruptured_elements = {}
        
        rupture_area = 0.0
        
        for sweep_id in self.esweep_ids:
            eid = self.sys.events.sweeps[sweep_id].seid
            try:
                tmp = ruptured_elements[eid]
            except KeyError:
                rupture_area += self.sys.geometry.elements[eid].sideLength() ** 2.0
                ruptured_elements[eid] = True
        
        return rupture_area
    
    def slipMap(self):
        ruptured_elements = {}
        
        for sweep_id in self.esweep_ids:
            # ele_info is:
            # {ele_id, total_slip, shear_init, normal_init, shear_final, normal_final, number_of_slips}
            #total_slip += VCEvents.sweeps[sweep_id].slip
            eid = self.sys.events.sweeps[sweep_id].seid
            try:
                ruptured_elements[eid]['total_slip'] += self.sys.events.sweeps[sweep_id].slip
                ruptured_elements[eid]['shear_final'] = self.sys.events.sweeps[sweep_id].shear_final
                ruptured_elements[eid]['normal_final'] = self.sys.events.sweeps[sweep_id].normal_final
                ruptured_elements[eid]['number_of_slips'] += 1
            except KeyError:
                ruptured_elements[eid] = {  'ele_id'                :eid,
                                            'total_slip'            :self.sys.events.sweeps[sweep_id].slip,
                                            'shear_init'            :self.sys.events.sweeps[sweep_id].shear_init,
                                            'normal_init'           :self.sys.events.sweeps[sweep_id].normal_init,
                                            'shear_final'           :self.sys.events.sweeps[sweep_id].shear_final,
                                            'normal_final'          :self.sys.events.sweeps[sweep_id].normal_final,
                                            'number_of_slips'       :1,
                                            'sweep_of_first_slip'   :self.sys.events.sweeps[sweep_id].snum
                                            }
        return ruptured_elements
                
    def averageSlip(self):
        ruptured_elements = {}
        
        total_slip = 0.0
        
        for sweep_id in self.esweep_ids:
            total_slip += self.sys.events.sweeps[sweep_id].slip
            eid = self.sys.events.sweeps[sweep_id].seid
            try:
                tmp = ruptured_elements[eid]
            except KeyError:
                ruptured_elements[eid] = True
        
        return total_slip/float(len(ruptured_elements))
    
    def minMaxSlip(self):
        slips = sorted([self.sys.events.sweeps[sweep_id].slip for sweep_id in self.esweep_ids])
        
        return (slips[0], slips[-1])
            
    def elementInEvent(self, eid):
        for sweep_id in self.esweep_ids:
            if self.sys.events.sweeps[sweep_id].seid == eid:
                return True
        
        return False
    
    def exportKML(self,f,events,min_year,max_year,fade_steps):
        
        allowed_year_range = self.sys.events.max_allowed_year - self.sys.events.start_year
        full_year_range = max_year - min_year
        if full_year_range > allowed_year_range:
            modified_full_year_range = math.modf(full_year_range/allowed_year_range)[1] * allowed_year_range
        else:
            modified_full_year_range = full_year_range
        curr_year = self.year - min_year + self.sys.events.start_year
        
        if curr_year < self.sys.events.start_year + modified_full_year_range:
            
            time_step = modified_full_year_range/500.0
            
            if curr_year > self.sys.events.max_allowed_year:
                mod_curr_year = curr_year
                while mod_curr_year > self.sys.events.max_allowed_year:
                    mod_curr_year -= allowed_year_range
            else:
                mod_curr_year = curr_year
            
            f.write('<Folder>\n')
            f.write('       <name>Event %i on Year %f (%f)</name>\n'%(self.eid,curr_year,self.year))

            f.write('       <Style>\n')
            f.write('               <ListStyle>\n')
            f.write('                       <listItemType>checkHideChildren</listItemType>\n')
            f.write('               </ListStyle>\n')
            f.write('       </Style>\n')
            
            involved_elements = {}
            for sid in self.esweep_ids:
                try:
                    involved_elements[self.sys.events.sweeps[sid].seid] += self.sys.events.sweeps[sid].slip
                except KeyError:
                    involved_elements[self.sys.events.sweeps[sid].seid] = self.sys.events.sweeps[sid].slip
                    
            year, month, day = self.sys.geometry.converter.yearDecimalToYearMonthDay(mod_curr_year)
            year = math.modf(curr_year)[1]
            
            for eid,eslip in involved_elements.iteritems():
                    for num in range(0,fade_steps + 1):
                            
                        year_end, month_end, day_end = self.sys.geometry.converter.yearDecimalToYearMonthDay(mod_curr_year + time_step*float(num+1))
                        year_end = math.modf(curr_year + time_step*float(num+1))[1]
                        
                        self.sys.geometry.elements[eid].exportKML(f,'elementBreakStyle%i_%i'%(int(round(self.magnitude)),num),(year, month, day),(year_end, month_end, day_end))
                        
                        year, month, day = year_end, month_end, day_end
            
            f.write('</Folder>\n')
    
    def exportSweepKML(self, f, fade_steps):
        curr_year = 2011.0 + self.year 
        
        f.write('<Folder>\n')
        f.write('       <name>Event %i on Year %f (%f)</name>\n'%(self.eid,curr_year,self.year))

        f.write('       <Style>\n')
        f.write('               <ListStyle>\n')
        f.write('                       <listItemType>checkHideChildren</listItemType>\n')
        f.write('               </ListStyle>\n')
        f.write('       </Style>\n')
        
        year, month, day = self.sys.geometry.converter.yearDecimalToYearMonthDay(curr_year)
        
        min_slip, max_slip = self.minMaxSlip()
        
        end_year = year + len(self.esweep_ids) - 1
            
        for sweep in self.sweeps.itervalues():
            sweep.exportKML(f, events.geometry.m_elements, fade_steps, year, max_slip, min_slip, end_year)
        
        f.write('</Folder>\n')

        
    def addSweep(self, sid):
        if self.esweep_ids is None:
            self.esweep_ids = []
        self.esweep_ids.append(sid)
        self.sys.events.addSweep(sid)
  
    def __str__(self):
        pass

class VCSweep(object):
    def __init__(self,sys):
        self.sid            = None # this is an internal unique id
        self.snum           = None # this is the sweep num from the simulation. there can be multiple sids per snum
        self.seid           = None
        self.slip           = None
        self.shear_init     = None
        self.normal_init    = None
        self.shear_final    = None
        self.normal_final   = None
        self.sys            = sys
        #super(VCSweep,self).__init__()

    def __str__(self):
        pass
'''