#!/usr/bin/env python
from . import VCSys

class VCGeometry(VCSys):
    def __init__(self,sim_file):
        super(VCGeometry, self).__init__(sim_file)
        
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