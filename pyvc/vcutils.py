#!/usr/bin/env python
from mpl_toolkits.basemap import Basemap

import numpy as np

import sys
import math
import Queue
import multiprocessing
from operator import itemgetter
import time
import itertools
import os
import re
import subprocess

import tables

import quakelib

'''
from geographiclib.geodesic import Geodesic
import calendar
'''



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
# A class to handle parallel field calculations
#-------------------------------------------------------------------------------
class VCFieldProcessor(multiprocessing.Process):
    def __init__(self, work_queue, result_queue, field_1d, event_element_data, event_element_slips, lat_size, lon_size, cutoff, type='displacement'):#, min_lat, min_lon, max_lat, max_lon):
        
        # base class initialization
        #multiprocessing.Process.__init__(self)
        super(VCFieldProcessor,self).__init__()
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.type = type
        self.field_1d = field_1d
        self.event_element_data = event_element_data
        self.event_element_slips = event_element_slips
        self.lat_size = lat_size
        self.lon_size = lon_size
        self.cutoff = cutoff
        #self.event_center = event_center
        #self.event_radius = event_radius
    
        #self.counter = counter
        #self.total_tasks = total_tasks
	

    
    def run(self):
        while not self.kill_received:
            # get a task
            try:
                start, end = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            sys.stdout.write('{} - {}, '.format(start, end))
            sys.stdout.flush()
            #print 'processing elements {} - {}'.format(start, end)
            
            # create a element list
            elements = quakelib.EventElementList()
            
            # create elements and add them to the element list
            for element in self.event_element_data[start:end]:
                #print element
                ele = quakelib.EventElement4()
                ele.set_rake(element['rake_rad'])
                
                ele.set_slip(self.event_element_slips[element['block_id']])
                ele.set_vert(0, element['m_x_pt1'], element['m_y_pt1'], element['m_z_pt1'])
                ele.set_vert(1, element['m_x_pt2'], element['m_y_pt2'], element['m_z_pt2'])
                ele.set_vert(2, element['m_x_pt3'], element['m_y_pt3'], element['m_z_pt3'])
                ele.set_vert(3, element['m_x_pt4'], element['m_y_pt4'], element['m_z_pt4'])
                elements.append(ele)
            
            # create an event
            event = quakelib.Event()
            
            # add the elements to the event
            event.add_elements(elements)
            
            # lame params
            lame_lambda = 3.2e10
            lame_mu = 3.0e10
            
            if self.type == 'displacement':
                # calculate the displacements
                if self.cutoff is None:
                    disp_1d = event.event_displacements(self.field_1d, lame_lambda, lame_lambda)
                else:
                    disp_1d = event.event_displacements(self.field_1d, lame_lambda, lame_lambda, self.cutoff)
                disp = np.array(disp_1d).reshape((self.lat_size,self.lon_size))
                
                # empty arrays to store the results
                dX = np.empty((self.lat_size, self.lon_size))
                dY = np.empty((self.lat_size, self.lon_size))
                dZ = np.empty((self.lat_size, self.lon_size))
                
                it = np.nditer(dX, flags=['multi_index'])
                while not it.finished:
                    dX[it.multi_index] = disp[it.multi_index][0]
                    dY[it.multi_index] = disp[it.multi_index][1]
                    dZ[it.multi_index] = disp[it.multi_index][2]
                    it.iternext()
                
                # store the result
                processed_displacements = {}
                processed_displacements['dX'] = dX
                processed_displacements['dY'] = dY
                processed_displacements['dZ'] = dZ
                self.result_queue.put(processed_displacements)
            elif self.type == 'gravity':
                # calculate the gravity changes
                if self.cutoff is None:
                    dGrav_1d = event.event_gravity_changes(self.field_1d, lame_lambda, lame_mu)
                else:
                    dGrav_1d = event.event_gravity_changes(self.field_1d, lame_lambda, lame_mu, self.cutoff)
                dGrav = np.array(dGrav_1d).reshape((self.lat_size,self.lon_size))
                
                # store the result
                self.result_queue.put(dGrav)
            elif self.type == 'dilat_gravity':
                # calculate the gravity changes
                if self.cutoff is None:
                    dGrav_1d = event.event_dilat_gravity_changes(self.field_1d, lame_lambda, lame_mu)
                else:
                    dGrav_1d = event.event_dilat_gravity_changes(self.field_1d, lame_lambda, lame_mu, self.cutoff)
                dGrav = np.array(dGrav_1d).reshape((self.lat_size,self.lon_size))
                
                # store the result
                self.result_queue.put(dGrav)
                

#-------------------------------------------------------------------------------
# Parent class for all of the field calculations
#-------------------------------------------------------------------------------
class VCField(object):
    def __init__(self, min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding, map_res, map_proj):
        # These are constrained this way so we can plot on 1024x780 for the
        # animations
        max_map_width = 690.0
        max_map_height = 658.0
        #max_map_width = 309.0
        #max_map_height = 309.0
        
        # A conversion instance for doing the lat-lon to x-y conversions
        self.convert = quakelib.Conversion(base_lat, base_lon)
        
        # Calculate the lat-lon range based on the min-max and the padding
        lon_range = max_lon - min_lon
        lat_range = max_lat - min_lat
        max_range = max((lon_range, lat_range))
        self.min_lon = min_lon - lon_range*padding
        self.min_lat = min_lat - lat_range*padding
        self.max_lon = max_lon + lon_range*padding
        self.max_lat = max_lat + lat_range*padding
        
        # We need a map instance to calculate the aspect ratio
        map = Basemap(
            llcrnrlon=self.min_lon,
            llcrnrlat=self.min_lat,
            urcrnrlon=self.max_lon,
            urcrnrlat=self.max_lat,
            lat_0=(self.max_lat+self.min_lat)/2.0,
            lon_0=(self.max_lon+self.min_lon)/2.0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )

        # Using the aspect ratio (h/w) to find the actual map width and height
        # in pixels
        if map.aspect > max_map_height/max_map_width:
            map_height = max_map_height
            map_width = max_map_height/map.aspect
        else:
            map_width = max_map_width
            map_height = max_map_width*map.aspect
        
        #print map.aspect, map_width, map_height, max_map_height/max_map_width
        
        self.lons_1d = np.linspace(self.min_lon,self.max_lon,int(map_width))
        self.lats_1d = np.linspace(self.min_lat,self.max_lat,int(map_height))
        
        _lons_1d = quakelib.FloatList()
        _lats_1d = quakelib.FloatList()
        
        for lon in self.lons_1d:
            _lons_1d.append(lon)
        
        for lat in self.lats_1d:
            _lats_1d.append(lat)
        
        self.field_1d = self.convert.convertArray2xyz(_lats_1d,_lons_1d)

    def calculate_field_values(self, event_element_data, event_element_slips, cutoff, type='displacement'):
        #-----------------------------------------------------------------------
        # Break up the work and start the workers
        #-----------------------------------------------------------------------
        
        # How many elements in the event
        event_size = float(len(event_element_slips))
        
        # How many seperate CPUs do we have
        num_processes = available_cpu_count()
        
        sys.stdout.write('{} processors : '.format(num_processes))
        
        # Figure out how many segments we will break the task up into
        seg = int(round(event_size/float(num_processes)))
        if seg < 1:
            seg = 1
        
        # Break up the job.
        segmented_elements_indexes = []
        if event_size < num_processes:
            segments = int(event_size)
        else:
            segments = int(num_processes)
        
        for i in range(segments):
            if i == num_processes - 1:
                end_index = len(event_element_slips)
            else:
                end_index = seg*int(i + 1)
            start_index = int(i) * seg
            if start_index != end_index:
                segmented_elements_indexes.append((start_index, end_index))
        
        # Add all of the jobs to a work queue
        work_queue = multiprocessing.Queue()
        for job in segmented_elements_indexes:
            work_queue.put(job)
        
        # Create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()
        
        # Spawn workers
        for i in range(len(segmented_elements_indexes)):
            #print available_cpu_count()
            worker = VCFieldProcessor(work_queue,
                result_queue,
                self.field_1d,
                event_element_data,
                event_element_slips,
                self.lats_1d.size,
                self.lons_1d.size,
                cutoff,
                type=type                
            )
            worker.start()
        
        # Collect the results off the queue
        self.results = []
        for i in range(len(segmented_elements_indexes)):
            self.results.append(result_queue.get())

        #sys.stdout.write('\033[30C\r')
        #sys.stdout.flush()

#-------------------------------------------------------------------------------
# A class to handle calculating event gravity changes
#-------------------------------------------------------------------------------
class VCGravityField(VCField):
    def __init__(self, min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=0.01, map_res='i', map_proj='cyl',dilat=True):
        
        super(VCGravityField,self).__init__(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding, map_res, map_proj)
        
        # Define how the cutoff value scales if it is not explitly set
        self.cutoff_min_size = 20.0
        self.cutoff_min = 20.0
        self.cutoff_p2_size = 65.0
        self.cutoff_p2 = 90.0
        
        self.dG = None
        self.dG_min = sys.float_info.max
        self.dilat = dilat


    #---------------------------------------------------------------------------
    # Sets up the gravity change calculation and then passes it to the
    # VCFieldProcessor class.
    #---------------------------------------------------------------------------
    def calculate_field_values(self, event_element_data, event_element_slips, cutoff=None, save_file_prefix=None):
        
        #-----------------------------------------------------------------------
        # If the cutoff is none (ie not explicitly set) calculate the cutoff for
        # this event.
        #-----------------------------------------------------------------------
        event_size = float(len(event_element_slips))
        if cutoff is None:
            if  event_size >= self.cutoff_min_size:
                cutoff = linear_interp(
                    event_size,
                    self.cutoff_min_size,
                    self.cutoff_p2_size,
                    self.cutoff_min,
                    self.cutoff_p2
                    )
            else:
                cutoff = self.cutoff_min
    
        sys.stdout.write('{:0.2f} cutoff : '.format(cutoff))
        sys.stdout.flush()
        #-----------------------------------------------------------------------
        # Run the field calculation. The results are stored in self.results
        #-----------------------------------------------------------------------
        if self.dilat:
            super(VCGravityField,self).calculate_field_values(event_element_data, event_element_slips, cutoff, type='dilat_gravity')
        else:
            super(VCGravityField,self).calculate_field_values(event_element_data, event_element_slips, cutoff, type='gravity')
        
        #-----------------------------------------------------------------------
        # Combine the results
        #-----------------------------------------------------------------------
        for result in self.results:
            '''
            min = np.amin(np.fabs(result[result.nonzero()]))
            if min < self.dG_min:
                self.dG_min = min
            '''
            if self.dG is None:
                self.dG = result
            else:
                self.dG += result
        
        #-----------------------------------------------------------------------
        # If the save file is set then we need to combine the results to be
        # saved. This is done seperately from above because self.dG is a
        # cumulative result that could include contributions from multiple
        # events. We only want to save the results of a single calculation.
        #-----------------------------------------------------------------------
        if save_file_prefix is not None:
            dG = None
            for result in self.results:
                if dG is None:
                    dG = result
                else:
                    dG += result

                # New, not thoroughly tested
                if self.free_air:
                    np.save('{}dG_total.npy'.format(save_file_prefix), dG)
                else:
                    np.save('{}dG_dilat.npy'.format(save_file_prefix), dG)



    def init_field(self, value):
        self.dG = np.empty((self.lats_1d.size, self.lons_1d.size))
        self.dG.fill(value)

    def load_field_values(self, file_prefix, factor=1.0, subtract=False):
        if self.dG is None:
            self.init_field(0.0)
        
        try:
            # New, not thoroughly tested
            if self.free_air:
                dG = np.load('{}dG_total.npy'.format(file_prefix))
            else:
                dG = np.load('{}dG_dilat.npy'.format(file_prefix))

            '''
            min = np.amin(np.fabs(dG[dG.nonzero()]))
            if min < self.dG_min:
                self.dG_min = min
            '''
            
            if not subtract:
                self.dG += factor*dG
            else:
                self.dG -= factor*dG
            return True
        except IOError:
            return False
            
            
            
    def get_field_value(self,lat,lon):
        # returns gravity field value (in microgals)
        #              at nearest lat/lon grid point 
        delta_lat = self.lats_1d[7] - self.lats_1d[6]
        delta_lon = self.lons_1d[7] - self.lons_1d[6]
    
        i_lat = np.where(np.logical_and(lat-delta_lat*0.5<self.lats_1d,self.lats_1d<lat+delta_lat*0.5))[0][0]
        i_lon = np.where(np.logical_and(lon-delta_lon*0.5<self.lons_1d,self.lons_1d<lon+delta_lon*0.5))[0][0]
    
        #print "\nlat_in: {}, matched: {}".format(lat,self.lats_1d[i_lat])
        #print "\nlon_in: {}, matched: {}".format(lon,self.lons_1d[i_lon])
    
        return self.dG[i_lat][i_lon]*pow(10,8)
        
        
            
    def save_field_values(self,file_prefix):
        np.save('{}dG.npy'.format(save_file_prefix), self.dG)
        
        
        
    def save_lat_lon_values(self,file_prefix):
        np.save('{}lats.npy'.format(file_prefix), self.lats_1d)
        np.save('{}lons.npy'.format(file_prefix), self.lons_1d)



    def shrink_field(self, percentage):
        self.dG *= percentage
        #zeros = np.zeros(self.dG.shape)
        #self.dG = np.where(self.dG >= self.dG_min, self.dG, zeros)
    '''
    def __imul__(self, value):
        if self.dG is None:
            self.init_field(0.0)
        self.dG *= value

        return self
    '''

#-------------------------------------------------------------------------------
# A class to handle calculating event displacements
#-------------------------------------------------------------------------------
class VCDisplacementField(VCField):
    def __init__(self, min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=0.01, map_res='i', map_proj='cyl'):
        
        super(VCDisplacementField,self).__init__(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding, map_res, map_proj)
        
        # Define how the cutoff value scales if it is not explitly set
        self.cutoff_min_size = 20.0
        self.cutoff_min = 46.5
        self.cutoff_p2_size = 65.0
        self.cutoff_p2 = 90.0
    
        self.dX = None
        self.dY = None
        self.dZ = None
    
        self.dX_min = sys.float_info.max
        self.dY_min = sys.float_info.max
        self.dZ_min = sys.float_info.max
    
    #---------------------------------------------------------------------------
    # Sets up the displacement calculation and then passes it to the
    # VCFieldProcessor class.
    #---------------------------------------------------------------------------
    def calculate_field_values(self, event_element_data, event_element_slips, cutoff=None, save_file_prefix=None):
    
        #-----------------------------------------------------------------------
        # If the cutoff is none (ie not explicitly set) calculate the cutoff for
        # this event.
        #-----------------------------------------------------------------------
        event_size = float(len(event_element_slips))
        if cutoff is None:
            if  event_size >= self.cutoff_min_size:
                cutoff = linear_interp(
                    event_size,
                    self.cutoff_min_size,
                    self.cutoff_p2_size,
                    self.cutoff_min,
                    self.cutoff_p2
                    )
            else:
                cutoff = self.cutoff_min
    
        sys.stdout.write('{:0.2f} cutoff : '.format(cutoff))
        sys.stdout.flush()
        #-----------------------------------------------------------------------
        # Run the field calculation. The results are stored in self.results
        #-----------------------------------------------------------------------
        super(VCDisplacementField,self).calculate_field_values(event_element_data, event_element_slips, cutoff, type='displacement')
        
        #-----------------------------------------------------------------------
        # Combine the results
        #-----------------------------------------------------------------------
        for result in self.results:
            '''
            min_x = np.amin(np.fabs(result['dX'][result['dX'].nonzero()]))
            if min_x < self.dX_min:
                self.dX_min = min_x
            min_y = np.amin(np.fabs(result['dY'][result['dY'].nonzero()]))
            if min_y < self.dY_min:
                self.dY_min = min_y
            min_z = np.amin(np.fabs(result['dZ'][result['dZ'].nonzero()]))
            if min_z < self.dZ_min:
                self.dZ_min = min_z
            '''
            if self.dX is None:
                self.dX = result['dX']
            else:
                self.dX += result['dX']
            
            if self.dY is None:
                self.dY = result['dY']
            else:
                self.dY += result['dY']
                
            if self.dZ is None:
                self.dZ = result['dZ']
            else:
                self.dZ += result['dZ']
        
        #-----------------------------------------------------------------------
        # If the save file is set then we need to combine the results to be
        # saved. This is done seperately from above because self.dX etc. is a
        # cumulative result that could include contributions from multiple
        # events. We only want to save the results of a single calculation.
        #-----------------------------------------------------------------------
        if save_file_prefix is not None:
            dX = None
            dY = None
            dZ = None
            for result in self.results:
                if dX is None:
                    dX = result['dX']
                else:
                    dX += result['dX']
                
                if dY is None:
                    dY = result['dY']
                else:
                    dY += result['dY']
                    
                if dZ is None:
                    dZ = result['dZ']
                else:
                    dZ += result['dZ']

            np.save('{}dX.npy'.format(save_file_prefix), dX)
            np.save('{}dY.npy'.format(save_file_prefix), dY)
            np.save('{}dZ.npy'.format(save_file_prefix), dZ)
            


    def init_field(self, value):
        self.dX = np.empty((self.lats_1d.size, self.lons_1d.size))
        self.dY = np.empty((self.lats_1d.size, self.lons_1d.size))
        self.dZ = np.empty((self.lats_1d.size, self.lons_1d.size))
        
        self.dX.fill(value)
        self.dY.fill(value)
        self.dZ.fill(value)
        
        
    def save_lat_lon_values(self,file_prefix):
        np.save('{}lats.npy'.format(save_file_prefix), self.lats_1d)
        np.save('{}lons.npy'.format(save_file_prefix), self.lons_1d)
        
    
    def load_field_values(self, file_prefix):
        if self.dX is None or self.dY is None or self.dZ is None:
            self.init_field(0.0)
        
        try:
            dX = np.load('{}dX.npy'.format(file_prefix))
            dY = np.load('{}dY.npy'.format(file_prefix))
            dZ = np.load('{}dZ.npy'.format(file_prefix))
            '''
            min_x = np.amin(np.fabs(dX[dX.nonzero()]))
            if min_x < self.dX_min:
                self.dX_min = min_x
            min_y = np.amin(np.fabs(dY[dY.nonzero()]))
            if min_y < self.dY_min:
                self.dY_min = min_y
            min_z = np.amin(np.fabs(dZ[dZ.nonzero()]))
            if min_z < self.dZ_min:
                self.dZ_min = min_z
            '''
            self.dX += dX
            self.dY += dY
            self.dZ += dZ
            
            return True
        except IOError:
            return False
    
    def shrink_field(self, percentage):
        if self.dX is None or self.dY is None or self.dZ is None:
            self.init_field(0.0)
        
        self.dX *= percentage
        self.dY *= percentage
        self.dZ *= percentage
        
        #print percentage, self.dX_min, self.dY_min, self.dZ_min
        '''
        zeros = np.zeros(self.dX.shape)
        self.dX = np.where(self.dX >= self.dX_min, self.dX, zeros)
        self.dY = np.where(self.dY >= self.dY_min, self.dY, zeros)
        self.dZ = np.where(self.dZ >= self.dZ_min, self.dZ, zeros)
        '''
    
    '''
    def __imul__(self, value):
        if self.dX is None or self.dY is None or self.dZ is None:
            self.init_field(0.0)
        self.dX *= value
        self.dY *= value
        self.dZ *= value

        return self
    '''

        
