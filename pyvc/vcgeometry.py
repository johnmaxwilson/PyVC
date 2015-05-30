#!/usr/bin/env python
from . import VCSys
from . import vcexceptions
from operator import itemgetter
import sys

#-------------------------------------------------------------------------------
# A class representing the geometry of a Virtual California simulation.
#-------------------------------------------------------------------------------
class VCGeometry(VCSys):
    def __init__(self,sim_data):
        super(VCGeometry, self).__init__(sim_data)
        
        # the table we need from the simulation data
        self.geometry_data = self.sim_data.file.root.block_info_table
        self.model_extents = self.sim_data.file.root.model_extents
        self.base_lat_lon = self.sim_data.file.root.base_lat_lon
    
    @property
    def base_lat(self):
        return self.base_lat_lon[0]
    
    @property
    def base_lon(self):
        return self.base_lat_lon[1]
    
    @property
    def min_lat(self):
        return self.model_extents.col('min_lat')[0]
    
    @property
    def max_lat(self):
        return self.model_extents.col('max_lat')[0]
    
    @property
    def min_lon(self):
        return self.model_extents.col('min_lon')[0]
    
    @property
    def max_lon(self):
        return self.model_extents.col('max_lon')[0]
    
    @property
    def min_x(self):
        return self.model_extents.col('min_x')[0]
    
    @property
    def max_x(self):
        return self.model_extents.col('max_x')[0]
    
    @property
    def min_y(self):
        return self.model_extents.col('min_y')[0]
    
    @property
    def max_y(self):
        return self.model_extents.col('max_y')[0]
    
    @property
    def min_z(self):
        return self.model_extents.col('min_z')[0]
    
    @property
    def max_z(self):
        return self.model_extents.col('max_z')[0]
    
    def get_fault_traces(self):
        traces = {}
        for block in self.geometry_data.where('(m_trace_flag_pt1 != 0) & (m_trace_flag_pt4 != 0)'):
            try:
                traces[block['section_id']].append( (block['m_x_pt1'], block['m_y_pt1'], block['m_z_pt1']) )
                traces[block['section_id']].append( (block['m_x_pt4'], block['m_y_pt4'], block['m_z_pt4']) )
            except KeyError:
                traces[block['section_id']] = [ (block['m_x_pt1'], block['m_y_pt1'], block['m_z_pt1']) ]
                traces[block['section_id']].append( (block['m_x_pt4'], block['m_y_pt4'], block['m_z_pt4']) )
        return traces
    
    def sections_with_elements(self, elements):
        if len(elements) > 1:
            ele_getter = itemgetter(*elements)
            return set(ele_getter(self.geometry_data.read(field='section_id')))
        else:
            return [self.geometry_data[elements[0]]['section_id']]
    
    
#    def bounds_xyz_from_sections(self, sections):
#        min_x = min(element_data[0]['m_x_pt1'], element_data[0]['m_x_pt4'])
#        max_x = max(element_data[0]['m_x_pt1'], element_data[0]['m_x_pt4'])
#        min_y = min(element_data[0]['m_y_pt1'], element_data[0]['m_y_pt4'])
#        max_y = max(element_data[0]['m_y_pt1'], element_data[0]['m_y_pt4'])
#        if len(element_data) > 1:
#            for block in element_data:
#                if block['m_x_pt1'] < min_x: min_x = block['m_x_pt1']
#                if block['m_x_pt4'] < min_x: min_x = block['m_x_pt4']
#                if block['m_x_pt1'] > max_x: max_x = block['m_x_pt1']
#                if block['m_x_pt4'] > max_x: max_x = block['m_x_pt4']
#                if block['m_y_pt1'] < min_y: min_y = block['m_y_pt1']
#                if block['m_y_pt4'] < min_y: min_y = block['m_y_pt4']
#                if block['m_y_pt1'] > max_y: max_y = block['m_y_pt1']
#                if block['m_y_pt4'] > max_y: max_y = block['m_y_pt4']
#        if len(element_data) < 10:
#            min_x = min_x - (max_x-min_x)*.5
#            max_x = max_x + (max_x-min_x)*.5
#            min_y = min_y - (max_y-min_y)*.5
#            max_y = max_y + (max_y-min_y)*.5
#        return [min_x,max_x,min_y,max_y]
            
    
    
    def get_slip_rates(self,section_filter=None,per_year=True):
        rates = {}
        # Convert slip rates from meters/second to meters/(decimal year)
        CONVERSION = 3.15576*pow(10,7)
        
        if section_filter is not None:
            for secid in section_filter['filter']:
                for block in self.geometry_data.read_where('{type} == {value}'.format(type='section_id', value=secid)):
                    if per_year:
                        rates[int(block['block_id'])] = block['slip_velocity']*CONVERSION
                    else:
                        rates[int(block['block_id'])] = block['slip_velocity']
        else:
            for block in self.geometry_data:
                if per_year:
                    rates[int(block['block_id'])] = block['slip_velocity']*CONVERSION
                else:
                    rates[int(block['block_id'])] = block['slip_velocity']

        return rates
        
        
    def get_element_slip_rate(self,block_id,per_year=True):
        rates = {}
        # Convert slip rates from meters/second to meters/(decimal year)
        CONVERSION = 3.15576*pow(10,7)
        
        if per_year:
            return self.geometry_data[block_id]['slip_velocity']*CONVERSION
        else:
            return self.geometry_data[block_id]['slip_velocity']
        
        
        
    def get_slip_time_series(self,events_in_range,event_element_slips,DT=0.1,start_year=0.0,duration=100.0,section_filter=None):
        from numpy import arange
        # event_element_slips = dictionary indexed by event_id with entries being dictionaries of slips indexed by block_id
        # slip_time_series    = dictionary indexed by block_id with entries being arrays of absolute slip at each time step
       
       
        # Convert slip rates from meters/second to meters/(decimal year)
        CONVERSION = 3.15576*pow(10,7)
        # DT = 0.1yr evaluates field every 36.5 days
        
        # Get slip rates for the elements
        slip_rates = self.get_slip_rates(section_filter=section_filter)
        

        #Initialize blocks with 0.0 slip at time t=0.0
        slip_time_series  = {block_id:[0.0] for block_id in slip_rates.keys()}
        #event_time_series = [None for each in slip_time_series[slip_time_series.keys()[0]]]
    

        #Initialize time steps to evaluate slip    
        time_values = arange(start_year+DT,start_year+duration+DT,DT)

        if len(events_in_range['event_number']) == 0:
            raise vcexceptions.NoEventsFound(event_range={'type':'year', 'filter':(start_year,start_year+duration)}, magnitude_filter=None, section_filter=section_filter)
        else:
            for k in range(len(time_values)):
                if k>0:
                    # current time in simulation
                    right_now = time_values[k]
            
                    # back slip all elements by subtracting the slip_rate*dt
                    for block_id in slip_time_series.keys():
                        last_slip = slip_time_series[block_id][k-1]
                        this_slip = slip_rates[block_id]*CONVERSION*DT
                        slip_time_series[block_id].append(last_slip-this_slip)

                    # check if any elements slip as part of simulated event in the window of simulation time
                    # between (current time - DT, current time), add event slips to the slip at current time 
                    # for elements involved
                    for j in range(len(events_in_range['event_number'])):
                        evid    = events_in_range['event_number'][j]
                        ev_year = events_in_range['event_year'][j]
                        if right_now-DT < ev_year <= right_now:
                            for block_id in event_element_slips[evid].keys():
                                try:
                                    slip_time_series[block_id][k] += event_element_slips[evid][block_id]
                                except KeyError:
                                    # When faults within the specified section_filter have events that produce slip on elements
                                    #   outside the section_filter, the slip history for the external elements
                                    #   is added to the slip_time_series and only includes backslip plus any
                                    #   co-seismic slips from events selected with the current 
                                    #   section/magnitude_filter/event_range.
                                    slip_time_series[block_id] = [0.0]
                                    slip_rates[block_id] = self.get_element_slip_rate(block_id)
                                    
                                    for l in range(k):
                                        last_slip = slip_time_series[block_id][l]
                                        this_slip = slip_rates[block_id]*CONVERSION*DT
                                        slip_time_series[block_id].append(last_slip-this_slip)
                                
                                #if len(event_time_series[k]) == 1:
                                #    event_time_series[k]       = events_in_range['magnitude'][evid]
                                #else:
                                #    event_time_series[k].append(events_in_range['magnitude'][evid])
                            
            #return slip_time_series,event_time_series
            return slip_time_series



    def get_average_slip_time_series(self,events_in_range,event_element_slips,dt=0.1,start_year=0.0,duration=100.0,section_filter=None):
        from numpy import zeros

        slip_time_series  = self.get_slip_time_series(events_in_range,event_element_slips,DT=dt,start_year=start_year,duration=duration,section_filter=section_filter)

        BLOCK_IDS       = slip_time_series.keys()
        num_steps       = len(slip_time_series[BLOCK_IDS[0]])
        average_series  = zeros(num_steps)

        for k in range(num_steps):
            average_series[k]  = sum([slip_time_series[b][k] for b in BLOCK_IDS])
            average_series[k] /= float(len(BLOCK_IDS))


        return average_series

    
    
    def events_on_section(self, secid):
        try:
            return self.sim_data.file.root.events_by_section._v_children['section_{}'.format(secid)].read()
        except KeyError:
            raise vcexceptions.BadSectionID(secid)

    def get_section_name(self, secid):
        bis = self.geometry_data.read_where('{type} == {value}'.format(type='section_id', value=secid))
        return bis[0]['fault_name']

    def get_section_info(self, section_filter=None, section_id=None):
        if section_filter is not None:
            section_info = {}
            for section in section_filter['filter']:
                # get all of the blocks in the section
                bis = self.geometry_data.read_where('{type} == {value}'.format(type=section_filter['type'], value=section))
                # get the number of blocks along the strike and dip (add 1 cause
                # its zero based)
                bas = max(bis[:]['das_id']) + 1
                bad = max(bis[:]['depth_id']) + 1
                if section_filter['type'] == 'section_id':
                    sn = bis[0]['fault_name']
                    sid = section
                else:
                    sn = section
                    sid = bis[0]['section_id']
                section_info[sid] = {'name':sn, 'blocks_along_strike':bas, 'blocks_along_dip':bad}
            return section_info
        elif section_id is not None:
            return {}
        else:
            section_info = {}
            curr_das = None
            curr_sec = None
            for block in self.geometry_data:
                sec = block['section_id']
                das_id = block['das_id']
                depth_id = block['depth_id']
                if sec != curr_sec:
                    #new sec
                    curr_sec = sec
                    max_bad = 0
                    section_info[curr_sec] = {'name':block['fault_name'], 'blocks_along_strike':0, 'blocks_along_dip':max_bad}
                if das_id != curr_das:
                    #new das
                    section_info[curr_sec]['blocks_along_strike'] += 1
                    curr_das = das_id
                if depth_id > max_bad:
                    max_bad = depth_id
                    section_info[curr_sec]['blocks_along_dip'] = max_bad
            return section_info
    
    def __getitem__(self, bid):
        return self.geometry_data[bid]



            
    