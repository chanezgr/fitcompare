"""
FITCOMPARE by Grégory Chanez / nakan.ch
This program is intended to run in a Docker container

This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
# #############################
# IMPORT section

import fitparse
import argparse
import configparser
import re
import os
import sys
import math
import datetime
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib 
from scipy.signal import savgol_filter
import numpy as np
import csv
import json

# Import advanced HR analysis functions
from fitcompare_advanced import *

# #############################
# INIT section

# Init the seaborn/matplotlib for graphs
sns.set()

# Define CONST
SCRIPT_VER = "2.7.0"
# TODO: 
# - Clean the filtering method of HRV
# - Create a configuration line on the project.yaml to remove the gray dotted line on HR chart
# 
# CHANGELOG:
# 2.7.0: Add support for 5hz GPS on Garmin FIT files, import HRV from Suunto JSON files and skip specific timestamps
# 2.6.1: Add values in charts titles / Clean error if no HRV data / Remove "half last point" for SIGMA devices
# 2.6.0: Add hrvCsv option for fit files in order to provide a separate HRV CSV file (specific for Polar watches)
# 2.5.4: Add nktool_battery as standard charge field
# 2.5.3: Add "Track" as a GNSS/GPS mode
# 2.5.2: Add a command line parameter to change the project configuration file
# 2.5.1: Add the HR monitor name in HRV chart
# 2.5.0: Force re-align of data arrays in order to get good HR values (peak). Add a vertical line at max HR on heart_rate graph
# 2.4.4: Fix the data loading for heart_rate with values None (will be set to 0)
# 2.4.3: Add a function to remove abnormal HRV values and setting the percentage
# 2.4.2: Slight adjust of the HR score for some extreme cases 
# 2.4.1: Add the ability to set custom threshold for HRV peaks removal
# 2.4.0: Add a beginning of basic support of HRV data analyzis
# 2.3.0: Improve reliability of HR score with a rewrite of reduce_latency function
# 2.2.1: Add the ability to enter manually start and end battery level in yaml project file
# 2.2.0: Improve mapbox map output file. Split private functions in separate file to publish fitcompare core
# 2.1.10: If map is disabled in the config file, disable all the referecnes to the map in the summary
# 2.1.9: Battery burn rate is displayed in HH h MM instead of decimal

# Get the global configuration of the script
config = configparser.ConfigParser()
config.read('config.ini')
MAPBOX_API_KEY = config['map']['mapbox_api_key']

# This script is run in a container. Define the working directory (mounted dir)
APP_PATH = "/project/"

# Output the running dialog with version number at the very beginning
print("Running fitcompare v%s" % (SCRIPT_VER))

# Default configurations, can be overriden by project.yaml file:
project_conf_file_exists = False
project_conf_altitude_gap = 1
delta_values = {}
hrvCsv_values = {}
hrvSuunto_values = {}
hrvDelta_values = {}
project_conf_zoom = False
project_conf_zoom_range = [0, 0]
project_conf_ignore = []
values_to_compare = ['heart_rate', 'altitude', 'distance']
project_conf_map = True
project_conf_map_style = 'satellite-streets-v12'
project_conf_align = True
project_conf_remove_hrv_abnormal = False
project_conf_remove_hrv_abnormal_threshold = 20
custom_graphs_values = []
charge = {}
conf_has_custom_graphs = False
project_conf_inc_smoothed_alt = False

# #############################
# ARGS section

# Handle the arguments and basic help
parser = argparse.ArgumentParser(description='Compare two or more FIT files')
parser.add_argument('fitfilesarg', metavar='FITFILE', nargs='+', help='Fit Files to compare')
parser.add_argument('--reference-file', '-r', dest='reference_file', help='Set the reference FIT File')
parser.add_argument('--prefix', '-p', dest='project_prefix', help='Set the project prefix for output files')
parser.add_argument('--debug', '-d', action='store_true', help='Enable debug')
parser.add_argument('--export', '-e', action='store_true', help='Export graphs values also as CSV')
parser.add_argument('--config', '-c', dest='project_config', help='Use an alternative configuration YAML file')
parser.add_argument('--listfields', '-l', action='store_true', help='List all fields for FITFILE')
args = parser.parse_args()

# List fields:
config_list_fields = args.listfields
# If debug mode
if (args.debug): print("[debug] Enable debug mode")

# Define the configuration file:
project_conf_file = APP_PATH + '/project.yaml'
if args.project_config is not None:
  project_conf_file = APP_PATH + '/' + args.project_config
  if (args.debug): print("[debug] Project configuration file set to " + args.project_config)

# Define if there is a reference file, define vars accordingly
with_reference_file = False
if args.reference_file is not None:
  with_reference_file = True
  reference_file = args.reference_file
  if (args.debug): print("[debug] Reference file set to " + args.reference_file)

# If there is a prefix to the project
if args.project_prefix is not None:
  project_prefix = args.project_prefix
  if (args.debug): print("[debug] Project prefix set to " + args.project_prefix)
else: 
  project_prefix = ''
  if (args.debug): print("[debug] Project prefix is empty")
# Build the fit files list. If a reference file is given, the first element of the list is the reference
fitfiles = []
if (with_reference_file == True):
  fitfiles.append(args.reference_file)
for fitfile in args.fitfilesarg:
  fitfiles.append(fitfile)

# If fit files number is only 1, then no align of data
if (len(fitfiles) == 1):
  project_conf_align = False
  if (args.debug): print("[debug] Project contains only one file, align is disabled")
  
# #############################
# PROJECT CONFIG section
  
# Read the project configuration, if there is one, and override defaults
if (args.debug): print("[debug] Now, trying to open the configuration file" + project_conf_file)
if (os.path.isfile(project_conf_file)):
  project_conf_file_exists = True
  if (args.debug): print("[debug] A project configuration file is present")
  with open(project_conf_file, 'r') as conf_file:
    project_conf = yaml.safe_load(conf_file)
  # Align data (True or False)
  if ("align" in project_conf['project']):
    project_conf_align = project_conf['project']['align']
    if (args.debug): print("[debug] Read configuration file: 'align' value set to " + str(project_conf['project']['align']))
  # Include smoothed alt on graph plot (True or False)
  if ("includeSmoothedAlt" in project_conf['project']):
    project_conf_inc_smoothed_alt = project_conf['project']['includeSmoothedAlt']
    if (args.debug): print("[debug] Read configuration file: 'includeSmoothedAlt' value set to " + str(project_conf['project']['includeSmoothedAlt']))
  # Zoom to a certain part of the project (list of begining and end relative time)
  if ("zoom" in project_conf['project']):
    project_conf_zoom = True
    project_conf_zoom_range = project_conf['project']['zoom']
    if (args.debug): print("[debug] Read configuration file: 'zoom' value set to [%i, %i]" % (project_conf_zoom_range[0], project_conf_zoom_range[1]))
  # Ignore certains timestamps (relative, first datapoint is 0)
  if ("ignore" in project_conf['project']):
    project_conf_ignore = project_conf['project']['ignore']
    if (args.debug): print("[debug] Read configuration file: 'ignore' value contains points")
  # Gap altitude data (wait for x seconds to plot altitude)
  if ("altitudeGap" in project_conf['project']):
    project_conf_altitude_gap = project_conf['project']['altitudeGap']
    if (args.debug): print("[debug] Read configuration file: 'altitudeGap' value set to %i" % (project_conf['project']['altitudeGap']))
  # Generate or not map (True or False)
  if ("map" in project_conf['project']):
    project_conf_map = project_conf['project']['map']
    if (args.debug): print("[debug] Read configuration file: 'map' value set to " + str(project_conf['project']['map']))
  # Change map style (mapbox map style)
  if ("mapStyle" in project_conf['project']):
    project_conf_map_style = project_conf['project']['mapStyle']
    if (args.debug): print("[debug] Read configuration file: 'mapStyle' value set to " + project_conf['project']['mapStyle'])
  # List of values to graph
  if ("graphs" in project_conf['project']):
    values_to_compare = project_conf['project']['graphs']
    if (args.debug): print("[debug] Read configuration file: 'graphs' value set to " + str(project_conf['project']['graphs']))
  # Remove aberrant values for HRV (more than 3 seconds)
  if ("removeAbnormalHrv" in project_conf['project']):
    project_conf_remove_hrv_abnormal = project_conf['project']['removeAbnormalHrv']
    if (args.debug): print("[debug] Read configuration file: 'removeAbnormalHrv' value set to " + str(project_conf['project']['removeAbnormalHrv']))
  # Remove aberrant values for HRV (more than 3 seconds)
  if ("removeAbnormalHrvThreshold" in project_conf['project']):
    project_conf_remove_hrv_abnormal_threshold = project_conf['project']['removeAbnormalHrvThreshold']
    if (args.debug): print("[debug] Read configuration file: 'removeAbnormalHrvThreshold' value set to " + str(project_conf['project']['removeAbnormalHrvThreshold']))     
      
  # Generate a list of custom graphs fields:
  if (("customGraphs" in project_conf) and (len(project_conf['customGraphs']) > 0)):
    if (args.debug): print("[debug] Read configuration file: 'customGraphs' is present")
    conf_has_custom_graphs = True
    for custom_graph in project_conf['customGraphs']:
      for custom_value in custom_graph['values']:
        custom_graphs_values.append(custom_value['field'])
  # Configuration for each fit file
  for ffile in fitfiles:
    # If we have the file
    if ffile in project_conf:
      # Delta allow adjustement in alignement of data
      if "delta" in project_conf[ffile]:
        if (args.debug): print("[debug] Read configuration file: 'delta' value set to %i for file %s" % (project_conf[ffile]['delta'], ffile))
        delta_values[ffile] = project_conf[ffile]['delta']
      if "hrvCsv" in project_conf[ffile]:
        if (args.debug): print("[debug] Read configuration file: 'hrvCsv' value set to %s for file %s" % (project_conf[ffile]['hrvCsv'], ffile))
        hrvCsv_values[ffile] = project_conf[ffile]['hrvCsv']
      if "hrvSuunto" in project_conf[ffile]:
        if (args.debug): print("[debug] Read configuration file: 'hrvSuunto' value set to %s for file %s" % (project_conf[ffile]['hrvSuunto'], ffile))
        hrvSuunto_values[ffile] = project_conf[ffile]['hrvSuunto']
      if "hrvDelta" in project_conf[ffile]:
        if (args.debug): print("[debug] Read configuration file: 'hrvDelta' value set to %i for file %s" % (project_conf[ffile]['hrvDelta'], ffile))
        hrvDelta_values[ffile] = project_conf[ffile]['hrvDelta']
      if "charge" in project_conf[ffile]:
        if (args.debug): print("[debug] Read configuration file: 'charge' value set to %i/%1 for file %s" % (project_conf[ffile]['charge'][0], project_conf[ffile]['charge'][1], ffile))
        charge[ffile] = project_conf[ffile]['charge']
else:
  if (args.debug): print("[debug] No configuration file to open, continuing with defaults")
  
# #############################
# FUNCTIONS section

# This function transform a FIT timestamp to Python datetime object
def fit_ts_to_dt(timestamp):
  # Specific FIT epoch
  fit_epoch = datetime.datetime(1989, 12, 31, 0, 0, 0)
  # Convert to usable datetime
  datetime_fit = fit_epoch + datetime.timedelta(seconds=timestamp)
  return datetime_fit

# This function load fit sessions (one for single activity, multiple for multisport)
# Input: 
# - fitname (fit file name)
# - summary (array of summary data)
# Output: 
# - Array of sessions (0 -> sport / 1 -> start_time / 2 -> total_elapsed_time)
def loadFitSession(fitname, summary):
  global APP_PATH
  r_sessions = []
  data = fitparse.FitFile(APP_PATH + fitname)
  i = 0
  for session in data.get_messages("session"):
    r_session = []
    r_session.append(session.get_value('sport'))
    r_session.append(session.get_value('start_time'))
    r_session.append(session.get_value('total_elapsed_time'))
    duration_seconds = (session.get_value('start_time')-summary[4]).total_seconds()
    r_session.append(duration_seconds)
    # Add everything to the main return value
    r_sessions.append(r_session)
  return r_sessions

# This function load fit data in an array
# Input:
# - fitname (fit file name)
# - summary (array of summary data)
# - array of fields to add to the fit data (graphs and custom graphs) - IMPORTANT: prive a copy of array [:]
# Output: 
# - Array of data for each point, a dict for each fields
def loadFitData(fitname, summary, fields):
  global delta_values, project_conf_zoom, project_conf_zoom_range, project_conf_map, custom_graphs_values, APP_PATH, config_list_fields
  # By default we include the timestamp in the data collected
  fields.append('timestamp')
  
  # And we also add the custom graphs fields:
  if (len(custom_graphs_values) > 0):
    for field in custom_graphs_values:
      fields.append(field)
  
  delta = 0
  start_point = summary[4] + datetime.timedelta(0,project_conf_zoom_range[0])
  end_point = summary[4] + datetime.timedelta(0,project_conf_zoom_range[1])
  
  # Order list for special fields:
  priority_fields = {}
  priority_fields['altitude'] = ['enhanced_altitude', 'altitude']
  priority_fields['speed'] = ['enhanced_speed', 'speed']
  priority_fields['charge'] = ['nktool_battery']
  
  if fitname in delta_values:
    delta = delta_values[fitname]
    if (args.debug): print("[debug] [loadFitData] Delta value to apply for file %s: %i" % (fitname, delta))
  # Load in an array the following data:
  # We include the position if map is enabled
  if (project_conf_map):
    fields.append('position')
    gps5hz_data = load5hzGPS(fitname, delta)
    if (args.debug) and (gps5hz_data): print("[debug] [loadFitData] Fitfile %s has 5hz GPS points" % (fitname))
    if (args.debug) and (not gps5hz_data): print("[debug] [loadFitData] Fitfile %s has NO 5hz GPS points" % (fitname))
  # timestamp, heart_rate, altitude, distance, power
  data = fitparse.FitFile(APP_PATH + fitname)
  # For each point
  all_values = []
  altitude = 0
  
  # Keep previous value for rare cases where heart_rate contains None for one point
  hr_previous_value = 0
  
  i = 0
  all_file_values = []
  for record in data.get_messages('record'):
    # If we have no zoom, or in range
    record_timestamp = record.get_value('timestamp') + datetime.timedelta(0,delta)
    if ((project_conf_zoom == False) or ((record_timestamp >= start_point) and (record_timestamp <= end_point))):
      # New point 
      this_value = {}
      if ((i == 20) and (config_list_fields)):
        print("*********************************************************")
        print("Fields for file %s:" % (fitname))
        for record_data in record:
          print(" - %s" % (record_data.name))
        print("*********************************************************")
      for value in fields:
        # If value is timestamp, then add the delta
        if (value == "timestamp"):
          this_value[value] = record.get_value(value) + datetime.timedelta(0,delta)
        # If in prority, run trough all the possible fields
        elif (value in priority_fields):
          this_single_value = None
          for possible_field in priority_fields[value]:
            if (record.get_value(possible_field) != None):
              this_single_value = record.get_value(possible_field)
              break
          this_value[value] = this_single_value
        # If field is position, put an array instead of a single value:
        elif (value == 'position'):
          this_value[value] = []
          # Search if we have 5hz GPS position for this point:
          multiple_point = False
          for g5hz_pt in gps5hz_data:
            if (g5hz_pt['ts'] == record_timestamp):
              this_value[value].append({'lat': g5hz_pt['lat'], 'long': g5hz_pt['long']})
              multiple_point = True
              break
          if (multiple_point == False):
            this_value[value].append({'lat': record.get_value('position_lat'), 'long': record.get_value('position_long')})
        # No priority list, just put the value if not none
        else:
          if (value == 'heart_rate'):
            if (record.get_value(value) == None):
              this_value[value] = hr_previous_value
              if (args.debug): print("[debug] [loadFitData] NOTICE: A value 'None' was found in heart_rate loading data from file %s " % (fitname))
            else:
              this_value[value] = record.get_value(value)
              hr_previous_value = record.get_value(value)
          else:
            this_value[value] = record.get_value(value)
      all_values.append(this_value)
      i = i+1
  
  return all_values

# This function loads a hrv array from a CSV
def loadCsvHrv(csv_file, hrvDelta):
  global APP_PATH, project_conf_remove_hrv_abnormal, project_conf_remove_hrv_abnormal_threshold

  rrintervals = []
  last_value = 0
  with open(APP_PATH + csv_file, mode='r') as file:
    csv_reader = csv.reader(file)
    i = 0
    for row in csv_reader:
      i = i+1
      if ((i > 1) and (i > hrvDelta)):
        this_value = int(row[0])
        if (last_value == 0):
          hrv_percentage = 0
        else:
          hrv_percentage = abs(100-(this_value*100/last_value))
                  
        # The soft and percentage filter
        if ((last_value != 0 and project_conf_remove_hrv_abnormal) and (hrv_percentage > project_conf_remove_hrv_abnormal_threshold)):
          rrintervals.append(last_value)
        else:
          rrintervals.append(this_value)
          last_value = this_value
          
  return rrintervals

# This function loads a hrv array from a Suunto JSON
def loadSuuntoHrv(json_file, hrvDelta):
  global APP_PATH, project_conf_remove_hrv_abnormal, project_conf_remove_hrv_abnormal_threshold
  with open(APP_PATH + json_file, 'r', encoding='utf-8') as file:
    data = json.load(file)
    rr_values = data['DeviceLog']['R-R']['Data']

  rrintervals = []
  i = 0
  last_value = 0
  for rr in rr_values:
    i += 1
    if (i > hrvDelta):
      this_value = rr
      if (last_value == 0):
        hrv_percentage = 0
      else:
        hrv_percentage = abs(100-(this_value*100/last_value))

      # The soft and percentage filter
      if ((last_value != 0 and project_conf_remove_hrv_abnormal) and (hrv_percentage > project_conf_remove_hrv_abnormal_threshold)):
        rrintervals.append(last_value)
      else:
        rrintervals.append(this_value)
        last_value = this_value

  return rrintervals

# This function loads a hrv array from a FIT
def loadFitHrv(fitname, hrvDelta):
  global APP_PATH, project_conf_remove_hrv_abnormal, project_conf_remove_hrv_abnormal_threshold
  
  data = fitparse.FitFile(APP_PATH + fitname)
  rrintervals = []
  last_value = 0
  i = 0
  for record in data.get_messages('hrv'):
    for record_data in record:
      for RR_interval in record_data.value:
        if RR_interval is not None:
          i += 1
          if (i > hrvDelta):
            this_value = RR_interval*1000
            if (last_value == 0):
              hrv_percentage = 0
            else:
              hrv_percentage = abs(100-(this_value*100/last_value))
            
            # The soft and percentage filter
            if ((last_value != 0 and project_conf_remove_hrv_abnormal) and (hrv_percentage > project_conf_remove_hrv_abnormal_threshold)):
                rrintervals.append(last_value)
            else:
              rrintervals.append(this_value)
              last_value = this_value
            
  return rrintervals

# This function take a fit file name and "decode" all the values
# Fit file name is something like: 
#    MakeModel_HRSource_GNSSConfig_DistanceSensor.fit
# where:
#    MakeModel:  Describe the make and model of the watch. 
#    HRSource:   Describe the belt or sensor model. OHR means optical heart rate sensor of the watch
#    GNSSConfig: Describe the GNSS configuration for the recoding. This can be:
#                - SatIQ:    for Garmin SatIQ automatic
#                - Ultra:    for UltraTrack or equivalent
#                - GPS:      for GPS only
#                - GNSS:     for all GNSS systems
#                - GNSSDual: for all GNSS systems + GPS Multiband
#                - Track:    for running mode specific for running track
#    DistanceSensor: Describe the distance sensor brand and model used (if any)
# Return an array with labels for each part as follow:
# [MakeModel, HRSource, GNSSConfig, DistanceSensor]
def decodeFitName(fitname):
  fitnameparts = os.path.splitext(os.path.basename(fitname))[0].split("_")
  # Add space before each capital
  makemodel = fitnameparts[0]
  # HR Source
  if (fitnameparts[1] == "OHR"):
    hrsource = "Cardio optique au poignet"
  else:
    hrsource = fitnameparts[1]
  # GNSS Config:
  if (fitnameparts[2] == "GPS"):
    gnsssource = "GPS Seul"
  elif (fitnameparts[2] == "GNSS"):
    gnsssource = "Tous les systèmes GNSS"
  elif (fitnameparts[2] == "GNSSDual"):
    gnsssource = "Tous les systèmes GNSS et GPS Multibande"
  elif (fitnameparts[2] == "SatIQ"):
    gnsssource = "Mode GNSS et GPS Multibande automatique"
  elif (fitnameparts[2] == "Track"):
    gnsssource = "Profil spécifique pour la piste"
  elif (fitnameparts[2] == "NONE"):
    gnsssource = "Non spécifié"
  # Distance sensor
  if (len(fitnameparts) == 4):
    distancesource = fitnameparts[3]
  else: 
    distancesource = None
  return [makemodel, hrsource, gnsssource, distancesource]
  
# This function read the summary informations of a fit file
# Input:
# - fitname (fit file name)
# Output:
# - Complete array of summary data of the fit file
def fitSummary(fitname):
  global APP_PATH
  # Load fit file
  data = fitparse.FitFile(APP_PATH + fitname)
  # Get profile version
  profile_ver = data.profile_version
  # Get protocol version
  protocol_ver = data.protocol_version
  # Get manufacturer (create an exception for health fit export)
  try:
    manufacturer = data.messages[0].get_value('manufacturer')
  except:
    manufacturer = None
  # Get time created (create an exception for health fit export)
  try:
    time_created = data.messages[0].get_value('time_created')
  except:
    time_created = None
  
  # Defaults
  start_position_lat = 0
  start_position_long = 0
  total_ascent = 0
  total_descent = 0
  total_elapsed_time = 0
  total_moving_time = 0
  total_distance = 0
  sport = ''
  sub_sport = ''
  avg_lat = []
  avg_long = []
  
  # Get all the details from the session
  for session in data.get_messages("session"):
    for sess in session:
      if (sess.name == 'start_position_lat'):
        start_position_lat = sess.value
      if (sess.name == 'start_position_long'):
        start_position_long = sess.value
      if (sess.name == 'sport'):
        sport = sess.value
      if (sess.name == 'sub_sport'):
        sub_sport = sess.value
      if (sess.name == 'total_ascent'):
        total_ascent = sess.value
      if (sess.name == 'total_descent'):
        total_descent = sess.value
      if (sess.name == 'total_elapsed_time'):
        total_elapsed_time = sess.value
      if (sess.name == 'total_moving_time'):
        total_moving_time = sess.value
      if (sess.name == 'total_distance'):
        total_distance = sess.value
  
  # Iterate trough all record points
  i=0
  dist=0
  start_battery = 0
  end_battery = 0
  start_alt = 0
  alt = 0
  avg_lat_final = 0
  avg_long_final = 0
  
  for record in data.get_messages('record'):
    i+=1
    if i <= 1: 
      start_battery = record.get_value('nktool_battery')
      timestamp = record.get_value('timestamp')
    if ((i >= project_conf_altitude_gap) and (start_alt == 0)):
      if (record.get_value('enhanced_altitude') != None):
        start_alt = record.get_value('enhanced_altitude')
      elif (record.get_value('altitude') != None):
        start_alt = record.get_value('altitude')
    if (record.get_value('enhanced_altitude') != None):
      alt = record.get_value('enhanced_altitude')
    elif (record.get_value('altitude') != None):
      alt = record.get_value('altitude')
    dist = record.get_value('distance')
    end_battery = record.get_value('nktool_battery')
    if ((record.get_value('position_lat') != None) and (record.get_value('position_long') != None) and (project_conf_map == True)):
      avg_lat.append(record.get_value('position_lat'))
      avg_long.append(record.get_value('position_long'))
  if (project_conf_map): 
    avg_lat_final = (sum(avg_lat) / len(avg_lat)) * (180/pow(2,31))
    avg_long_final = (sum(avg_long) / len(avg_long)) * (180/pow(2,31))
  
  if (dist == None):
    dist = total_distance
  
  # If battery detail is manually entered into the project YAML file
  if (start_battery == None):
    try:
      start_battery = charge[fitname][0]
      end_battery = charge[fitname][1]
    except:
      start_battery = None
      end_battery = None  

  return_val = []
  return_val.append(profile_ver)        # 0 -> ANT/FIT profil version
  return_val.append(protocol_ver)       # 1 -> ANT/FIT protocol version
  return_val.append(manufacturer)       # 2 -> manufacturer
  return_val.append(time_created)       # 3 -> time created (of the session data)
  return_val.append(timestamp)          # 4 -> first timestamp in record data
  return_val.append(dist)               # 5 -> Total distance
  return_val.append(i)                  # 6 -> Number of records
  return_val.append(total_elapsed_time) # 7 -> Total time in session data
  return_val.append(total_moving_time)  # 8 -> Total moving time in session data
  return_val.append(sport)              # 9 -> Sport (major)
  return_val.append(sub_sport)          # 10 -> Sport (sub)
  return_val.append(total_ascent)       # 11 -> Total ascent (value of session)
  return_val.append(total_descent)      # 12 -> Total descent (value of session)
  return_val.append(start_battery)      # 13 -> Battery value at beginning of activity
  return_val.append(end_battery)        # 14 -> Battery value at end of activity
  return_val.append(avg_lat_final)      # 15 -> Average latitude to center the map
  return_val.append(avg_long_final)     # 16 -> Average longitude to center the map
  return_val.append(start_alt)          # 17 -> Altitude at beginning of activity
  return_val.append(alt)                # 18 -> Altitude at end of activity
  
  return return_val
    
# This function returns the smoothed data of altitude for a file
# Input:
# - file_data: array of value for a fit file
# Output:
# - array of smoothed data of the altitude values
def smoothAltitude(file_data):
  global project_conf_zoom_range, project_conf_altitude_gap
  # If the zoom window is smaller than 70 (regular smoothed data)
  if ((project_conf_zoom) and ((project_conf_zoom_range[1] - project_conf_zoom_range[0]) <= 70)):
    smooth_window = project_conf_zoom_range[1] - project_conf_zoom_range[0]
  else:
    smooth_window = 70
  a_alt = []
  i = 0
  for point in file_data:
    i += 1
    skipThisPoint = False
    # If we are still in the skip window for altitude
    if (i <= project_conf_altitude_gap):
      skipThisPoint = True
    if (skipThisPoint == False):
      a_alt.append(point['altitude'])
  smoothed_altitude = savgol_filter(a_alt, smooth_window, 3) # window size 51, polynomial order 3
  return smoothed_altitude.tolist()
  
# This function commpute D+ from smoothed altitude data
# Input: 
# - smoothed altitude data array
# Output:
# - a computed value (float) of gain of altitude
def normalizedAltGain(smoothed_alt_data):
  norm_alt_pos = 0
  i_alt = 9000
  for apoint in smoothed_alt_data:
    if (i_alt < apoint):
      norm_alt_pos += apoint - i_alt
    i_alt = apoint
  return norm_alt_pos

# This function commpute D+ from smoothed altitude data
# Input: 
# - smoothed altitude data array
# Output:
# - a computed value (float) of loss of altitude
def normalizedAltLoss(smoothed_alt_data):
  norm_alt_loss = 0
  i_alt = -800
  for apoint in smoothed_alt_data:
    if (i_alt > apoint):
      norm_alt_loss += i_alt - apoint
    i_alt = apoint
  return norm_alt_loss

# This function fill a data array to the given value
# Input: 
# - data: array of value for a fit file
# - lenght: lenght to fill the data array to
# - fields: array of fields to add to the fit data (graphs and custom graphs) - IMPORTANT: prive a copy of array [:]
def fillDataArray(data, lenght, fields):
  global project_conf_map, custom_graphs_values
  # Fill value
  fill_value = None
  # We include the position if maap is enabled
  if (project_conf_map):
    fields.append('position_lat')
    fields.append('position_long')
  
  # And we also add the custom graphs fields:
  if (len(custom_graphs_values) > 0):
    for field in custom_graphs_values:
      fields.append(field)
      
  this_lenght = len(data)
  to_target_lenght = lenght - this_lenght
  start_timestamp = data[this_lenght -1]['timestamp']
  for i in range(to_target_lenght):
    this_point = {}
    # Fill each point with null data
    this_timestamp = start_timestamp + datetime.timedelta(0,i+1)
    this_point['timestamp'] = this_timestamp
    for field in fields:
      this_point[field] = fill_value
    data.append(this_point)
  return data

# This function will read and load additionnal positions of 5hz record 
def load5hzGPS(fitname, delta):
  global APP_PATH
  # Load fit file
  fitfile = fitparse.FitFile(APP_PATH + fitname)
  
  gps5hz_data = []
  # Get all the "unknown_467" messages (GPS 5hz)
  for message in fitfile.get_messages('unknown_467'):
    values = message.get_values()
    timestamp = fit_ts_to_dt(values['unknown_253']) + datetime.timedelta(0,delta)
    latitudes = values['unknown_1']
    longitudes = values['unknown_2']
    if ((timestamp is not None) and (latitudes is not None) and (longitudes is not None)):
      gps5hz_data.append({'ts': timestamp, 'lat': latitudes, 'long': longitudes})
    
  return gps5hz_data


# #############################
# PROCESS section

# Iterate through the fit files, to store all the relevant informations into an array
i=0
ff_summary = {}
ff_data = {}
textOutput = []
max_nb_points = 0
for ffile in fitfiles:
  if (args.debug): print("[debug] Processing file %s" % (ffile))
  i+=1

  # Get the relevant details from the fit file content
  ff_summary[ffile] = fitSummary(ffile)
  summary = ff_summary[ffile]
  
  # Get the max number of points
  if (args.debug): print("[debug] File %s has %i points" % (ffile, summary[6]))
  if (max_nb_points < summary[6]):
    if (args.debug): print("[debug] File %s has yet the higher number of points" % (ffile))
    max_nb_points = summary[6]
  
  # Load data of fit file in array
  if (args.debug): print("[debug] Call loadFitData for file %s" % (ffile))
  ff_data[ffile] = loadFitData(ffile, summary, values_to_compare[:])

  if (args.debug): print("[debug] Call loadFitSession for file %s" % (ffile))
  ff_session = loadFitSession(ffile, summary)
  
  # If altitude in the graphs list, we compute smoothed alt for all devices as well ad normalized alt gain/loss
  if ("altitude" in values_to_compare):
    if (args.debug): print("[debug] Altitude is in the field list, so compute smoothed altitude for file %s" % (ffile))
    smoothed_altitude = smoothAltitude(ff_data[ffile])
    normalized_alt_gain = normalizedAltGain(smoothed_altitude)
    normalized_alt_loss = normalizedAltLoss(smoothed_altitude)
  
  # Is this the reference file?  
  thisIsReferenceFile = ""
  if ((with_reference_file) and (i==1)):
    thisIsReferenceFile = " (reference file)"
  # Get the relevant details from FIT file name
  if (args.debug): print("[debug] Call decodeFitName for file %s" % (ffile))
  fitfiletags = decodeFitName(ffile)

  # Compute battery rate
  if (args.debug): print("[debug] Check for battery values on file %s" % (ffile))
  battery_rate = None
  battery_projection = None
  if ((summary[13] != None) and (summary[14] != None) and (summary[7] != None)):
    battery_rate = ((summary[13] - summary[14])/(summary[7]/3600))
    if (battery_rate > 0):
      battery_projection = 100/battery_rate
      # Compute HH:MM for battery projection
      battery_projection_hours = int(battery_projection)
      battery_projection_minutes = (battery_projection*60) % 60

  # Display the summary of the fit file
  if (args.debug): print("[debug] Start summary output for file %s" % (ffile))
  textOutput.append("=========================================================================\n")
  textOutput.append("FIT FILE: " + os.path.basename(ffile) + thisIsReferenceFile + "\n")
  textOutput.append("-------------------------------------------------------------------------\n")
  textOutput.append(" Device used:                  " + fitfiletags[0] + "\n")
  textOutput.append(" HR Measurement:               " + fitfiletags[1] + "\n")
  textOutput.append(" GNSS Mode:                    " + fitfiletags[2] + "\n")
  if (fitfiletags[3] != None):
    textOutput.append(" Distance Measurement:         " + fitfiletags[3] + "\n")
  textOutput.append("-------------------------------------------------------------------------\n")
  textOutput.append(" FIT Profile version:          %.2f\n" % (summary[0]))
  textOutput.append(" FIT Protocol version:         %.2f\n" % (summary[1]))
  if (summary[2] != None):
    textOutput.append(" FIT Manufacturer:             %s\n" % (summary[2]))
  if (summary[3] != None):
    textOutput.append(" Creation timestamp   :        " + summary[3].strftime("%m/%d/%Y, %H:%M:%S") + "\n")
  textOutput.append(" First point timestamp:        " + summary[4].strftime("%m/%d/%Y, %H:%M:%S") + "\n")
  if (len(ff_session) == 1):
    if (summary[5] != None):
      textOutput.append(" Total distance:               %i\n" % (summary[5]))
    textOutput.append(" Total number of points:       %i  (%s)\n" % (summary[6], datetime.timedelta(seconds=summary[6])))
    textOutput.append(" Total elapsed time:           %.2f (%s)\n" % (summary[7], datetime.timedelta(seconds=int(summary[7]))))
    if (summary[8] != None):
      textOutput.append(" Total moving time:            %.2f (%s)\n" % (summary[8], datetime.timedelta(seconds=int(summary[8]))))
    textOutput.append(" Sport / Sub Sport:            %s / %s\n" % (summary[9], summary[10]))
    if ((summary[11] != None) and (summary[12] != None) and ("altitude" in values_to_compare)):
      textOutput.append(" Total ascent / descent:       %.2f / %.2f\n" % (summary[11], summary[12]))
      textOutput.append(" Normalized ascent / descent:  %.2f / %.2f\n" % (normalized_alt_gain, normalized_alt_loss))
  elif (len(ff_session) > 1):
    textOutput.append(" Multisession activity:\n")
    for sess_details in ff_session:
      textOutput.append(" --> Session type %s\n" % (sess_details[0]))
      textOutput.append("     Session start:            %s (%i)\n" % (sess_details[1].strftime("%m/%d/%Y, %H:%M:%S"), sess_details[3]))
      textOutput.append("     Session duration:         %.2f (%s)\n" % (sess_details[2], datetime.timedelta(seconds=int(sess_details[2]))))
  if ((summary[13] != None) and (summary[14] != None)):
    textOutput.append(" Battery level start / end:    %.2f / %.2f\n" % (summary[13], summary[14]))
    if (battery_projection != None):
      textOutput.append(" Battery burn rate:            %.2f%%/hr (projection: %02dh%02d)\n" % (battery_rate, battery_projection_hours, battery_projection_minutes))
  textOutput.append("=========================================================================\n\n")

if (args.debug): print("[debug] End of files processing and summary output")

# Display project values:
if (args.debug): print("[debug] Starting project summary output")
# Store current date/time
now = datetime.datetime.now()
now_str = now.strftime("%d/%m/%Y %H:%M:%S")

textOutput.append("=========================================================================\n")
textOutput.append(" PROJECT VALUES\n")
textOutput.append("-------------------------------------------------------------------------\n")
textOutput.append(" Script version:                     %s\n" % (SCRIPT_VER))
textOutput.append(" Python version:                     %i.%i.%i\n" % (sys.version_info[0], sys.version_info[1], sys.version_info[2]))
textOutput.append(" Date/time of execution:             %s\n" % (now_str))
textOutput.append(" Project file configuration exists:  %s\n" % (project_conf_file_exists))
textOutput.append(" Zoom on certain points:             %s\n" % (project_conf_zoom))
if (project_conf_zoom):
  textOutput.append(" Zoom from / to:                     %i / %i\n" % (project_conf_zoom_range[0], project_conf_zoom_range[1]))
textOutput.append("-------------------------------------------------------------------------\n")
for ffile in fitfiles:
  textOutput.append(" Configuration values for %s\n" % (ffile))
  if ffile in delta_values:
    textOutput.append("  Delta: %i\n" % (delta_values[ffile]))
  if ffile in charge:
    textOutput.append("  Charge value: %i -> %i\n" % (charge[ffile][0], charge[ffile][1]))
  if ffile in hrvCsv_values:
    textOutput.append("  HRV CSV File: %s\n" % (hrvCsv_values[ffile]))
textOutput.append("=========================================================================\n")

# Display the project output and write it to project_logfile.txt
if (args.debug): print("[debug] Writing complete output to logfile")
print("".join(textOutput))
if (project_prefix != ''):
  flog_file = APP_PATH + project_prefix + "_" + 'logfile.txt'
else:
  flog_file = APP_PATH + "logfile.txt"
flog = open(flog_file, "w")
flog.write("".join(textOutput))
flog.close()

# Build an array with all the timestamps of all fit files
# This is only needed when more than one fit file is analyzed
if (project_conf_align):
  if (args.debug): print("[debug] Align values configured: build an array of all common timestamps")
  all_timestamp = []
  i=0
  for ffile in fitfiles:
    this_timestamp = []
    # For each record
    for record in ff_data[ffile]:
      this_timestamp.append(record['timestamp'])
    all_timestamp.append(this_timestamp)

  # Then build a common_timestamps array
  common_timestamp = []
  i = 0
  for ts in all_timestamp[0]:
    # By default, the timestamp is considered as OK
    thisPoint = True
    # If this relative point is in ignore list
    if (i in project_conf_ignore):
      thisPoint = False
    
    # Loop over all files
    for filearray in all_timestamp:
      # Loop over all timestamps
      if ts not in filearray:
        # The ts is not common to all files
        thisPoint = False
    # If point were found in each file, add to the common array
    if (thisPoint == True):
      common_timestamp.append(ts)
    i += 1
  print(" Common timestamps:                  %i" % (len(common_timestamp)))
  # Now, remove all the timestamps in the fffiles arrays which are not in the common 
  if (args.debug): print("[debug] Align: removing all timestamps points not in the common list")
  for ffile in fitfiles:
    for record in ff_data[ffile][:]:
      if record['timestamp'] not in common_timestamp:
        ff_data[ffile].remove(record)

else:
  # If we don't align, we have to fill the shortest dataset to have the same amount of points
  # First get all the file timestamps array lengh:
  if (args.debug): print("[debug] Align values disabled")
  longest_ts_array = 0
  # If we have a zoom, then the longest is the window of the zoom:
  if project_conf_zoom:
    longest_ts_array = project_conf_zoom_range[1] - project_conf_zoom_range[0]
  # If no zoom, measure all file lenght:
  else: 
    for ffile in fitfiles:
      this_ffile_timestamps = []
      for record in ff_data[ffile]:
        this_ffile_timestamps.append(record['timestamp'])
      this_ffile_lenght = len(this_ffile_timestamps)
      if (args.debug): print("[debug] Before filling, %s file has %i points" % (ffile, this_ffile_lenght))
      if (this_ffile_lenght > longest_ts_array):
        longest_ts_array = this_ffile_lenght
  print(" Longest timestamps:                  %i" % (longest_ts_array))
  
  # Put all the files at the same lenght:
  for ffile in fitfiles:
    ff_data[ffile] = fillDataArray(ff_data[ffile], longest_ts_array, values_to_compare[:])
    if (args.debug): print("[debug] Filling file %s to %i points" % (ffile, longest_ts_array))
    this_ffile_timestamps = []
    for record in ff_data[ffile]:
      this_ffile_timestamps.append(record['timestamp'])
    this_ffile_lenght = len(this_ffile_timestamps)
    if (args.debug): print("[debug] After filling, %s file has %i points" % (ffile, this_ffile_lenght))
print("=========================================================================")

# ==============
# GRAPHS
def generateGraph(APP_PATH, project_prefix, compare_value, chartData, args, project_conf_align, chartTitle, hr_max_pos): 
  # Generate the graph
  pathlib.Path(APP_PATH + "pnggraphs").mkdir(exist_ok=True)
  if (project_prefix != ''):
    graph_file = APP_PATH + "pnggraphs/" + project_prefix + "_" + re.sub( '(?<!^)(?=[A-Z])', '_', compare_value ).lower()
  else:
    graph_file = APP_PATH + "pnggraphs/" + re.sub( '([A-Z])', r'+\1', compare_value ).lower()
  
  # Create a Pandas DataSet for this graph
  chartDataFrame = pd.DataFrame(chartData)
  # If the CSV export is enabled
  if (args.export):
    chartDataFrame.to_csv(graph_file + '.csv', sep=',', decimal='.')
  sns.set_theme(font='Montserrat')
  sns.set(rc = {'figure.figsize':(20, 10)})
  # Max number of points
  if (project_conf_align):
    global common_timestamp
    max_nb_points = len(common_timestamp)
  else:
    global longest_ts_array
    max_nb_points = longest_ts_array
  thisPlot, thisAx = sns.lineplot(x=None, y=None, data=chartDataFrame, linewidth=1, dashes=False).set(title=chartTitle, xlim=(-5,max_nb_points+5))
  plt.grid(True)
  
  # If heart_rate graph and analyzis, generate a vertical line at max heart rate
  if ((compare_value == "heart_rate") and project_conf_align):
    for mhp in hr_max_pos:
      plt.axvline(x=mhp, color='gray', linewidth=1, linestyle='dotted')
  
  plt.savefig(graph_file + '.png', bbox_inches='tight', pad_inches=0.3)
  plt.clf()

# Start with the generation of the comparaison data sets
shortest_hrv = 0
for compare_value in values_to_compare:
  
  if (args.debug): print("[debug] Configuring output for field %s" % (compare_value))
  # We print what we are doing
  print("Generating data for %s" % (compare_value))
  # Build a complete dataset
  if (compare_value == 'heart_rate'):
    chartTitle = "Analyse de la fréquence cardiaque (bpm)"
  elif (compare_value == 'altitude'):
    chartTitle = "Analyse de l'altitude (m)"
  elif (compare_value == 'distance'):
    chartTitle = "Analyse de l'accumulation de distance (m)"
  elif (compare_value == 'power'):
    chartTitle = "Analyse des données de puissance (W)"
  elif (compare_value == 'hrv'):
    chartTitle = "Analyse des données R-R (ms)"
  else:
    chartTitle = "Analyse du champ de données \"%s\"" % (compare_value)
 
  # Loop over the fitfiles
  #thisData = {}
  #thisLabel = {}
  chartData = {}
  chartData.clear()
  
  # ##############################
  # Generate all "standard" graphs
  # ##############################
  hr_max_pos = []
  for ffile in fitfiles:
    # Parse the fitfile
    # Build a dataset for this file
    #data = []
    thisPoint = 0
    a_values = []
    # Loop over all points
    i = 0
    
    # For the HR stats data
    average_hr_gap = {}
    average_hr_gap['average'] = []
    average_hr_gap['max'] = 0
    
    hr_analyze = False
    
    # For altitude stats
    start_alt = None
    end_alt = None

    for record in ff_data[ffile]:
      skipThisPoint = False
      # If we have a sync timestamps, check that this timestamp is in the list
      if (((project_conf_align) and (record['timestamp'] in common_timestamp)) or (project_conf_align == False)):
        i += 1
        
        # Special process

        # If the current field is HR
        if ((compare_value == 'heart_rate') and (record['heart_rate'] != None)):
          thisPoint = record['heart_rate']
        # If the current field is heart_rate, and we have a reference file, and we have at least two files:
        if ((compare_value == 'heart_rate') and (len(fitfiles) >= 2) and (with_reference_file)):
          # If this file is not the reference file
          if (ffile != reference_file):
            # Than we can compute the HR score, starting after one minute
            hr_analyze = True
            if ((i >= 60) and (record['heart_rate'] != None)):
              # Current bpm
              cur_bpm = record['heart_rate']
              # Get the HR value of the reference file for this timestamp
              average_hr_gap = bpm_new_point(cur_bpm, record['timestamp'], average_hr_gap, ff_data, reference_file, i)
        
        # If the current field is altitude
        if (compare_value == 'altitude'):
          if ((i == 1) or (start_alt == None)):
            if (record['altitude'] != None):
              start_alt = record['altitude']
          if (i <= project_conf_altitude_gap):
            skipThisPoint = True
          elif (record['altitude'] != None):
            thisPoint = record['altitude']
          if (record['altitude'] != None):
            end_alt = record['altitude']

        # All other fields
        if (record[compare_value] != None):
          thisPoint = record[compare_value]
        
        # Record this point
        if (skipThisPoint == False):
          a_values.append(thisPoint)
    # Get the ffile decode
    legend = decodeFitName(ffile)
    # Get the summary
    summary = ff_summary[ffile]
    
    # If we have altitude data, get the smoothed and normalized values
    if ("altitude" in values_to_compare):
      smoothed_altitude = smoothAltitude(ff_data[ffile])
      normalized_alt_gain = normalizedAltGain(smoothed_altitude)
      normalized_alt_loss = normalizedAltLoss(smoothed_altitude)
    
    # Compute the HR score
    if (hr_analyze):
      hr_adv_data = adv_hr_sum(average_hr_gap)  
      hr_max_pos.append(hr_adv_data['max_gap_position'])
    
    if (compare_value == 'heart_rate'):
      hr_summary = ''
      if (hr_analyze):
        hr_summary = " Ecart moyen: %.2f - Ecart max: %.2f - Score: %.1f%%" % (hr_adv_data['average_gap'], hr_adv_data['max_gap'], hr_adv_data['hr_score'])
      chart_legend = legend[0] + " (mesure cardio: " + legend[1] + ")" + hr_summary
    elif (compare_value == 'altitude'):
      chart_legend = "%s (D+: %.1f / D-: %.1f / Altitude de départ: %.1f / Altitude d'arrivée: %.1f" % (legend[0], normalized_alt_gain, normalized_alt_loss, summary[17], summary[18])
    elif (compare_value == 'distance'):
      if (summary[5] == None):
        summary[5] = 0
      if (legend[3] != None):
        chart_legend = "%s (Distance mesurée par: %s): %.2f m" % (legend[0], legend[3], summary[5])
      else:
        chart_legend = "%s: %.2f m" % (legend[0], summary[5])
    
    # This part for the HRV graph
    elif (compare_value == "hrv"):
      hrvDelta = 0
      chart_legend = "%s (%s)" % (legend[0], legend[1])
      if ffile in hrvDelta_values: 
        hrvDelta = hrvDelta_values[ffile]
      # If the current file has a HRV parameter
      if ffile in hrvCsv_values:
        a_values = loadCsvHrv(hrvCsv_values[ffile], hrvDelta)
      elif ffile in hrvSuunto_values:
        a_values = loadSuuntoHrv(hrvSuunto_values[ffile], hrvDelta)
      else:
        a_values = loadFitHrv(ffile, hrvDelta)
      if (args.debug): print("[debug] Number of HRV points for %s: %i" % (ffile, len(a_values)))
      # If we have 0 points, then rise error, it's not possible to go ahead with HRV...
      if (len(a_values) == 0): 
      	print("ERROR: No valid HRV data. Add a CSV or ensure HRV is correctly set in FIT file")
      	break
      if (shortest_hrv == 0) or (shortest_hrv > len(a_values)):
        shortest_hrv = len(a_values)
    else:
      chart_legend = "%s" % (legend[0])
    # Get the lengh of datatable to align next values
    graph_lengh = len(a_values)
    # Add the dataset to chart
    chartData[chart_legend] = a_values
    # If altitude graph AND include smoothed altitude
    if ((compare_value == 'altitude') and (project_conf_inc_smoothed_alt)):
      lengh_diff = len(smoothed_altitude) - graph_lengh
      smoothed_altitude_aligned = smoothed_altitude[lengh_diff:]
      # Add smoothed alt data to chart
      chartData['%s (smoothed altitude)' % (legend[0])] = smoothed_altitude_aligned
      
  # Check for HRV data lenght
  if (compare_value == "hrv"):
    i = 0
    for line in chartData:
      key, value = list(chartData.items())[i]
      if len(value) > shortest_hrv: 
        line = value[:shortest_hrv]
        chartData[key] = line
      i = i+1
  
  # Generate graph
  generateGraph(APP_PATH, project_prefix, compare_value, chartData, args, project_conf_align, chartTitle, hr_max_pos)
  
# ###############################
# Generate all the customs graphs
# ###############################
if (conf_has_custom_graphs):
  for cust_graph in project_conf['customGraphs']:
    chartData = {}
    chartData.clear()
    graph_name = cust_graph['name']
    chartTitle = graph_name
    print("Generating custom graph: %s" % (graph_name))
    for cg_value in cust_graph['values']:
      a_values = []
      for record in ff_data[cg_value['file']]:
        # If we have a sync timestamps, check that this timestamp is in the list
        if ((project_conf_align) and (record['timestamp'] in common_timestamp) or (project_conf_align == False)):
          i += 1
          a_values.append(record[cg_value['field']])
      legend = decodeFitName(cg_value['file'])
      chart_legend = "%s - %s" % (legend[0], cg_value['label'])
      chartData[chart_legend] = a_values

    # Generate the graph
    pathlib.Path(APP_PATH + "pnggraphs").mkdir(exist_ok=True)
    if (project_prefix != ''):
      graph_file = APP_PATH + "pnggraphs/" + project_prefix + "_" + graph_name.lower().replace(" ", "")
    else:
      graph_file = APP_PATH + "pnggraphs/" + graph_name.lower().replace(" ", "")
    # Create a Pandas DataSet for this graph
    chartDataFrame = pd.DataFrame(chartData)
    # If the CSV export is enabled
    if (args.export):
      chartDataFrame.to_csv(graph_file + '.csv', sep=',', decimal='.')
    sns.set_theme(font='Montserrat')
    sns.set(rc = {'figure.figsize':(20, 10)})
    # Max number of points
    if (project_conf_align):
      max_nb_points = len(common_timestamp)
    else:
      max_nb_points = longest_ts_array
    thisPlot, thisAx = sns.lineplot(x=None, y=None, data=chartDataFrame, linewidth=1, dashes=False).set(title=chartTitle, xlim=(-5,max_nb_points+5))
    plt.grid(True)
    
    plt.savefig(graph_file + '.png', bbox_inches='tight', pad_inches=0.3)
    plt.clf()
    
    
# Generate a GPS MAP
def generateMapboxMap(fitfiles, ff_data, project_prefix, MAPBOX_API_KEY, project_conf_map_style, ff_summary, APP_PATH):

  print("Generating map")
  gpx_data = {}
  gpx_colors = ['#0000ff', '#ff0000', '#00ff00', '#bf00ff', '#6e6e6e', '#D7DF01', '#A9BCF5', '#A9F5A9', '#F5A9A9', '#000000', '#01DFD7', '#F5A9E1', '#FF8000', '#08088A']
  for ffile in fitfiles:
    gpx_data[ffile] = []
    for point in ff_data[ffile]:
      # If it's a 5hz GPS point
      if (isinstance(point['position'][0]['long'], (tuple, list))):
          i = 0
          for gps_point in point['position'][0]['long']:
            if (gps_point is not None and point['position'][0]['lat'][i] is not None):
              gpx_data[ffile].append([gps_point * (180/pow(2,31)), point['position'][0]['lat'][i] * (180/pow(2,31))])
            i += 1
      # It's a standard point
      elif ((point['position'][0]['long'] != None) and (point['position'][0]['lat'] != None)):
        gpx_data[ffile].append([point['position'][0]['long'] * (180/pow(2,31)), point['position'][0]['lat'] * (180/pow(2,31))])
    start_lat = ff_summary[ffile][15]
    start_long = ff_summary[ffile][16]

  pathlib.Path(APP_PATH + "map").mkdir(exist_ok=True)
  if (project_prefix != ''):
    map_file = APP_PATH + "map/" + project_prefix + "_" + 'map.html'
  else:
    map_file = APP_PATH + "map/" + 'map.html'

  fmap = open(map_file, "w")
  fmap.write('<html lang="en">\n')
  fmap.write('<head>\n')
  fmap.write('<meta charset="utf-8">\n')
  fmap.write('<script src="https://unpkg.com/leaflet@1.7.1/dist/leaflet.js" integrity="sha512-XQoYMqMTK8LvdxXYG3nZ448hOEQiglfqkJs1NOQV44cWnUrBc8PkAOcXy20w0vlaXaVUearIOBhiXZ5V3ynxwA==" crossorigin=""></script>\n')
  fmap.write('<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet-gpx/1.3.1/gpx.min.js"></script>\n')
  fmap.write('<link href="https://api.mapbox.com/mapbox-gl-js/v2.2.0/mapbox-gl.css" rel="stylesheet">\n')
  fmap.write('<script src="https://api.mapbox.com/mapbox-gl-js/v2.2.0/mapbox-gl.js"></script>\n')
  fmap.write('<script src="https://api.mapbox.com/mapbox-gl-js/plugins/mapbox-gl-compare/v0.4.0/mapbox-gl-compare.js"></script>\n')
  fmap.write('<link href="https://fonts.googleapis.com/css?family=Montserrat" rel="stylesheet\n">')
  fmap.write('<style>\n')
  fmap.write('body {\n')
  fmap.write('    font-family: \'Montserrat\'; font-size: 16px;\n')
  fmap.write('}\n')
  fmap.write('.mapLegend {\n')
  fmap.write('    font-family: \'Montserrat\'; font-size: 16px;\n')
  fmap.write('}\n')
  fmap.write('</style>\n')
  fmap.write('<link rel="stylesheet" href="https://api.mapbox.com/mapbox-gl-js/plugins/mapbox-gl-compare/v0.4.0/mapbox-gl-compare.css" type="text/css"></head><body>\n')
  fmap.write('<br><br><div align="center" style="width: 1200px; height: 800px; padding-left: 30px;">\n')
  fmap.write('<div align="center" id="mapid" style="width: 100%; height: 680px;"></div>\n')
  fmap.write('<div align="left">\n')
  i=0
  for ffile in fitfiles:
    details = decodeFitName(ffile)
    fmap.write('<div class="mapLegend" align="left" style="margin-right: 8px; padding-left: 70px;"><font color="%s">&#9679;</font>%s (Mode GNSS: %s)</div>\n' % (gpx_colors[i], details[0], details[2]))
    i+=1
  fmap.write('</div>\n')
  fmap.write('<script>\n')
  fmap.write('mapboxgl.accessToken = \'%s\';\n' % (MAPBOX_API_KEY))
  fmap.write('var map = new mapboxgl.Map({\n')
  fmap.write('        container: \'mapid\',\n')
  fmap.write('        style: \'mapbox://styles/mapbox/%s\',\n' % (project_conf_map_style))
  fmap.write('        center: [%f, %f],\n' % (start_long, start_lat))
  fmap.write('        zoom: 13\n')
  fmap.write('});')

  fmap.write('    map.on(\'style.load\', () => {\n')
  fmap.write('        map.addSource(\'mapbox-dem\', {\n')
  fmap.write('            \'type\': \'raster-dem\',\n')
  fmap.write('            \'url\': \'mapbox://mapbox.mapbox-terrain-dem-v1\',\n')
  fmap.write('            \'tileSize\': 512,\n')
  fmap.write('            \'maxzoom\': 14\n')
  fmap.write('        });\n')
  fmap.write('        // add the DEM source as a terrain layer with exaggerated height\n')
  fmap.write('        map.setTerrain({ \'source\': \'mapbox-dem\', \'exaggeration\': 1.5 });\n')
  fmap.write('    });\n')
  fmap.write('map.addControl(new mapboxgl.FullscreenControl());')

  fmap.write('map.on(\'load\', function () {\n')

  i = 0
  for ffile in fitfiles:
    fmap.write('map.addSource(\'route%i\', {\n' % (i))
    fmap.write('\'type\': \'geojson\',\n')
    fmap.write('\'data\': {\n')
    fmap.write('\'type\': \'Feature\',\n')
    fmap.write('\'properties\': {},\n')
    fmap.write('\'geometry\': {\n')
    fmap.write('\'type\': \'LineString\',\n')
    fmap.write('\'coordinates\': [')

    for point in gpx_data[ffile]:
      fmap.write('[%f,%f],' % (point[0], point[1]))
    fmap.write(']}}});\n')
    fmap.write('map.addLayer({\n')
    fmap.write('\'id\': \'route%i\',\n' % (i))
    fmap.write('\'type\': \'line\',\n'),
    fmap.write('\'source\': \'route%i\',\n' % (i))
    fmap.write('\'layout\': {\n')
    fmap.write('\'line-join\': \'round\',\n')
    fmap.write('\'line-cap\': \'round\'\n')
    fmap.write('},\n')
    fmap.write('\'paint\': {\n')
    fmap.write('\'line-color\': \'%s\',\n' % (gpx_colors[i]))
    fmap.write('\'line-opacity\': 0.8,\n')
    fmap.write('\'line-width\': 4\n')
    fmap.write('}});\n')
    i += 1
  fmap.write('});\n')

  fmap.write('</script>')
  fmap.close()

# Generate Mapbox map if map is enabled
if (project_conf_map):
  generateMapboxMap(fitfiles, ff_data, project_prefix, MAPBOX_API_KEY, project_conf_map_style, ff_summary, APP_PATH)

# Generate an example config file
# Get all sessions if more than 1
sessions = []
times = {}
all_sessions = loadFitSession(fitfiles[0], ff_summary[fitfiles[0]])
if (len(all_sessions) > 1):
  i = 0
  all_sessions_start = []
  for session in all_sessions:
    session_start = []
    for ffile in fitfiles:
      f_sessions = loadFitSession(ffile, ff_summary[ffile])
      session_start.append(f_sessions[i][3])
    i += 1
    all_sessions_start.append(session_start)
    
fexconf = open(APP_PATH + "project.yaml.example", "w")
fexconf.write('project:\n')
fexconf.write('  align: False\n')
# If multisession:
if (len(all_sessions) > 1):
  i = 0
  for session in all_sessions_start:   
    this_start = np.amax(session)
    if (i >= 1):
      fexconf.write('  zoom: [%i, %i] # %s\n' % (old_start, this_start, all_sessions[i-1][0]))
    old_start = this_start
    i += 1
  fexconf.write('  zoom: [%i, %i] # %s\n' % (old_start, len(common_timestamp), all_sessions[i-1][0]))
else:
  # Not multisession, an arbitrary example
  fexconf.write('  zoom: [90, 120]\n')
fexconf.write('  altitudeGap: 8  # Seconds\n')
fexconf.write('  map: false\n')
fexconf.write('  mapStyle: outdoors-v12\n')
fexconf.write('  graphs: [\'heart_rate\', \'altitude\', \'distance\']\n')
fexconf.write('  includeSmoothedAlt: false\n')
fexconf.write('  removeAbnormalHrv: false\n')
fexconf.write('  removeAbnormalHrvThreshold: 20 # percentage of the previous value\n')
fexconf.write('customGraphs:\n')
fexconf.write('  - name: Altitude baro vs GPS\n')
fexconf.write('    values:\n')
fexconf.write('      - file: %s\n' % (fitfiles[0]))
fexconf.write('        field: enhanced_altitude\n')
fexconf.write('        label: Altitude baro\n')
fexconf.write('      - file: %s\n' % (fitfiles[0]))
fexconf.write('        field: GPS altitude\n')
fexconf.write('        label: Altitude GPS\n')
for ffile in fitfiles:
  fexconf.write('%s:\n' % (ffile))
  fexconf.write('  delta: 0\n')
  fexconf.write('#  charge: [99, 87]\n')
  fexconf.write('#  hrvCsv: polar_hrv.csv\n')
  fexconf.write('#  hrvCsv: polar_hrv.csv\n')
fexconf.close()
