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

--------------------------------------------------------------------------
Architecture (v3.0.0 refactor):
Every FIT file is parsed exactly once into a pandas DataFrame indexed by
timestamp (see FitFileData). All downstream steps - alignment, smoothing,
scores, chart data and CSV export - operate on those DataFrames instead of
lists of dicts and positional summary arrays.
Migration to version 3.X.X from older versions should not alter the results 
in any ways, but I might not have yet tested all the scenarios.
"""

import argparse
import configparser
import csv
import datetime
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import fitparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
import seaborn as sns
import yaml
from scipy.signal import savgol_filter

from fitcompare_advanced import compute_hr_score, compute_gps_score

# #############################
# CONSTANTS

SCRIPT_VER = "3.1.0"
# CHANGELOG:
# 3.1.0: Add --gpx-ref / -G: an absolute GPX trace used as GPS reference for the
#        FIT files. Each FIT GPS point is matched to the nearest GPX point, a
#        1 m tolerance is applied (retained gap = max(0, dist_cm - 100)), and a
#        GPS score (mirroring the HR score) plus per-point cm deviation and the
#        max deviation are reported and plotted.
# 3.0.0: Entire cide refactoring, replacing spaghetti code
#        Improvement: single FIT parse, full pandas-backed data model
#        Improvement: vectorized alignment/scores, structured functions. 
#        Improvement: D+/D- is now computed once on the aligned data (previously the log used pre-alignment values while the chart legend used post-alignment ones). 
#        Improvement: The map HTML is now built from a template + json.dumps (route coords rounded to 6 decimals / ~11cm) instead of ~90 write() calls; the rendered map is unchanged.
#        Improvement: example config file is not generated on each run anymore, but using --gen-config
#        Improvement: general code readability
# 2.7.2: Console output for hr score
# 2.7.1: Fix a bug if some altitude value are "None"
# 2.7.0: Add support for 5hz GPS on Garmin FIT files, import HRV from Suunto JSON files and skip specific timestamps

# FIT timestamps are seconds since this specific epoch
FIT_EPOCH = datetime.datetime(1989, 12, 31, 0, 0, 0)
# Semicircle -> degrees conversion factor
SEMICIRCLE_TO_DEG = 180 / pow(2, 31)
# This script runs in a container. Working directory (mounted dir).
APP_PATH = "/project/"
# Colors used for the map routes
GPX_COLORS = ['#0000ff', '#ff0000', '#00ff00', '#bf00ff', '#6e6e6e', '#D7DF01',
              '#A9BCF5', '#A9F5A9', '#F5A9A9', '#000000', '#01DFD7', '#F5A9E1',
              '#FF8000', '#08088A']
# Fields whose graph values should be rendered as integers (as the legacy code did)
INT_FIELDS = {'heart_rate', 'cadence'}
# Priority order for fields that have an "enhanced" variant
PRIORITY_FIELDS = {
  'altitude': ['enhanced_altitude', 'altitude'],
  'speed': ['enhanced_speed', 'speed'],
  'charge': ['nktool_battery'],
}


# #############################
# CONFIGURATION

@dataclass
class Config:
  """All runtime configuration: CLI arguments + project.yaml, with defaults."""
  # From CLI
  fitfiles: list = field(default_factory=list)
  reference_file: str = None
  with_reference_file: bool = False
  gpx_ref: str = None
  project_prefix: str = ''
  debug: bool = False
  export: bool = False
  list_fields: bool = False
  gen_config: bool = False
  # From project.yaml (with defaults)
  align: bool = True
  include_smoothed_alt: bool = False
  zoom: bool = False
  zoom_range: tuple = (0, 0)
  ignore: list = field(default_factory=list)
  altitude_gap: int = 1
  map: bool = True
  map_style: str = 'satellite-streets-v12'
  draw_gpx_ref: bool = False
  values_to_compare: list = field(default_factory=lambda: ['heart_rate', 'altitude', 'distance'])
  remove_hrv_abnormal: bool = False
  remove_hrv_abnormal_threshold: int = 20
  # Per-file options
  delta: dict = field(default_factory=dict)
  hrv_csv: dict = field(default_factory=dict)
  hrv_suunto: dict = field(default_factory=dict)
  hrv_delta: dict = field(default_factory=dict)
  charge: dict = field(default_factory=dict)
  # Custom graphs
  custom_graphs: list = field(default_factory=list)
  custom_graphs_values: list = field(default_factory=list)
  # Raw project config (for custom graph definitions) and existence flag
  raw: dict = field(default_factory=dict)
  conf_exists: bool = False

  def dbg(self, message):
    if self.debug:
      print("[debug] " + message)


def parse_args():
  parser = argparse.ArgumentParser(description='Compare two or more FIT files')
  parser.add_argument('fitfilesarg', metavar='FITFILE', nargs='+', help='Fit Files to compare')
  parser.add_argument('--reference-file', '-r', dest='reference_file', help='Set the reference FIT File')
  parser.add_argument('--gpx-ref', '-G', dest='gpx_ref',
                      help='GPX file used as absolute GPS reference for the FIT files')
  parser.add_argument('--prefix', '-p', dest='project_prefix', help='Set the project prefix for output files')
  parser.add_argument('--debug', '-d', action='store_true', help='Enable debug')
  parser.add_argument('--export', '-e', action='store_true', help='Export graphs values also as CSV')
  parser.add_argument('--config', '-c', dest='project_config', help='Use an alternative configuration YAML file')
  parser.add_argument('--listfields', '-l', action='store_true', help='List all fields for FITFILE')
  parser.add_argument('--gen-config', '-g', dest='gen_config', action='store_true',
                      help='Generate an example project.yaml.example from the FIT files, then exit')
  return parser.parse_args()


def load_config(args):
  """Build a Config from the CLI arguments and the optional project.yaml file."""
  cfg = Config()
  cfg.debug = args.debug
  cfg.export = args.export
  cfg.list_fields = args.listfields
  cfg.gen_config = args.gen_config
  cfg.reference_file = args.reference_file
  cfg.with_reference_file = args.reference_file is not None
  cfg.gpx_ref = args.gpx_ref
  if cfg.gpx_ref:
    cfg.dbg("GPX reference trace set to " + cfg.gpx_ref)
  cfg.project_prefix = args.project_prefix if args.project_prefix is not None else ''

  if cfg.debug:
    print("[debug] Enable debug mode")
  cfg.dbg("Project prefix set to '%s'" % cfg.project_prefix)

  # Build the fit files list (reference file first, if any)
  if cfg.with_reference_file:
    cfg.fitfiles.append(args.reference_file)
    cfg.dbg("Reference file set to " + args.reference_file)
  cfg.fitfiles.extend(args.fitfilesarg)

  # A single file cannot be aligned against anything
  if len(cfg.fitfiles) == 1:
    cfg.align = False
    cfg.dbg("Project contains only one file, align is disabled")

  # Locate the project configuration file
  conf_file = APP_PATH + '/project.yaml'
  if args.project_config is not None:
    conf_file = APP_PATH + '/' + args.project_config
    cfg.dbg("Project configuration file set to " + args.project_config)

  cfg.dbg("Now, trying to open the configuration file " + conf_file)
  if not os.path.isfile(conf_file):
    cfg.dbg("No configuration file to open, continuing with defaults")
    return cfg

  cfg.conf_exists = True
  cfg.dbg("A project configuration file is present")
  with open(conf_file, 'r') as handle:
    cfg.raw = yaml.safe_load(handle)
  _apply_project_section(cfg)
  _apply_custom_graphs(cfg)
  _apply_per_file_options(cfg)
  return cfg


def _apply_project_section(cfg):
  """Override defaults from the 'project' section of project.yaml."""
  project = cfg.raw.get('project', {})
  simple = {
    'align': 'align',
    'includeSmoothedAlt': 'include_smoothed_alt',
    'altitudeGap': 'altitude_gap',
    'map': 'map',
    'mapStyle': 'map_style',
    'drawGpxRef': 'draw_gpx_ref',
    'graphs': 'values_to_compare',
    'removeAbnormalHrv': 'remove_hrv_abnormal',
    'removeAbnormalHrvThreshold': 'remove_hrv_abnormal_threshold',
    'ignore': 'ignore',
  }
  for yaml_key, attr in simple.items():
    if yaml_key in project:
      setattr(cfg, attr, project[yaml_key])
      cfg.dbg("Read configuration file: '%s' value set to %s" % (yaml_key, project[yaml_key]))
  if 'zoom' in project:
    cfg.zoom = True
    cfg.zoom_range = project['zoom']
    cfg.dbg("Read configuration file: 'zoom' value set to [%i, %i]" % (cfg.zoom_range[0], cfg.zoom_range[1]))


def _apply_custom_graphs(cfg):
  custom = cfg.raw.get('customGraphs')
  if not custom:
    return
  cfg.dbg("Read configuration file: 'customGraphs' is present")
  cfg.custom_graphs = custom
  for custom_graph in custom:
    for custom_value in custom_graph['values']:
      cfg.custom_graphs_values.append(custom_value['field'])


def _apply_per_file_options(cfg):
  mapping = {
    'delta': cfg.delta,
    'hrvCsv': cfg.hrv_csv,
    'hrvSuunto': cfg.hrv_suunto,
    'hrvDelta': cfg.hrv_delta,
    'charge': cfg.charge,
  }
  for ffile in cfg.fitfiles:
    file_conf = cfg.raw.get(ffile)
    if not file_conf:
      continue
    for yaml_key, target in mapping.items():
      if yaml_key in file_conf:
        target[ffile] = file_conf[yaml_key]
        cfg.dbg("Read configuration file: '%s' set to %s for file %s" % (yaml_key, file_conf[yaml_key], ffile))


# #############################
# FIT FILE PARSING (single pass)

@dataclass
class FitSummary:
  """Named, self-documenting replacement for the old positional summary list."""
  profile_ver: float
  protocol_ver: float
  manufacturer: object
  time_created: object
  first_timestamp: object       # first record timestamp (without delta)
  total_distance: float
  n_records: int
  total_elapsed_time: float
  total_moving_time: float
  sport: str
  sub_sport: str
  total_ascent: float
  total_descent: float
  start_battery: float
  end_battery: float
  avg_lat: float
  avg_long: float
  start_alt: float
  end_alt: float


@dataclass
class FitFileData:
  """Everything extracted from one FIT file, parsed a single time."""
  name: str
  df: pd.DataFrame            # working records (zoom applied), indexed by timestamp
  summary: FitSummary
  sessions: list             # list of [sport, start_time, total_elapsed_time, duration_seconds]
  gps5hz: dict               # {timestamp: (lat_tuple, long_tuple)}


def fit_ts_to_dt(timestamp):
  """Convert a raw FIT timestamp (seconds since FIT epoch) to a datetime."""
  return FIT_EPOCH + datetime.timedelta(seconds=timestamp)


def _records_dataframe(fitfile, delta):
  """Build the per-record DataFrame with the derived priority columns."""
  rows = [msg.get_values() for msg in fitfile.get_messages('record')]
  df = pd.DataFrame(rows)
  if 'timestamp' not in df:
    return df, None
  raw_first_ts = df['timestamp'].iloc[0]
  # Apply the per-file delta to every timestamp, then index by it
  df['timestamp'] = df['timestamp'] + datetime.timedelta(0, delta)

  # Collapse the "enhanced" variants into the canonical field name
  for target, candidates in PRIORITY_FIELDS.items():
    present = [c for c in candidates if c in df.columns]
    if present:
      combined = df[present[0]]
      for extra in present[1:]:
        combined = combined.combine_first(df[extra])
      df[target] = combined
    elif target not in df.columns:
      df[target] = np.nan

  df = df.set_index('timestamp', drop=False).sort_index()
  return df, raw_first_ts


def _load_5hz_gps(fitfile, delta):
  """Read the optional 5hz GPS records (Garmin 'unknown_467' messages)."""
  gps5hz = {}
  for message in fitfile.get_messages('unknown_467'):
    values = message.get_values()
    timestamp = fit_ts_to_dt(values['unknown_253']) + datetime.timedelta(0, delta)
    latitudes = values['unknown_1']
    longitudes = values['unknown_2']
    if timestamp is not None and latitudes is not None and longitudes is not None:
      gps5hz[timestamp] = (latitudes, longitudes)
  return gps5hz


def _load_sessions(fitfile, first_timestamp):
  """Sessions: [sport, start_time, total_elapsed_time, seconds_since_first_point]."""
  sessions = []
  for session in fitfile.get_messages('session'):
    start_time = session.get_value('start_time')
    duration = (start_time - first_timestamp).total_seconds()
    sessions.append([
      session.get_value('sport'),
      start_time,
      session.get_value('total_elapsed_time'),
      duration,
    ])
  return sessions


def _build_summary(fitfile, df, raw_first_ts, cfg, ffile):
  """Compute the FitSummary from the full (pre-zoom) record DataFrame + sessions."""
  # Session-level aggregates (legacy kept the values of the last session message)
  sess = {'sport': '', 'sub_sport': '', 'total_ascent': 0, 'total_descent': 0,
          'total_elapsed_time': 0, 'total_moving_time': 0, 'total_distance': 0}
  for session in fitfile.get_messages('session'):
    for name in list(sess.keys()):
      value = session.get_value(name)
      if value is not None:
        sess[name] = value

  try:
    manufacturer = fitfile.messages[0].get_value('manufacturer')
  except (IndexError, KeyError):
    manufacturer = None
  try:
    time_created = fitfile.messages[0].get_value('time_created')
  except (IndexError, KeyError):
    time_created = None

  altitude = df['altitude'] if 'altitude' in df else pd.Series(dtype=float)
  # Start altitude: first available altitude at/after the altitude gap
  start_alt = 0
  after_gap = altitude.iloc[cfg.altitude_gap - 1:].dropna() if len(altitude) else altitude
  if len(after_gap):
    start_alt = after_gap.iloc[0]
  end_alt = altitude.dropna().iloc[-1] if altitude.dropna().size else 0

  # Total distance: last record value, or session total if absent
  total_distance = df['distance'].iloc[-1] if 'distance' in df and len(df) else None
  if total_distance is None or (isinstance(total_distance, float) and pd.isna(total_distance)):
    total_distance = sess['total_distance']

  # Battery: first / last record value, or manual charge config
  start_battery = df['charge'].iloc[0] if 'charge' in df and len(df) else None
  end_battery = df['charge'].iloc[-1] if 'charge' in df and len(df) else None
  if start_battery is not None and pd.isna(start_battery):
    start_battery = None
  if end_battery is not None and pd.isna(end_battery):
    end_battery = None
  if start_battery is None:
    try:
      start_battery = cfg.charge[ffile][0]
      end_battery = cfg.charge[ffile][1]
    except (KeyError, IndexError):
      start_battery = None
      end_battery = None

  # Map center: average of all GPS points
  avg_lat = avg_long = 0
  if cfg.map and 'position_lat' in df and 'position_long' in df:
    lat = df['position_lat'].dropna()
    lon = df['position_long'].dropna()
    if len(lat) and len(lon):
      avg_lat = lat.mean() * SEMICIRCLE_TO_DEG
      avg_long = lon.mean() * SEMICIRCLE_TO_DEG

  return FitSummary(
    profile_ver=fitfile.profile_version,
    protocol_ver=fitfile.protocol_version,
    manufacturer=manufacturer,
    time_created=time_created,
    first_timestamp=raw_first_ts,
    total_distance=total_distance,
    n_records=len(df),
    total_elapsed_time=sess['total_elapsed_time'],
    total_moving_time=sess['total_moving_time'],
    sport=sess['sport'],
    sub_sport=sess['sub_sport'],
    total_ascent=sess['total_ascent'],
    total_descent=sess['total_descent'],
    start_battery=start_battery,
    end_battery=end_battery,
    avg_lat=avg_lat,
    avg_long=avg_long,
    start_alt=start_alt,
    end_alt=end_alt,
  )


def load_fit(ffile, cfg):
  """Parse a FIT file a single time and return a FitFileData."""
  delta = cfg.delta.get(ffile, 0)
  if delta:
    cfg.dbg("[loadFitData] Delta value to apply for file %s: %i" % (ffile, delta))

  fitfile = fitparse.FitFile(APP_PATH + ffile)
  df, raw_first_ts = _records_dataframe(fitfile, delta)

  if cfg.list_fields and len(df) > 20:
    _print_fields(ffile, fitfile)

  summary = _build_summary(fitfile, df, raw_first_ts, cfg, ffile)
  sessions = _load_sessions(fitfile, raw_first_ts)
  gps5hz = _load_5hz_gps(fitfile, delta) if cfg.map else {}
  if cfg.map:
    cfg.dbg("[loadFitData] Fitfile %s has %s 5hz GPS points" % (ffile, "some" if gps5hz else "NO"))

  # Zoom: keep only the records inside the [start, end] window
  if cfg.zoom:
    start = raw_first_ts + datetime.timedelta(0, cfg.zoom_range[0])
    end = raw_first_ts + datetime.timedelta(0, cfg.zoom_range[1])
    df = df[(df['timestamp'] >= start) & (df['timestamp'] <= end)]

  return FitFileData(name=ffile, df=df, summary=summary, sessions=sessions, gps5hz=gps5hz)


def _print_fields(ffile, fitfile):
  """Reproduce the --listfields output (fields of the 21st record)."""
  records = list(fitfile.get_messages('record'))
  print("*********************************************************")
  print("Fields for file %s:" % ffile)
  for record_data in records[20]:
    print(" - %s" % record_data.name)
  print("*********************************************************")


# #############################
# FIT FILE NAME DECODING

def decode_fit_name(fitname):
  """Decode MakeModel_HRSource_GNSSConfig[_DistanceSensor].fit into readable labels."""
  parts = os.path.splitext(os.path.basename(fitname))[0].split("_")
  makemodel = parts[0]
  hrsource = "Cardio optique au poignet" if parts[1] == "OHR" else parts[1]
  gnss_labels = {
    "GPS": "GPS Seul",
    "GNSS": "Tous les systèmes GNSS",
    "GNSSDual": "Tous les systèmes GNSS et GPS Multibande",
    "SatIQ": "Mode GNSS et GPS Multibande automatique",
    "Track": "Profil spécifique pour la piste",
    "NONE": "Non spécifié",
  }
  gnsssource = gnss_labels.get(parts[2])
  distancesource = parts[3] if len(parts) == 4 else None
  return [makemodel, hrsource, gnsssource, distancesource]


# #############################
# ALTITUDE HELPERS

def smooth_altitude(alt_values, cfg):
  """Savitzky-Golay smoothing of an altitude series (already gap-trimmed)."""
  if cfg.zoom and (cfg.zoom_range[1] - cfg.zoom_range[0]) <= 70:
    window = cfg.zoom_range[1] - cfg.zoom_range[0]
  else:
    window = 70
  return savgol_filter(alt_values, window, 3).tolist()


def normalized_alt_gain(smoothed):
  """Cumulative positive elevation change (D+) of a smoothed altitude series."""
  diff = np.diff(np.asarray(smoothed, dtype=float))
  return float(diff[diff > 0].sum())


def normalized_alt_loss(smoothed):
  """Cumulative negative elevation change (D-) of a smoothed altitude series."""
  diff = np.diff(np.asarray(smoothed, dtype=float))
  return float(-diff[diff < 0].sum())


# #############################
# ALIGNMENT

@dataclass
class AlignedData:
  """Per-file working DataFrames plus the point count used for the x axis."""
  frames: dict               # {ffile: DataFrame ordered on the common points}
  max_points: int
  common_timestamps: list    # only meaningful when aligned


def align_files(fitdatas, cfg):
  """Align the files on their common timestamps, or pad them to equal length."""
  if cfg.align:
    return _align_common(fitdatas, cfg)
  return _align_padded(fitdatas, cfg)


def _align_common(fitdatas, cfg):
  """Keep only the timestamps present in every file (minus the ignored ones)."""
  index_sets = [set(fd.df.index) for fd in fitdatas.values()]
  common_all = set.intersection(*index_sets)
  first_index = list(next(iter(fitdatas.values())).df.index)
  common = [ts for pos, ts in enumerate(first_index)
            if ts in common_all and pos not in cfg.ignore]

  frames = {name: fd.df.reindex(common) for name, fd in fitdatas.items()}
  return AlignedData(frames=frames, max_points=len(common), common_timestamps=common)


def _align_padded(fitdatas, cfg):
  """Non-aligned mode: pad every file to the longest one (or the zoom window)."""
  if cfg.zoom:
    longest = cfg.zoom_range[1] - cfg.zoom_range[0]
  else:
    longest = max(len(fd.df) for fd in fitdatas.values())
  # Frames are returned untrimmed; padding to `longest` happens per-field at
  # value-extraction time (a field is padded by repeating its last value).
  frames = {name: fd.df for name, fd in fitdatas.items()}
  return AlignedData(frames=frames, max_points=longest, common_timestamps=[])


def field_values(frame, field_name, cfg, max_points):
  """Return the graph values of one field for one file (fill + carry + pad)."""
  if field_name in frame.columns:
    series = frame[field_name]
  else:
    series = pd.Series(np.nan, index=frame.index)
  # Carry the last known value forward, leading gaps become 0 (legacy thisPoint)
  filled = series.ffill().fillna(0)
  if field_name in INT_FIELDS:
    filled = filled.astype(int)
  values = filled.tolist()
  # Non-aligned mode pads to the longest file by repeating the last value
  if not cfg.align and len(values) < max_points:
    pad = values[-1] if values else 0
    values = values + [pad] * (max_points - len(values))
  return values


# #############################
# HRV LOADERS

def _filter_hrv(values, cfg, delta):
  """Shared HRV cleaning: skip the first `delta` points and drop abnormal spikes."""
  rrintervals = []
  last_value = 0
  for i, value in enumerate(values, start=1):
    if i <= delta:
      continue
    percentage = 0 if last_value == 0 else abs(100 - (value * 100 / last_value))
    if last_value != 0 and cfg.remove_hrv_abnormal and percentage > cfg.remove_hrv_abnormal_threshold:
      rrintervals.append(last_value)
    else:
      rrintervals.append(value)
      last_value = value
  return rrintervals


def load_csv_hrv(csv_file, cfg, delta):
  with open(APP_PATH + csv_file, mode='r') as handle:
    rows = list(csv.reader(handle))
  # Legacy skipped the header row (i > 1) as well as the first `delta` points
  values = [int(row[0]) for row in rows[1:]]
  return _filter_hrv(values, cfg, delta)


def load_suunto_hrv(json_file, cfg, delta):
  with open(APP_PATH + json_file, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
  values = data['DeviceLog']['R-R']['Data']
  return _filter_hrv(values, cfg, delta)


def load_fit_hrv(ffile, cfg, delta):
  fitfile = fitparse.FitFile(APP_PATH + ffile)
  values = []
  for record in fitfile.get_messages('hrv'):
    for record_data in record:
      for rr_interval in record_data.value:
        if rr_interval is not None:
          values.append(rr_interval * 1000)
  return _filter_hrv(values, cfg, delta)


# #############################
# GPS REFERENCE (GPX)

def load_gpx(gpx_file):
  """Parse a GPX file and return its track points as a list of (lat_deg, lon_deg)."""
  tree = ET.parse(APP_PATH + gpx_file)
  root = tree.getroot()
  # GPX elements are namespaced; iterate namespace-agnostically on 'trkpt'.
  ns = ''
  if root.tag.startswith('{'):
    ns = root.tag.split('}')[0] + '}'
  coords = []
  for pt in root.iter(ns + 'trkpt'):
    lat = pt.get('lat')
    lon = pt.get('lon')
    if lat is not None and lon is not None:
      coords.append((float(lat), float(lon)))
  return coords


def _fit_gps_coords(frame):
  """Return the per-record (lat_deg, lon_deg) list of one aligned file's GPS fix.

  Points without a fix are (nan, nan) so the per-point gap array stays aligned
  with the other graphs' x-axis. Only the 1 Hz record fixes are used here (the
  5 Hz Garmin bursts are a map-only expansion).
  """
  coords = []
  if 'position_lat' in frame.columns and 'position_long' in frame.columns:
    for lat, lon in zip(frame['position_lat'].tolist(), frame['position_long'].tolist()):
      if pd.notna(lat) and pd.notna(lon):
        coords.append((lat * SEMICIRCLE_TO_DEG, lon * SEMICIRCLE_TO_DEG))
      else:
        coords.append((float('nan'), float('nan')))
  return coords


def compute_gps_scores(fitdatas, aligned, cfg, gpx_coords):
  """Compute the GPS score of every FIT file against the GPX reference trace."""
  scores = {}
  for ffile in fitdatas:
    fit_coords = _fit_gps_coords(aligned.frames[ffile])
    score = compute_gps_score(fit_coords, gpx_coords)
    if score is None:
      print("GPS Score: %s has no usable GPS fix (skipped)" % ffile)
      continue
    scores[ffile] = score
    print("GPS Score: %s - Ecart moyen: %.1f cm - Ecart max: %.1f cm - Score: %.1f%%"
          % (ffile, score['average_gap'], score['max_gap'], score['gps_score']))
  return scores


def generate_gps_graph(gps_scores, aligned, fitdatas, cfg):
  """Plot the per-point retained GPS deviation (cm) of each file vs the GPX ref."""
  print("Generating data for gps")
  chart_data = {}
  vlines = []
  for ffile in fitdatas:
    score = gps_scores.get(ffile)
    if score is None:
      continue
    gaps = list(score['gaps'])
    if len(gaps) < aligned.max_points:
      gaps = gaps + [float('nan')] * (aligned.max_points - len(gaps))
    tags = decode_fit_name(ffile)
    legend = "%s (Ecart moyen: %.1f cm - Ecart max: %.1f cm - Score GPS: %.1f%%)" % (
      tags[0], score['average_gap'], score['max_gap'], score['gps_score'])
    chart_data[legend] = gaps
    vlines.append(score['max_gap_position'])
  _render_chart(chart_data, "Analyse de l'ecart GPS vs trace de reference (cm)",
                _graph_path('gps', cfg), cfg, aligned.max_points, vlines=vlines)


# #############################
# TEXT REPORT

def _battery_projection(summary):
  """Return (rate, hours, minutes) battery burn projection, or (None, None, None)."""
  if summary.start_battery is None or summary.end_battery is None or summary.total_elapsed_time is None:
    return None, None, None
  rate = (summary.start_battery - summary.end_battery) / (summary.total_elapsed_time / 3600)
  if rate <= 0:
    return rate, None, None
  projection = 100 / rate
  return rate, int(projection), (projection * 60) % 60


def build_report(fitdatas, cfg, alt_norms, gps_scores=None):
  """Build the human-readable text report (identical layout to the legacy output)."""
  out = []
  for idx, (ffile, fd) in enumerate(fitdatas.items(), start=1):
    s = fd.summary
    tags = decode_fit_name(ffile)
    is_ref = " (reference file)" if (cfg.with_reference_file and idx == 1) else ""
    rate, proj_h, proj_m = _battery_projection(s)

    out.append("=========================================================================\n")
    out.append("FIT FILE: " + os.path.basename(ffile) + is_ref + "\n")
    out.append("-------------------------------------------------------------------------\n")
    out.append(" Device used:                  " + tags[0] + "\n")
    out.append(" HR Measurement:               " + tags[1] + "\n")
    out.append(" GNSS Mode:                    " + tags[2] + "\n")
    if tags[3] is not None:
      out.append(" Distance Measurement:         " + tags[3] + "\n")
    out.append("-------------------------------------------------------------------------\n")
    out.append(" FIT Profile version:          %.2f\n" % s.profile_ver)
    out.append(" FIT Protocol version:         %.2f\n" % s.protocol_ver)
    if s.manufacturer is not None:
      out.append(" FIT Manufacturer:             %s\n" % s.manufacturer)
    if s.time_created is not None:
      out.append(" Creation timestamp   :        " + s.time_created.strftime("%m/%d/%Y, %H:%M:%S") + "\n")
    out.append(" First point timestamp:        " + s.first_timestamp.strftime("%m/%d/%Y, %H:%M:%S") + "\n")
    if len(fd.sessions) == 1:
      if s.total_distance is not None:
        out.append(" Total distance:               %i\n" % s.total_distance)
      out.append(" Total number of points:       %i  (%s)\n" % (s.n_records, datetime.timedelta(seconds=s.n_records)))
      out.append(" Total elapsed time:           %.2f (%s)\n" % (s.total_elapsed_time, datetime.timedelta(seconds=int(s.total_elapsed_time))))
      if s.total_moving_time is not None:
        out.append(" Total moving time:            %.2f (%s)\n" % (s.total_moving_time, datetime.timedelta(seconds=int(s.total_moving_time))))
      out.append(" Sport / Sub Sport:            %s / %s\n" % (s.sport, s.sub_sport))
      if s.total_ascent is not None and s.total_descent is not None and "altitude" in cfg.values_to_compare:
        out.append(" Total ascent / descent:       %.2f / %.2f\n" % (s.total_ascent, s.total_descent))
        gain, loss = alt_norms.get(ffile, (0, 0))
        out.append(" Normalized ascent / descent:  %.2f / %.2f\n" % (gain, loss))
    elif len(fd.sessions) > 1:
      out.append(" Multisession activity:\n")
      for sess in fd.sessions:
        out.append(" --> Session type %s\n" % sess[0])
        out.append("     Session start:            %s (%i)\n" % (sess[1].strftime("%m/%d/%Y, %H:%M:%S"), sess[3]))
        out.append("     Session duration:         %.2f (%s)\n" % (sess[2], datetime.timedelta(seconds=int(sess[2]))))
    if s.start_battery is not None and s.end_battery is not None:
      out.append(" Battery level start / end:    %.2f / %.2f\n" % (s.start_battery, s.end_battery))
      if rate is not None and proj_h is not None:
        out.append(" Battery burn rate:            %.2f%%/hr (projection: %02dh%02d)\n" % (rate, proj_h, proj_m))
    if gps_scores and ffile in gps_scores:
      gs = gps_scores[ffile]
      out.append(" GPS score (vs GPX ref):       %.1f%%\n" % gs['gps_score'])
      out.append(" GPS ecart moyen / max:        %.1f cm / %.1f cm @ point %i\n"
                 % (gs['average_gap'], gs['max_gap'], gs['max_gap_position']))
    out.append("=========================================================================\n\n")

  out.extend(_build_project_report(fitdatas, cfg))
  return out


def _build_project_report(fitdatas, cfg):
  now = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
  out = []
  out.append("=========================================================================\n")
  out.append(" PROJECT VALUES\n")
  out.append("-------------------------------------------------------------------------\n")
  out.append(" Script version:                     %s\n" % SCRIPT_VER)
  out.append(" Python version:                     %i.%i.%i\n" % (sys.version_info[0], sys.version_info[1], sys.version_info[2]))
  out.append(" Date/time of execution:             %s\n" % now)
  out.append(" Project file configuration exists:  %s\n" % cfg.conf_exists)
  out.append(" Zoom on certain points:             %s\n" % cfg.zoom)
  if cfg.zoom:
    out.append(" Zoom from / to:                     %i / %i\n" % (cfg.zoom_range[0], cfg.zoom_range[1]))
  out.append(" GPX reference trace:                %s\n" % (cfg.gpx_ref if cfg.gpx_ref else "none"))
  out.append("-------------------------------------------------------------------------\n")
  for ffile in fitdatas:
    out.append(" Configuration values for %s\n" % ffile)
    if ffile in cfg.delta:
      out.append("  Delta: %i\n" % cfg.delta[ffile])
    if ffile in cfg.charge:
      out.append("  Charge value: %i -> %i\n" % (cfg.charge[ffile][0], cfg.charge[ffile][1]))
    if ffile in cfg.hrv_csv:
      out.append("  HRV CSV File: %s\n" % cfg.hrv_csv[ffile])
  out.append("=========================================================================\n")
  return out


def write_logfile(text, cfg):
  print("".join(text))
  name = (cfg.project_prefix + "_logfile.txt") if cfg.project_prefix else "logfile.txt"
  with open(APP_PATH + name, "w") as handle:
    handle.write("".join(text))


# #############################
# GRAPHS

def _graph_path(compare_value, cfg):
  pathlib.Path(APP_PATH + "pnggraphs").mkdir(exist_ok=True)
  if cfg.project_prefix:
    snake = re.sub('(?<!^)(?=[A-Z])', '_', compare_value).lower()
    return APP_PATH + "pnggraphs/" + cfg.project_prefix + "_" + snake
  return APP_PATH + "pnggraphs/" + re.sub('([A-Z])', r'+\1', compare_value).lower()


def _render_chart(chart_data, chart_title, graph_file, cfg, max_points, vlines=None):
  """Turn a {legend: values} dict into a PNG (and optionally a CSV)."""
  frame = pd.DataFrame(chart_data)
  if cfg.export:
    frame.to_csv(graph_file + '.csv', sep=',', decimal='.')
  sns.set_theme(font='Montserrat')
  sns.set(rc={'figure.figsize': (20, 10)})
  sns.lineplot(x=None, y=None, data=frame, linewidth=1, dashes=False).set(
    title=chart_title, xlim=(-5, max_points + 5))
  plt.grid(True)
  for x in (vlines or []):
    plt.axvline(x=x, color='gray', linewidth=1, linestyle='dotted')
  plt.savefig(graph_file + '.png', bbox_inches='tight', pad_inches=0.3)
  plt.clf()


CHART_TITLES = {
  'heart_rate': "Analyse de la fréquence cardiaque (bpm)",
  'altitude': "Analyse de l'altitude (m)",
  'distance': "Analyse de l'accumulation de distance (m)",
  'power': "Analyse des données de puissance (W)",
  'hrv': "Analyse des données R-R (ms)",
}


def _altitude_series(values, cfg):
  """Gap-trim, smooth and compute D+/D- for one altitude value list."""
  trimmed = values[cfg.altitude_gap:]
  smoothed = smooth_altitude(trimmed, cfg)
  return trimmed, smoothed, normalized_alt_gain(smoothed), normalized_alt_loss(smoothed)


def compute_alt_norms(aligned, fitdatas, cfg):
  """Pre-compute per-file normalized D+/D- on the aligned data (used by the report)."""
  norms = {}
  if "altitude" not in cfg.values_to_compare:
    return norms
  for ffile in fitdatas:
    values = field_values(aligned.frames[ffile], 'altitude', cfg, aligned.max_points)
    _, _, gain, loss = _altitude_series(values, cfg)
    norms[ffile] = (gain, loss)
  return norms


def _legend_for(compare_value, ffile, fitdatas, cfg, values, aligned, hr_max_pos, alt_norms):
  """Build the chart legend string for a file/field, mirroring the legacy text."""
  tags = decode_fit_name(ffile)
  s = fitdatas[ffile].summary

  if compare_value == 'heart_rate':
    hr_summary = ''
    if cfg.with_reference_file and len(cfg.fitfiles) >= 2 and ffile != cfg.reference_file:
      ref_hr = field_values(aligned.frames[cfg.reference_file], 'heart_rate', cfg, aligned.max_points)
      score = compute_hr_score(ref_hr, values)
      if score is not None:
        hr_summary = " Ecart moyen: %.2f - Ecart max: %.2f - Score: %.1f%%" % (
          score['average_gap'], score['max_gap'], score['hr_score'])
        print(f"HR Score: {hr_summary}")
        if score['max_gap_position'] is not None:
          hr_max_pos.append(score['max_gap_position'])
    return tags[0] + " (mesure cardio: " + tags[1] + ")" + hr_summary

  if compare_value == 'altitude':
    gain, loss = alt_norms.get(ffile, (0, 0))
    return "%s (D+: %.1f / D-: %.1f / Altitude de départ: %.1f / Altitude d'arrivée: %.1f" % (
      tags[0], gain, loss, s.start_alt, s.end_alt)

  if compare_value == 'distance':
    dist = s.total_distance if s.total_distance is not None else 0
    if tags[3] is not None:
      return "%s (Distance mesurée par: %s): %.2f m" % (tags[0], tags[3], dist)
    return "%s: %.2f m" % (tags[0], dist)

  if compare_value == 'hrv':
    return "%s (%s)" % (tags[0], tags[1])

  return "%s" % tags[0]


def generate_standard_graph(compare_value, aligned, fitdatas, cfg, alt_norms):
  """Build and render one of the configured comparison graphs."""
  print("Generating data for %s" % compare_value)
  chart_title = CHART_TITLES.get(compare_value, 'Analyse du champ de données "%s"' % compare_value)
  chart_data = {}
  hr_max_pos = []
  shortest_hrv = 0

  for ffile in fitdatas:
    if compare_value == 'hrv':
      values = _hrv_values(ffile, cfg)
      if not values:
        print("ERROR: No valid HRV data. Add a CSV or ensure HRV is correctly set in FIT file")
        return
      shortest_hrv = len(values) if shortest_hrv == 0 else min(shortest_hrv, len(values))
    else:
      values = field_values(aligned.frames[ffile], compare_value, cfg, aligned.max_points)
      if compare_value == 'altitude':
        values, smoothed, _, _ = _altitude_series(values, cfg)

    legend = _legend_for(compare_value, ffile, fitdatas, cfg, values, aligned, hr_max_pos, alt_norms)
    chart_data[legend] = values

    if compare_value == 'altitude' and cfg.include_smoothed_alt:
      tags = decode_fit_name(ffile)
      chart_data['%s (smoothed altitude)' % tags[0]] = smoothed[len(smoothed) - len(values):]

  if compare_value == 'hrv':
    chart_data = {legend: vals[:shortest_hrv] for legend, vals in chart_data.items()}

  _render_chart(chart_data, chart_title, _graph_path(compare_value, cfg), cfg,
                aligned.max_points, vlines=hr_max_pos if compare_value == 'heart_rate' and cfg.align else None)


def _hrv_values(ffile, cfg):
  delta = cfg.hrv_delta.get(ffile, 0)
  if ffile in cfg.hrv_csv:
    values = load_csv_hrv(cfg.hrv_csv[ffile], cfg, delta)
  elif ffile in cfg.hrv_suunto:
    values = load_suunto_hrv(cfg.hrv_suunto[ffile], cfg, delta)
  else:
    values = load_fit_hrv(ffile, cfg, delta)
  cfg.dbg("Number of HRV points for %s: %i" % (ffile, len(values)))
  return values


def generate_custom_graphs(aligned, fitdatas, cfg):
  for custom_graph in cfg.custom_graphs:
    graph_name = custom_graph['name']
    print("Generating custom graph: %s" % graph_name)
    chart_data = {}
    for cg_value in custom_graph['values']:
      frame = aligned.frames[cg_value['file']]
      values = field_values(frame, cg_value['field'], cfg, aligned.max_points)
      tags = decode_fit_name(cg_value['file'])
      chart_data["%s - %s" % (tags[0], cg_value['label'])] = values

    if cfg.project_prefix:
      graph_file = APP_PATH + "pnggraphs/" + cfg.project_prefix + "_" + graph_name.lower().replace(" ", "")
    else:
      graph_file = APP_PATH + "pnggraphs/" + graph_name.lower().replace(" ", "")
    pathlib.Path(APP_PATH + "pnggraphs").mkdir(exist_ok=True)
    _render_chart(chart_data, graph_name, graph_file, cfg, aligned.max_points)


# #############################
# MAP

def _route_coordinates(fd, aligned, cfg):
  """Return the list of [long_deg, lat_deg] coordinates for a file's route."""
  frame = aligned.frames[fd.name]
  coords = []
  for ts, row in frame.iterrows():
    lat = row.get('position_lat')
    lon = row.get('position_long')
    # 5hz GPS burst for this timestamp: expand every sub-point
    if ts in fd.gps5hz:
      lats, lons = fd.gps5hz[ts]
      for i, sub_lon in enumerate(lons):
        if sub_lon is not None and lats[i] is not None:
          coords.append(_to_deg(sub_lon, lats[i]))
    elif pd.notna(lon) and pd.notna(lat):
      coords.append(_to_deg(lon, lat))
  return coords


def _to_deg(lon_semi, lat_semi):
  """Semicircle -> degrees, rounded to ~11cm (6 decimals, as the legacy %f output)."""
  return [round(lon_semi * SEMICIRCLE_TO_DEG, 6), round(lat_semi * SEMICIRCLE_TO_DEG, 6)]


# Static <head> + page skeleton up to the (dynamic) legend list
_MAP_HEAD = """<html lang="en">
<head>
<meta charset="utf-8">
<script src="https://unpkg.com/leaflet@1.7.1/dist/leaflet.js" integrity="sha512-XQoYMqMTK8LvdxXYG3nZ448hOEQiglfqkJs1NOQV44cWnUrBc8PkAOcXy20w0vlaXaVUearIOBhiXZ5V3ynxwA==" crossorigin=""></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet-gpx/1.3.1/gpx.min.js"></script>
<link href="https://api.mapbox.com/mapbox-gl-js/v2.2.0/mapbox-gl.css" rel="stylesheet">
<script src="https://api.mapbox.com/mapbox-gl-js/v2.2.0/mapbox-gl.js"></script>
<script src="https://api.mapbox.com/mapbox-gl-js/plugins/mapbox-gl-compare/v0.4.0/mapbox-gl-compare.js"></script>
<link href="https://fonts.googleapis.com/css?family=Montserrat" rel="stylesheet">
<style>
body { font-family: 'Montserrat'; font-size: 16px; }
.mapLegend { font-family: 'Montserrat'; font-size: 16px; }
</style>
<link rel="stylesheet" href="https://api.mapbox.com/mapbox-gl-js/plugins/mapbox-gl-compare/v0.4.0/mapbox-gl-compare.css" type="text/css"></head><body>
<br><br><div align="center" style="width: 1200px; height: 800px; padding-left: 30px;">
<div align="center" id="mapid" style="width: 100%; height: 680px;"></div>
<div align="left">
"""

# Static terrain/DEM setup, injected right after the map is created
_MAP_TERRAIN = """map.on('style.load', () => {
    map.addSource('mapbox-dem', {
        'type': 'raster-dem',
        'url': 'mapbox://mapbox.mapbox-terrain-dem-v1',
        'tileSize': 512,
        'maxzoom': 14
    });
    // add the DEM source as a terrain layer with exaggerated height
    map.setTerrain({ 'source': 'mapbox-dem', 'exaggeration': 1.5 });
});
map.addControl(new mapboxgl.FullscreenControl());
"""


def _map_legend(index, ffile):
  tags = decode_fit_name(ffile)
  return ('<div class="mapLegend" align="left" style="margin-right: 8px; padding-left: 70px;">'
          '<font color="%s">&#9679;</font>%s (Mode GNSS: %s)</div>\n'
          % (GPX_COLORS[index], tags[0], tags[2]))


def _map_route_js(index, coords):
  """GeoJSON source + line layer for one route, serialized with json.dumps."""
  source = {'type': 'geojson', 'data': {'type': 'Feature', 'properties': {},
            'geometry': {'type': 'LineString', 'coordinates': coords}}}
  layer = {'id': 'route%i' % index, 'type': 'line', 'source': 'route%i' % index,
           'layout': {'line-join': 'round', 'line-cap': 'round'},
           'paint': {'line-color': GPX_COLORS[index], 'line-opacity': 0.8, 'line-width': 4}}
  return ("map.addSource('route%i', %s);\n" % (index, json.dumps(source))
          + "map.addLayer(%s);\n" % json.dumps(layer))


def _map_gpx_ref_js(coords):
  """GeoJSON source + dashed black line layer for the GPX reference trace."""
  source = {'type': 'geojson', 'data': {'type': 'Feature', 'properties': {},
            'geometry': {'type': 'LineString', 'coordinates': coords}}}
  layer = {'id': 'gpxRef', 'type': 'line', 'source': 'gpxRef',
           'layout': {'line-join': 'round', 'line-cap': 'round'},
           'paint': {'line-color': '#000000', 'line-opacity': 0.9, 'line-width': 3,
                     'line-dasharray': [2, 2]}}
  return ("map.addSource('gpxRef', %s);\n" % json.dumps(source)
          + "map.addLayer(%s);\n" % json.dumps(layer))


def _map_gpx_ref_legend():
  return ('<div class="mapLegend" align="left" style="margin-right: 8px; padding-left: 70px;">'
          '<font color="#000000">&#9679;</font>GPX (trace de reference)</div>\n')


def generate_map(fitdatas, aligned, cfg, mapbox_key, gpx_coords=None):
  print("Generating map")
  files = list(fitdatas)
  routes = [_route_coordinates(fd, aligned, cfg) for fd in fitdatas.values()]
  # Legacy centered the map on the last file's average position
  last_summary = fitdatas[files[-1]].summary
  map_init = {
    'container': 'mapid',
    'style': 'mapbox://styles/mapbox/%s' % cfg.map_style,
    'center': [last_summary.avg_long, last_summary.avg_lat],
    'zoom': 13,
  }

  html = [_MAP_HEAD]
  html += [_map_legend(i, f) for i, f in enumerate(files)]
  if cfg.draw_gpx_ref and gpx_coords:
    html.append(_map_gpx_ref_legend())
  html.append('</div>\n<script>\n')
  html.append("mapboxgl.accessToken = '%s';\n" % mapbox_key)
  html.append("var map = new mapboxgl.Map(%s);\n" % json.dumps(map_init))
  html.append(_MAP_TERRAIN)
  html.append("map.on('load', function () {\n")
  html += [_map_route_js(i, coords) for i, coords in enumerate(routes)]
  if cfg.draw_gpx_ref and gpx_coords:
    # GPX points are stored as (lat, lon); the map expects [lon, lat]
    ref_coords = [[round(lon, 6), round(lat, 6)] for lat, lon in gpx_coords]
    html.append(_map_gpx_ref_js(ref_coords))
  html.append('});\n</script>\n</body>\n</html>\n')

  pathlib.Path(APP_PATH + "map").mkdir(exist_ok=True)
  map_file = APP_PATH + "map/" + ((cfg.project_prefix + "_map.html") if cfg.project_prefix else "map.html")
  with open(map_file, "w") as fmap:
    fmap.write("".join(html))


# #############################
# EXAMPLE CONFIG

def generate_example_config(fitdatas, cfg, common_count):
  """Write project.yaml.example, pre-filled with sensible/zoom values."""
  first = list(fitdatas.values())[0]
  all_sessions = first.sessions
  multisession_starts = []
  if len(all_sessions) > 1:
    for i in range(len(all_sessions)):
      multisession_starts.append([fd.sessions[i][3] for fd in fitdatas.values()])

  with open(APP_PATH + "project.yaml.example", "w") as out:
    out.write('project:\n')
    out.write('  align: False\n')
    if len(all_sessions) > 1:
      old_start = None
      for i, session in enumerate(multisession_starts):
        this_start = np.amax(session)
        if i >= 1:
          out.write('  zoom: [%i, %i] # %s\n' % (old_start, this_start, all_sessions[i - 1][0]))
        old_start = this_start
      out.write('  zoom: [%i, %i] # %s\n' % (old_start, common_count, all_sessions[len(multisession_starts) - 1][0]))
    else:
      out.write('  zoom: [90, 120]\n')
    out.write('  altitudeGap: 8  # Seconds\n')
    out.write('  map: false\n')
    out.write('  mapStyle: outdoors-v12\n')
    out.write('  drawGpxRef: false # Draw the --gpx-ref trace on the map (dashed black)\n')
    out.write('  graphs: [\'heart_rate\', \'altitude\', \'distance\']\n')
    out.write('  includeSmoothedAlt: false\n')
    out.write('  removeAbnormalHrv: false\n')
    out.write('  removeAbnormalHrvThreshold: 20 # percentage of the previous value\n')
    out.write('customGraphs:\n')
    out.write('  - name: Altitude baro vs GPS\n')
    out.write('    values:\n')
    out.write('      - file: %s\n' % first.name)
    out.write('        field: enhanced_altitude\n')
    out.write('        label: Altitude baro\n')
    out.write('      - file: %s\n' % first.name)
    out.write('        field: GPS altitude\n')
    out.write('        label: Altitude GPS\n')
    for ffile in fitdatas:
      out.write('%s:\n' % ffile)
      out.write('  delta: 0\n')
      out.write('#  charge: [99, 87]\n')
      out.write('#  hrvCsv: polar_hrv.csv\n')


# #############################
# MAIN

def main():
  print("Running fitcompare v%s" % SCRIPT_VER)

  global_conf = configparser.ConfigParser()
  global_conf.read('config.ini')
  mapbox_key = global_conf['map']['mapbox_api_key']

  cfg = load_config(parse_args())

  # Parse every FIT file exactly once
  fitdatas = {}
  for ffile in cfg.fitfiles:
    cfg.dbg("Processing file %s" % ffile)
    fitdatas[ffile] = load_fit(ffile, cfg)
    cfg.dbg("File %s has %i points" % (ffile, fitdatas[ffile].summary.n_records))

  # Align (or pad) the working data
  aligned = align_files(fitdatas, cfg)

  # Generate an example project configuration, then stop (--gen-config)
  if cfg.gen_config:
    generate_example_config(fitdatas, cfg, aligned.max_points)
    print(" Example configuration written to %sproject.yaml.example" % APP_PATH)
    return

  # D+/D- computed once, on the aligned data, and reused by report + legends
  alt_norms = compute_alt_norms(aligned, fitdatas, cfg)

  # GPS reference (GPX): score every FIT file against the absolute trace
  gps_scores = None
  gpx_coords = None
  if cfg.gpx_ref:
    gpx_coords = load_gpx(cfg.gpx_ref)
    cfg.dbg("GPX reference trace loaded with %i points" % len(gpx_coords))
    if gpx_coords:
      gps_scores = compute_gps_scores(fitdatas, aligned, cfg, gpx_coords)
    else:
      print("WARNING: GPX reference file %s contains no track points" % cfg.gpx_ref)

  # Text report
  report = build_report(fitdatas, cfg, alt_norms, gps_scores=gps_scores)
  write_logfile(report, cfg)
  if cfg.align:
    print(" Common timestamps:                  %i" % aligned.max_points)
  else:
    print(" Longest timestamps:                  %i" % aligned.max_points)
  print("=========================================================================")

  # Graphs
  for compare_value in cfg.values_to_compare:
    cfg.dbg("Configuring output for field %s" % compare_value)
    generate_standard_graph(compare_value, aligned, fitdatas, cfg, alt_norms)
  generate_custom_graphs(aligned, fitdatas, cfg)
  if gps_scores:
    generate_gps_graph(gps_scores, aligned, fitdatas, cfg)

  # Map
  if cfg.map:
    generate_map(fitdatas, aligned, cfg, mapbox_key, gpx_coords=gpx_coords)


if __name__ == "__main__":
  main()
