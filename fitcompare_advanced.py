"""
Advanced heart-rate and GPS analysis for fitcompare.

The public entry points are ``compute_hr_score`` and ``compute_gps_score``:

* ``compute_hr_score``: given the reference and candidate heart-rate series
  *already aligned on the common timestamps*, it reproduces the legacy
  latency-compensated gap scoring in a single O(n) pass (the previous
  implementation was O(n^2) because it re-scanned the whole reference list for
  every point).

* ``compute_gps_score``: given the candidate FIT GPS points and an absolute
  GPX reference trace, it finds the nearest GPX point for every FIT point,
  applies a 1 m tolerance (retained gap = max(0, dist_cm - 100)) and derives a
  score mirroring the HR one. The nearest-point search uses a KD-tree built on
  an equirectangular projection (meters), which is accurate enough for the
  cm-level, local-area comparisons this tool is designed for.
"""

import math

import numpy as np
from scipy.spatial import cKDTree

# Mean Earth radius used for the equirectangular projection (meters)
EARTH_RADIUS_M = 6371000.0
# Free GPS tolerance: the first 500cm of every point-to-point gap is ignored
GPS_MARGIN_CM = 500.0


def _nearest_value(window, value):
  """Return the element of ``window`` closest to ``value`` (legacy find_nearest_value)."""
  window = np.asarray(window, dtype=float)
  return float(window[np.abs(window - value).argmin()])


def compute_hr_score(ref_hr, cand_hr, start_index=59, latency=5):
  """Compute the heart-rate comparison score between two aligned series.

  Parameters
  ----------
  ref_hr, cand_hr : sequence of float
      Heart-rate values of the reference and candidate files, aligned point by
      point on the common timestamps (same length).
  start_index : int
      0-based index of the first point taken into account (legacy started the
      analysis at the 60th point, i.e. 0-based index 59).
  latency : int
      Size of the backward window used to compensate HR measurement latency
      (legacy looked at the last 5 seconds).

  Returns
  -------
  dict or None
      ``None`` if no point could be scored, otherwise a dict with keys
      ``average_gap``, ``max_gap``, ``max_gap_position`` and ``hr_score``.
  """
  gaps = []
  max_gap = 0
  max_gap_position = None

  for pos in range(len(cand_hr)):
    if pos < start_index:
      continue
    ref_bpm = ref_hr[pos]
    # Skip points where the reference has no usable value (legacy get_bpm_ts)
    if ref_bpm is None or ref_bpm == 0:
      continue
    # Closest reference value within the latency window (legacy reduce_latency:
    # ff_data[ref][pos-latency+1 : pos+1], i.e. the current point and the
    # previous `latency-1` ones).
    window = ref_hr[max(0, pos - latency + 1):pos + 1]
    gap = abs(cand_hr[pos] - _nearest_value(window, cand_hr[pos]))
    gaps.append(gap)
    if gap > max_gap:
      max_gap = gap
      # Legacy stored a_position == i == pos + 1 (1-based)
      max_gap_position = pos + 1

  if not gaps:
    return None

  return _summarize(gaps, max_gap, max_gap_position)


def _summarize(gaps, max_gap, max_gap_position):
  """Turn the collected gaps into the final score (legacy adv_hr_sum)."""
  average_gap = sum(gaps) / len(gaps)

  # Average-gap penalty: nothing below 0.5 bpm, then 10 points per bpm, capped at 60
  avg_coef = 0 if average_gap <= 0.5 else abs(average_gap) - 0.5
  avg_gap_score = min(avg_coef * 10, 60)

  # Max-gap penalty: proportional up to 80 bpm, then capped at 50
  max_gap_score = abs(max_gap / 1.7) if max_gap <= 80 else 50

  return {
    'average_gap': average_gap,
    'max_gap': max_gap,
    'max_gap_position': max_gap_position,
    'hr_score': 100 - (avg_gap_score + max_gap_score),
  }


def _to_xy(coords, lat0_rad):
  """Equirectangular projection of (lat_deg, lon_deg) pairs to planar meters."""
  lat = np.radians(coords[:, 0])
  lon = np.radians(coords[:, 1])
  x = lon * EARTH_RADIUS_M * math.cos(lat0_rad)
  y = lat * EARTH_RADIUS_M
  return np.column_stack([x, y])


def compute_gps_score(fit_coords, gpx_coords, margin_cm=GPS_MARGIN_CM, measured_dist=None, real_dist=None, laps=1):
  """Compute the GPS comparison score of a FIT trace against a GPX reference.
  
  Parameters
  ----------
  fit_coords : sequence of (lat_deg, lon_deg)
      GPS points of one FIT file, in the same order as the records (so the
      per-point gap array stays aligned with the other graphs' x-axis). Points
      with no fix are represented by (nan, nan) and are skipped.
  gpx_coords : sequence of (lat_deg, lon_deg)
      Every track point of the GPX reference trace.
    margin_cm : float
       Free tolerance applied per point: retained_gap = max(0, dist_cm - margin).
       Defaults to 500cm.
    measured_dist : float
       Total distance measured by the device (m).
    real_dist : float
       Total distance of one loop of the reference GPX (m).
    laps : int
       Number of loops performed.
  
  Returns
  -------
  dict or None
      ``None`` if the FIT file has no usable fix or the GPX trace is empty,
      otherwise a dict with keys:
  
      * ``average_gap``  - mean of the retained gaps (cm)
      * ``max_gap``      - largest retained gap (cm)
      * ``max_gap_position`` - 1-based index of the largest gap in ``gaps``
      * ``gps_score``    - score on 100 (100 = perfect)
      * ``gaps``         - per-point retained gaps (cm), NaN where no fix
      * ``raw_dists``    - per-point raw distances (cm) to the nearest GPX pt
      * ``measured_dist`` - measured distance (m)
      * ``real_dist``     - actual total reference distance (m)
  """
  n = len(fit_coords)
  gaps = [float('nan')] * n
  raw_dists = [float('nan')] * n
  
  valid_idx = [i for i, c in enumerate(fit_coords)
               if c is not None and len(c) == 2
               and c[0] == c[0] and c[1] == c[1]   # not NaN
               and c[0] != 0.0 and c[1] != 0.0]
  if not valid_idx or len(gpx_coords) == 0:
    return None
  
  gpx = np.asarray(gpx_coords, dtype=float)
  fit_valid = np.asarray([fit_coords[i] for i in valid_idx], dtype=float)
  lat0 = math.radians(float(np.mean(gpx[:, 0])))
  tree = cKDTree(_to_xy(gpx, lat0))
  dists_m, _ = tree.query(_to_xy(fit_valid, lat0))
  dists_cm = dists_m * 100.0
  
  retained = []
  for k, i in enumerate(valid_idx):
    cm = float(dists_cm[k])
    raw_dists[i] = cm
    gap = cm - margin_cm
    if gap < 0:
      gap = 0.0
    gaps[i] = gap
    retained.append(gap)
  
  average_gap = float(sum(retained) / len(retained))
  max_gap = float(max(retained))
  max_gap_position = valid_idx[int(np.argmax(retained))] + 1  # 1-based, like HR
  
  # GPS score (cm units after the margin). The average gap drives most of the
  # penalty; the max gap is penalized gently so a single spike does not tank the
  # whole score.
  avg_penalty = min(average_gap * 0.05, 45)
  max_penalty = min(max_gap * 0.008, 45)
  
  # Distance penalty
  dist_penalty = 0
  actual_total_dist = None
  if real_dist is not None and measured_dist is not None:
    actual_total_dist = real_dist * laps
    diff_pct = abs(measured_dist - actual_total_dist) / actual_total_dist
    if diff_pct > 0.01:
      dist_penalty = min((diff_pct - 0.01) * 100 * 2, 20)
  
  return {
    'average_gap': average_gap,
    'max_gap': max_gap,
    'max_gap_position': max_gap_position,
    'gps_score': 100 - (avg_penalty + max_penalty + dist_penalty),
    'gaps': gaps,
    'raw_dists': raw_dists,
    'measured_dist': measured_dist,
    'real_dist': actual_total_dist,
  }
