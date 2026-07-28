"""
Advanced heart-rate analysis for fitcompare.

The public entry point is ``compute_hr_score``: given the reference and
candidate heart-rate series *already aligned on the common timestamps*, it
reproduces the legacy latency-compensated gap scoring in a single O(n) pass
(the previous implementation was O(n^2) because it re-scanned the whole
reference list for every point).
"""

import numpy as np


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
