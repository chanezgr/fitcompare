import numpy as np

# This function find the closest value to "value" in "array"
def find_nearest_value(array, value):
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return array[idx]

# This function uses the find_nearest_value in order to compensate latency of HR measurement or slight misalignments
def reduce_latency(ff_data, reference_file, a_position, value):
  """
  This function search in the last 5 seconds the closest value to "value"
  """
  slice_hr = []
  
  if (a_position < 5):
    slice = ff_data[reference_file][:a_position]
    for record in slice:
      slice_hr.append(record['heart_rate'])
  else:
    slice = ff_data[reference_file][a_position-5:a_position]
    for record in slice:
      slice_hr.append(record['heart_rate'])
  
  nv = find_nearest_value(slice_hr, value)
  return nv
  
# This function get the HR value at a timestamp, and for the next 4 seconds
# Input:
# - file_data: array of value for a fit file
# - the specific timestamp to compare
# Output:
# - The HR for the specific timestamp
def get_bpm_ts(file_data, timestamp):
  # Loop over all data
  for point in file_data:
    cur_ts = point['timestamp']
    if (cur_ts == timestamp):
      if (point['heart_rate'] != None):
        return point['heart_rate']
      elif (prev_bpm != None):
        return prev_bpm
      else:
        return 0
    prev_bpm = point['heart_rate']

# This function handles a new HR comparison point
def bpm_new_point(cur_bpm, ts, average_hr_gap, ff_data, reference_file, a_position):
  ref_bpm = get_bpm_ts(ff_data[reference_file], ts)
  if ((ref_bpm != 0) and (ref_bpm != None)):
    closest_value = reduce_latency(ff_data, reference_file, a_position, cur_bpm)
    # Compute this point gap
    bpm_gap = abs(cur_bpm - closest_value)
    # Add this point gap to array
    average_hr_gap['average'].append(bpm_gap)
    # If the gap is greater than previous gap, it's a new bigger one
    if (bpm_gap > average_hr_gap['max']):
      average_hr_gap['max'] = bpm_gap
      average_hr_gap['max_position'] = a_position
      #print("%i - %i" % (a_position, bpm_gap))
  return average_hr_gap

def adv_hr_sum(average_hr_gap):
  # Average BPM diff:
  avg_bpm_gap_final = sum(average_hr_gap['average']) / len(average_hr_gap['average']) 
  # Score of average bpm
  if (avg_bpm_gap_final <= 0.5):
    avg_bpm_coef = 0
  else:
    avg_bpm_coef = abs(avg_bpm_gap_final) - 0.5
  hr_gap_score = avg_bpm_coef*10
  if (hr_gap_score >= 60):
    hr_gap_score = 60
  if (average_hr_gap['max'] <= 80):
    max_bpm_gap_score = abs(average_hr_gap['max']/1.7)
  else:
    max_bpm_gap_score = 50

  hr_score = 100 - (hr_gap_score + max_bpm_gap_score)
  
  returnValues = {}
  returnValues['average_gap'] = avg_bpm_gap_final
  returnValues['max_gap'] = average_hr_gap['max']
  returnValues['max_gap_position'] = average_hr_gap['max_position']
  returnValues['hr_score'] = hr_score
  return returnValues