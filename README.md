# fitcompare
Python tool to compare two FIT files data

This script is only intended to run in a docker image. Any other use is not recommended.

Note that this is the fitcompare tool used for nakan.ch compare graphs and data.

## Build the Docker image
1. Get the file `fitcompare.py`, `fitcompare_advanced.py` and `Dockerfile`, place all the files in a directory
2. Create a Mapbox API key and put it in a file named `config.ini` with the following format: 

```
[map]
mapbox_api_key=YOUR_PERSONAL_MAPBOX_KEY_HERE
```

3. Build the Docker image using: `docker build . -t fitcompare`

## Run the container=
The long command to run the container is:
`docker run -v .:/project fitcompare OPTIONS`

I use to place someting like this in my .bashrc/.zshrc:

`alias fitcompare="docker run -v .:/project fitcompare"`

This way, I can run only "fitcompare OPTIONS" in the directory where FIT files are 

## Command line usage

```
usage: fitcompare.py [-h] [--reference-file REFERENCE_FILE]
                     [--prefix PROJECT_PREFIX] [--debug] [--export]
                     [--config PROJECT_CONFIG]
                     FITFILE [FITFILE ...]

Compare two or more FIT files

positional arguments:
  FITFILE               Fit Files to compare

options:
  -h, --help            show this help message and exit
  --reference-file, -r REFERENCE_FILE
                        Set the reference FIT File
  --prefix, -p PROJECT_PREFIX
                        Set the project prefix for output files
  --debug, -d           Enable debug
  --export, -e          Export graphs values also as CSV
  --config, -c PROJECT_CONFIG
                        Use an alternative configuration YAML file
```

## Name of the FIT files
  
FIT File name have to be correctly formatted:

`<DEVICE_NAME>_<HR_DEVICE>_<GNSS_CONFIG>[_<OPTIONAL_DISTANCE_MEASUREMENT_DEVICE>].fit`

Where:
`<DEVICE_NAME>` = Name without spaces
`<HR_DEVICE>` = Name without spaces
`<GNSS_CONFIG>` =  - `SatIQ`:    for Garmin SatIQ automatic
                 - `Ultra`:    for UltraTrack or equivalent
                 - `GPS`:      for GPS only
                 - `GNSS`:     for all GNSS systems
                 - `GNSSDual`: for all GNSS systems + GPS Multiband
                 - `SatIQ`:    for automatic mode on Garmin watches
                 - `NONE`:     for non GPS records
Example:

`fitcompare -r SuuntoRace_PolarH10_GNSSDual_Stryd.fit GarminFenix9_OHR_GNSSDual.fit GarmninInstinct4_OHR_GNSSDual.fit -p compare_prototypes`

## Configuration YAML file

In addition to command line, options can be configured in the project.yaml in the same directory.
After each successful run, a new project.yaml.example is generated with example values

```
project: # This is a section that configure global aspects of the project
  align: True/False # True is the default value, the one to compare two files generated at the same time. False is to compare differetn files from different time. For information only, cannot be a valid data processing, because fake data are added to the shorter file.
  zoom: [90, 120] # Zoom between two timestamps (relative seconds of activity)
  altitudeGap: 8  # Number of seconds ignored at the beginning of activity, if some files start at 0 altitude and then put the correct one.
  map: True/False # If a map should be generated. It is mandatory to put a map: False for activities without GPS data
  mapStyle: outdoors-v12 # Should be a valid Mapbox style https://docs.mapbox.com/api/maps/styles/
  graphs: ['heart_rate', 'altitude', 'distance'] # Fields for which a graph shoud be generated. Usual values are: heart_rate, distance, speed, altitude, cadence, power, hrv
  includeSmoothedAlt: False # Should the data of smoothed altitude be included into the elevation graph
  removeAbnormalHrv: false # If HRV values are plotted, this will remove abnormal spikes in HRV values
  removeAbnormalHrvThreshold: 20 # percentage of the previous value a HRV point will be considered abnormal
  
customGraphs: # In this section, we can configure custom graphs
  - name: Altitude baro vs GPS # Name of the custom graph
    values: # Values
      - file: GarminFenix8_WahooTRACKR_GNSSDual.fit # Fit file with the data field
        field: enhanced_altitude # data field to include in the customer graph
        label: Altitude baro # label of the field fir this file
      - file: GarminFenix8_WahooTRACKR_GNSSDual.fit 
        field: GPS altitude
        label: Altitude GPS
GarminFenix8_WahooTRACKR_GNSSDual.fit: # options for each fit files
  delta: 0 # delta in second ti apply to the timestamps of this file (can be positive or negative int value)
  charge: [99, 87] # Battery level at beginning / end of activity to generate battery life estimation, if not integrated in FIT
  hrvCsv: polar_hrv.csv # Source for HRV data in CSV file (for Polar R-R records in CSV format)
AppleWatchSeries10_OHR_GNSS.fit:
  delta: 0
```
