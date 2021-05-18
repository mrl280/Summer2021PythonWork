import math
import os
import time
import glob
import h5py

import pandas as pd
import numpy as np

from wd_beam_num import wd_beam_num


def PickleLongPulse(station, date):
    """
    Take a RISR_HDF5 long pulse file and pickle it
    Works for both RISR_HDF5-N and RISR_HDF5-C, although some RISR_HDF5-C files seem not to have altitude data
    
    Parameters of Interest:
        - Epoch time (seconds since 1970.01.01)
        - Date time (human readable date time)
        - Year
        - Month
        - Day
        - Hour
        - Minute
        - Second (included microseconds)
        - Range
    
        1D Parameters (parameters that vary only with time)
            Beam id 'bmid' (this is usually a 5 digit number that uniquely identifies a RISR_HDF5 beam)
            World-Day Beam Number 'bmnum' (1 through 11)
            Beam azimuth 'bmazm'
            Elevation Angle 'elv' (0=horizontal,90=vert)
            Transmission Frequency 'tfreq'
            Geodetic Latitude and Longitude 'gdlat', 'glon'
            Code baud count per pulse 'cbadl'

        2D Parameters (parameters that vary with range and time)
            Log10 of Density 'log10_Ne' with error 'Errolog10_Ne_err'
            Ion Temperature 'ion_temp' with error 'ion_temp_err'
            Electron Temperature 'e_temp' with error 'e_temp_err'
            Line of Sight Ion Velocity 'los_ion_vel' with error 'los_ion_vel_err'
            Corrected Geomagnetic Latitude 'cgm_lat'
            Corrected Geomagnetic Longitude 'cgm_lon'
            Composition [O+]/Ne 'comp' with error 'comp_err' 
            Altitude 'alt'
    """

    pattern = '%Y-%m-%d %H:%M:%S'

    # Loop through all the files for this station/date and pickle all Long Pulse files
    in_dir = "data/" + station + "/" + station + date
    for in_file in glob.iglob(in_dir + "/*.h5"):

        try:
            file = h5py.File(in_file, "r")
        except:
            continue

        # Look at the file type, we are only interested in Long Pulse Files here
        exp_notes = file['/Metadata/Experiment Notes']
        if station == "ran":
            file_type = str(exp_notes[25][0])
        elif station == "ras":
            file_type = str(exp_notes[39][0])
        else:
            raise Exception("Station " + station + " not recognized")
        file_type = file_type[10:len(file_type) - 1].rstrip()  # Strip out everything but the file type name
        if file_type != "Long Pulse":
            print("     " + in_file + " is a " + file_type + " file, skipping it...")
        else:
            # The file is a long pulse file, we are good to go ahead
            print("     " + in_file + " is a " + file_type + " file, reading in data...")

            keys = [key for key in file['/Data/Array Layout/'].keys()]  # These will be the list of beams
            try:
                data_time = file['/Data/Array Layout/' + keys[0] + '/timestamps']
            except:
                raise Exception("Error: no keys found, unable to obtain data resolution")
            resolution = int((data_time[1] - data_time[0]) / 60)
            print("         Data time resolution: " + str(resolution) + " minute data")
            # Strip the "x.h5" and replace with data resolution and "LongPulse.pkl"
            out_file = in_file[:len(in_file) - 4] + str(resolution) + "min.LongPulse.pkl"

            # Pre-allocate list for 1D parameters (only vary with time)
            epoch = []
            year, month, day = [], [], []
            hour, minute, second = [], [], []
            date_time = []
            range = []
            bmid, bmnum, bmazm = [], [], []
            cbadl = []  # Code baud count per pulse  # Data Quality Parameter
            elv = []
            tfreq = []
            gdlat, gdlon = [], []  # Geodetic

            # Pre-allocate lists for 2D parameters (vary with range and time)
            log10_Ne, log10_Ne_err = [], []
            ion_temp, ion_temp_err = [], []
            e_temp, e_temp_err = [], []
            los_ion_vel, los_ion_vel_err = [], []
            cgm_lat, cgm_lon = [], []
            comp, comp_err = [], []
            alt = []  # Geographic coordinates

            # Loop through all the beams
            for beam in file['/Data/Array Layout/'].keys():
                print("         Looping through beam "
                      + str(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][0]))

                # Loop through all the times and ranges here
                for epoch_here_idx, epoch_here in enumerate(file['/Data/Array Layout/' + beam + '/timestamps']):

                    # Compute required parameters so they don't need to be re-computed at every range
                    date_time_here = time.strftime(pattern, time.gmtime(epoch_here))
                    year_here = int(date_time_here[:4])
                    month_here = int(date_time_here[5:7])
                    day_here = int(date_time_here[8:10])
                    hour_here = int(date_time_here[11:13])
                    minute_here = int(date_time_here[14:16])
                    second_here = int(date_time_here[17:19])

                    bmid_here = file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][epoch_here_idx]
                    bmnum_here = wd_beam_num(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][epoch_here_idx])
                    bmazm_here = file['/Data/Array Layout/' + beam + '/1D Parameters/azm'][epoch_here_idx]
                    cbadl_here = file['/Data/Array Layout/' + beam + '/1D Parameters/cbadl'][epoch_here_idx]
                    elv_here = file['/Data/Array Layout/' + beam + '/1D Parameters/elm'][epoch_here_idx]
                    tfreq_here = file['/Data/Array Layout/' + beam + '/1D Parameters/tfreq'][epoch_here_idx]
                    gdlat_here = file['/Data/Array Layout/' + beam + '/1D Parameters/gdlat'][epoch_here_idx]
                    gdlon_here = file['/Data/Array Layout/' + beam + '/1D Parameters/glon'][epoch_here_idx]

                    for range_here_idx, range_here in enumerate(file['/Data/Array Layout/' + beam + '/range']):
                        epoch.append(epoch_here)
                        range.append(range_here)

                        date_time.append(date_time_here)
                        year.append(year_here)
                        month.append(month_here)
                        day.append(day_here)
                        hour.append(hour_here)
                        minute.append(minute_here)
                        second.append(second_here)

                        # The 1D parameters are just accessed in time
                        # These Could be moved out of this range loop but it is simpler to leave them here)
                        bmid.append(bmid_here)
                        bmnum.append(bmnum_here)
                        bmazm.append(bmazm_here)
                        cbadl.append(cbadl_here)
                        elv.append(elv_here)
                        tfreq.append(tfreq_here)
                        gdlat.append(gdlat_here)
                        gdlon.append(gdlon_here)

                        # The 2D parameters are accessed in time and range
                        log10_Ne.append(file['/Data/Array Layout/' + beam + '/2D Parameters/nel'][range_here_idx][epoch_here_idx])
                        log10_Ne_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dnel'][range_here_idx][epoch_here_idx])
                        ion_temp.append(file['/Data/Array Layout/' + beam + '/2D Parameters/ti'][range_here_idx][epoch_here_idx])
                        ion_temp_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dti'][range_here_idx][epoch_here_idx])
                        e_temp.append(file['/Data/Array Layout/' + beam + '/2D Parameters/te'][range_here_idx][epoch_here_idx])
                        e_temp_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dte'][range_here_idx][epoch_here_idx])
                        los_ion_vel.append(file['/Data/Array Layout/' + beam + '/2D Parameters/vo'][range_here_idx][epoch_here_idx])
                        los_ion_vel_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dvo'][range_here_idx][epoch_here_idx])
                        cgm_lat.append(file['/Data/Array Layout/' + beam + '/2D Parameters/cgm_lat'][range_here_idx][epoch_here_idx])
                        cgm_lon.append(file['/Data/Array Layout/' + beam + '/2D Parameters/cgm_long'][range_here_idx][epoch_here_idx])
                        comp.append(file['/Data/Array Layout/' + beam + '/2D Parameters/po+'][range_here_idx][epoch_here_idx])
                        comp_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dpo+'][range_here_idx][epoch_here_idx])
                        try:
                            alt.append(file['/Data/Array Layout/' + beam + '/2D Parameters/gdalt'][range_here_idx][epoch_here_idx])
                        except:
                            # This parameter might not exist
                            alt.append(math.nan)

            # Put the data into a dataframe
            print("     Building data frame...")
            df = pd.DataFrame(
                {'stationId': [station] * len(epoch),
                 'dateTime': date_time,
                 'epoch': epoch,
                 'decimalTime': np.asarray(hour) + np.asarray(minute) / 60.0 + np.asarray(second) / 3600.0,
                 'year': year, 'month': month, 'day': day,
                 'hour': hour, 'minute': minute, 'second': second,
                 'range': range,
                 'bmId': bmid, 'wdBmnum': bmnum, 'bmazm': bmazm,
                 'cbadl': cbadl,
                 'elv': elv,
                 'transFreq': tfreq,
                 'gdlat': gdlat, 'gdlon': gdlon,
                 'log10Ne': log10_Ne,       'log10NeErr': log10_Ne_err,
                 'ionTemp': ion_temp,       'ionTempErr': ion_temp_err,
                 'eTemp': e_temp,           'eTempErr': e_temp_err,
                 'losIonVel': los_ion_vel,  'losIonVelErr': los_ion_vel_err,
                 'cgmLat': cgm_lat,         'cgmLon': cgm_lon,
                 'comp': comp,              'compErr': comp_err,
                 'alt': alt,
                 })

            # Save to file
            print("     Pickling as " + out_file + "...")
            df.to_pickle(out_file)

        file.close()


if __name__ == '__main__':
    """
    Handler to call PickleFITACF on RISR_HDF5 Long Pulse data files
    """
    PICKLE_ALL = False  # To prevent accidentally pickling all data

    if PICKLE_ALL:
        print("Pickling all downloaded RISR_HDF5 data...")
        for station in os.listdir("data/"):
            for in_dir in os.listdir("data/" + station):
                print("\nStarting " + in_dir)
                PickleLongPulse(station, in_dir[3:])
    else:
        station = "ran"
        date = "20161012"
        print("Pickling " + station + date + "...")
        PickleLongPulse(station, date)
