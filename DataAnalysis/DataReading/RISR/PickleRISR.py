import calendar
import os
import time
import glob

import pandas as pd
import numpy as np

from DataAnalysis.DataReading.RISR_HDF5.wd_beam_num import wd_beam_num


def PickleRISR(station, date):
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
            Aspect angle 'aspect'  The angle between the magnetic field line and the radar beam.
            Transmission Frequency 'tfreq'
            Geodetic Latitude and Longitude 'gdlat', 'glon'
            Code baud count per pulse 'cbadl'
            Pulse length in seconds 'pl'

        2D Parameters (parameters that vary with range and time)
           Density 'Ne' with error 'Ne_err'
            Ion Temperature 'ion_temp' with error 'ion_temp_err'
            Electron Temperature 'e_temp' with error 'e_temp_err'
            Line of Sight Ion Velocity 'los_ion_vel' with error 'los_ion_vel_err'
            Corrected Geomagnetic Latitude 'cgm_lat'
            Corrected Geomagnetic Longitude 'cgm_lon'
            Composition [O+]/Ne 'comp' with error 'comp_err'
    """

    # Loop through all the files for this station/date and pickle all txt files
    in_dir = "data/" + station + "/" + station + date

    # Build a temp file all the unnecessary data stripped out
    for in_file in glob.iglob(in_dir + "/*_raw.txt"):

        # First just obtain the header
        with open(in_file, 'r') as infile:
            for line in infile:
                if line.strip():
                    if line.strip()[0] == "Y":
                        header = line

        # Strip everything that is not data
        temp_file = in_dir + "/temp_" + station + date + ".txt"
        with open(in_file, 'r') as infile, open(temp_file, 'w+') as outfile:
            outfile.write(header.lstrip())  # Add in the file header
            for line in infile:
                stripped_line = line.strip()
                if not stripped_line or stripped_line[0] != "2":
                    continue

                outfile.write(line.lstrip())  # non-empty line. Write it to output

        # Read in the data from the adjusted temp file
        df = pd.read_csv(temp_file, delim_whitespace=True, lineterminator='\r')
        df.drop(df.tail(1).index, inplace=True)  # drop last row, it is just white space

        df['MONTH'] = df['MONTH'].astype(int)
        df['DAY'] = df['DAY'].astype(int)
        df['HOUR'] = df['HOUR'].astype(int)
        df['MIN'] = df['MIN'].astype(int)
        df['SEC'] = df['SEC'].astype(int)

        pattern = '%Y-%m-%d %H:%M:%S'
        year = [int(df['YEAR'][0])] * len(df['YEAR'])  # The years near the bottom of the file can get messed up
        epoch = []
        wdBmnum = []
        date_time = []
        # Loop through all the rows are compute things that need to be computed
        for row in range(len(year)):
            wdBmnum.append(wd_beam_num(df['BEAMID'][row]))

            date_time_here = str(year[row]) + "-" + str(df['MONTH'][row]) + "-" + str(df['DAY'][row]) + " " + \
                             str(df['HOUR'][row]) + ":" + str(df['MIN'][row]) + ":" + str(df['SEC'][row])
            date_time.append(date_time_here)
            epoch.append(calendar.timegm(time.strptime(date_time_here, pattern)))

        adjusted_df = pd.DataFrame(
            {'stationId': [station] * len(epoch),
             'dateTime': date_time,
             'epoch': epoch,
             'decimalTime': np.asarray(df['HOUR']) + np.asarray(df['MIN']) / 60.0 + np.asarray(df['SEC']) / 3600.0,
             'year': year, 'month': df['MONTH'], 'day': df['DAY'],
             'hour': df['HOUR'], 'minute': df['MIN'], 'second': df['SEC'],
             'range': df['RANGE'],
             'bmId': df['BEAMID'], 'wdBmnum': wdBmnum, 'bmazm': df['AZM'],
             'cbadl': df['CBADL'],
             'elv': df['ELM'],
             'transFreq': df['TFREQ'], 'reciveFreq': df['RFREQ'],
             'gdlat': df['GDLAT'], 'gdlon': df['GLON'],
             'Ne': df['NE'], 'NeErr': df['DNE'],
             'ionTemp': df['TI'], 'ionTempErr': df['DTI'],
             'eTemp': df['TE'], 'eTempErr': df['DTE'],
             'losIonVel': df['VO'], 'losIonVelErr': df['DVO'],
             'cgmLat': df['CGM_LAT'], 'cgmLon': df['CGM_LONG'],
             'comp': df['PO+'], 'compErr': df['DPO+'],
             'aspect': df['ASPECT'],
             'pulseLength': df['PL']
             })

        # Compute Data Resolution in minutes
        temp_df = adjusted_df.drop_duplicates(subset=['dateTime'])
        resolution = int((temp_df['minute'].iloc[1] + temp_df['second'].iloc[1] / 60.0)
                         - (temp_df['minute'].iloc[0] + temp_df['second'].iloc[0] / 60.0))

        # Save to file
        out_file = in_file[:len(in_file) - 8] + "." + str(resolution) + "min.pkl"
        print("     Pickling as " + out_file + "...")
        adjusted_df.to_pickle(out_file)

        # Remove the temporary temp file
        for file in glob.iglob(in_dir + "/temp*.txt"):
            os.remove(file)


if __name__ == '__main__':
    """
    Handler to call PickleRISR on RISR txt data files
    """
    PICKLE_ALL = True  # To prevent accidentally pickling all data

    if PICKLE_ALL:
        print("Pickling all downloaded RISR data...")
        for station in os.listdir("data/"):
            for in_dir in os.listdir("data/" + station):
                print("\nStarting " + in_dir)
                PickleRISR(station, in_dir[3:])
    else:
        station = "ran"
        date = "20091107"
        print("Pickling " + station + date + "...")
        PickleRISR(station, date)