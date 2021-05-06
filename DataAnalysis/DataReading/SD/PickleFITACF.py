import calendar
import pydarn
import bz2
import pandas as pd
import time
import glob

from DataAnalysis.DataReading.SD.radar_code_from_station_id import radar_code_from_station_id

if __name__ == '__main__':
    """
    Take a fitACF file, turn it into a reduced Pandas DataFrame, and pickle it for later
    More on fitACF here: https://radar-software-toolkit-rst.readthedocs.io/en/latest/references/general/fitacf/
    Parameters of Interest:
        Scalar Parameters
            - Time (Saved as epoch - seconds since 1970.01.01)
                    This is a combination of year 'time.yr', month 'time.mo', day 'time.dy', 
                    hour 'time.hr', minute 'time.mt', second 'time.sc', and microsecond 'time.us'
            - Beam Number 'bmnum' (starting at Beam 0)
            - Transmitted Frequency 'tfreq' (in KHz)
            - Station Identifier 'stid'
        
        Vector Parameters
            - Range Gate 'slist'
            - Ground Scatter Flag 'qflg'
            - Quality Flag 'gflg'
            - Velocity 'v' with Error 'v_e'
            - Power (actually SNR) from lambda fit 'p_l' with Error 'p_l_e'
            - Spectral width from lambda fit (w_l) 'w_l' with Error 'w_l_e'
            - Elevation angle 'elv' with low 'elv_low' and high 'elv_high' estimates
    """

    station = "rkn"
    date = "20161012"

    pattern = '%Y.%m.%d %H:%M:%S.%f'  # This is the pattern we will use to convert time info to epoch
    # Loop through all the files for this station/date and pickle
    in_dir = "data/" + station + date
    for in_file in glob.iglob(in_dir + "/*.fitacf.bz2"):

        # Unpack and open the file
        print("\nReading " + in_file + "...")
        with bz2.open(in_file) as fp:
            fitacf_stream = fp.read()
        sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
        fitacf_data = sdarn_read.read_fitacf()

        # print("List of available parameters: " + str(fitacf_data[0].keys()))

        # Pre-allocate lists for the scalar parameters
        times = []  # time is the name of a library used herein
        bmnum = []
        tfreq = []
        # stid will be the same for the whole event, so we will build it at the end

        # Pre-allocate lists for the vector parameters
        slist = []
        qflg = []
        gflg = []
        v = []
        v_e = []
        p_l = []
        p_l_e = []
        w_l = []
        w_l_e = []
        elv = []
        elv_low = []
        elv_high = []

        # Loop through every record and build parallel arrays
        print("Building Parallel Arrays...")
        for record in range(len(fitacf_data)):

            # Convert date and time to epoch
            # Note the conversion of us to s (the 0. are stripped and only the after decimal part is used)
            date_time = str(fitacf_data[record]['time.yr']) + "." + str(fitacf_data[record]['time.mo']) + "." \
                        + str(fitacf_data[record]['time.dy']) + " " + str(fitacf_data[record]['time.hr']) + ":" \
                        + str(fitacf_data[record]['time.mt']) + ":" + str(fitacf_data[record]['time.sc']) + "." \
                        + str((fitacf_data[record]['time.us'] / 1e6))[2:]
            epoch = calendar.timegm(time.strptime(date_time, pattern))  # Doesn't keep microseconds

            # Loop through all the gates
            try:
                num_gates_reporting = len(fitacf_data[record]['slist'])
            except:
                # Sometimes there won't be any gates reporting, this is legacy behaviour and results in partial records
                num_gates_reporting = 0

            for element in range(num_gates_reporting):
                # Build up the vectored data
                slist.append(fitacf_data[record]['slist'][element])
                qflg.append(fitacf_data[record]['qflg'][element])
                gflg.append(fitacf_data[record]['gflg'][element])
                v.append(fitacf_data[record]['v'][element])
                v_e.append(fitacf_data[record]['v_e'][element])
                p_l.append(fitacf_data[record]['p_l'][element])
                p_l_e.append(fitacf_data[record]['p_l_e'][element])
                w_l.append(fitacf_data[record]['w_l'][element])
                w_l_e.append(fitacf_data[record]['w_l_e'][element])
                elv.append(fitacf_data[record]['elv'][element])
                elv_low.append(fitacf_data[record]['elv_low'][element])
                elv_high.append(fitacf_data[record]['elv_high'][element])

                # Build up the scalar data, these values are the same for every gate
                bmnum.append(fitacf_data[record]['bmnum'])
                tfreq.append(fitacf_data[record]['tfreq'])
                times.append(epoch)

        # Put the data into a dataframe
        print("Building and filtering the data frame...")
        stid = [radar_code_from_station_id(fitacf_data[0]['stid'])] * len(times)
        df = pd.DataFrame({'stid': stid,
                           'time': times,
                           'bmnum': bmnum,
                           'gate': slist,
                           'tfreq': tfreq,
                           'qflg': qflg,
                           'gflg': gflg,
                           'vel': v,
                           'vel_err': v_e,
                           'pwr': p_l,
                           'pwr_err': p_l_e,
                           'wdt': w_l,
                           'wdt_err': w_l_e,
                           'elv': elv,
                           'elv_low': elv_low,
                           'elv_high': elv_high})

        # Filter with the ground and quality flags
        df = df.loc[(df['gflg'] == 0) & (df['qflg'] == 1) & (df['pwr'] >= 3.0)]

        # Drop the ground scatter and quality flags and re-index
        df.drop(['gflg', 'qflg'], inplace=True, axis=1)
        df.reset_index(drop=True)

        # Save to file
        # Strip the "xx.xxx.fitacf.bz2" and replace with the station code and pickle file extension
        out_file = in_file[:len(in_file) - 17]
        out_file += radar_code_from_station_id(fitacf_data[0]['stid']) + ".pkl"
        print("Pickling as " + out_file + "...")
        df.to_pickle(out_file)
