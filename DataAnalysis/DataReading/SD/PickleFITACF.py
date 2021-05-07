import calendar
import pydarn
import bz2
import pandas as pd
import time
import glob

if __name__ == '__main__':
    """
    Take a fitACF file, turn it into a reduced Pandas DataFrame, and pickle it for later
    More on fitACF here: https://radar-software-toolkit-rst.readthedocs.io/en/latest/references/general/fitacf/
    
    Takes a days worth of fitACF data and combine it all into one pickled data frame.
    
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

    # Pre-allocate lists for the scalar parameters
    # Station id will be the same for the whole event, so we will build it at the end
    times = []  # time is the name of a library used herein
    bmnum = []
    tfreq = []

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

    # Loop through all the files for this station/date
    in_dir = "data/" + station + date
    for in_file in glob.iglob(in_dir + "/*.fitacf.bz2"):

        # Unpack and open the file
        print("Reading " + in_file + "...")
        with bz2.open(in_file) as fp:
            fitacf_stream = fp.read()
        sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
        fitacf_data = sdarn_read.read_fitacf()

        # print("List of available parameters: " + str(fitacf_data[0].keys()))

        # Loop through every record in this file and add to the parallel arrays
        print("Adding to Parallel Arrays...")
        for record in range(len(fitacf_data)):

            # Convert date and time to epoch
            # Note the conversion of us to s:
            #   - the 0. are stripped and only the after decimal part is used)
            #   - we have to ensure it doesn't flip to scientific notation otherwise time.strptime can't handle it
            #   - not really necessary because calendar.timegm doesn't keep microseconds anyway
            date_time = str(fitacf_data[record]['time.yr']) + "." + str(fitacf_data[record]['time.mo']) + "." \
                        + str(fitacf_data[record]['time.dy']) + " " + str(fitacf_data[record]['time.hr']) + ":" \
                        + str(fitacf_data[record]['time.mt']) + ":" + str(fitacf_data[record]['time.sc']) + "." \
                        + str(f"{fitacf_data[record]['time.us'] / 1e6:.6f}")[2:]
            epoch = calendar.timegm(time.strptime(date_time, pattern))  # Doesn't keep microseconds

            try:
                num_gates_reporting = len(fitacf_data[record]['slist'])
            except:
                # Sometimes there won't be any gates reporting, this is legacy behaviour and results in partial records
                num_gates_reporting = 0

            # Loop through all the gates reporting
            for gate_idx in range(num_gates_reporting):
                # Build up the vectored data
                slist.append(fitacf_data[record]['slist'][gate_idx])
                qflg.append(fitacf_data[record]['qflg'][gate_idx])
                gflg.append(fitacf_data[record]['gflg'][gate_idx])
                v.append(fitacf_data[record]['v'][gate_idx])
                v_e.append(fitacf_data[record]['v_e'][gate_idx])
                p_l.append(fitacf_data[record]['p_l'][gate_idx])
                p_l_e.append(fitacf_data[record]['p_l_e'][gate_idx])
                w_l.append(fitacf_data[record]['w_l'][gate_idx])
                w_l_e.append(fitacf_data[record]['w_l_e'][gate_idx])
                elv.append(fitacf_data[record]['elv'][gate_idx])
                elv_low.append(fitacf_data[record]['elv_low'][gate_idx])
                elv_high.append(fitacf_data[record]['elv_high'][gate_idx])

                # Build up the scalar data, these values are the same for every gate
                # These could be pulled out of this loop but it is simpler to leave them here
                bmnum.append(fitacf_data[record]['bmnum'])
                tfreq.append(fitacf_data[record]['tfreq'])
                times.append(epoch)

    # Put the data into a dataframe
    print("Building and filtering the data frame...")
    stid = [station] * len(times)
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
    df.reset_index(drop=True, inplace=True)

    # Save to file
    out_file = in_dir + "/" + station + date + ".pkl"
    print("Pickling as " + out_file + "...")
    df.to_pickle(out_file)
