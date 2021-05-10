import calendar
import pydarn
import bz2
import pandas as pd
import time
import glob


def PickleFITACF(station, date):
    """
    Take a fitACF file, turn it into a reduced Pandas DataFrame, and pickle it for later
    More on fitACF here: https://radar-software-toolkit-rst.readthedocs.io/en/latest/references/general/fitacf/
    
    Takes a days worth of fitACF data and combine it all into one pickled data frame.
    
    Parameters of Interest:
        Scalar Parameters
            - Epoch time (seconds since 1970.01.01)
                    This is a combination of year 'time.yr', month 'time.mo', day 'time.dy', 
                    hour 'time.hr', minute 'time.mt', second 'time.sc', and microsecond 'time.us'
            - Year
            - Month
            - Day
            - Hour
            - Minute
            - Second (included microseconds)
            - Beam azimuth 'bmazm'
            - Integration Seconds (include microseconds) 'intt' 
            - Beam Number 'bmnum' (starting at Beam 0)
            - Transmitted Frequency 'tfreq' (in KHz)
            - Station Identifier 'stid'
            - Distance to the first range gate 'frang'
            - Range separation 'rsep' 
            - Control program name or command line, or comment buffer 'combf' 
            - Major version number of the FITACF algorithm 'fitacf.revision.major'
            - Minor version number of the FITACF algorithm 'fitacf.revision.minor'
        
        Vector Parameters
            - Range Gate 'slist'
            - Ground Scatter Flag 'qflg'
            - Quality Flag 'gflg'
            - Velocity 'v' with Error 'v_e'
            - Power (actually SNR) from lambda fit 'p_l' with Error 'p_l_e'
            - Spectral width from lambda fit (w_l) 'w_l' with Error 'w_l_e'
            - Standard deviation of lambda fit 'sd_l' and phase fit 'sd_phi'
            - Elevation angle 'elv' with low 'elv_low' and high 'elv_high' estimates
    """

    pattern = '%Y.%m.%d %H:%M:%S'  # This is the pattern we will use to convert time info to epoch

    # Pre-allocate lists for the scalar parameters
    # Station id will be the same for the whole event, so we will build it at the end
    epoch = []
    year = []
    month = []
    day = []
    hour = []
    minute = []
    second = []
    bmnum = []
    bmazm = []
    intt = []
    tfreq = []
    frang = []
    rsep = []
    combf = []
    fitACF_rev_major = []
    fitACF_rev_minor = []

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
    sd_l = []
    sd_phi = []
    elv = []
    elv_low = []
    elv_high = []

    # Loop through all the files for this station/date
    in_dir = "data/" + station + "/" + station + date
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
            date_time = str(fitacf_data[record]['time.yr']) + "." + str(fitacf_data[record]['time.mo']) + "." \
                        + str(fitacf_data[record]['time.dy']) + " " + str(fitacf_data[record]['time.hr']) + ":" \
                        + str(fitacf_data[record]['time.mt']) + ":" + str(fitacf_data[record]['time.sc'])
            epoch_here = calendar.timegm(time.strptime(date_time, pattern))  # Doesn't keep microseconds

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
                sd_l.append(fitacf_data[record]['sd_l'][gate_idx])
                sd_phi.append(fitacf_data[record]['sd_phi'][gate_idx])
                elv.append(fitacf_data[record]['elv'][gate_idx])
                elv_low.append(fitacf_data[record]['elv_low'][gate_idx])
                elv_high.append(fitacf_data[record]['elv_high'][gate_idx])

                # Build up the scalar data, these values are the same for every gate
                # These could be pulled out of this loop but it is simpler to leave them here
                epoch.append(epoch_here)
                year.append(fitacf_data[record]['time.yr'])
                month.append(fitacf_data[record]['time.mo'])
                day.append(fitacf_data[record]['time.dy'])
                hour.append(fitacf_data[record]['time.hr'])
                minute.append(fitacf_data[record]['time.mt'])
                second.append(fitacf_data[record]['time.sc'] + fitacf_data[record]['time.us'] / 1e6)
                bmnum.append(fitacf_data[record]['bmnum'])
                bmazm.append(fitacf_data[record]['bmazm'])
                intt.append(fitacf_data[record]['intt.sc'] + fitacf_data[record]['intt.us'] / 1e6)
                tfreq.append(fitacf_data[record]['tfreq'])
                frang.append(fitacf_data[record]['frang'])
                rsep.append(fitacf_data[record]['rsep'])
                fitACF_rev_major.append(fitacf_data[record]['fitacf.revision.major'])
                fitACF_rev_minor.append(fitacf_data[record]['fitacf.revision.minor'])
                combf.append(fitacf_data[record]['combf'])

    # Put the data into a dataframe
    print("Building the data frame...")
    stid = [station] * len(epoch)
    df = pd.DataFrame({'stationId': stid,
                       'epoch': epoch,
                       'year': year,
                       'month': month,
                       'day': day,
                       'hour': hour,
                       'minute': minute,
                       'second': second,
                       'beammNumber': bmnum,
                       'beamAzimuth': bmazm,
                       'intgrtnTime': intt,
                       'gate': slist,
                       'transFreq': tfreq,
                       'firstRang': frang,
                       'rangeSep': rsep,
                       'fitACFMajorRev': fitACF_rev_major,
                       'fitACFMinorRev': fitACF_rev_minor,
                       'CtrlPrgrm': combf,

                       'qflg': qflg,
                       'gflg': gflg,
                       'vel': v,
                       'vel_err': v_e,
                       'pwr': p_l,
                       'pwr_err': p_l_e,
                       'wdt': w_l,
                       'wdt_err': w_l_e,
                       'stdDevLambda': sd_l,
                       'stdDevPhase': sd_phi,
                       'elv': elv,
                       'elv_low': elv_low,
                       'elv_high': elv_high
                       })

    # Save to file
    out_file = in_dir + "/" + station + date + ".pkl"
    print("Pickling as " + out_file + "...")
    df.to_pickle(out_file)


if __name__ == '__main__':
    station = "rkn"
    date = "20161012"
    PickleFITACF(station, date)

    # Filter with the ground and quality flags
    # df = df.loc[(df['gflg'] == 0) & (df['qflg'] == 1) & (df['pwr'] >= 3.0)]

    # Drop the ground scatter and quality flags and re-index
    # df.drop(['gflg', 'qflg'], inplace=True, axis=1)
    # df.reset_index(drop=True, inplace=True)

