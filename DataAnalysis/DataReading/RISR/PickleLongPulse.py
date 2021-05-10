import pandas as pd
import glob
import h5py

from DataAnalysis.DataReading.RISR.wd_beam_num import wd_beam_num

if __name__ == '__main__':
    """
    Take a RISR-N long pulse file and pickle it
    Note: Not tested with RISR-C data
    
    Parameters of Interest:
        Epoch time
        Range
    
        1D Parameters (parameters that vary only with time)
            Beam id 'bmid' (this is usually a 5 digit number that uniquely identifies a RISR beam)
            World-Day Beam Number 'bmnum' (1 through 11)
            Elevation Angle 'elv' (0=horizontal,90=vert)
            Transmission Frequency 'tfreq'

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

    station = "ran"
    date = "201432"

    # Loop through all the files for this station/date and pickle all Long Pulse files
    in_dir = "data/" + station + "/" + station + date
    for in_file in glob.iglob(in_dir + "/*.h5"):

        try:
            file = h5py.File(in_file, "r")
        except:
            continue

        # Look at the file type, we are only interested in Long Pulse Files here
        exp_notes = file['/Metadata/Experiment Notes']
        file_type = str(exp_notes[25][0])
        file_type = file_type[10:len(file_type) - 1].rstrip()  # Strip out everything but the file type name
        if file_type != "Long Pulse":
            print("\n" + in_file + " is a " + file_type + " file, skipping it...")
        else:
            # The file is a long pulse file, we are good to go ahead
            print("\n" + in_file + " is a " + file_type + " file, reading in data...")

            # Pre-allocate list for 1D parameters (only vary with time)
            epoch = []
            ranges = []
            bmid = []
            bmnum = []
            elv = []
            tfreq = []

            # Pre-allocate lists for 2D parameters (vary with range and time)
            log10_Ne = []
            log10_Ne_err = []
            ion_temp = []
            ion_temp_err = []
            e_temp = []
            e_temp_err = []
            los_ion_vel = []
            los_ion_vel_err = []
            cgm_lat = []
            cgm_lon = []
            comp = []
            comp_err = []
            alt = []

            # Loop through all the beams
            for beam in file['/Data/Array Layout/'].keys():
                print("Looping through beam " + str(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][0]))

                # Loop through all the times and ranges here
                for time_idx, time in enumerate(file['/Data/Array Layout/' + beam + '/timestamps']):
                    for range_idx, range in enumerate(file['/Data/Array Layout/' + beam + '/range']):
                        epoch.append(time)
                        ranges.append(range)

                        # The 1D parameters are just accessed in time
                        # These Could be moved out of this range loop but it is simpler to leave them here)
                        bmid.append(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][time_idx])
                        bmnum.append(
                            wd_beam_num(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][time_idx]))
                        elv.append(file['/Data/Array Layout/' + beam + '/1D Parameters/elm'][time_idx])
                        tfreq.append(file['/Data/Array Layout/' + beam + '/1D Parameters/tfreq'][time_idx])

                        # The 2D parameters are accessed in time and range
                        log10_Ne.append(file['/Data/Array Layout/' + beam + '/2D Parameters/nel'][range_idx][time_idx])
                        log10_Ne_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dnel'][range_idx][time_idx])
                        ion_temp.append(file['/Data/Array Layout/' + beam + '/2D Parameters/ti'][range_idx][time_idx])
                        ion_temp_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dti'][range_idx][time_idx])
                        e_temp.append(file['/Data/Array Layout/' + beam + '/2D Parameters/te'][range_idx][time_idx])
                        e_temp_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dte'][range_idx][time_idx])
                        los_ion_vel.append(file['/Data/Array Layout/' + beam + '/2D Parameters/vo'][range_idx][time_idx])
                        los_ion_vel_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dvo'][range_idx][time_idx])
                        cgm_lat.append(file['/Data/Array Layout/' + beam + '/2D Parameters/cgm_lat'][range_idx][time_idx])
                        cgm_lon.append(file['/Data/Array Layout/' + beam + '/2D Parameters/cgm_long'][range_idx][time_idx])
                        comp.append(file['/Data/Array Layout/' + beam + '/2D Parameters/po+'][range_idx][time_idx])
                        comp_err.append(file['/Data/Array Layout/' + beam + '/2D Parameters/dpo+'][range_idx][time_idx])
                        alt.append(file['/Data/Array Layout/' + beam + '/2D Parameters/gdalt'][range_idx][time_idx])

            # Put the data into a dataframe
            print("Building data frame...")
            df = pd.DataFrame(
                {'stationId': [station] * len(epoch),
                 'epoch': epoch,
                 'range': ranges,
                 'bmId': bmid,
                 'wdBmnum': bmnum,
                 'elv': elv,
                 'transFreq': tfreq,
                 'log10Ne': log10_Ne,
                 'log10NeErr': log10_Ne_err,
                 'ionTemp': ion_temp,
                 'ionTempErr': ion_temp_err,
                 'eTemp': e_temp,
                 'eTempErr': e_temp_err,
                 'losIonVel': los_ion_vel,
                 'losIonVelErr': los_ion_vel_err,
                 'cgmLat': cgm_lat,
                 'cgmLon': cgm_lon,
                 'comp': comp,
                 'compErr': comp_err,
                 'alt': alt,
                 })

            # Save to file
            out_file = in_file[:len(in_file) - 4] + "LongPulse.pkl"  # Strip the "x.h5" and replace with "LongPulse.pkl"
            print("Pickling as " + out_file + "...")
            df.to_pickle(out_file)

        file.close()
