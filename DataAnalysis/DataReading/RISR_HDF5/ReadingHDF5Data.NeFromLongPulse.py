import h5py

if __name__ == '__main__':
    """
    Read in and look at the parameters in an HDF5 file representing the World Day experiment
    Long Pulse Uncorrected Ne
    Doesn't produce anything, just for looking at the file structure
    """
    in_file = "ran2020728.3"

    station = in_file[0:3]
    in_dir = station + "/" + in_file[:-2]
    path = "data/" + in_dir + "/" + in_file
    file = h5py.File(path + ".h5", "r")

    # The file type is held in the CKINDAT parameter of experimental notes
    exp_notes = file['/Metadata/Experiment Notes']
    if station == "ran":
        file_type = str(exp_notes[25][0])
    elif station == "ras":
        file_type = str(exp_notes[39][0])
    else:
        raise Exception("Station " + station + " not recognized")
    file_type = file_type[10:len(file_type) - 1].rstrip()  # Strip out everything but the file type name
    print("File type: " + file_type)

    # Look at the data layout description
    data_layout_desc = file['/Data/Array Layout/Array with beamid=56954 /Layout Description']
    print(data_layout_desc[()])

    # The file has both data and metadata
    print("\nFile Keys: " + str(file.keys()))

    # Look at the metadata
    meta_data = file['/Metadata']
    print("\nMetadata keys: " + str(meta_data.keys()))
    data_params = file['/Metadata/Data Parameters']
    print(data_params)
    # print(data_params[()])
    exp_notes = file['/Metadata/Experiment Notes']
    print(exp_notes)
    # print(exp_notes[()])
    exp_params = file['/Metadata/Experiment Parameters']
    print(exp_params)
    # print(exp_params[()])

    # Look at the data
    data = file['/Data']
    print("\nData keys: " + str(data.keys()))

    data_array = file['/Data/Array Layout']
    print("Array Layout Keys: " + str(data_array.keys()))

    data_array_beam56954 = file['/Data/Array Layout/Array with beamid=56954 ']
    print("Beam 56954 Keys: " + str(data_array_beam56954.keys()))
    data_range = file['/Data/Array Layout/Array with beamid=56954 /range']
    print(data_range)
    data_time = file['/Data/Array Layout/Array with beamid=56954 /timestamps']
    print(str(data_time) + "\n")

    # Here are all the 1D parameters
    data_1d_params = file['/Data/Array Layout/Array with beamid=56954 /1D Parameters']
    print("1D parameter keys: " + str(data_1d_params.keys()))
    data_1Dparam_params = file['/Data/Array Layout/Array with beamid=56954 /1D Parameters/Data Parameters']
    print(data_1Dparam_params[()])
    data_1Dparam_beamid = file['/Data/Array Layout/Array with beamid=56954 /1D Parameters/beamid']
    print(data_1Dparam_beamid)

    # Here are all the 2d parameters
    data_2d_params = file['/Data/Array Layout/Array with beamid=56954 /2D Parameters']
    print("\n2D parameter keys: " + str(data_2d_params.keys()))
    data_2Dparam_params = file['/Data/Array Layout/Array with beamid=56954 /2D Parameters/Data Parameters']
    print(data_2Dparam_params[()])
    data_2Dparam_dpopl = file['/Data/Array Layout/Array with beamid=56954 /2D Parameters/dpopl']
    print(data_2Dparam_dpopl)

    file.close()

