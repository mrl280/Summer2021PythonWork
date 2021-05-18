import h5py

if __name__ == '__main__':
    """
    Read in and look at the parameters in an HDF5 file representing the World Day experiment
    Vector velocity from Long Pulse
    Doesn't produce anything, just for looking at the file structure
    """
    # in_file = "ran2018612.4"
    in_file = "ran2017515.6"
    # in_file = "ran2020728.3"

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
    data_layout_desc = file['/Data/Array Layout/Layout Description']
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

    data_table = file['/Data/Table Layout']  # The data doesn't seem to be as well organized in the table layout
    # print(data_table[()])

    data_array = file['/Data/Array Layout']
    print("Array Layout Keys: " + str(data_array.keys()))

    data_cgm_lat = file['/Data/Array Layout/cgm_lat']
    print(data_cgm_lat)
    data_time = file['/Data/Array Layout/timestamps']
    print(str(data_time) + "\n")

    # Here are all the 1D parameters
    data_1d_params = file['/Data/Array Layout/1D Parameters']
    print("1D parameter keys: " + str(data_1d_params.keys()))
    data_1Dparam_params = file['/Data/Array Layout/1D Parameters/Data Parameters']
    print(data_1Dparam_params[()])
    data_1Dparam_altb = file['/Data/Array Layout/1D Parameters/altb']
    print(data_1Dparam_altb)
    data_1Dparam_alte = file['/Data/Array Layout/1D Parameters/alte']
    print(data_1Dparam_alte)
    data_1Dparam_inttms = file['/Data/Array Layout/1D Parameters/inttms']
    print(data_1Dparam_inttms)

    # Here are all the 2d parameters
    data_2d_params = file['/Data/Array Layout/2D Parameters']
    print("\n2D parameter keys: " + str(data_2d_params.keys()))
    data_2d_params_params = file['/Data/Array Layout/2D Parameters/Data Parameters']
    print(data_2d_params_params[()])
    data_2Dparam_vi1 = file['/Data/Array Layout/2D Parameters/vi1']
    print(data_2Dparam_vi1)

    file.close()
