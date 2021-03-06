import h5py

if __name__ == '__main__':
    """
    Read in and look at the parameters in an HDF5 file representing the World Day experiment
    Long Pulse Uncorrected Ne
    Doesn't produce anything, just for looking at the file structure
    """
    # in_file = "ran2020728.3"
    in_file = "ran2009128.8"

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

    # The file has both data and metadata
    print("\nFile Keys: " + str(file.keys()))

    # Look at the metadata
    meta_data = file['/Metadata']
    print("\nMetadata keys: " + str(meta_data.keys()))
    data_params = file['/Metadata/Data Parameters']
    print(data_params)
    print(data_params[()])
    exp_notes = file['/Metadata/Experiment Notes']
    print(exp_notes)
    # print(exp_notes[()])
    exp_params = file['/Metadata/Experiment Parameters']
    print(exp_params)
    # print(exp_params[()])

    # Look at the data
    data = file['/Data']
    print("\nData keys: " + str(data.keys()))

    data_table = file['/Data/Table Layout']
    print("Table Layout Keys: " + str(data_table))
    print(data_table[0])

