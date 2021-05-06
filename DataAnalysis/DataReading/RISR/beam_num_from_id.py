import glob
import math

import h5py


def beam_num_from_id(id):
    """
    :param id: The RISR-N beam id
    :return: The RISR-N beam number (1 - 11) from the world day mode.  May return Nan if the beam id is not recognized.
    """
    if id == 55748:
        return 1
    elif id == 57656:
        return 2
    elif id == 56954:
        return 3
    elif id == 57782:
        return 4
    elif id == 60551:
        return 5
    elif id == 60617:
        return 6
    elif id == 60683:
        return 7
    elif id == 63443:
        return 8
    elif id == 64280:
        return 9
    elif id == 63587:
        return 10
    elif id == 65486:
        return 11
    else:
        return math.nan


if __name__ == '__main__':
    station = "ran"
    date = "20161012"

    # Loop through all long pulse files and test the beam number conversion
    in_dir = "data/" + station + date
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
            for beam in file['/Data/Array Layout/'].keys():
                a = 1
                print("Beam id: " + str(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][0]))
                print("Beam number: " + str(beam_num_from_id(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][0])))

    # Test something that is not valid to make sure we get NaN back
    print("\n")
    print(beam_num_from_id(56))
    print(math.isnan(beam_num_from_id(56)))
