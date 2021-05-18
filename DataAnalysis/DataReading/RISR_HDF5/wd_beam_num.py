import glob
import math

import h5py


def wd_beam_num(beam_id):
    """
    Note: Applicable only when RISR_HDF5 is operating in World Day mode
    These beam numbers ar the same for both RISR_HDF5-N and RISR_HDF5-C
    :param beam_id: The RISR_HDF5 beam id
    :return: The RISR_HDF5 beam number (1 - 11) from the world-day mode.  May return Nan if the beam id is not recognized.
    """
    if beam_id == 55748:
        return 1
    elif beam_id == 57656:   # This WD RISR_HDF5-N beam overlaps with rkn beam 5 (gates 30-40)
        return 2
    elif beam_id == 56954:
        return 3
    elif beam_id == 57782:
        return 4
    elif beam_id == 60551:
        return 5
    elif beam_id == 60617:
        return 6
    elif beam_id == 60683:
        return 7
    elif beam_id == 63443:
        return 8
    elif beam_id == 64280:
        return 9
    elif beam_id == 63587:   # This WD RISR_HDF5-N beam overlaps with inv beam 12 (gates 31-38)
        return 10
    elif beam_id == 65486:
        return 11
    else:
        return math.nan


if __name__ == '__main__':
    """
    Testing
    """
    station = "ran"
    date = "20161012"

    # Loop through all long pulse files for this station/date and test the beam number conversion
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
                print("Beam number: " + str(
                    wd_beam_num(file['/Data/Array Layout/' + beam + '/1D Parameters/beamid'][0])))

    # Test something that is not valid to make sure we get NaN back
    if wd_beam_num(56) is not math.nan:
        print("Error: WD beam id 56 should not have a beam number")
