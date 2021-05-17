# Michael Luciuk
# Sept 14, 2019

# Read in columns of data from text file


def getData(in_file):
    """
    Purpose:
        Read in columned data from text file
    Pre-conditions:
        :param: in_file: must be valid file in Data folder
    Post-conditions:
        none
    Return:
        An list of lists containing gate, time, freq, and velocity data
    """

    # Read in data from FILE
    file = open("Data/" + in_file + ".txt", "r")
    lines = file.readlines()
    gate = []
    time = []
    freq = []
    vel = []

    # Go through each line and add data to arrays:
    for line in lines:
        gate.append(float(line.split()[0]))
        time.append(float(line.split()[1]))
        freq.append(float(line.split()[2]))
        vel.append(float(line.split()[3]))

    file.close()

    return [gate, time, freq, vel]


if __name__ == "__main__":
    result = getData("20160905rkn_1.5-8.0 UT")
    # print(len(result[0]))
