import pydarn
import bz2


if __name__ == '__main__':
    """
    Read in and look at a SuperDARN data file
    """
    fitacf_file = "data/rkn/rkn20190318/20190318.0001.00.rkn.fitacf.bz2"
    with bz2.open(fitacf_file) as fp:
        fitacf_stream = fp.read()

    sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
    fitacf_data = sdarn_read.read_fitacf()  # fitacf_data is a list of keyed dictionaries

    print(fitacf_data[0].keys())

    # print(fitacf_data[:]['slist'])

    slist = []  # an example vector parameter
    year = []  # an example scalar parameter
    for i in range(5):
        num_gates_reporting = len(fitacf_data[i]['slist'])
        year.extend([fitacf_data[i]['time.yr']] * num_gates_reporting)
        slist.extend(fitacf_data[i]['slist'])
        # slist = slist + fitacf_data[i]['slist']

    print(year)
    print(slist)

    print(len(year))
    print(len(slist))

    # print(len(fitacf_data[31]['slist']))
    # # print(len(fitacf_data[32]['slist']))
    # print(fitacf_data[33]['slist'])
    #
    # for record in range(0):
    #     print("\nTime: " + str(fitacf_data[record]['time.mt']) + " " + str(fitacf_data[record]['time.sc']))
    #     print("Beam: " + str(fitacf_data[record]['bmnum']))
    #     print("Vel: " + str(fitacf_data[record]['v']))
    #     print("Ground scatter flag: " + str(fitacf_data[record]['gflg']))
    #     print("Gate: " + str(fitacf_data[record]['slist']))
    #     print("Power: " + str(fitacf_data[record]['p_l']))
    #     print("Quality Flag Array: " + str(fitacf_data[record]['qflg']))
