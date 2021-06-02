import pydarn
import bz2
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters

if __name__ == '__main__':
    fitacf_file = "data/rkn/rkn20160926/20160926.1001.00.rkn.fitacf.bz2"
    with bz2.open(fitacf_file) as fp:
        fitacf_stream = fp.read()

    sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
    fitacf_data = sdarn_read.read_fitacf()
    register_matplotlib_converters()
    print(fitacf_data[0].keys())

    pydarn.Fan.plot_fan(fitacf_data)
    # plt.title("Radar {:d}, Beam {:d}".format(fitacf_data[0]['stid'], fitacf_data[0]['bmnum']))

    plt.show()
