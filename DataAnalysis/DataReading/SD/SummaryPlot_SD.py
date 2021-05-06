import pydarn
import bz2
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters

if __name__ == '__main__':
    fitacf_file = "data/20190318.1001.00.rkn.fitacf.bz2"
    with bz2.open(fitacf_file) as fp:
        fitacf_stream = fp.read()

    sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
    fitacf_data = sdarn_read.read_fitacf()
    register_matplotlib_converters()

    pydarn.RTP.plot_summary(fitacf_data, beam_num=2, slant=False)

    plt.show()
