import matplotlib
import _pickle

print('matplotlib: {}'.format(matplotlib.__version__))

# Check to see what pickle protocols are available
print('_pickle highest protocol: {}'.format(_pickle.HIGHEST_PROTOCOL))