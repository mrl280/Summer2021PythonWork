import numpy as np
from matplotlib import pyplot as plt

""" 
Test using numpy's least squares function (numpy.linalg.lstsq)

Following the example here: https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
"""

x_data = np.asarray(
    [18.4085, 15.50605, 12.05065, 6.745565, 9.06777, 18.49445, 16.8882, 6.3907, 12.89735, 6.448525, 12.4599,
     20.53535, 0.544503, 10.03475, 7.65597, 11.84405, 10.05775, 10.03465, 0.792501, 9.51731, 2.33615, 1.67687,
     5.79429])

y_data = np.asarray(
    [12.55, 9.01, 16., 3.595, 4.61, 11.45, 12.95, 4.72, 7.64, 4.565, 8., 10.085, 0.418, 6.985, 4.33, 9.06, 6.86,
     4.75, 1.1015, 6.32, 2.335, 1.25, 5.74])

fig, ax = plt.subplots(figsize=(8, 9), sharex='col', dpi=300, nrows=1, ncols=1)
ax.set_ylim((0, 30))
ax.set_xlim((0, 30))

ax.plot(x_data, y_data, 'ro')
ax.plot([ax.get_ylim()[0], ax.get_ylim()[1]], [ax.get_xlim()[0], ax.get_xlim()[1]],
        linestyle='--', linewidth=2, color="purple")

A = np.vstack([x_data, np.ones(len(x_data))]).T
m, c = np.linalg.lstsq(A, y_data, rcond=None)[0]
ax.plot(x_data, m * x_data + c, 'g', label='Fitted line')

ax.legend(loc='upper left')

plt.show()
