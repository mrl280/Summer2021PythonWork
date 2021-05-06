import matplotlib.pyplot as plt
import numpy as np
import idlsave
import os


def test1():
    myPlot = plt.figure()
    x = np.array([1, 7, 20, 50, 80])
    y = np.array([10, 19, 30, 35, 51])
    plt.scatter(x, y)
    fit = np.polyfit(np.log(x), y, 1)
    dummy_x = np.arange(1, 100, 1)
    plt.plot(dummy_x, fit[0]*np.log(dummy_x)+fit[1])

    plt.show()


def test2():
    myPlot = plt.figure()
    x = np.array([-1, -7, -20, -50, -80])
    y = np.array([10, 19, 30, 35, 51])
    plt.scatter(x, y)
    fit = np.polyfit(np.log(np.negative(x)), y, 1)
    dummy_x = np.arange(-100, 0, 1)
    plt.plot(dummy_x, fit[0]*np.log(np.negative(dummy_x))+fit[1])

    plt.show()


def test3():
    myPlot = plt.figure()
    x = np.array([-1, -7, -20, -50, -80, -150])
    y = np.array([-10, -19, -30, -35, -51, -50])
    plt.scatter(x, y)
    fit = np.polyfit(np.power(x, 2), y, 1)
    dummy_x = np.arange(-100, 0, 1)
    plt.plot(dummy_x, fit[0]*np.power(dummy_x, 2)+fit[1])

    plt.show()


def test4():
    myPlot = plt.figure()
    x = np.array([-1, -7, -20, -50, -80, -150])
    y = np.array([-10, -19, -30, -35, -51, -50])
    plt.scatter(x, y)
    fit = np.polyfit(1/x, y, 1)
    dummy_x = np.arange(-150, 0, 1)
    plt.plot(dummy_x, fit[0]/dummy_x+fit[1])

    plt.show()



def test5():
    # When testing on the Python console use:
    # import idlsave
    # s = idlsave.read("Data/20160306rkn.sav")

    # try to read in an IDL .sav file
    cur_path = os.path.dirname(__file__)  # where we are

    s = idlsave.read("Data/20160306rkn.sav")
    vel = s.fit['v_e']
    beam = s.prm['bmnum']
    freq = s.prm['tfreq']
    print(len(vel))
    print(len(beam))
    print(freq[1:5])


if __name__ == "__main__":
    test5()
