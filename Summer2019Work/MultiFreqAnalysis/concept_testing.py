import matplotlib.pyplot as plt
import numpy as np
import idlsave
import os


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


def test6():
    for i in range(5):
        print(i)


if __name__ == "__main__":
    test6()
    print(1/60)

