# Read in and show a figure

from PIL import Image as Im
import os


def ReadInFig():
    VelScatFile = "06032016_Velocity_Comparison_19.1-19.19UT"
    cur_path = os.path.dirname(__file__)  # where we are
    with Im.open(cur_path + '/Processed_data/' + VelScatFile + '.eps') as image:
        width, height = image.size
        image.show()


if __name__ == "__main__":
    ReadInFig()
