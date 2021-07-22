import pathlib
import bz2

import pandas as pd
import pydarn
import _pickle as cPickle

from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.basic_SD_df_filter import basic_SD_df_filter
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch


def handler(df, time_units, area):
    """

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units
    :param area:

    :return df: pandas.DataFrame:
        The restricted dataframe for the area of interest
    """
    if area == 1:
        return pick_out_area_of_interest_1(df=df.copy(), time_units=time_units)
    elif area == 2:
        return pick_out_area_of_interest_2(df=df.copy(), time_units=time_units)
    elif area == 3:
        return pick_out_area_of_interest_3(df=df.copy(), time_units=time_units)
    elif area == 4:
        return pick_out_area_of_interest_4(df=df.copy(), time_units=time_units)
    elif area == 5:
        return pick_out_area_of_interest_5(df=df.copy(), time_units=time_units)

    # After the first 5 areas were cut out and considered individually, it was decided that the third area should be
    #  further split into 3 subareas.
    elif area == "3a":
        return pick_out_area_of_interest_3a(df=df.copy(), time_units=time_units)
    elif area == "3b":
        return pick_out_area_of_interest_3b(df=df.copy(), time_units=time_units)
    elif area == "3c":
        return pick_out_area_of_interest_3c(df=df.copy(), time_units=time_units)

    else:
        raise Exception("Area not recognized, the handler doesn't know who to call")


def pick_out_area_of_interest_1(df, time_units):
    """

    Pick out area 1 from the Sept 26, 2106 event.

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 1
    """

    # Break it into frequencies and handle them each individually
    df10 = df[(df['transFreq_MHz'] == 10)].copy()
    df10 = df10.loc[(df10['gate'] >= 40)]
    df10 = df10.loc[(df10[time_units] >= 0.6) & (df10[time_units] <= 1.5)]

    df12 = df[(df['transFreq_MHz'] == 12)].copy()
    df12 = df12.loc[(df12['gate'] >= 50)]
    df12 = df12.loc[(df12[time_units] >= 0.5) & (df12[time_units] <= 1.25)]

    df13 = df[(df['transFreq_MHz'] == 13)].copy()
    df13 = df13.loc[(df13['gate'] >= 60)]
    df13 = df13.loc[(df13[time_units] >= 0.75) & (df13[time_units] <= 1.25)]

    # This area is not seen in 14 MHz data
    df14 = df[(df['transFreq_MHz'] == 14)].copy()
    df14 = df14.loc[(df14['gate'] >= 60)]
    df14 = df14.loc[(df14[time_units] >= 0.75) & (df14[time_units] <= 1.25)]

    """"""

    df = pd.concat([df10, df12, df13, df14])

    return df


def pick_out_area_of_interest_2(df, time_units):
    """

    Pick out area 2 from the Sept 26, 2106 event.

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 2
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1['gate'] >= 30) & (df10_1['gate'] <= 42)]
    df10_1 = df10_1.loc[(df10_1[time_units] >= 1.15) & (df10_1[time_units] <= 1.4)]

    df10_2 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_2 = df10_2.loc[(df10_2['gate'] >= 26) & (df10_2['gate'] < 30)]
    df10_2 = df10_2.loc[(df10_2[time_units] >= 1.25) & (df10_2[time_units] <= 1.45)]

    df10_3 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_3 = df10_3.loc[(df10_3['gate'] >= 18) & (df10_3['gate'] < 26)]
    df10_3 = df10_3.loc[(df10_3[time_units] >= 1.35) & (df10_3[time_units] <= 1.5)]

    df10 = pd.concat([df10_1, df10_2, df10_3])

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1['gate'] >= 30) & (df12_1['gate'] <= 51)]
    df12_1 = df12_1.loc[(df12_1[time_units] >= 0.95) & (df12_1[time_units] <= 1.4)]

    df12_2 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_2 = df12_2.loc[(df12_2['gate'] >= 25) & (df12_2['gate'] < 30)]
    df12_2 = df12_2.loc[(df12_2[time_units] >= 1.2) & (df12_2[time_units] <= 1.474)]

    df12_3 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_3 = df12_3.loc[(df12_3['gate'] >= 21) & (df12_3['gate'] < 25)]
    df12_3 = df12_3.loc[(df12_3[time_units] >= 1.35) & (df12_3[time_units] <= 1.474)]

    df12 = pd.concat([df12_1, df12_2, df12_3])

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1['gate'] >= 42) & (df13_1['gate'] <= 52)]
    df13_1 = df13_1.loc[(df13_1[time_units] >= 0.9) & (df13_1[time_units] <= 1.25)]

    df13_2 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_2 = df13_2.loc[(df13_2['gate'] >= 35) & (df13_2['gate'] < 42)]
    df13_2 = df13_2.loc[(df13_2[time_units] >= 1.0) & (df13_2[time_units] <= 1.35)]

    df13_3 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_3 = df13_3.loc[(df13_3['gate'] >= 31) & (df13_3['gate'] < 35)]
    df13_3 = df13_3.loc[(df13_3[time_units] >= 1.1) & (df13_3[time_units] <= 1.4)]

    df13_4 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_4 = df13_4.loc[(df13_4['gate'] >= 20) & (df13_4['gate'] < 31)]
    df13_4 = df13_4.loc[(df13_4[time_units] >= 1.2) & (df13_4[time_units] <= 1.5)]

    df13 = pd.concat([df13_1, df13_2, df13_3, df13_4])

    """"""

    df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_1 = df14_1.loc[(df14_1['gate'] >= 42) & (df14_1['gate'] <= 52)]
    df14_1 = df14_1.loc[(df14_1[time_units] >= 0.9) & (df14_1[time_units] <= 1.25)]

    df14_2 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_2 = df14_2.loc[(df14_2['gate'] >= 35) & (df14_2['gate'] < 42)]
    df14_2 = df14_2.loc[(df14_2[time_units] >= 1.0) & (df14_2[time_units] <= 1.35)]

    df14_3 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_3 = df14_3.loc[(df14_3['gate'] >= 31) & (df14_3['gate'] < 35)]
    df14_3 = df14_3.loc[(df14_3[time_units] >= 1.1) & (df14_3[time_units] <= 1.4)]

    df14_4 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_4 = df14_4.loc[(df14_4['gate'] >= 20) & (df14_4['gate'] < 31)]
    df14_4 = df14_4.loc[(df14_4[time_units] >= 1.2) & (df14_4[time_units] <= 1.5)]

    df14 = pd.concat([df14_1, df14_2, df14_3, df14_4])

    """"""

    df = pd.concat([df10, df12, df13, df14])

    return df


def pick_out_area_of_interest_3(df, time_units):
    """

    Pick out area 3 from the Sept 26, 2106 event.

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 3
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1[time_units] >= 1.5) & (df10_1[time_units] <= 1.6)]
    df10_1 = df10_1.loc[(df10_1['gate'] >= 18) & (df10_1['gate'] <= 33)]

    df10_2 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_2 = df10_2.loc[(df10_2[time_units] > 1.6) & (df10_2[time_units] <= 1.7)]
    df10_2 = df10_2.loc[(df10_2['gate'] >= 15) & (df10_2['gate'] <= 30)]

    df10_3 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_3 = df10_3.loc[(df10_3[time_units] > 1.7) & (df10_3[time_units] <= 3.5)]
    df10_3 = df10_3.loc[(df10_3['gate'] >= 5) & (df10_3['gate'] <= 30)]

    df10 = pd.concat([df10_1, df10_2, df10_3])

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1[time_units] >= 1.5) & (df12_1[time_units] <= 1.6)]
    df12_1 = df12_1.loc[(df12_1['gate'] >= 17) & (df12_1['gate'] <= 35)]

    df12_2 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_2 = df12_2.loc[(df12_2[time_units] > 1.6) & (df12_2[time_units] <= 3.35)]
    df12_2 = df12_2.loc[(df12_2['gate'] >= 5) & (df12_2['gate'] <= 33)]

    df12_3 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_3 = df12_3.loc[(df12_3[time_units] > 3.35) & (df12_3[time_units] <= 3.5)]
    df12_3 = df12_3.loc[(df12_3['gate'] >= 5) & (df12_3['gate'] <= 25)]

    df12 = pd.concat([df12_1, df12_2, df12_3])

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1[time_units] >= 1.5) & (df13_1[time_units] <= 1.6)]
    df13_1 = df13_1.loc[(df13_1['gate'] >= 17) & (df13_1['gate'] <= 40)]

    df13_2 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_2 = df13_2.loc[(df13_2[time_units] > 1.6) & (df13_2[time_units] <= 3.35)]
    df13_2 = df13_2.loc[(df13_2['gate'] >= 5) & (df13_2['gate'] <= 33)]

    df13_3 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_3 = df13_3.loc[(df13_3[time_units] > 3.35) & (df13_3[time_units] <= 3.5)]
    df13_3 = df13_3.loc[(df13_3['gate'] >= 5) & (df13_3['gate'] <= 25)]

    df13 = pd.concat([df13_1, df13_2, df13_3])

    """"""

    # df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    # df14_1 = df14_1.loc[(df14_1[time_units] >= 1.4) & (df14_1[time_units] <= 1.5)]
    # df14_1 = df14_1.loc[(df14_1['gate'] >= 32) & (df14_1['gate'] <= 40)]

    df14_2 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_2 = df14_2.loc[(df14_2[time_units] > 1.5) & (df14_2[time_units] <= 3.0)]
    df14_2 = df14_2.loc[(df14_2['gate'] >= 5) & (df14_2['gate'] <= 39)]

    df14_3 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_3 = df14_3.loc[(df14_3[time_units] > 3.0) & (df14_3[time_units] <= 3.5)]
    df14_3 = df14_3.loc[(df14_3['gate'] >= 5) & (df14_3['gate'] <= 30)]

    # df14 = pd.concat([df14_1, df14_2, df14_3])
    df14 = pd.concat([df14_2, df14_3])

    """"""

    df = pd.concat([df10, df12, df13, df14])

    return df


def pick_out_area_of_interest_4(df, time_units):
    """

    Pick out area 4 from the Sept 26, 2106 event.

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 4
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1['gate'] >= 40) & (df10_1['gate'] <= 53)]
    df10_1 = df10_1.loc[(df10_1[time_units] >= 3.1) & (df10_1[time_units] <= 3.5)]

    df10_2 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_2 = df10_2.loc[(df10_2['gate'] >= 30) & (df10_2['gate'] < 40)]
    df10_2 = df10_2.loc[(df10_2[time_units] >= 3.2) & (df10_2[time_units] <= 3.7)]

    df10_3 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_3 = df10_3.loc[(df10_3['gate'] >= 17) & (df10_3['gate'] < 30)]
    df10_3 = df10_3.loc[(df10_3[time_units] >= 3.45) & (df10_3[time_units] <= 3.9)]

    df10 = pd.concat([df10_1, df10_2, df10_3])

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1['gate'] >= 35) & (df12_1['gate'] <= 50)]
    df12_1 = df12_1.loc[(df12_1[time_units] >= 3.1) & (df12_1[time_units] <= 3.7)]

    df12_2 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_2 = df12_2.loc[(df12_2['gate'] >= 25) & (df12_2['gate'] < 35)]
    df12_2 = df12_2.loc[(df12_2[time_units] >= 3.35) & (df12_2[time_units] <= 3.8)]

    df12_3 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_3 = df12_3.loc[(df12_3['gate'] >= 17) & (df12_3['gate'] < 25)]
    df12_3 = df12_3.loc[(df12_3[time_units] >= 3.5) & (df12_3[time_units] <= 3.9)]

    df12 = pd.concat([df12_1, df12_2, df12_3])

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1['gate'] >= 32) & (df13_1['gate'] <= 50)]
    df13_1 = df13_1.loc[(df13_1[time_units] >= 3.1) & (df13_1[time_units] <= 3.7)]

    df13_2 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_2 = df13_2.loc[(df13_2['gate'] >= 25) & (df13_2['gate'] < 32)]
    df13_2 = df13_2.loc[(df13_2[time_units] >= 3.3) & (df13_2[time_units] <= 3.8)]

    df13_3 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_3 = df13_3.loc[(df13_3['gate'] >= 17) & (df13_3['gate'] < 25)]
    df13_3 = df13_3.loc[(df13_3[time_units] >= 3.5) & (df13_3[time_units] <= 3.9)]

    df13 = pd.concat([df13_1, df13_2, df13_3])

    """"""

    df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_1 = df14_1.loc[(df14_1['gate'] >= 39) & (df14_1['gate'] <= 50)]
    df14_1 = df14_1.loc[(df14_1[time_units] >= 3.1) & (df14_1[time_units] <= 3.55)]

    df14_2 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_2 = df14_2.loc[(df14_2['gate'] >= 25) & (df14_2['gate'] < 39)]
    df14_2 = df14_2.loc[(df14_2[time_units] >= 3.2) & (df14_2[time_units] <= 3.8)]

    df14_3 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_3 = df14_3.loc[(df14_3['gate'] >= 17) & (df14_3['gate'] < 25)]
    df14_3 = df14_3.loc[(df14_3[time_units] >= 3.55) & (df14_3[time_units] <= 3.9)]

    df14 = pd.concat([df14_1, df14_2, df14_3])

    df = pd.concat([df10, df12, df13, df14])

    return df


def pick_out_area_of_interest_5(df, time_units):
    """

    Pick out area 5 from the Sept 26, 2106 event.
    This area was kind of an after thought

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 5
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1[time_units] >= 0.5) & (df10_1[time_units] <= 0.6)]
    df10_1 = df10_1.loc[(df10_1['gate'] >= 10) & (df10_1['gate'] <= 22)]

    df10_2 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_2 = df10_2.loc[(df10_2[time_units] > 0.6) & (df10_2[time_units] <= 1.0)]
    df10_2 = df10_2.loc[(df10_2['gate'] >= 5) & (df10_2['gate'] <= 23)]

    df10_3 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_3 = df10_3.loc[(df10_3[time_units] > 1.0) & (df10_3[time_units] <= 1.1)]
    df10_3 = df10_3.loc[(df10_3['gate'] >= 5) & (df10_3['gate'] < 12)]

    df10 = pd.concat([df10_1, df10_2, df10_3])

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1[time_units] >= 0.5) & (df12_1[time_units] <= 0.6)]
    df12_1 = df12_1.loc[(df12_1['gate'] >= 11) & (df12_1['gate'] <= 21)]

    df12_2 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_2 = df12_2.loc[(df12_2[time_units] > 0.6) & (df12_2[time_units] <= 1.0)]
    df12_2 = df12_2.loc[(df12_2['gate'] >= 5) & (df12_2['gate'] <= 20)]

    df12_3 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_3 = df12_3.loc[(df12_3[time_units] > 1.0) & (df12_3[time_units] <= 1.1)]
    df12_3 = df12_3.loc[(df12_3['gate'] >= 5) & (df12_3['gate'] < 13)]

    df12 = pd.concat([df12_1, df12_2, df12_3])

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1[time_units] >= 0.5) & (df13_1[time_units] <= 0.6)]
    df13_1 = df13_1.loc[(df13_1['gate'] >= 10) & (df13_1['gate'] <= 21)]

    df13_2 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_2 = df13_2.loc[(df13_2[time_units] > 0.6) & (df13_2[time_units] <= 0.95)]
    df13_2 = df13_2.loc[(df13_2['gate'] >= 5) & (df13_2['gate'] <= 21)]

    df13_3 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_3 = df13_3.loc[(df13_3[time_units] > 0.95) & (df13_3[time_units] <= 1.1)]
    df13_3 = df13_3.loc[(df13_3['gate'] >= 5) & (df13_3['gate'] < 13)]

    df13 = pd.concat([df13_1, df13_2, df13_3])

    """"""

    df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_1 = df14_1.loc[(df14_1[time_units] >= 0.45) & (df14_1[time_units] <= 0.6)]
    df14_1 = df14_1.loc[(df14_1['gate'] >= 10) & (df14_1['gate'] <= 23)]

    df14_2 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_2 = df14_2.loc[(df14_2[time_units] > 0.6) & (df14_2[time_units] <= 0.95)]
    df14_2 = df14_2.loc[(df14_2['gate'] >= 5) & (df14_2['gate'] <= 22)]

    df14_3 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_3 = df14_3.loc[(df14_3[time_units] > 0.95) & (df14_3[time_units] <= 1.1)]
    df14_3 = df14_3.loc[(df14_3['gate'] >= 5) & (df14_3['gate'] <= 13)]

    df14 = pd.concat([df14_1, df14_2, df14_3])

    """"""

    df = pd.concat([df10, df12, df13, df14])

    return df


def pick_out_area_of_interest_3a(df, time_units):
    """

    Pick out area 3a from the Sept 26, 2106 event.
    This area is characterized by high velocity echoes

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 3a
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1[time_units] >= 1.5) & (df10_1[time_units] <= 1.65)]
    df10_1 = df10_1.loc[(df10_1['gate'] >= 19) & (df10_1['gate'] <= 34)]

    df10_2 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_2 = df10_2.loc[(df10_2[time_units] > 1.65) & (df10_2[time_units] <= 2.0)]
    df10_2 = df10_2.loc[(df10_2['gate'] >= 7) & (df10_2['gate'] <= 30)]

    df10_3 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_3 = df10_3.loc[(df10_3[time_units] > 2.0) & (df10_3[time_units] <= 2.4)]
    df10_3 = df10_3.loc[(df10_3['gate'] >= 7) & (df10_3['gate'] <= 22)]

    df10 = pd.concat([df10_1, df10_2, df10_3])

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1[time_units] >= 1.45) & (df12_1[time_units] <= 1.6)]
    df12_1 = df12_1.loc[(df12_1['gate'] >= 20) & (df12_1['gate'] <= 33)]

    df12_2 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_2 = df12_2.loc[(df12_2[time_units] > 1.6) & (df12_2[time_units] <= 2.4)]
    df12_2 = df12_2.loc[(df12_2['gate'] >= 8) & (df12_2['gate'] <= 27)]

    df12 = pd.concat([df12_1, df12_2])

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1[time_units] >= 1.45) & (df13_1[time_units] <= 1.6)]
    df13_1 = df13_1.loc[(df13_1['gate'] >= 19) & (df13_1['gate'] <= 33)]

    df13_2 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_2 = df13_2.loc[(df13_2[time_units] > 1.6) & (df13_2[time_units] <= 2.4)]
    df13_2 = df13_2.loc[(df13_2['gate'] >= 8) & (df13_2['gate'] <= 28)]

    df13 = pd.concat([df13_1, df13_2])

    """"""

    df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_1 = df14_1.loc[(df14_1[time_units] >= 1.45) & (df14_1[time_units] <= 1.6)]
    df14_1 = df14_1.loc[(df14_1['gate'] >= 19) & (df14_1['gate'] <= 36)]

    df14_2 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_2 = df14_2.loc[(df14_2[time_units] > 1.6) & (df14_2[time_units] <= 2.4)]
    df14_2 = df14_2.loc[(df14_2['gate'] >= 8) & (df14_2['gate'] <= 28)]

    df14 = pd.concat([df14_1, df14_2])

    """"""

    df = pd.concat([df10, df12, df13, df14])

    return df


def pick_out_area_of_interest_3b(df, time_units):
    """

    Pick out area 3b from the Sept 26, 2106 event.
    This area is characterized by lower velocity echoes

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 3b
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1[time_units] >= 2.5) & (df10_1[time_units] < 2.85)]
    df10_1 = df10_1.loc[(df10_1['gate'] >= 5) & (df10_1['gate'] <= 25)]

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1[time_units] >= 2.5) & (df12_1[time_units] < 2.85)]
    df12_1 = df12_1.loc[(df12_1['gate'] >= 5) & (df12_1['gate'] <= 30)]

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1[time_units] >= 2.5) & (df13_1[time_units] < 2.8)]
    df13_1 = df13_1.loc[(df13_1['gate'] >= 7) & (df13_1['gate'] <= 35)]

    """"""

    df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_1 = df14_1.loc[(df14_1[time_units] >= 2.45) & (df14_1[time_units] < 2.8)]
    df14_1 = df14_1.loc[(df14_1['gate'] >= 7) & (df14_1['gate'] <= 35)]

    """"""

    df = pd.concat([df10_1, df12_1, df13_1, df14_1])

    return df


def pick_out_area_of_interest_3c(df, time_units):
    """

    Pick out area 3c from the Sept 26, 2106 event.
    This area is characterized by high velocity echoes

    :param df: pandas.DataFrame:
        Data for the Sept 26, 2106 event
    :param time_units: str:
        Time units

    :return df: pandas.DataFrame:
        Just the data for area 3c
    """

    # Break it into frequencies and handle them each individually
    df10_1 = df[(df['transFreq_MHz'] == 10)].copy()
    df10_1 = df10_1.loc[(df10_1[time_units] >= 3.0) & (df10_1[time_units] < 3.5)]
    df10_1 = df10_1.loc[(df10_1['gate'] >= 6) & (df10_1['gate'] <= 29)]

    """"""

    df12_1 = df[(df['transFreq_MHz'] == 12)].copy()
    df12_1 = df12_1.loc[(df12_1[time_units] >= 2.95) & (df12_1[time_units] < 3.5)]
    df12_1 = df12_1.loc[(df12_1['gate'] >= 7) & (df12_1['gate'] <= 29)]

    """"""

    df13_1 = df[(df['transFreq_MHz'] == 13)].copy()
    df13_1 = df13_1.loc[(df13_1[time_units] >= 2.95) & (df13_1[time_units] < 3.5)]
    df13_1 = df13_1.loc[(df13_1['gate'] >= 7) & (df13_1['gate'] <= 30)]

    """"""

    df14_1 = df[(df['transFreq_MHz'] == 14)].copy()
    df14_1 = df14_1.loc[(df14_1[time_units] >= 2.9) & (df14_1[time_units] < 3.5)]
    df14_1 = df14_1.loc[(df14_1['gate'] >= 7) & (df14_1['gate'] <= 30)]

    """"""

    df = pd.concat([df10_1, df12_1, df13_1, df14_1])

    return df


if __name__ == "__main__":
    """ Testing/Running """

    station = "rkn"
    year = "2016"
    month = "09"
    day = "26"

    time_units = 'ut'
    start_hour = 0
    end_hour = 4

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + ".pbz2"
    data_stream = bz2.BZ2File(in_file, "rb")
    df = cPickle.load(data_stream)

    # Restrict data to within the desired hour range
    _, start_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=start_hour)
    _, end_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=end_hour)
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    # Add in decimal time
    radar_id = pydarn.read_hdw_file(station).stid
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id,
                                date_time_est=(df['datetime'].iat[0]).to_pydatetime())

    # Filter for beam and gate ranges of interest
    # All events use data within the same gate and beam range
    beam_range = (7, 7)

    # Near range gates can sometimes be funky.  And because of the way elevation angle is resolved,
    #  the nearest gates can't have valid elevation angles in 15 km mode
    gate_range = (5, 74)
    df = df.loc[((df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1])) &
                ((df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1]))]

    # Make sure we only use 15 km data
    df = df.loc[(df['firstRang'] == 90) & (df['rangeSep'] == 15)]
    df.reset_index(drop=True, inplace=True)

    # Filter out ground, low quality, and low power
    df = basic_SD_df_filter(df)

    # Put frequencies in MHz and round, this makes them easier to compare
    df['transFreq_MHz'] = round(df['transFreq'] * 1e-3, 0)

    loc_root = str((pathlib.Path().parent.absolute()))
    out_dir = loc_root + "/data"

    # Go through each area, cutting them out and saving to file
    areas = [1, 2, 3, 4, 5, "3a", "3b", "3c"]
    dfs = dict()
    for area in areas:
        limited_df = handler(df=df, time_units=time_units, area=area)

        out_file = out_dir + "/" + station + year + month + day + "_area" + str(area) + ".pbz2"

        print("Saving dataframe as " + out_file)
        with bz2.BZ2File(out_file, "w") as file:
            cPickle.dump(limited_df, file)
