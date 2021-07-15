

def basic_SD_df_filter(df):
    """
    Performs the basic standard filtering one usually needs to perform on SuperDARN data.
    Operation is completed inplace.

    Only keeps SuperDARN data where:
        - Ground scatter flag is 0
        - Quality flag is 1
        - Power is >= 3 dB

    :param df: pandas.DataFrame:
        SuperDARN data frame
    :return df: pandas.DataFrame:
        SuperDARN data frame filtered for ground scatter, bad quality data, and low power data. Data is re-indexed.
    """
    try:
        df = df.loc[(df['gflg'] == 0) & (df['qflg'] == 1) & (df['pwr'] >= 3.0)]
        df.reset_index(drop=True, inplace=True)
        return df
    except Exception as e:
        print("Error in basic_SD_df_filter(): ")
        raise e
