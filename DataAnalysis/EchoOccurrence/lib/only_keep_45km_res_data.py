import warnings


def only_keep_45km_res_data(df):
    """
    Restrict to 45 km mode data.
    Warns if non-45 km data is being removed

    :param df: pandas.DataFrame:
            SuperDARN dataframe with 'frang' and 'rsep' parameters.  If these columns don't exist nothing will happen
            and a warning will be printing.
    :return: pandas.DataFrame:
            The provided dataframe, except with only 45 km resolution data
    """

    number_of_records_before_spatial_resolution_check = df.shape[0]

    try:
        df = df.loc[(df['frang'] == 180) & (df['rsep'] == 45)]
        df.reset_index(drop=True, inplace=True)
    except KeyError:
        warnings.warn("only_keep_45km_res_data() was not able to restrict your data because the columns 'frang' "
                      "and 'rsep' could not be found.", category=Warning)
    except BaseException as e:
        raise e

    number_of_records_after_spatial_resolution_check = df.shape[0]
    if number_of_records_before_spatial_resolution_check > number_of_records_after_spatial_resolution_check:
        warnings.warn("Not all data within the specified range is 45 km resolution data.  "
                      "All non-45 km data was discarded.", category=Warning)
    return df
