from lib.get_data_handler import get_data_handler


if __name__ == '__main__':
    """ 
    Check to see if the data we are looking for in a pickle read is actually in the file or not.
    
    Pickle files don't always have all the data for month, day, gate, or beam, so this program exists to investigate 
     the file, and find out if you have everything you need
    """

    station = "dcn"
    year_range = (2019, 2021)
    local_testing = False

    df = get_data_handler(station, year_range=year_range, month_range=None, day_range=None,
                          gate_range=None, beam_range=None, freq_range=None, occ_data=True, local_testing=local_testing)

    print(df.keys())

    print(df.head())

    starting_datetime = df['datetime'].iat[0]
    print("The first datetime in the df is:")
    print(starting_datetime)

    ending_datetime = df['datetime'].iat[-1]
    print("The final datetime in the df is:")
    print(ending_datetime)

    year = []
    month = []
    for i in range(len(df)):
        datetime_obj = df['datetime'].iat[i]
        year.append(datetime_obj.year)
        month.append(datetime_obj.month)

    df['year'] = year
    df['month'] = month

    years = [2019, 2020, 2021]
    months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    for y in years:
        print("Now looking at year=" + str(y))
        df_yy = df.loc[(df['year'] == y)]

        for m in months:
            df_mm = df_yy.loc[(df_yy['month'] == m)]
            print("     For month " + str(m) + " len(df) =" + str(len(df_mm)))




