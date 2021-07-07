#  This code is to be inserted around line 70 of plot_solar_data.py

# df = df.loc[(df['decimal_year'] >= year_range[0]) & (df['decimal_year'] <= year_range[1])]
# df.reset_index(drop=True, inplace=True)
#
# # print(len(df))  # 731
#
# # Remove some elements
# flux = df['observed_flux'].copy()  # len 5900
# flux[100:110] = np.asarray([math.nan] * len(flux[100:110]))
# df['observed_flux'] = flux
#
# # Plot the raw data
# ax.plot(df['decimal_year'], df['observed_flux'], color="cornflowerblue", linestyle='-', linewidth=0.5)
#
# df = df.loc[df['observed_flux'].notna()]
#
# print(len(df))
# print(df.head())