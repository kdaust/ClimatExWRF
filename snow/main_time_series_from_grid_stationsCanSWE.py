"""
Script to extract time series based on lat lon for WRF, CaSR and ERA5Land.
"""
import pandas as pd

from climatex import process_all_stations, process_all_stations_neighbors

if __name__ == '__main__':

    base_path = '/climatex/'
    #station_timeseries = process_stations(summary_df, base_path)
    summary_df = pd.read_csv(base_path+'stations_summary.csv', index_col=0)
    source = 'era5'

    # Process all stations
    summary_df['station_id'] = summary_df.index
    # snow_data = process_all_stations(summary_df, source='wrf', base_path=base_path)
    # snow_data = process_all_stations(summary_df, source=source,
    #                                   base_path='climatex/casr_dataset/year/')

    filepath='climatex/era5_dataset/swe_daily_1999-2023.nc'
    snow_data = process_all_stations(summary_df, source=source, filepath=filepath)

    # Access individual station data
    for station_id, df in snow_data.items():
        print(f"\n{station_id} summary:")
        print(df.describe())

        # Optionally save to CSV
        df.to_csv(base_path+source+'_data/'+f'{station_id}_snow_timeseries.csv')
