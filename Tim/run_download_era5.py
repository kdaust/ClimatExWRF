#this script includes the variables needed for wrf
import cdsapi
import sys
import os.path

c = cdsapi.Client()

# takes an input date for the last day of the week you want to include
if len(sys.argv) == 5:
    download_directory = sys.argv[1]
    download_year = sys.argv[2]   
    download_month = sys.argv[3]
    download_land = sys.argv[4]
   
else:
    raise Exception("Invalid input entries. Needs DOWNLOAD_DIRECTORY DOWNLOAD_YEAR DOWNLOAD_MONTH DOWNLOAD_LAND=Y,N")


variables_single = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                    '2m_temperature', 'mean_sea_level_pressure', 'sea_surface_temperature','skin_temperature', 'surface_pressure',
                    'snow_depth','snow_density',
                    'soil_temperature_level_1','soil_temperature_level_2','soil_temperature_level_3','soil_temperature_level_4',
                    'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3','volumetric_soil_water_layer_4',
                    'sea_ice_cover','land_sea_mask','geopotential','total_precipitation']

variables_pressure = ['geopotential', 'relative_humidity',
                      'u_component_of_wind', 'v_component_of_wind', 'temperature', 
                      'vertical_velocity']


#variables_land = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
#                    '2m_temperature', 'skin_temperature', 'surface_pressure',
#                    'snow_depth_water_equivalent']

variables_land = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                    '2m_temperature', 'skin_temperature', 'surface_pressure',
                    'snow_depth_water_equivalent',
                    'soil_temperature_level_1','soil_temperature_level_2','soil_temperature_level_3','soil_temperature_level_4',
                    'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3','volumetric_soil_water_layer_4']

months_year = ['01','02','03','04','05','06','07','08','09','10','11','12']
months_yearm1 = ['09','10','11','12']
months_yearp1 = ['01']


days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
        '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']

#days = ['10','11','12','13','14','15','16','17','18','19','20','21','22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
#days = ['27', '28', '29', '30', '31']


times_3hr = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
times_1hr = ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00',
             '13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00']
    
    
def retrieve_pressure(years,months,days,savename,times,download_directory):
    year_name = str(years)
    month_name = str(months)
    day_name = str(days)
    hour_name = str(times)
    if not os.path.exists(download_directory + '/' + 'ERA5_pressure_' + year_name + "_" + month_name + "_" + day_name +"_"+hour_name  + '.grib'):
        c.retrieve('reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'variable': variables_pressure,
            'year': years,
            'month': months,
            'day': days,
            'time': times,
            'pressure_level': ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', 
                               '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000'],
            'area': [90, -180, 25, -85],
            'format': 'grib',
        },
        download_directory + '/' + 'ERA5_pressure_' + year_name + "_" + month_name + "_" + day_name +"_"+hour_name  + '.grib')
    else:
        print("File exists; continuing")

def retrieve_single(years,months,days,savename,times,download_directory):
    year_name = str(years)
    month_name = str(months)
    day_name = str(days)
    hour_name = str(times)
    if not os.path.exists(download_directory + '/' + 'ERA5_single_' + year_name + "_" + month_name + "_" + day_name +"_"+hour_name  + '.grib'):
        c.retrieve('reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': variables_single,
            'year': years,
            'month': months,
            'day': days,
            'time': times,
            'area': [90, -180, 25, -85],
            'format': 'grib',
        },
        download_directory + '/' + 'ERA5_single_' + year_name + "_" + month_name + "_" + day_name +"_"+hour_name  + '.grib')
    else:
        print("File exists; continuing")

def retrieve_land(years,months,days,savename,times,download_directory):
    year_name = str(years)
    month_name = str(months)
    day_name = str(days)
    hour_name = str(times)
    if not os.path.exists(download_directory + '/' + 'ERA5_land_' + year_name + "_" + month_name + "_" + day_name +"_"+hour_name  + '.grib'):
        c.retrieve('reanalysis-era5-land', # changed from land to single levels
        {
            'product_type': 'reanalysis',
            'variable': variables_land,
            'year': years,
            'month': months,
            'day': days,
            'time': times,
            'area': [90, -180, 25, -85],
            'format': 'grib',
        },
        download_directory + '/' + 'ERA5_land_' + year_name + "_" + month_name + "_" + day_name +"_"+hour_name  + '.grib')
    else:
        print("File exists; continuing")


########## Main ###############

# Main simulation year

year = download_year
month = download_month

for day in days:
    try:
        for time in times_3hr:
            print("Getting pressure data for {}-{}-{}_{} to {}".format(year,month,day,time,download_directory))
            retrieve_pressure(year,month,day,"00",time,download_directory)
            print("Getting single data for {}-{}-{}_{} to {}".format(year,month,day,time,download_directory))
            retrieve_single(year,month,day,"00",time,download_directory)
            if download_land == "Y":
                print("Getting land data for {}-{}-{}_{} to {}".format(year,month,day,time,download_directory))
                retrieve_land(year,month,day,"00",time,download_directory)

    except:
        break


"""
year = startyear_input
for month in months_year:
    for day in days:
        for time in times_3hr:
            print("Getting pressure data for {}-{}-{}_{}".format(year,month,day,time))
            retrieve_pressure(year,month,day,"00",time,download_directory))
            print("Getting single data for {}-{}-{}_{}".format(year,month,day,time))
            retrieve_single(year,month,day,"00",time,download_directory)

# January of following year
year = str(int(startyear_input)+1)
for month in months_yearp1:
    for day in days:
        for time in times_3hr:
            print("Getting pressure data for yearp1 {}-{}-{}_{}".format(year,month,day,time))
            retrieve_pressure(year,month,day,"00",time,download_directory))
            print("Getting single data for yearp1 {}-{}-{}_{}".format(year,month,day,time))
            retrieve_single(year,month,day,"00",time,download_directory)

# Simulation start year
year = str(int(startyear_input)-1)
for month in months_yearm1:
    for day in days:
        for time in times_3hr:
            print("Getting pressure data for yearm1 {}-{}-{}_{}".format(year,month,day,time))
            retrieve_pressure(year,month,day,"00",time,download_directory))
            print("Getting single data for yearm1 {}-{}-{}_{}".format(year,month,day,time))
            retrieve_single(year,month,day,"00",time,download_directory)

            if month == months_yearm1[0]:
                print("Getting land data for yearm1 {}-{}-{}_{}".format(year,month,day,time))
                retrieve_land(year,month,day,"00",time,download_directory)
                
"""
