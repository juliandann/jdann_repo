import arcpy,sys, os, string, glob 
import os
import numpy as np
import pandas as pd
from os.path import isfile, join


# Check out any necessary licenses  
arcpy.CheckOutExtension("spatial")
# Overwrite pre-existing files  
arcpy.env.overwriteOutput = True  
#Get the raster datasets in the input workspace and loop through them from the start

def zonal_stats(shapefiles,raster_path):
    out_raster = arcpy.Raster(raster_path)
    #shapefiles_teller = [s for s in shapefiles if 'TL' in s]
    #shapefiles_teller_lidar = [s for s in shapefiles_teller if 'TL_SAR_41' not in s]
    #print shapefiles_teller_lidar

    for shape in shapefiles:
        
        
        if float(shape_str) < 63335:       
            #table output
            out_tbl = shape_str+'_Zstats.dbf'
            print out_tbl
            z = arcpy.gp.ZonalStatisticsAsTable_sa(shape, "Counter" ,out_raster , out_tbl, "DATA","ALL")

def dbf_to_csv():
    #print arcpy.ListFiles("*.dbf")
    
    #find dbf files
    for dbf_file in arcpy.ListFiles("*Zstats.dbf"):
        dbf_str = float(dbf_file.replace('_Zstats.dbf',''))
        if (dbf_str >31251) & (dbf_str<65001):
            outLocation = 'Z:/AKSeward/Data/GIS/Teller/2019_Snow/zonal_stats_onbuffers/'
            outName = dbf_file.replace('.dbf','.csv')
            arcpy.TableToTable_conversion(dbf_file, outLocation, outName)
            print outName
def concat_csv(path):
    #getting filenames of csvs
    files = []
    for file in os.listdir(path):      
        if file.endswith(".csv"):  
            files.append(join(path,file))
    df_final = pd.DataFrame()
    for f in files:
        df = pd.read_csv(f)
        df_final = pd.concat([df,df_final])
    print df_final
    df_final.to_csv(path+'all_baptiste_snoD_vs_insitu_3m_Zstats_nans.csv')
    
#setting working environment
#arcpy.env.workspace = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/Buffers_may_2017_hydrosense_pts/individual_3m_buffer_polygons/'
arcpy.env.workspace = 'Z:/AKSeward/Data/GIS/Teller/2019_Snow/individual_2019_snow_buffer_shapefiles/'
arcpy.env.extent = 'MINOF'

shapefiles = arcpy.ListFeatureClasses()
#print shapefiles
raster = 'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/imagesnoD_nans.tif'

#zonal_stats(shapefiles,raster)
dbf_to_csv()
#path = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/STD_Slope/'
path = 'Z:/AKSeward/Data/GIS/Teller/2019_Snow/zonal_stats_onbuffers/'
concat_csv(path)

#shapefilepath = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/individual_3m_buffers/'


'''
files = []
for file in os.listdir(shapefilepath):      
    if file.endswith(".shp"):  
        files.append(join(shapefilepath,file))

print files
'''
