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
    shapefiles_teller = [s for s in shapefiles if 'TL' in s]
    shapefiles_teller_lidar = [s for s in shapefiles_teller if 'TL_SAR_41' not in s]
    print shapefiles_teller_lidar
    for shape in shapefiles_teller_lidar:
        shape_str = shape.replace('.shp','')
        
        #table output
        out_tbl = shape_str+'_Zstats.dbf'
        print out_tbl
        z = arcpy.gp.ZonalStatisticsAsTable_sa(shape, "name" ,out_raster , out_tbl, "DATA","ALL")

def dbf_to_csv():
    #print arcpy.ListFiles("*.dbf")
    
    #find dbf files
    for dbf_file in arcpy.ListFiles("*Zstats.dbf"):
        outLocation = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/STD_Slope/'
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
    df_final.to_csv(path+'all_may_2017_slope_0.5mDTM_3mbuffer_Zstats.csv')
    
#setting working environment
arcpy.env.workspace = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/Buffers_may_2017_hydrosense_pts/individual_3m_buffer_polygons/'
shapefiles = arcpy.ListFeatureClasses()
raster = 'Z:/AKSeward/DEM/Lidar/2017_UAS/DEMs/JBD_Products/DTM_JBD/Focal_statistics/TL27_0.5m_DTM_slope_focal_std_3m_window.tif'

#zonal_stats(shapefiles,raster)
#dbf_to_csv()
path = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/STD_Slope/'
concat_csv(path)

#shapefilepath = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/individual_3m_buffers/'


'''
files = []
for file in os.listdir(shapefilepath):      
    if file.endswith(".shp"):  
        files.append(join(shapefilepath,file))

print files
'''
