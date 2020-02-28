import pandas as pd
import numpy as np
import matplotlib.pyplot
import h5py
from osgeo import gdal
import osr
from SAR_gen_functions import *
from scipy import stats
import matplotlib.colors as mcolors

def mat_to_dataframe(mat_file,x='',y=''):
    #converts mat_files to
    f = h5py.File(mat_file,'r')
    keys = [key for key in f.keys()]
    if (len(x) == 0) and (len(y)==0):
        for key in keys:
            #find x or y matrices
            if ('x' or 'lon') in key:
                y = pd.Series(f.get(key).value.flatten())
                print(y)
            if ('y' or 'lat') in key:
                x = pd.Series(f.get(key).value.flatten())
                print(x)

    places = ['x','y','lat','lon']
    for key in keys:

        if all(x not in key for x in places):
            print(key)
            try:
                var1
            except NameError:
                print(y,x)
                var1 = pd.DataFrame(f.get(key).value.T,columns=y,index=x)
                var1.name = key
            else:
                var2 = pd.DataFrame(f.get(key).value.T,columns=y,index=x)
                var2.name = key
    return x,y

def write_geotiff(raster, gt, wkt, outputpath, dtype=gdal.GDT_Float32, options=['COMPRESS=LZW'], color_table=0, nbands=1, nodata=False):
    '''
    function from Jon Schwenk to convert to a geotif
    '''
    width = np.shape(raster)[1]
    height = np.shape(raster)[0]
    print(width,height)


    # Prepare destination file
    driver = gdal.GetDriverByName("GTiff")
    if options != 0:
        dest = driver.Create(outputpath, width, height, nbands, dtype, options)
    else:
        dest = driver.Create(outputpath, width, height, nbands, dtype)

    # Write output raster
    if color_table != 0:
        dest.GetRasterBand(1).SetColorTable(color_table)

    dest.GetRasterBand(1).WriteArray(raster)

    if nodata is not False:
        dest.GetRasterBand(1).SetNoDataValue(nodata)

    # Set transform and projection
    dest.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    #srs.ImportFromWkt(wkt)
    srs.ImportFromEPSG(32603)
    dest.SetProjection(srs.ExportToWkt())

    # Close output raster dataset
    dest = None
def get_geotransform(raster):
    y = raster.index.values
    y_min = y.min()
    y_max = y.max()

    x = raster.columns.values
    x_min = x.min()
    x_max = x.max()

    pixelwidth = x[1] - x[0]
    pixelheight = y[1] - y[0]
    print(pixelwidth,pixelheight)

    gt = ([x_min,round(pixelwidth,1),0,y_max,0,round(pixelheight,1)])
    print(gt)
    return gt
def rsquare_snow(df):
    #plot a one to one line against values and get an R_squared
    r2 = r_square_indi(df['Corrected_Depth'],df['UAS_MEAN'],df.index,0.0)
    return r2

def clean_mag(df,x,y):
    '''Program meant to clean up dataframe'''
    #get rid of 0s and negatives in either x,y column
    df = df[(df[x] > 0) | (df[y] > 0)]
    #drop na values in x,y
    #df = df.dropna(axis=0,subset=[x,y])
    return df

def raster_correlation(raster1,raster2,bins=np.arange(0.0,115,5),xlabel='',ylabel='',title=''):
    '''
    Function meant to plot the pixel values for random sample against eachother and get a linear correlation
    '''
    r1 = gdal.Open(raster1)
    r2 =gdal.Open(raster2)

    #reading rasters as array
    r1_arr = r1.ReadAsArray()
    r2_arr = r2.ReadAsArray()

    valid_index = np.where((r1_arr >= 0) & (r2_arr>=0))

    #binning by the second array
    r2_vals = r2_arr[valid_index]
    r1_vals = r1_arr[valid_index]
    bin_ind = np.digitize(r2_arr[valid_index],bins)
    #print(bin_ind)
    for bin in range(0,len(bins)):
        print(bins[bin])
        ind = np.where(bin_ind == bin)
        #print(bins[bin],ind,r2_vals[ind])
        #print(np.mean(r1_vals[ind]),len(r2_vals[ind].tolist()))
    binned_means,bin_edges,_ = stats.binned_statistic(r2_arr[valid_index],r1_arr[valid_index],'mean',bins=bins)
    binned_count,bin_edges,_ = stats.binned_statistic(r2_arr[valid_index],r1_arr[valid_index],'count',bins=bins)
    binned_std,bin_edges,_ = stats.binned_statistic(r2_arr[valid_index],r1_arr[valid_index],'std',bins=bins)

    #print(binned_means,bin_edges)
    #setting color ramp for visualization
    clist = [(0, "wheat"), (0.125, "lightgreen"), (0.25, "lightgreen"), (0.5, "green"),
         (0.7, "forestgreen"), (0.75, "forestgreen"), (1, "darkgreen")]
    rvb = mcolors.LinearSegmentedColormap.from_list("", clist)

    N = len(binned_means.tolist())
    print(N)
    x = np.arange(N).astype(float)

    fig,ax = plt.subplots(figsize=(10,7))
    bars = ax.bar(bin_edges[0:len(binned_means)],binned_means,yerr=binned_std,align='edge',linewidth=0.5,capsize=5,width=(bins[1]- bins[0]),edgecolor='k',color=rvb(np.arange(0,len(binned_means))/N))
    ax.set_title(title,fontsize=16)
    ax.set_xlabel(xlabel,fontsize=16)
    ax.set_ylabel(ylabel,fontsize=16)
    ax.set_ylim(0,3)

    #create dataframes of x,y,label
    df = pd.DataFrame(list(zip(bin_edges[0:len(binned_means)].tolist(),binned_means.tolist(),binned_count.tolist(),binned_std.tolist())), columns =['center', 'height','counts','std'])
    df["counts"] = df["counts"].astype(int)
    df['center'] = df['center']#+(bins[1]-bins[0])
    print(df)

    label_point(df['center'],df['height']+(ax.get_ylim()[1]-ax.get_ylim()[0])/15.0,df['counts'],ax)
    #plt.scatter(r2_arr[valid_index],r1_arr[valid_index])
    plt.savefig('Z:/AKSeward/Data/GIS/Teller/2019_Snow/Correlations/density_vs_depth.png',dpi=500)
    plt.close()
    #plt.show()
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x'], point['y'], 'n = '+str(point['val'].astype(int)),rotation=90)

def binned_columns(x,y,bins=np.arange(0.0,115,5),xlabel='',ylabel='',title=''):
    '''
    Binning by the x variable and plotting a bar plot
    '''
    #binning
    binned_means,bin_edges,_ = stats.binned_statistic(x,y,'mean',bins=bins)
    binned_count,bin_edges,_ = stats.binned_statistic(x,y,'count',bins=bins)
    binned_std,bin_edges,_ = stats.binned_statistic(x,y,'std',bins=bins)

    N = len(binned_means.tolist())
    clist = [(0, "wheat"), (0.125, "lightgreen"), (0.25, "lightgreen"), (0.5, "green"),
         (0.7, "forestgreen"), (0.75, "forestgreen"), (1, "darkgreen")]
    rvb = mcolors.LinearSegmentedColormap.from_list("", clist)

    fig,ax = plt.subplots(figsize=(10,7))
    bars = ax.bar(bin_edges[0:len(binned_means)],binned_means,yerr=binned_std,align='edge',linewidth=0.5,capsize=5,width=(bins[1]- bins[0]),edgecolor='k',color=rvb(np.arange(0,len(binned_means))/N))
    ax.set_title(title,fontsize=16)
    ax.set_xlabel(xlabel,fontsize=16)
    ax.set_ylabel(ylabel,fontsize=16)
    ax.set_ylim(0,2.5)

    #create dataframes of x,y,label
    df = pd.DataFrame(list(zip(bin_edges[0:len(binned_means)].tolist(),binned_means.tolist(),binned_count.tolist(),binned_std.tolist())), columns =['center', 'height','counts','std'])
    df["counts"] = df["counts"].astype(int)
    df['center'] = df['center']#+(bins[1]-bins[0])
    print(df)

    label_point(df['center'],df['height']+(ax.get_ylim()[1]-ax.get_ylim()[0])/15.0,df['counts'],ax)
    #plt.scatter(r2_arr[valid_index],r1_arr[valid_index])
    plt.savefig('Z:/AKSeward/Data/GIS/Teller/2019_Snow/Correlations/shrubheight_vs_depth_in_situ.png',dpi=500)
    plt.close()
    #plt.show()
def main():
    #raster_subsample_correlation
    snow = "Z:/AKSeward/Data/GIS/Teller/2019_Snow/Correlations/imagesnoD_clipped.tif"
    shrub_height = "Z:/AKSeward/Data/GIS/Teller/2019_Snow/Correlations/shrub_height_clipped.tif"
    shrub_density ="Z:/AKSeward/Data/GIS/Teller/2019_Snow/Correlations/shrub_density_clipped.tif"
    shrubs = pd.read_csv('Z:/AKSeward/Data/GIS/Teller/2019_Snow/snow_shrub_height_merge.csv')
    #raster_correlation(snow,shrub_density,xlabel='Shrub Density (counts/pixel)',ylabel='Average Snow Depth (m)',title='Impact of Shrub Density on Snow Depth')
    print(shrubs['shrub_height'])
    #binned_columns(shrubs['shrub_height'],shrubs['Corrected_Depth']/100.0,bins=np.arange(0.0,5,0.25),xlabel='Shrub Height (m)',ylabel='In-Situ Snow Depth (m)',title='Impact of Shrub Height on Snow Depth')
    x= pd.read_csv('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/xp.csv')
    y= pd.read_csv('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/yp.csv')
    filepath = 'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/snow2019_assess.mat'
    x,y = mat_to_dataframe(filepath)
    x.to_csv('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/xp.csv')
    y.to_csv('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/yp.csv')
    #gt_imagesno = get_geotransform(imagesno)
    #print(gt_imagesno)
    wkt = 'PROJCS["WGS 84 / UTM zone 3N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32603"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
    #write_geotiff(imagesno.values,gt_imagesno,wkt,'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/imagesno_nans.tif')
    write_geotiff(imagesnoD.values,gt_imagesno,wkt,'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/imagesnoD_nans.tif')
    '''

    snow_path = 'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/2019_snow_depth_comp.csv'
    snow_buffer = 'Z:/AKSeward/Data/GIS/Teller/2019_Snow/valid_only/2019_snow_comp.csv'
    snow = pd.read_csv(snow_buffer)

    #grab only the values that do not have shrub heights
    ind = pd.merge(snow['Counter'],shrubs['Counter'],on='Counter')
    false_ind = snow['Counter'][~snow.Counter.isin(ind.values.flatten())]
    snow = snow[snow.Counter.isin(false_ind.values.flatten())]

    snow = clean_mag(snow,'Corrected_Depth','UAS_MEAN')
    snow_inrange = snow[(snow['Corrected_Depth']<=snow['UAS_MAX']) & (snow['Corrected_Depth']>=snow['UAS_MIN'])]
    snow_outrange = snow[(snow['Corrected_Depth']>snow['UAS_MAX']) | (snow['Corrected_Depth']<snow['UAS_MIN'])]
    print(snow['Corrected_Depth'],snow['UAS_MEAN'])
    r2 = rsquare_snow(snow)
    r2_1 = r_squared_1_1_line(snow['Corrected_Depth'],snow['UAS_MEAN'])
    print(r2_1)
    plt.scatter(snow_inrange['Corrected_Depth'],snow_inrange['UAS_MEAN'],c='b',alpha=0.8,label='In Range')
    plt.scatter(snow_outrange['Corrected_Depth'],snow_outrange['UAS_MEAN'],c='r',alpha=0.2,label='Out of Range')
    plt.title('2019 Snow Comparison')
    plt.xlabel('Snow Depth from UAS (3m buffer avg) (cm)')
    plt.ylabel('In-Situ Snow Depth (cm)')
    plt.text(0.15,0.85,r'R $^2$: '+str(round(r2['r_sq'],3)),horizontalalignment='left',transform=plt.gcf().transFigure)
    plt.text(0.15,0.82,'Slope: '+str(round(r2['slope'],3)),horizontalalignment='left',transform=plt.gcf().transFigure)
    plt.text(0.85,0.85,r'R $^2$ 1:1 : '+str(round(r2_1,3)),horizontalalignment='right',transform=plt.gcf().transFigure)
    xy = np.arange(0,250,1)
    plt.plot(xy,xy,'r--')
    plt.xlim(0,2.50)
    plt.ylim(0,2.50)
    plt.legend(loc=4)

    plt.savefig('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/Figures/2019_comp_3m_buffer_range_noshrub.png',dpi=500)
    plt.close()
    #plt.show()
    #ax.text(0.05,0.9,r'R $^2$ : '+str(round(r_sq,3)),horizontalalignment='left',transform=ax.transAxes,fontsize=font_size
    '''
if __name__ == "__main__":
    main()
