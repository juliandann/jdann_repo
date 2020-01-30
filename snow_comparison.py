import pandas as pd
import numpy as np
import matplotlib.pyplot
import h5py
import gdal
import osr
def mat_to_dataframe(mat_file):
    #converts mat_files to
    f = h5py.File(mat_file,'r')
    keys = [key for key in f.keys()]

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
                var1 = pd.DataFrame(f.get(key).value.T,columns=y,index=x).fillna(0.0)
                var1.name = key
            else:
                var2 = pd.DataFrame(f.get(key).value.T,columns=y,index=x).fillna(0.0)
                var2.name = key
    return var1,var2

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
def main():
    filepath = 'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/snow2019_assess.mat'
    imagesno,imagesnoD = mat_to_dataframe(filepath)

    gt_imagesno = get_geotransform(imagesno)
    print(gt_imagesno)
    wkt = 'PROJCS["WGS 84 / UTM zone 3N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32603"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
    write_geotiff(imagesno.values,gt_imagesno,wkt,'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/imagesno.tif')
    write_geotiff(imagesnoD.values,gt_imagesno,wkt,'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/imagesnoD.tif')
    '''
    f = h5py.File(filepath,'r')
    [key for key in f.keys()]

    x = pd.Series(f.get('xp').value.flatten())
    y = pd.Series(f.get('yp').value.flatten())
    image2019sno = pd.DataFrame(f.get('image2019sno').value,columns=y,index=x).stack()
    image2019snoD = pd.DataFrame(f.get('image2019snoD').value,columns=y,index=x).stack()
    image2019sno.name='image2019sno'
    image2019snoD.name='image2019snoD'
    image2019snoD.to_csv('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/snoD.csv')
    image2019sno.to_csv('Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/sno.csv')

    print(image2019sno)
    print(image2019snoD)
    df = pd.merge(image2019sno, image2019snoD,how='outer')
    print(df)
    '''
if __name__ == "__main__":
    main()
