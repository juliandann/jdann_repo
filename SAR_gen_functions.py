import numpy as np
import pandas as pd
#import geopy.distance
import xarray as xr
import time
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib
import matplotlib.dates as mdates
import pdb
from scipy.optimize import curve_fit
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.pyplot import cm
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn import linear_model
import statsmodels.api as sm
# Import function to automatically create polynomial features!
from sklearn.preprocessing import PolynomialFeatures
# Import Linear Regression and a regularized regression function
from sklearn.linear_model import LassoCV
# Finally, import function to make a machine learning pipeline
from sklearn.pipeline import make_pipeline

from os import listdir
from os.path import isfile, join
from pygam import LinearGAM
from sklearn.metrics import mean_squared_error


class Paths:
    def __init__(self, figures, data_save,data_load):
        self.figures = figures
        self.data_save = data_save
        self.data_load = data_load
def path_file_reader(path_class,file):
    df = pd.read_csv(path_class.data_load + file)
    return df

#work comp class
work_paths = Paths('Z:/JDann/Documents/Documents/Julian_Python/SAR_programs_20181003/Figures/','Z:/AKSeward/2017_SAR/ABoVE_Soil_Moisture_Products/JBD_Products/','Z:/AKSeward/2017_SAR/ABoVE_Soil_Moisture_Products/JBD_Products/')

#personal comp paths
personal_paths = Paths('/Users/juliandann/Documents/LANL/SAR_DATA_AND_Programs/jdann_repo/Figures/','/Users/juliandann/Documents/LANL/qgis/CSV/','/Users/juliandann/Documents/LANL/qgis/CSV/')
def powspace(start, stop, power, num):
    start = np.power(start, 1/float(power))
    stop = np.power(stop, 1/float(power))
    return np.power( np.linspace(start, stop, num=num), power)

def dist(lat1, long1, lat2, long2):
    """Planar geometric distance function

    Parameters
    ----------
    lat1 : float
        latitude
    long1 : float
        longitude
    lat2 : float
        latitude coordinate two.
    long2 : float
        Longitude coordinate two.

    Returns
    -------
    type
        Description of returned object.

    """

    return np.sqrt((lat1-lat2)**2.0+ (long1-long2)**2.0)

def find_SAR_pixel(df,lat1,long1, lat2, long2,export_columns):
    """function that grabs the closest dataframe measurement using dist() function

    Parameters
    ----------
    df : Pandas Dataframe
        Description of parameter `df`.
    lat1 : string
        column within df for latitude
    long1 : string
        column within df for longitude
    lat2 : Float
        Value for latitude
    long2 : float
        Value for longitude
    export_columns : string
        Names of columns to put into the original dataframe that correspond to the nearest value.

    Returns
    -------
    Dataframe
        Updated dataframe with columns from the nearest value.

    """
    distances = df.apply(
        lambda row: dist(lat, long, row[lat1], row[lat2]),
        axis=1)
    return df.loc[distances.idxmin(), export_columns]

def closest_point_index_orig(point, points):
    """ Find closest point from a list of points. """
    return cdist([point], points).argmin()

def closest_point_index(df1,df2,lat1,lon1,lat2,lon2,*export_columns,prefix='df2'):
    """Grab the closest geodesic point from a list using scipy

    Parameters
    ----------
    point : float
        Can be a list of point pairs or a single pair
    points : float
        array of point pairs

    Returns
    -------
    index value(s)
        index values from the points array

    """

    print('Making tuples')
    #make tuple of lat lon from each dataframe
    df1['point'] = [(x, y) for x,y in zip(df1[lat1], df1[lon1])]
    df2['point'] = [(x, y) for x,y in zip(df2[lat2], df2[lon2])]
    print('Calculating nearest point')
    #calculate the closest point index value between the two arrays coordinates
    index = [cdist([x], list(df2['point'])).argmin() for x in df1['point']]

    #adding new columns to df1
    for column in export_columns:
        df1[prefix+'_'+column] = df2[column].iloc[index].values

    return df1

def csv_combined_file_maker(files,combined_filepath='combined_file.csv'):
    """Function that combines a list of csv files into one csv file that can be imported easily

    Parameters
    ----------
    files : list of strings
        The filepaths of the files you wish to combine. They must all have the same column names
    combined_filepath : string
        Filepath for the combined file.

    Returns
    -------
    csv
        Creates a CSV in the directory python is running in unless combined filepath is specified

    """

    list = []
    for i in files:
        df = pd.read_csv(i,sep=',')
        list.append(df)

    df = pd.concat(list,ignore_index=True)
    df.to_csv(combined_filepath,sep=',')

def distancer(lat1,lon1,lat2,lon2):
    """Function meant to be used in df.apply() in order to compute distance in WGS84 coordinate system

    Parameters
    ----------
    lat1 : float
        latitude for pair one.
    lon1 : float
        Longitude for pair one.
    lat2 : float
        Latitude for pair two
    lon2 : float
        Longitude for pair two.

    Returns
    -------
    float
        Distance in meters

    """
    coords_1 = (lat1, lon1)
    coords_2 = (lat2, lon2)
    return geopy.distance.VincentyDistance(coords_1, coords_2).m

def distance_calc(df,lat1,lon1,lat2,lon2,new_col_name):
    """Function used in conjunction with distancer to find distance in meters for columns of pandas dataframes that are in WGS84.

    Parameters
    ----------
    df : Pandas Dataframe
        Contains pairs of lat/long
    lat1 : float
        latitude for pair one.
    lon1 : float
        Longitude for pair one.
    lat2 : float
        Latitude for pair two
    lon2 : float
        Longitude for pair two.
    new_col_name : string
        The column name that will be added to the original dataframe

    Returns
    -------
    Pandas Dataframe
        Dataframe with a new column with 'new_col_name' that has distance in meters between the two coordinates.




    """

    df[new_col_name] = df.apply(lambda x: distancer(x[lat1],x[lon1],x[lat2],x[lon2]), axis=1)
    return df

def average_SM_at_pixel(df,pixel_index,savename='pixel_hydro_vs_above.png'):

    #get rid of barrow
    df.Locale=df.Locale.str.replace(' ','')
    df = df[df.Locale != 'Barrow']
    print(df.Locale.unique())
    #groupby pixel index and measurement depth
    pixel_group = df.groupby([pixel_index,'VWC_Measurement_Depth','SAR_Plot','Locale'])

    #get average of VWC values
    avg_pix = pixel_group['VWC','above_0.06','above_0.12','above_0.2','lbc_0.06','lbc_0.12','lbc_0.2'].agg(np.mean)

    #getting error of VWC
    std_vwc = pixel_group['VWC'].agg(np.std)

    avg_pix['VWC_std'] = std_vwc
    avg_pix[['VWC','VWC_std']] = avg_pix[['VWC','VWC_std']]/100.0

    avg_pix.reset_index(level=['VWC_Measurement_Depth','SAR_Plot','Locale'], inplace=True)

    #print(avg_pix)
    #plotting 6cm vwc vs. above 6cm closest pixel
    df_6cm = avg_pix[avg_pix['VWC_Measurement_Depth']==6]
    df_12cm = avg_pix[avg_pix['VWC_Measurement_Depth']==12]
    df_20cm = avg_pix[avg_pix['VWC_Measurement_Depth']==20]

    dfs = [df_6cm,df_12cm,df_20cm]
    y_col = ['above_0.06','above_0.12','above_0.2']

    #plot settings
    symbol =['^','s','o']
    title = ' In-Situ Soil Moisture vs. derived ABoVE SAR Soil Moisture Per Pixel'
    alpha = [1,0.6,0.4,0.2]
    fill_style = ['top','full','bottom',]
    edge_colors =['yellow','purple','orange']
    size = [10,10,10]
    save_name = savename
    fig,ax=plt.subplots(figsize=(15,10))
    col_name = 'VWC'
    y_col =['above_0.06','above_0.12','above_0.2']
    depth = ['0.06cm','0.12cm','0.2cm']
    set = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999"]

    colors = {'Teller':set[2],'Kougarok':set[1],'Council':set[0],'Barrow':set[3]}

    #fit a linear regression to x,y
    j=0

    for df in dfs:
        x= np.array(df['VWC'].values).reshape((-1, 1))
        y = df[y_col[j]].values
        model = LinearRegression().fit(x,y)
        r_sq = model.score(x,y)
        print(depth[j])
        print('r_sq for linear regression:',r_sq)
        print('intercept:', model.intercept_)
        print('slope:', model.coef_)
        j=j+1

    r2_6cm = r2_score(df_6cm['above_0.06'],df_6cm['VWC'])
    r2_12cm = r2_score(df_12cm['above_0.12'],df_12cm['VWC'])
    r2_20cm = r2_score(df_20cm['above_0.2'],df_20cm['VWC'])
    print('R-Squared values for 1:1 line: ')
    print('r^2 6cm: '+str(r2_6cm))
    print('r^2 12cm: '+str(r2_12cm))
    print('r^2 20cm: '+str(r2_20cm))
    i=0

    #calculating an overarching r squared and line of best fit

    for df1 in dfs:
        #print(list(df1))
        grouped = df1.groupby(pixel_index)
        for key,group in grouped:
            #print('Key: '+str(key))

            #print('Group: '+str(group))
            key = (group['Locale'].iloc[0]).strip()
            depth=group['VWC_Measurement_Depth'].iloc[0]
            #print('Key: '+str(key))
            #print(group[['Locale','SAR_Plot']])
            #print('Here!'+str(col_name)+str(y_col))


            #average groupby
            x_data = group[col_name]
            y_data = group[y_col[i]]
            xerr = group['VWC_std']


            #group.plot(ax=ax,x=x_data,y=y_data,xerr=xerr,yerr=yerr,kind='scatter',s=size[i],color=colors[key],edgecolors='k',linewidths=1.5 ,marker=MarkerStyle(symbol[i]))#,fillstyle=fill_style[i]),edgecolors=edge_colors[i])#,edgecolors=edge_colors[i],linewidths=2)
            ax.errorbar(x_data,y_data,xerr=xerr,markersize=size[i],color=colors[key],fmt=symbol[i])
            plt.rcParams['lines.linewidth'] = 0.5
        i=i+1
    plt.xlim(0,1)
    plt.ylim(0,1)

    #creating custom legends
    legend_elements = [ Line2D([0], [0], marker='^', color='w', label='6cm',
                          markerfacecolor='grey', markersize=10,markeredgecolor='k',markeredgewidth=2),
                      Line2D([0], [0], marker='s', color='w', label='12cm',
                            markerfacecolor='grey', markersize=15,markeredgecolor='k',markeredgewidth=2),
                        Line2D([0], [0], marker='o', color='w', label='20cm',
                            markerfacecolor='grey', markersize=25,markeredgecolor='k',markeredgewidth=2)]

    legend_elements2 = [ Line2D([0], [0], marker='o', color='w', label='Teller',
                          markerfacecolor=set[2], markersize=15,markeredgecolor='k'),
                      Line2D([0], [0], marker='o', color='w', label='Kougarok',
                            markerfacecolor=set[1], markersize=15,markeredgecolor='k'),
                        Line2D([0], [0], marker='o', color='w', label='Council',
                            markerfacecolor=set[0], markersize=15,markeredgecolor='k')]

    leg = plt.legend(handles=legend_elements, loc='upper right',prop={'size':20})
    ax.add_artist(leg)

    ax.legend(handles = legend_elements2, loc='upper left',prop={'size':16})
    #ax.add_artist(leg2)

    #create 1:1 line
    y =np.linspace(0,1,100)
    x = np.linspace(0,1,100)
    plt.plot(x,y,ls='--',c='k')
    plt.xlabel('Volumetric Water Content (%/100)',fontsize=18)
    plt.ylabel('ABoVE SAR P-Band Flight (%/100)',fontsize=18)
    plt.title(title,fontsize=22)
    plt.savefig(save_name,dpi=500)
    #plt.show()
    plt.close()




    '''ax= df_6cm.plot(x='VWC',y='above_0.06',c='b',label='6 cm',kind='scatter')
    df_12cm.plot( x='VWC',y='above_0.12',color='g',ax = ax,label='12 cm',kind='scatter')
    df_20cm.plot( x='VWC',y='above_0.2',color='r',ax = ax,label='20 cm',kind='scatter')

    x = np.arange(-0.2,1.2,0.1)
    plt.plot(x,x,'k--')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()'''

def adjusting_SAR_Plot_names(df,column,prefix,replacement):
    """Short summary.

    Parameters
    ----------
    df : Pandas Dataframe
        Data
    column : string
        Column within dataframe
    prefix : string
        single element or list of strings that are the start of a value in the column you wish to change
    replacement : string
        Replacement value for the prefix

    Returns
    -------
    Pandas Dataframe
        Dataframe with changed column.

    """
    i = 0
    for pre in prefix:
        df[column][df[column].str.startswith(pre)] = replacement[i]
        i=i+1
    return df

def engstrom_SM(df,dielectric):
    """SAR soil moisture calibration from Engstrom et al. (2005)

    Parameters
    ----------
    df : dataframe
        Dataframe of ABoVE
    dielectric : string
        Column name for dielectric constant

    Returns
    -------
    dataframe
        Dataframe with extra column

    """
    df['test_VWC'] = -2.50 + (2.508*df[dielectric]) - (0.03634*(df[dielectric]**2.0)) + (0.0002394*(df[dielectric]**3.0))
    df.to_csv('Z:/AKSeward/2017_SAR/ABoVE_Soil_Moisture_Products/JBD_Products/discrete_all_sites_engstromtest.csv')

def td_comparison(df,pixel_index,var1,var2,var2_error,save_folder):
    """Comparison of in-situ thaw depth with ABoVE data using pixel binning.

    Parameters
    ----------
    df : Pandas dataframe
        Dataframe from csv from closest_point_index()
    pixel_index : string
        Column name of thing to use to bin same values
    save_folder : string
        folder to save plots into
    var1 : string
        Column name of first variable to compare with var2
    var2 : string
        Column name of second variable to compare with var1
    var2_error : string
        Column name of error associated with var2

    Returns
    -------
        Plots can be found in save_folder

    """
    '''
    #get rid of -9999 values
    df = df.drop(df[df.Thaw_Depth == -9999].index)

    pixel_group = df.groupby([pixel_index])

    avg_td = pixel_group[var1,var2,var2_error].agg(['mean','count'])
    avg_td.to_csv(save_folder+'avg_td_per_pixel.csv')
    '''

    fig,ax=plt.subplots(figsize=(15,10))

    legend_elements = [ Line2D([0], [0], marker='o', color='red', label='0 - 5',
                          markerfacecolor='red', markersize=10),
                          Line2D([0], [0], marker='o', color='orange', label='6 - 10', markerfacecolor='orange', markersize=10),
                          Line2D([0], [0], marker='o', color='green', label='11 - 15',markerfacecolor='green', markersize=10),
                          Line2D([0], [0], marker='o', color='blue', label='16 -20',
                                                markerfacecolor='blue', markersize=10),]
    #plt.errorbar(df[var1],df[var2]*100.0,yerr=df[var2_error]*100.0,fmt='o',mfc=df['Count'])
    #df.plot(var1,var2,yerr=var2_error,c='Count',kind='scatter',colormap='viridis',ax = ax)
    colors = df.groupby('Color')
    r_score = r2_score(df[var1],df[var2])
    print(r_score)
    for name, group in colors:
        print(name,group)
        #group.plot(var1,var2,yerr=var2_error,c=str(name), kind='scatter')
        plt.errorbar(group[var1],group[var2],yerr=group[var2_error],ecolor=str(name),c=str(name),fmt='o')

    ax.legend(handles=legend_elements,title= 'Counts',fontsize=18)
    plt.text(78,25,'R-Squared ='+str(round(r_score,3)),horizontalalignment='right',fontsize=20)
    x = np.linspace(0,100,100)
    y = x

    plt.plot(x,y,'--k')
    plt.ylim(20,80)
    plt.xlim(20,80)
    plt.xlabel('In-Situ Thaw Depth (cm)',fontsize=18)
    plt.ylabel('ABoVE Pixel ALT (cm)',fontsize=18)
    plt.title('ABoVE ALT vs. In-Situ ALT',fontsize=24)
    plt.savefig(save_folder+'ALT_overall_comparison.png')
    plt.close()

def r_squared_1_1_line(x,y):
    """Compares x,y values with one to one line

    Parameters
    ----------
    x : Pandas column
        x values (in-situ hydrosense)
    y : Pandas Column
        y values (ABoVE values)

    Returns
    -------
    float
        R squared value compared to 1:1 line

    """
    #set up 1:1 line
    y_1 = x
    avg_y = np.mean(y)
    residual_1 = np.sum((y-y_1)**2.0)
    residual_null = np.sum((y-avg_y)**2.0)
    r2 = (residual_null - residual_1)/residual_null
    print('residual 1:1: '+str(residual_1))
    print('residual null: '+str(residual_null))
    return r2

def landcover_boxplots(df,class_cat,classes,mean,std,count):
    """Function written to analyze landcover type by h_value

    Parameters
    ----------
    df : dataframe
        Main data
    class_cat : string
        Overarching category
    classes : string
        Subclasses
    mean : string
        average value
    std : string
        standard deviation
    count : string
        counts

    Returns
    -------
    type
        Description of returned object.

    """
    print(list(df))
    grouped = df.groupby(classes)
    type = []
    fweighted_mean = []
    fweighted_std = []
    i=0
    ax = plt.figure(figsize=(24,12))
    for key,group in grouped:

        type.append(key)
        weighted_mean = []
        weighted_std = []
        j=0
        for row_index,row in group.iterrows():
            if j == 0:
                weighted_mean.append( (row[mean])*(row[count] / group[count].sum()))
                weighted_std.append( (row[std])*(row[count] / group[count].sum()))
            else:
                weighted_mean.append( weighted_mean[-1]+ (row[mean])*(row[count] / group[count].sum()))
                weighted_std.append(weighted_std[-1]+ (row[std])*(row[count] / group[count].sum()))
            j=j+1
        fweighted_mean.append(weighted_mean[-1])
        fweighted_std.append(weighted_std[-1])


        plt.errorbar(x=group[class_cat],y=group[mean],yerr=group[std],capsize=3,ls='',marker='o',label=type[-1])

        for row_index,row in group.iterrows():
            plt.text(i,0.5,str(round(row[count],2)),ha='center',rotation='vertical')
            i=i+1

        #print(group[['Ecosystem_LU','above_h_count','ABoVE_h_mean']])
    zippedlist = list(zip(type,fweighted_mean,fweighted_std))
    y_data = 'lbc_6_mean'
    y_err = 'lbc_6_std'
    df2 = pd.DataFrame(zippedlist, columns = ['Class' , y_data,y_err])

    #plt.tight_layout()


    bottom = df2[y_data]-df2[y_err]
    top = df2[y_data]+df2[y_err]
    df2['height'] = top-bottom
    df2['colors'] = ['blue','orange','green','red','purple','brown']
    #plt.bar(x=df2['Class'],height=df2['height'],bottom=bottom,linewidth=2.0,ecolor='black',edgecolor='black')
    #ax= df2.plot(x='Class',y='height',kind='bar',bottom=bottom,legend=False,ecolor='black',edgecolor='black')
    #df2.plot(x='Class',y=y_data,yerr=y_err,color=df2['colors'],kind='bar',legend=False)
    #plt.errorbar(x=df2['Class'],y=df2['h_mean'],yerr=df2['h_std'],capsize=3,ls='',color=df2['colors'])

    plt.legend(fontsize=15,loc=1)
    plt.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.48)
    plt.xticks(ha='center',rotation='vertical')
    plt.ylim(0.1,0.75)
    plt.xlabel('Vegetation Classification (GAP/LANDFIRE)',fontsize=18)
    plt.ylabel('ABoVE LBC 6cm (%/100)',fontsize=18)
    plt.title('LBC 6cm vs. Landcover Type',fontsize=24)
    #plt.show()
    plt.savefig('Z:/JDann/Documents/Documents/Julian_Python/SAR_programs_20181003/Figures/GapLandfire/lbc_6cm_discrete.png',dpi=500)
    plt.close()
def laura_bc_SM(df,e1,e2,z1,depth):
    #tau values from Hydrosense through R. Chen
    tau_12cm = (4.312 + np.sqrt(df[e1]))/5.354
    tau_20cm = (3.43 + np.sqrt(df[e2]))/3.158

    #Laura bc calibration equations
    mv1 = ((-24.28 * (tau_12cm**2.0)) + 134.55*tau_12cm - 110.245)/100.0
    mv2 = ((7.693 * (tau_20cm**2.0)) + 1.641*tau_20cm - 12.341)/100.0

    depth_tag = 'lbc_'+str(depth)

    #conditions
    df[depth_tag] = np.where(df[z1] > depth, mv1,(mv1 * df[z1] + mv2 *(depth - df[z1]))/depth )
    return df

def laura_bc_indi(e1,e2,z1,depth):
    tau_12cm = (4.312 + np.sqrt(e1))/5.354
    tau_20cm = (3.43 + np.sqrt(e2))/3.158

    #Laura bc calibration equations
    mv1 = ((-24.28 * (tau_12cm**2.0)) + 134.55*tau_12cm - 110.245)/100.0
    mv2 = ((7.693 * (tau_20cm**2.0)) + 1.641*tau_20cm - 12.341)/100.0

    if z1 > depth:
        mv=mv1
    else:
        mv = ((mv1 * z1) + mv2* (depth - z1))/depth

    return mv

def engstrom_SM_mv(df,mv1,mv2,z1,depth):
    depth_tag = 'engstrom_'+str(depth)
    df[depth_tag] = np.where(df[z1] > depth, df[mv1],(df[mv1] * df[z1] + df[mv2] *(depth - df[z1]))/depth )
    return df

def compare_above_vs_hydrosense(df_hydro,df_above,path,xlabel='',ylabel='',title='',save_name='test.png'):
    '''This function is meant to take the ABoVE SAR values with discrete measurements at 6,12,20 cm and compare with the hydrosense data'''

    #use ggplot
    matplotlib.style.use('ggplot')

    #make date into datetime object
    df_hydro['Date'] = pd.to_datetime(df_hydro['TimeStamp'],errors='coerce')

    #take only late dates
    VWC_df_late = df_hydro[(df_hydro['Date']>pd.Timestamp(2017,7,1))]

    symbol =['^','s','o']
    alpha = [1,0.6,0.4,0.2]

    #divide each column by 100.0
    VWC_df_late['VWC'] = VWC_df_late['VWC']/100.0

    size = [10,10,10]

    six = 'lbc_0.06_aug'
    twelve = 'lbc_0.12_aug'
    twenty = 'lbc_0.2_aug'
    '''
    six = 'engstrom_0.06'
    twelve = 'engstrom_0.12'
    twenty = 'engstrom_0.2'
    '''
    depths =[six,twelve,twenty]


    #grab SAR plot averages for each depth
    group_above = df_above.groupby('SAR_Plot')
    site_avg = []
    site_std = []
    site_name = []
    for key_ab,group_ab in group_above:
        if key_ab != np.nan:
            site_name.append(key_ab)
            avg_group = []
            std_group = []
            for depth in depths:
                avg_group.append(group_ab[depth].mean())
                std_group.append(group_ab[depth].std())
            site_avg.append(avg_group)
            site_std.append(std_group)
    site_avg = np.array(site_avg)
    site_std = np.array(site_std)

    df_ab = pd.DataFrame({'SAR_Plot':site_name,'six':site_avg[:,0],'twelve':site_avg[:,1],'twenty':site_avg[:,2],'six_err':site_std[:,0],'twelve_err':site_std[:,1],'twenty_err':site_std[:,2]})

    fig,ax=plt.subplots(figsize=(10,10))
    #colors = {'Teller':'royalblue','Kougarok':'salmon','Council':'seagreen'}
    set = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999"]

    colors = {'Teller':set[2],'Kougarok':set[1],'Council':set[0]}

    grouped =  VWC_df_late.groupby(['Locale','SAR_Plot'])
    VWC_depth = [6,12,20]
    x_array = []
    y_array = []
    key_array = []
    sar_array =[]
    depth_array = []
    xerr_array =[]
    model_std = []
    for i in range(0,len(VWC_depth)):
        for key,group in grouped:
            key2 = key
            key = group['Locale'].iloc[0]
            if key !='Barrow':
                #print(key)
                sarplot=group['SAR_Plot'].unique()[0]

                group_depth = group[group['VWC_Measurement_Depth'] == VWC_depth[i]]

                #average groupby
                x_data = group_depth['VWC'].mean()
                xerr = group_depth['VWC'].std()
                y = df_ab[df_ab['SAR_Plot'] == sarplot]
                if VWC_depth[i] == 6:
                    y_data = y['six']
                    yerr = y['six_err']
                if VWC_depth[i] == 12:
                    y_data = y['twelve']
                    yerr = y['twelve_err']
                if VWC_depth[i] == 20:
                    y_data = y['twenty']
                    yerr = y['twenty_err']

                x_array.append(x_data)
                y_array.append(float(y_data))
                key_array.append(key)
                sar_array.append(key2[1])
                depth_array.append(VWC_depth[i])
                xerr_array.append(xerr)
                model_std.append(float(yerr))
                ax.errorbar(x_data,y_data,xerr=xerr,yerr=yerr,markersize=size[i],color=colors[key],fmt=symbol[i])
                plt.rcParams['lines.linewidth'] = 0.5
    #print(x_array,'\n',y_array)
    r2 = r2_score(x_array,y_array)
    r2_other = r_square_indi(x_array,y_array,np.arange(0,len(x_array)),0.0)
    y_line = r2_other['slope']*np.array(x_array) + r2_other['intercept']
    print('Chi-Square',sum(((np.array(x_array)-np.array(y_line))**2.0)/(np.array(xerr_array)**2.0))/(len(x_array)-2.0))
    print('RMSE: ',np.sqrt(mean_squared_error(np.array(x_array),np.array(y_array))))
    plt.text(0.5,0.95,r'R$^2$ 1:1: '+str(round(r2,3))+ r'   R$^2$: '+str(round(r2_other['r_sq'],3)),horizontalalignment='center',fontsize=18)
    print(len(key_array),len(x_array),len(y_array))
    df2 = pd.DataFrame(list(zip(key_array,x_array,y_array,sar_array,depth_array,xerr_array,model_std)),columns=['Site','In_Situ','Model','SAR_Plot','Depth','In_Situ_std','Model_std'])
    df2.to_csv(path.data_load+'avg_plots.csv')
    plt.xlim(0,1)
    plt.ylim(0,1)

    #creating custom legends
    legend_elements = [ Line2D([0], [0], marker='^', color='w', label='6cm',
                          markerfacecolor='grey', markersize=10,markeredgecolor='k',markeredgewidth=2),
                      Line2D([0], [0], marker='s', color='w', label='12cm',
                            markerfacecolor='grey', markersize=15,markeredgecolor='k',markeredgewidth=2),
                        Line2D([0], [0], marker='o', color='w', label='20cm',
                            markerfacecolor='grey', markersize=25,markeredgecolor='k',markeredgewidth=2)]

    legend_elements2 = [ Line2D([0], [0], marker='o', color='w', label='Teller',
                          markerfacecolor=set[2], markersize=15,markeredgecolor='k'),
                      Line2D([0], [0], marker='o', color='w', label='Kougarok',
                            markerfacecolor=set[1], markersize=15,markeredgecolor='k'),
                        Line2D([0], [0], marker='o', color='w', label='Council',
                            markerfacecolor=set[0], markersize=15,markeredgecolor='k')]

    leg = plt.legend(handles=legend_elements, loc='upper right',prop={'size':20})
    ax.add_artist(leg)

    ax.legend(handles = legend_elements2, loc='upper left',prop={'size':16})
    #ax.add_artist(leg2)

    #create 1:1 line
    y =np.linspace(0,1,100)
    x = np.linspace(0,1,100)
    #plt.plot(x,y,ls='--',c='k')
    plt.xlabel(xlabel,fontsize=18)
    plt.ylabel(ylabel,fontsize=18)
    plt.title(title,fontsize=22)
    #plt.savefig(path.figures+save_name,dpi=500)
    print(path.figures+save_name)
    plt.show()
    #plt.close()

def weather(met_data,title='Temperature and Precipitation at Kougarok (mile 64)'):

    #make date into datetime object
    met_data['Date'] = pd.to_datetime(met_data['Date'],errors='coerce')
    met_data['Rain'][met_data['Rain'] == 6999] = np.nan
    met_data['Air Temp @ 1.5m'][met_data['Air Temp @ 1.5m'] == 6999] = np.nan
    fig, ax1 = plt.subplots(figsize=(15,10))
    ax1.plot(met_data['Date'],met_data['Air Temp @ 1.5m'])
    ax1.set_xlim(pd.Timestamp('2017-07-01'),pd.Timestamp('2017-09-01'))
    ax1.set_ylim(-20,30)
    ax1.set_ylabel('Air Temperature (Celsius)',color='b',fontsize=18)
    ax1.set_xlabel('Date',fontsize=18)
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=3))

    #initiate secondary y axis
    ax2 = ax1.twinx()
    ax2.plot(met_data['Date'],met_data['Rain'],color='r')
    ax2.set_ylabel('Rain (cm)',color='r',fontsize=16)
    ax2.tick_params(axis='y', labelcolor='r')

    above= plt.axvline(x=pd.Timestamp('2017-08-17'),color='k',linestyle='dashed',label='P-Band ABoVE SAR Flight')
    #plt.axvline(x=pd.Timestamp('2017-08-16'),color='k',linestyle='dashed',label='In-Situ Measurents')
    insitu = ax2.axvspan(pd.Timestamp('2017-08-17 00:00:00'),pd.Timestamp('2017-08-17 23:59:59'),label='In-Situ Measurents',color='crimson',alpha=0.3)
    plt.legend(handles=[above,insitu],fontsize=18,loc=1)
    plt.title(title,fontsize=24)

    plt.setp(ax1.get_xticklabels(), rotation=30)
    #plt.gcf().autofmt_xdate()
    #plt.show()
    plt.savefig('Z:/JDann/Documents/Documents/Julian_Python/SAR_programs_20181003/Figures/Weather/Kougarok_weather.png',dpi=500)
    plt.close()

def linear_regression_plots(df,pixel_index,savename='test.png'):

    #get rid of barrow
    df.Locale=df.Locale.str.replace(' ','')
    df = df[df.Locale != 'Barrow']

    #groupby pixel index and measurement depth
    pixel_group = df.groupby([pixel_index,'VWC_Measurement_Depth','Locale'])

    #get average of VWC values
    avg_pix = pixel_group['VWC','above_0.06','above_0.12','above_0.2','lbc_0.06','lbc_0.12','lbc_0.2'].agg(np.mean)


    #getting error of VWC
    std_vwc = pixel_group['VWC'].agg(np.std)

    avg_pix['VWC_std'] = std_vwc
    avg_pix.reset_index(level=['VWC_Measurement_Depth','Locale'], inplace=True)
    avg_pix[['VWC','VWC_std']] = avg_pix[['VWC','VWC_std']]/100.0

    grouped = avg_pix.groupby(['Locale','VWC_Measurement_Depth'])

    #set up figure
    fig,ax = plt.subplots(3,3,sharex='col',sharey='row',figsize=(15,12))
    ax = ax.flatten()

    #plot settings
    symbol ={'6':'^','12':'s','20':'o'}
    set = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999"]

    colors = {'Teller':set[2],'Kougarok':set[1],'Council':set[0],'Barrow':set[3]}


    j=0

    for key,group in grouped:
        site = key[0]
        depth = str(key[1])
        if depth == '6':
            y = 'lbc_0.06'
        if depth == '12':
            y = 'lbc_0.12'
        if depth == '20':
            y = 'lbc_0.2'


        #r-squared analysis
        x= np.array(group['VWC'].values).reshape((-1, 1))
        y = group[y].values
        model = LinearRegression().fit(x,y)
        r_sq = model.score(x,y)
        print('r_sq for linear regression:',r_sq)
        print('intercept:', model.intercept_)
        print('slope:', model.coef_)


        ax[j].errorbar(x,y,fmt=symbol[depth],color=colors[site])
        x_2 = np.linspace(0,1,20)
        ax[j].plot(x_2,model.coef_ *x_2 +model.intercept_,ls='--',c='k',lw=1.5)
        ax[j].set_xlim(0,1)
        ax[j].set_ylim(0,1)
        ax[j].set_title(site+' '+str(depth)+'cm')
        ax[j].text(0.05,0.9,r'r$^{2}$ ='+str(round(r_sq,2)))
        plt.rcParams['lines.linewidth'] = 0.5

        j=j+1

    plt.text(0.5,0.05,'In-Situ VWC (%/100)',horizontalalignment='center',transform=plt.gcf().transFigure,fontsize=18)
    plt.text(0.05,0.5,'P-Band SAR Laura BC Method (%/100)',verticalalignment='center',rotation='vertical',transform=plt.gcf().transFigure,fontsize=18)
    plt.text(0.5,0.95,r'R$^{2}$'+' Analyses',horizontalalignment='center',transform=plt.gcf().transFigure,fontsize=24)
    #creating custom legends
    legend_elements = [ Line2D([0], [0], marker='^', color='w', label='6cm',
                          markerfacecolor='grey', markersize=10,markeredgecolor='k',markeredgewidth=2),
                      Line2D([0], [0], marker='s', color='w', label='12cm',
                            markerfacecolor='grey', markersize=15,markeredgecolor='k',markeredgewidth=2),
                        Line2D([0], [0], marker='o', color='w', label='20cm',
                            markerfacecolor='grey', markersize=25,markeredgecolor='k',markeredgewidth=2)]

    legend_elements2 = [ Line2D([0], [0], marker='o', color='w', label='Teller',
                          markerfacecolor=set[2], markersize=15,markeredgecolor='k'),
                      Line2D([0], [0], marker='o', color='w', label='Kougarok',
                            markerfacecolor=set[1], markersize=15,markeredgecolor='k'),
                        Line2D([0], [0], marker='o', color='w', label='Council',
                            markerfacecolor=set[0], markersize=15,markeredgecolor='k')]

    leg = fig.legend(handles=legend_elements, loc='upper right',prop={'size':16})
    fig.add_artist(leg)

    fig.legend(handles = legend_elements2, loc='upper left',prop={'size':16})

    plt.savefig(savename,dpi=500)
    plt.close()
    plt.show()

def linear_regression_main(df,pixel_index,savename='test.png'):
    """Plot of all above VWC pixels compared wiith their vwc counterpoints and taking a linear regression and an R^2 value.

    Parameters
    ----------
    df : Pandas Dataframe
        closest ABoVE pixels to VWC
    pixel_index : string column name
        Description of parameter `pixel_index`.
    savename : string
        Save name for plot.

    Returns
    -------
    nothing

    """

    #get rid of barrow
    df.Locale=df.Locale.str.replace(' ','')
    df = df[(df.Locale != 'Barrow') ]

    print(df.Locale.unique())
    #groupby pixel index and measurement depth
    pixel_group = df.groupby([pixel_index,'VWC_Measurement_Depth','Locale'])

    #get average of VWC values
    avg_pix = pixel_group['VWC','above_0.06','above_0.12','above_0.2','lbc_0.06','lbc_0.12','lbc_0.2'].agg(np.mean)

    #getting error of VWC
    std_vwc = pixel_group['VWC'].agg(np.std)

    avg_pix['VWC_std'] = std_vwc
    avg_pix.reset_index(level=['VWC_Measurement_Depth','Locale'], inplace=True)
    avg_pix[['VWC','VWC_std']] = avg_pix[['VWC','VWC_std']]/100.0

    grouped = avg_pix.groupby(['VWC_Measurement_Depth'])

    #set up figure
    fig,ax = plt.subplots(3,1,figsize=(7,15),sharex=True)
    ax= ax.flatten()
    #plot settings
    symbol ={'6':'^','12':'s','20':'o'}
    set = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999"]

    colors = {'Teller':set[2],'Kougarok':set[1],'Council':set[0],'Barrow':set[3]}


    rsq = []
    intercepts = []
    coeff = []

    j=0
    for key,group in grouped:
        depth = str(key)
        if depth == '6':
            y = 'lbc_0.06'
        if depth == '12':
            y = 'lbc_0.12'
        if depth == '20':
            y = 'lbc_0.2'
        #r-squared analysis
        x= np.array(group['VWC'].values).reshape((-1, 1))
        y = group[y].values
        model = LinearRegression().fit(x,y)
        r_sq = model.score(x,y)
        print('r_sq for linear regression:',r_sq)
        print('intercept:', model.intercept_)
        print('slope:', model.coef_)
        rsq.append(r_sq)
        intercepts.append(model.intercept_)
        coeff.append(model.coef_)
        x_2 = np.linspace(0,1,20)
        ax[j].plot(x_2,model.coef_ *x_2 +model.intercept_,ls='--',c='k',lw=1.5)
        j=j+1

    grouped = avg_pix.groupby(['Locale','VWC_Measurement_Depth'])

    j=0
    for key,group in grouped:
        site = key[0]
        depth = str(key[1])
        if depth == '6':
            y = 'lbc_0.06'
            j=0
        if depth == '12':
            y = 'lbc_0.12'
            j=1
        if depth == '20':
            y = 'lbc_0.2'
            j=2

        ax[j].errorbar(group['VWC'],group[y].values,fmt=symbol[depth],color=colors[site])
        ax[j].set_xlim(0,1)
        ax[j].set_ylim(0,1)
        ax[j].set_title(str(depth)+'cm' )
        ax[j].text(0.05,0.9,r'r$^{2}$ ='+str(round(rsq[j],2)))
        plt.rcParams['lines.linewidth'] = 0.5



    plt.text(0.5,0.05,'In-Situ VWC (%/100)',horizontalalignment='center',transform=plt.gcf().transFigure,fontsize=18)
    plt.text(0.05,0.5,'P-Band SAR Laura BC Method (%/100)',verticalalignment='center',rotation='vertical',transform=plt.gcf().transFigure,fontsize=18)
    plt.text(0.5,0.95,r'R$^{2}$'+' Analyses',horizontalalignment='center',transform=plt.gcf().transFigure,fontsize=24)
    #creating custom legends
    legend_elements = [ Line2D([0], [0], marker='^', color='w', label='6cm',
                          markerfacecolor='grey', markersize=10,markeredgecolor='k',markeredgewidth=2),
                      Line2D([0], [0], marker='s', color='w', label='12cm',
                            markerfacecolor='grey', markersize=15,markeredgecolor='k',markeredgewidth=2),
                        Line2D([0], [0], marker='o', color='w', label='20cm',
                            markerfacecolor='grey', markersize=25,markeredgecolor='k',markeredgewidth=2)]

    legend_elements2 = [ Line2D([0], [0], marker='o', color='w', label='Kougarok',
                          markerfacecolor=set[1], markersize=15,markeredgecolor='k'),
                      Line2D([0], [0], marker='o', color='w', label='Kougarok',
                            markerfacecolor=set[1], markersize=15,markeredgecolor='k'),
                        Line2D([0], [0], marker='o', color='w', label='Council',
                            markerfacecolor=set[0], markersize=15,markeredgecolor='k')]

    leg = fig.legend(handles=legend_elements, loc='upper right',prop={'size':16})
    fig.add_artist(leg)

    fig.legend(handles = legend_elements2, loc='upper left',prop={'size':16})

    plt.savefig(savename,dpi=500)
    plt.close()
    #plt.show()
def sunfactor_func(df,aspect,sunfactor):
    df[sunfactor] = 0.5*(np.cos((df[aspect]-40.0)/360.*2.*np.pi)-np.sin((df[aspect]-40.)/360.*2.*np.pi))
    return df
def boxplots_macrotopology_vwc_plots(df,path):
    """Function that makes boxplots 3x1 for each lbc depth against aspect,slope,sunfactor, and curvature.

    Parameters
    ----------
    df : pandas dataframe
        Description of parameter `df`.

    Returns
    -------
    type
        Description of returned object.

    """

    #column names for easy changing
    aspect = 'aspect'
    slope = 'slope'
    sunfactor = 'sunfactor'
    curve = 'curvature'
    ndvi = 'NDVI_30m'
    fa = 'fa'
    pix_index = 'Index'
    six_cm = 'lbc_0.06_aug'
    twelve_cm = 'lbc_0.12_aug'
    twenty_cm = 'lbc_0.2_aug'
    h = 'h'
    h_err = 'h_uncertainty'
    alt = 'alt'
    alt_err = 'alt_uncertainty'
    vwc_depth = [six_cm,twelve_cm,twenty_cm]

    df = sunfactor_func(df,aspect,sunfactor)

    slope_bins = np.arange(0,95,5)
    slope_names= []
    list_slope = []

    aspect_names= []
    aspect_bins = np.arange(0,390,30)
    list_aspect = []

    sunfactor_names = []
    list_sunfactor = []
    sunfactor_bins = np.arange(-0.8,0.9,0.1)

    ndvi_bins = np.arange(-0.2,1,0.05)
    list_ndvi = []
    ndvi_names = []

    fa_bins =powspace(0,25000,4,10)
    list_fa = []
    fa_names = []

    curv_bins = np.arange(-1.0,1.1,0.1)
    list_curv =[]
    curv_names = []



    #separate pandas column into value 0-5 etc through slope
    for i in range(0,len(slope_bins)-1):
        list_slope.append( df[((df[slope] > float(slope_bins[i])) & (df[slope] <= float(slope_bins[i+1])))])
        slope_names.append(str(slope_bins[i])+' - '+str(slope_bins[i+1]))

    for i in range(0,len(aspect_bins)-1):
        list_aspect.append( df[((df[aspect] > float(aspect_bins[i])) & (df[aspect] <= float(aspect_bins[i+1])))])
        aspect_names.append(str(aspect_bins[i])+' - '+str(aspect_bins[i+1]))

    for i in range(0,len(sunfactor_bins)-1):
        list_sunfactor.append( df[((df[sunfactor] > float(sunfactor_bins[i])) & (df[sunfactor] <= float(sunfactor_bins[i+1])))])
        sunfactor_names.append(str(round(sunfactor_bins[i],2))+' - '+str(round(sunfactor_bins[i+1],2)))

    for i in range(0,len(curv_bins)-1):
        list_curv.append( df[((df[curve] > float(curv_bins[i])) & (df[curve] <= float(curv_bins[i+1])))])
        curv_names.append(str(round(curv_bins[i],2))+' - '+str(round(curv_bins[i+1],2)))

    for i in range(0,len(fa_bins)-1):
        list_fa.append( df[((df[fa] > float(fa_bins[i])) & (df[fa] <= float(fa_bins[i+1])))])
        fa_names.append(str(round(fa_bins[i],2))+' - '+str(round(fa_bins[i+1],2)))

    for i in range(0,len(ndvi_bins)-1):
        list_ndvi.append( df[((df[ndvi] > float(ndvi_bins[i])) & (df[ndvi] <= float(ndvi_bins[i+1])))])
        ndvi_names.append(str(round(ndvi_bins[i],2))+' - '+str(round(ndvi_bins[i+1],2)))


    #slope
    fig,ax = plt.subplots(1,3,figsize=(20,10))
    ax = ax.flatten()

    for j in range(0,len(vwc_depth)):
        for i in range(0,len(slope_bins)-1):
            ax[j].boxplot(list_slope[i][vwc_depth[j]].values,positions=[i],widths=0.6)
        ax[j].set_title('Slope vs. ' +str(vwc_depth[j]))
        ax[j].set_xticklabels(slope_names,rotation=45, horizontalalignment='right')
        ax[j].set_xlabel('Slope Range (deg)',fontsize=18)
        ax[j].set_ylim(0,1)
        ax[j].set_ylabel('VWC (%/100)',fontsize=18 )

    plt.savefig(path.figures+'ABoVE_vs_Macrotopology/slope.png',dpi=500)
    #plt.show()

    fig,ax = plt.subplots(1,3,figsize=(20,10))
    ax = ax.flatten()

    #aspect
    m = 0
    for k in range(0,len(vwc_depth)):
        for i in range(0,len(aspect_bins)-1):
            ax[k].boxplot(list_aspect[i][vwc_depth[m]].values,positions=[i],widths=0.6)
        ax[k].set_title('Aspect vs. '+str(vwc_depth[m]))
        ax[k].set_xticklabels(aspect_names,rotation=45, horizontalalignment='right')
        ax[k].set_xlabel('Aspect Range (deg)',fontsize=18)
        ax[k].set_ylim(0,1.0)
        ax[k].set_ylabel('VWC (%/100)',fontsize=18 )
        m=m+1
    plt.savefig(path.figures+'ABoVE_vs_Macrotopology/aspect.png',dpi=500)
    plt.close()


    fig,ax = plt.subplots(1,3,figsize=(20,10))
    ax = ax.flatten()

    #sunfactor
    n = 0
    for l in range(0,len(vwc_depth)):
        for i in range(0,len(sunfactor_bins)-1):
            ax[l].boxplot(list_sunfactor[i][vwc_depth[n]].values,positions=[i],widths=0.6)
        ax[l].set_xlabel('Sun Factor Range',fontsize=18)
        ax[l].set_ylim(0,1)
        ax[l].set_ylabel('VWC (%/100)',fontsize=18 )
        ax[l].set_title('Sun Factor vs. '+str(vwc_depth[n]))
        ax[l].set_xticklabels(sunfactor_names,rotation=45, horizontalalignment='right')
        n=n+1

    plt.savefig(path.figures+'ABoVE_vs_Macrotopology/sunfactor.png',dpi=500)
    plt.close()
    #plt.show()

    fig,ax = plt.subplots(1,3,figsize=(20,10))
    ax = ax.flatten()
    #curvature
    for j in range(0,len(vwc_depth)):
        for i in range(0,len(curv_bins)-1):
            ax[j].boxplot(list_curv[i][vwc_depth[j]].values,positions=[i],widths=0.6)
        ax[j].set_title('Curvature vs. ' +str(vwc_depth[j]))
        ax[j].set_xticklabels(curv_names,rotation=45, horizontalalignment='right')
        ax[j].set_xlabel('Curvature Range (deg)',fontsize=18)
        ax[j].set_ylim(0,1)
        ax[j].set_ylabel('VWC (%/100)',fontsize=18 )

    plt.savefig(path.figures+'ABoVE_vs_Macrotopology/curv.png',dpi=500)
    plt.close()

    fig,ax = plt.subplots(1,3,figsize=(20,10))
    ax = ax.flatten()
    #curvature
    for j in range(0,len(vwc_depth)):
        for i in range(0,len(ndvi_bins)-1):
            ax[j].boxplot(list_ndvi[i][vwc_depth[j]].values,positions=[i],widths=0.6)
        ax[j].set_title('NDVI vs. ' +str(vwc_depth[j]))
        ax[j].set_xticklabels(ndvi_names,rotation=45, horizontalalignment='right')
        ax[j].set_xlabel('NDVI Range',fontsize=18)
        ax[j].set_ylim(0,1)
        ax[j].set_ylabel('VWC (%/100)',fontsize=18 )

    plt.savefig(path.figures+'ABoVE_vs_Macrotopology/ndvi.png',dpi=500)
    plt.close()

    fig,ax = plt.subplots(1,3,figsize=(20,10))
    ax = ax.flatten()
    #curvature
    for j in range(0,len(vwc_depth)):
        for i in range(0,len(fa_bins)-1):
            ax[j].boxplot(list_fa[i][vwc_depth[j]].values,positions=[i],widths=0.6)
        ax[j].set_title('Flow Accumulation vs. ' +str(vwc_depth[j]))
        ax[j].set_xticklabels(fa_names,rotation=45, horizontalalignment='right')
        ax[j].set_xlabel('Flow Accumulation Range',fontsize=18)
        ax[j].set_ylim(0,1)
        ax[j].set_ylabel('VWC (%/100)',fontsize=18 )

    plt.savefig(path.figures+'ABoVE_vs_Macrotopology/fa.png',dpi=500)
    plt.close()

def alt_above_topo_comp(df,paths):
    slope = 'Slope'
    curvature = 'curvature'
    aspect = 'aspect'
    sunfactor = 'Sun_Factor'

    #sunfactor calculating from Min's email
    df = sunfactor_func(df,aspect,sunfactor)

    slope_bins = np.arange(0,95,5)
    aspect_bins = np.arange(0,390,30)
    curvature_bins = np. arange(-1,1,0.1)
    sunfactor_bins = np.arange(-1,1,0.1)

    #set up for loop lists
    bins = [slope_bins,aspect_bins,curvature_bins,sunfactor_bins]
    topo_name = [slope,aspect,curvature,sunfactor]


    fig,ax = plt.subplots(1,4,figsize=(25,10))
    ax = ax.flatten()
    #cycle through topo name
    for i in range(0,len(bins)):

        curv_names = []
        for j in range(0,len(bins[i])-1):
            ranged = df[((df[topo_name[i]] > float(bins[i][j])) & (df[topo_name[i]] <= float(bins[i][j+1])))]
            curv_names.append(str(round(bins[i][j],2))+' - '+str(round(bins[i][j+1],2)))
            ax[i].boxplot(ranged['alt'],widths=0.6,positions=[j])
        ax[i].set_title(topo_name[i].capitalize()+' vs. ' +' ALT ')
        ax[i].set_xticklabels(curv_names,rotation=45, horizontalalignment='right')
        ax[i].set_xlabel(topo_name[i].capitalize()+' Range (deg)',fontsize=18)
        ax[i].set_ylim(0,df['alt'].max()+(df['alt'].max()/10.0))
        ax[i].set_ylabel('ALT (m)',fontsize=18 )
    plt.savefig(paths.figures +'ALT/'+ 'alt_vs_macro_topo.png',dpi = 500)
    plt.close

def r_square_indi(x,y,ind,null_value):
    #combine x and y into same dataframe
    df = pd.DataFrame({'x':x, 'y':y,'Index':ind})
    #get rid of null value rows (like perfectly flat surface where the DEM has no values)
    flat = np.where(df['y'] == null_value)
    df_non_null = df.index.difference(df['Index'].iloc[flat])
    df = df.loc[df_non_null]

    #fill values with nan value with average
    y = df.y.fillna(df.y.mean())

    #r-squared analysis
    x= np.array(df.x.values).reshape((-1, 1))
    y = y.values
    model = LinearRegression().fit(x,y)
    r_sq = model.score(x,y)
    print('r_sq for linear regression:',r_sq)
    print('adj. r_sq for linear regression:',1 - (1-model.score(x, y))*(len(y)-1)/(len(y)-x.shape[1]-1))
    print('intercept:', model.intercept_)
    print('slope:', model.coef_)
    stat, p = stats.shapiro(df['y'].values)
    print('W Value: ',stat)
    xx = np.linspace(df.x.min(),df.x.max(),100)
    plt.plot(xx,model.intercept_ + model.coef_* xx,'k--')

    #plt.scatter(df.x,df.y,alpha=0.2)

    plt.xlim(df.x.min(),df.x.max())
    #plt.show()
    return {'r_sq':r_sq,'intercept':model.intercept_,'slope':model.coef_[0]}
def r_square_append(x,y,ind,null_value,ax,verbose=False,font_size=18):
    #combine x and y into same dataframe
    df = pd.DataFrame({'x':x, 'y':y,'Index':ind})
    #get rid of null value rows (like perfectly flat surface where the DEM has no values)
    flat = np.where(df['y'] == null_value)
    df_non_null = df.index.difference(df['Index'].iloc[flat])
    df = df.loc[df_non_null]

    #fill values with nan value with average
    y = df.y.fillna(df.y.mean())

    #r-squared analysis
    x= np.array(df.x.values).reshape((-1, 1))
    y = y.values
    model = LinearRegression().fit(x,y)
    r_sq = model.score(x,y)
    stat, p = stats.shapiro(df['y'].values)
    if verbose =='True':
        print('r_sq for linear regression:',r_sq)
        print('adj. r_sq for linear regression:',1 - (1-model.score(x, y))*(len(y)-1)/(len(y)-x.shape[1]-1))
        print('intercept:', model.intercept_)
        print('slope:', model.coef_)
        print('W Value: ',stat)
        print('p Value: ',p)
    xx = np.linspace(df.x.min(),df.x.max(),100)

    ax.plot(xx,model.intercept_ + (model.coef_* xx),'k--',linewidth=2)
    ax.text(0.05,0.9,r'R $^2$ : '+str(round(r_sq,3)),horizontalalignment='left',transform=ax.transAxes,fontsize=font_size)
    #ax.text(0.95,0.9,'W : '+str(round(stat,3)),horizontalalignment='right',transform=ax.transAxes,fontsize=font_size)
    #ax.set_xlim(df.x.min(),df.x.max())

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def no_nans(df,cols,null_values):
    orig_len = len(df)
    #number of columns
    num_col = len(cols)

    #drop all rows with nan values
    df = df.dropna()

    #scroll through columns to get rid of stuff
    for i in range(0,num_col):
        df = df.iloc[(np.where(df[cols[i]] != null_values[i]))]

    final_len = len(df)
    dropped = orig_len - final_len
    print('Total rows dropped: '+str(dropped))

    return df

def loadings(pca,features,components):

    loadings = pca.components_
    results = pd.DataFrame(loadings)
    results.columns = features
    results.index=components
    print(results)

def pca_test(df,features,target,path):
    #get rid of nans
    df = no_nans(df,['slope','curvature','aspect'],[0.0,0.0,0.0])

    # Separating out the features
    x = df.loc[:, features].values

    # Separating out the target
    y = pd.DataFrame(data = df.loc[:,[target]].values,columns =[target])

    # Standardizing the features
    x = StandardScaler().fit_transform(x)

    #determining number of dimensions to use to explain >80% of the variance
    pca2 = PCA()
    X_r = pca2.fit(x).transform(x)


    fig,ax = plt.subplots(figsize=(8,8))
    ax.plot(pca2.explained_variance_,linewidth=2.0)
    ax.set_title('Eigenvalues from Principal Component Analysis',fontsize=18)
    ax.set_ylabel('Eigenvalue',fontsize=16)
    ax.set_xlabel('Principal Components',fontsize=16)
    plt.axhline(y=1, linewidth=1, color='r', alpha=0.5,linestyle='--',label='Kaiser Criterion Threshold')
    plt.legend(fontsize=16)
    plt.ylim(0,max(pca2.explained_variance_+0.5))
    plt.savefig(path.figures+'PCA_analysis/eigenvalues.png',dpi=500)
    plt.close()

    #grabbing all eigenvalues above 1.0
    print(np.where(pca2.explained_variance_ >= 1.0))

    #minimizing dimensions
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])

    finalDf = pd.concat([principalDf, y], axis = 1)

    #downsampling
    df_sample = finalDf.sample(frac=0.05)

    #down sample for plotting
    target_bins = np.arange(0.0,0.95,0.1)
    target_names= []
    list_target = []

    #create column
    df_sample[target+'_cat'] =0.0

    #separate pandas column into value 0-5 etc through slope
    for i in range(0,len(target_bins)-1):
        df_sample[target+'_cat'] = np.where((df_sample[target] > float(target_bins[i])) & (df_sample[target] <= float(target_bins[i+1])), str(round(target_bins[i],2))+' - '+str(round(target_bins[i+1],2)),df_sample[target+'_cat'])
        list_target.append( df_sample[((df_sample[target] > float(target_bins[i])) & (df_sample[target] <= float(target_bins[i+1])))])
        target_names.append(str(round(target_bins[i],2))+' - '+str(round(target_bins[i+1],2)))


    #plotting the data
    fig,ax = plt.subplots(figsize = (12,8))
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('PCA 1 vs. 3', fontsize = 20)
    alpha_bins =np.linspace(1.0,0.3,len(list_target))
    for j in range(0,len(list_target)):
        ax.scatter(list_target[j].loc[:, 'principal component 1'], list_target[j].loc[:, 'principal component 2'], s = 50,label= 'VWC (%/100): '+target_names[j],alpha=alpha_bins[j])
    print('Explained Variance Ratio: ',pca.explained_variance_ratio_)
    print('Explained Variance: ',pca.explained_variance_)
    loadings(pca,features,['principal component 1', 'principal component 2','principal component 3'])
    plt.text(0.25,0.8,'Explained Variance Ratio: '+str(round(pca.explained_variance_ratio_[0],3))+', '+str(round(pca.explained_variance_ratio_[1],3)),horizontalalignment='left',transform=plt.gcf().transFigure)
    plt.legend(loc=1)
    plt.savefig(path.figures+'PCA_analysis/PCA1_vs_PCA2.png',dpi=500)
    plt.close()

    fig,ax = plt.subplots(figsize = (12,8))
    ax.set_xlabel('Principal Component 2', fontsize = 15)
    ax.set_ylabel('Principal Component 3', fontsize = 15)
    ax.set_title('PCA 2 vs. 3', fontsize = 20)
    for j in range(0,len(list_target)):
        ax.scatter(list_target[j].loc[:, 'principal component 2'], list_target[j].loc[:, 'principal component 3'], s = 50,label= 'VWC (%/100): '+target_names[j],alpha=alpha_bins[j])

    plt.text(0.25,0.8,'Explained Variance Ratio: '+str(round(pca.explained_variance_ratio_[0],3))+', '+str(round(pca.explained_variance_ratio_[1],3))+', '+str(round(pca.explained_variance_ratio_[2],3)),horizontalalignment='left',transform=plt.gcf().transFigure)
    plt.legend(loc=1)
    plt.savefig(path.figures+'PCA_analysis/PCA2_vs_PCA3.png',dpi=500)
    plt.close()

    fig,ax = plt.subplots(figsize = (12,8))
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 3', fontsize = 15)
    ax.set_title('PCA 1 vs. 3', fontsize = 20)
    for j in range(0,len(list_target)):
        ax.scatter(list_target[j].loc[:, 'principal component 1'], list_target[j].loc[:, 'principal component 3'], s = 50,label= 'VWC (%/100): '+target_names[j],alpha=alpha_bins[j])

    plt.text(0.25,0.8,'Explained Variance Ratio: '+str(round(pca.explained_variance_ratio_[0],3))+', '+str(round(pca.explained_variance_ratio_[1],3))+', '+str(round(pca.explained_variance_ratio_[2],3)),horizontalalignment='left',transform=plt.gcf().transFigure)
    plt.legend(loc=1)
    plt.savefig(path.figures+'PCA_analysis/PCA1_vs_PCA3.png',dpi=500)
    plt.close()
def curve_fitting(x,y,ind,null_value):

    #combine x and y into same dataframe
    df = pd.DataFrame({'x':x, 'y':y,'Index':ind})
    #get rid of null value rows (like perfectly flat surface where the DEM has no values)
    flat = np.where(df['y'] == null_value)
    df_non_null = df.index.difference(df['Index'].iloc[flat])
    df = df.loc[df_non_null]

    #fill values with nan value with average
    y = df.y.fillna(df.y.mean())

    #fit
    popt,pcov = curve_fit(func,df['x'].values,y)

    x_gen = np.linspace(0,1,100)

    df['y'] = y

    df_sample = df.sample(frac=0.05)

    plt.scatter(df_sample['x'],df_sample['y'],label='5% Sample')

    plt.plot(x_gen, func(x_gen, *popt), 'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.legend()
    plt.xlim(0,1)
    plt.show()

def multivariate_linear_correlation_matrix(df,features,target,path,fraction):

    #setting the number of columns/rows
    columns, rows = len(features)+1,len(features)+1

    plotting_x = np.append(features,target)
    plotting_y = np.append(features,target)

    #create subplots
    fig,ax = plt.subplots(columns,rows,figsize=(25,25))
    #plt.tight_layout()

    #take a subset of the data
    df_sample = df.sample(frac=fraction)

    #scroll through the subplots fitting each with an r^2 this assumes a linear relationship
    for y in range(0,columns):
        for x in range(0,rows):
            print(x,y)
            ax[y,x].scatter(df_sample[plotting_x[x]],df_sample[plotting_y[y]],alpha=0.1)
            if y==columns-1:
                ax[y,x].set_xlabel(plotting_x[x],fontsize=16)
            if x ==0:
                ax[y,x].set_ylabel(plotting_y[y],fontsize=16)
            print(plotting_x[x],plotting_y[y])
            r_square_append(df[plotting_x[x]],df[plotting_y[y]],df['Index'],0.0,ax[y,x])

    plt.savefig(path.figures+'Linear_Regression_Matrix/'+target+'.png',dpi=500)
    plt.close()
    plt.clf()

def linear_machine_learning_SM(df,variables,target):

    #separate the variable from the target variable
    X = df[variables]
    y = df[[target]]

    #separate test vs. training data for the model
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,random_state=1)





    print('Multiple Linear Regression Model: ')
    #initiate a linear regression function
    reg = LinearRegression()

    #set the model to train on specific variables
    reg.fit(X_train[variables], y_train)
    y_predicted = reg.predict(X_test[variables])
    print("Mean squared error: %.2f" % mean_squared_error(y_test, y_predicted))
    print('R²: %.2f' % r2_score(y_test, y_predicted))


    print('OLS Model: ')
    #try with OLS model too
    OLS_model = sm.OLS(y_train, X_train[variables]).fit()
    y_pred2 = OLS_model.predict(X_test[variables])
    print("Mean squared error: %.2f" % mean_squared_error(y_test, y_pred2))
    print('R²: %.2f' % r2_score(y_test, y_pred2))
    print(OLS_model.summary())

    print('Non-Linear regression: ')
    # Alpha (regularization strength) of LASSO regression
    lasso_eps = 0.0001
    lasso_nalpha=20
    lasso_iter=5000
    # Min and max degree of polynomials features to consider
    degree_min = 2
    degree_max = 8

    #grabbing values


    # Make a pipeline model with polynomial transformation and LASSO regression with cross-validation, run it for increasing degree of polynomial (complexity of the model)
    for degree in range(degree_min,degree_max+1):
        fig, ax = plt.subplots(figsize=(20,20))
        model = make_pipeline(PolynomialFeatures(degree, interaction_only=False), LassoCV(eps=lasso_eps,n_alphas=lasso_nalpha,max_iter=lasso_iter,normalize=True,cv=5))
        model.fit(X_train,y_train)
        test_pred = np.array(model.predict(X_test))
        #ax.scatter(y_test,test_pred)
        #RMSE=np.sqrt(np.sum(np.square(test_pred-y_test.values)))
        test_score = model.score(X_test,y_test)
        print(test_score)

    '''
    #plot the y_predictions vs. measured value
    fig, ax = plt.subplots(2,1)
    ax=ax.flatten()
    ax[0].scatter(y_test, y_predicted)
    ax[1].scatter(y_test,y_pred2,color='r')
    ax[0].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=4)
    ax[1].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=4)
    ax[0].set_xlabel('measured')
    ax[0].set_ylabel('predicted')
    ax[0].set_ylim(0,1.0)
    ax[1].set_ylim(0,1.0)
    plt.show()
    '''
    plt.show()

def comp_landcover(df,land_1,land_2,path,gapland = 'gaplandfire_legend.csv'):

    #getting gaplandfire table for comparison
    gapland = path_file_reader(path,gapland)

    grouped = df.groupby(land_1)

    #fig,ax = plt.subplots(4,4,figsize=(10,10))
    #ax = ax.flatten()
    nlcd_dic = {'11':'Open Water','12':'Perrenial Ice/Snow','21':'Developed, Open Space','22': 'Developed, Low Intensity','23':'Developed, Medium Intensity','24':'Developed, High Intensity','31':'Barren Land','41':'Deciduous Forest','42':'Evergreen Forest',
    '43':'Mixed Forest','51':'DwarfShrub','52':'Shrub/Scrub','71':'Grassland/Herbaceous','72':'Sedge/Herbaceous','74':'Moss','81':'Pasture/Hay','82':'Cultivated Crops','90':'Woody Wetlands','95':'Emergent Herbaceous Wetlands'}

    i=0
    for key,group in grouped:
        fig,ax = plt.subplots(figsize=(15,12))
        #ax = ax.flatten()
        #summing values for each unique gapland value
        counts = group[land_2].value_counts().sort_values(ascending=False)

        #max values and index
        max_frac = counts.max()/counts.sum()
        max_id = counts.idxmax()

        print('NLCD: '+ nlcd_dic[str(round(key))])
        #print('Gapland: '+str(max_id)+'  '+str(max_frac))

        counts = counts.to_frame()
        counts['above_h_value']=counts.index
        merged = pd.merge(counts,gapland,on='above_h_value')
        merged=merged.rename(columns={'gaplandfire_val':'counts'})
        merged['fraction']= merged['counts']/counts['gaplandfire_val'].sum()
        merged_2 = merged.groupby('Ecosystem_LU')['counts'].sum().sort_values(ascending=False).to_frame()
        merged_2['fraction'] = merged_2['counts']/merged_2['counts'].sum()
        #plotting this up


        x =(merged_2.index)[:3]
        y = merged_2['fraction'][:3]

        '''
        ax[i].bar(x,y)
        ax[i].set_title(nlcd_dic[str(round(key))])
        ax[i].set_xticklabels(x,rotation=90)
        ax[i].set_ylim(0,1)
        '''
        ax.bar(x,y)
        ax.set_title(nlcd_dic[str(round(key))])
        ax.set_xticklabels(x)
        ax.set_ylim(0,1)
        plt.savefig(work_paths.figures+'Comparing_NLCD_GaplandFire/'+str(round(key))+'.png')
        plt.close()
        i = i+1

    return counts
    #plt.show()
def gams(df,variables,target):
    var2 = variables

    var2.append(target)
    df_all = df[var2].dropna()
    print(df_all,len(df_all))
    X = df_all[variables].values
    y = df_all[target].values
    print(X,y)
    gam = LinearGAM(n_splines=10).gridsearch(X, y)
    print(gam.summary())

def five_year_avg_precip(path):
    #dfs = [pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()]
    #df_after = []

    #grab files in the path folder
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    print(onlyfiles)
    df_all = pd.DataFrame()
    for filepath in onlyfiles:
        #make dataframe of each file
        ds = xr.open_dataset(path+filepath)
        df = ds.to_dataframe()

        #crop to the Seward Peninsula
        lat_bounds = [64.161307,66.521784]
        lon_bounds = [-167.208134,-162.098092]
        df = df[(df['stn_lon'] > lon_bounds[0]) & (df['stn_lon'] < lon_bounds[1])]
        df = df[(df['stn_lat'] > lat_bounds[0]) & (df['stn_lat'] < lat_bounds[1])]

        #append to one datframe
        df_all = df_all.append(df)
        del df
    return df_all

def SAR_Plot_boxplots_depth(above_df,in_situ,soil_samples,add_samples=False):
    grouped = above_df.groupby('SAR_Plot')
    grouped_samples = soil_samples.groupby('SAR_Plot')

    depth_array = ['6','12','20']
    above_col_depth = ['lbc_0.06_aug','lbc_0.12_aug','lbc_0.2_aug']

    fig,ax = plt.subplots(figsize=(20,10))
    pos = 0
    key_array = []
    sample_range = []
    for key,group in grouped:
        key_array.append(key)
        in_situ_grouped = in_situ.groupby(['SAR_Plot','VWC_Measurement_Depth'])
        for key1,group1 in in_situ_grouped:
            if key1[0] == key:
                print(key1)
                if str(key1[1]) == '6':
                    above_col = above_col_depth[0]
                    colors = ['lightblue', 'salmon']
                    label = ['6cm','']
                if str(key1[1]) == '12':
                    above_col = above_col_depth[1]
                    colors = ['lightblue', 'salmon']
                    label = ['12cm','']
                if str(key1[1]) == '20':
                    above_col = above_col_depth[2]
                    colors = ['lightblue', 'salmon']
                    label = ['20cm','']



                #plotting boxplots
                box_insitu= ax.boxplot([group1['VWC'],group[above_col]*100.0],positions=[pos,pos+1],widths=0.5,patch_artist=True)

                #setting colors
                for patch, color in zip(box_insitu['boxes'], colors):
                    patch.set_facecolor(color)
                pos=pos+2
        if add_samples== True:
            #grab group of soil_samples
            sar_sample = grouped_samples.get_group(key)

            #setting labels for range of soil sample
            sample_range.append(str(sar_sample['Depth_min'].min())+' - '+str(sar_sample['Depth_max'].max()))

            box_sample =ax.boxplot(sar_sample['Volumetric_water_content']*100.0,positions=[pos],widths=0.5,patch_artist=True)

            #setting color
            colors=['palegreen']
            for patch, color in zip(box_sample['boxes'], colors):
                patch.set_facecolor(color)
            pos=pos+1
    if add_samples==True:
        x_ticks = []
        start = 0.5
        for i in range(0,len(key_array)):
            x_ticks.append((np.arange(start,start+5,2)).tolist())
            x_ticks.append([start+5.5])
            start=start+7
        x_ticks_flat = [val for sublist in x_ticks for val in sublist]

        labels = []
        for j in range(0,len(key_array)):
            labels.append(['6','12','20',sample_range[j]])
        labels_flat = [val for sublist in labels for val in sublist]
        section=7
        save_name = 'above_vs_hydro/boxplot_in_situ_vs_above_vs_sample.png'
        ax.legend([box_insitu["boxes"][0],box_insitu["boxes"][1],box_sample["boxes"][0]], ['In-Situ TDR', 'ABoVE','Soil Samples'], loc='lower right',fontsize=16)

    else:
        labels_flat = ['6','12','20']*len(key_array)
        x_ticks_flat =np.arange(0.5,pos+0.5,2)
        section =6
        save_name = 'above_vs_hydro/boxplot_in_situ_vs_above.png'
        ax.legend([box_insitu["boxes"][0],box_insitu["boxes"][1]], ['In-Situ TDR', 'ABoVE'], loc='lower right',fontsize=16)

    #plotting lines and sar plot text
    vert_lines =np.arange(section-0.5,pos+0.5,section)
    text_x = np.arange(section/2.0,pos+0.5,section)
    for j in range(0,len(vert_lines)):
        plt.axvline(x=vert_lines[j],linestyle='--')
        plt.text(text_x[j]-0.5,105,key_array[j],ha='center',fontsize=16)

    ax.set_xlim(-0.5, pos-0.5)

    plt.xticks(x_ticks_flat,labels_flat,ha='center',fontsize=14,rotation=90)
    plt.xlabel('Depth (cm)',fontsize=20)
    plt.ylabel('Volumetric Water Content (%)',fontsize=20)
    plt.yticks(fontsize=14)
    plt.savefig(work_paths.figures+save_name,dpi=500)
    plt.close()

def SAR_Plot_boxplots_plot(above_df,in_situ,soil_samples,add_samples=False):
    grouped = above_df.groupby('SAR_Plot')
    grouped_samples = soil_samples.groupby('SAR_Plot')

    depth_array = ['6','12','20']
    above_col_depth = ['lbc_0.06_aug','lbc_0.12_aug','lbc_0.2_aug']

    #average column of in_situ
    above_df['plot_avg'] = above_df[above_col_depth].mean(axis=1)

    fig,ax = plt.subplots(figsize=(20,10))
    pos = 0
    key_array = []
    sample_range = []
    in_situ_grouped = in_situ.groupby('SAR_Plot')


    #for loop through sar plots
    for key,group in grouped:
        if add_samples== True:
            #grab group of soil_samples
            sar_sample = grouped_samples.get_group(key)

            #setting labels for range of soil sample
            sample_range.append(str(sar_sample['Depth_min'].min())+' - '+str(sar_sample['Depth_max'].max()))

            box_sample =ax.boxplot(sar_sample['Volumetric_water_content']*100.0,positions=[pos],widths=0.5,patch_artist=True)

            #setting color
            colors=['palegreen']
            for patch, color in zip(box_sample['boxes'], colors):
                patch.set_facecolor(color)
            pos=pos+1

        key_array.append(key)
        in_situ_sample = in_situ_grouped.get_group(key)

        #plotting boxplots
        box_insitu= ax.boxplot([in_situ_sample['VWC'],group['plot_avg']*100.0],positions=[pos,pos+1],widths=0.5,patch_artist=True)

        colors = ['lightblue', 'salmon','lightgreen']
        #setting colors
        for patch, color in zip(box_insitu['boxes'], colors[:2]):
            patch.set_facecolor(color)
        pos=pos+2



    #setting file names
    if add_samples==True:
        save_name = 'above_vs_hydro/boxplot_in_situ_vs_above_vs_sample_plotavg.png'
        section = 3
        ax.legend([box_sample["boxes"][0],box_insitu["boxes"][0],box_insitu["boxes"][1]], ['Soil Samples','In-Situ TDR', 'ABoVE'], loc='lower right',fontsize=16)

    else:
        save_name = 'above_vs_hydro/boxplot_in_situ_vs_above_plot_avg.png'
        section = 2
        ax.legend([box_insitu["boxes"][0],box_insitu["boxes"][1]], ['In-Situ TDR', 'ABoVE'], loc='lower right',fontsize=16)
    #plotting lines and sar plot text
    vert_lines =np.arange(section-0.5,pos+0.5,section)
    text_x = np.arange(section/2.0,pos+0.5,section)
    for j in range(0,len(vert_lines)):
        plt.axvline(x=vert_lines[j],linestyle='--')
        plt.text(text_x[j]-0.5,105,key_array[j],ha='center',fontsize=16)

    ax.set_xlim(-0.5, pos-0.5)
    ax.set_xticks([])
    ax.set_xticklabels(['']*pos,rotation=90)
    #plt.xlabel('Depth (cm)',fontsize=20)
    plt.ylabel('Volumetric Water Content (%)',fontsize=20)
    plt.yticks(fontsize=14)
    plt.savefig(work_paths.figures+save_name,dpi=500)
    plt.close()

def soil_samples_vs_TDR(df,x,y,path):
    #use ggplot
    matplotlib.style.use('ggplot')

    #set up figure
    fig,ax = plt.subplots(figsize=(10,7))


    #add the data
    plt.scatter(df[x],df[y],color='firebrick',s=75)
    r_sq = r_square_append(df[x],df[y],df.index,0.0,ax,verbose=True)
    plt.xlabel('In-Situ VWC(%)',fontsize=18)
    plt.ylabel('Soil Sample VWC(%)',fontsize=18)
    plt.title('In-Situ vs. Soil Sample at SAR Plots',fontsize=24)
    plt.xlim(30,75)
    plt.ylim(30,75)
    plt.savefig(path.figures+'above_vs_hydro/lbc_vs_samples_rsq.png',dpi=500)
    plt.close()

def SM_change(df,in_situ,vegtype,depth_column,path):
    """Function intended to differentiate soil moisture in landcover types this assumes oct-aug.

    Parameters
    ----------
    df : Pandas DataFrame
        Includes all non-NAN values
    vegtype : string
        pandas column name
    depth_column : string
        pandas column name

    Returns
    -------
    plot
        Description of returned object.

    """

    nlcd_dic = {'11':'Open Water','12':'Perrenial Ice/Snow','21':'Developed, Open Space','22': 'Developed, Low Intensity','23':'Developed, Medium Intensity','24':'Developed, High Intensity','31':'Barren Land','41':'Deciduous Forest','42':'Evergreen Forest',
    '43':'Mixed Forest','51':'DwarfShrub','52':'Shrub/Scrub','71':'Grassland/Herbaceous','72':'Sedge/Herbaceous','74':'Moss','81':'Pasture/Hay','82':'Cultivated Crops','90':'Woody Wetlands','95':'Emergent Herbaceous Wetlands'}

    grouped = df.groupby(vegtype)

    key_array =[]
    if len(depth_column) == 1:


        fig,ax = plt.subplots(figsize=(15,15))
        pos= 0


        for key,group in grouped:
            plt.bar(x=pos,height=group[depth_column].mean()*100.0,yerr=(group[depth_column].std()*100.0)/np.sqrt(len(group[depth_column])),label=nlcd_dic[str(round(key))])
            print(nlcd_dic[str(round(key))],len(group[depth_column]))
            key_array.append(nlcd_dic[str(round(key))])
            pos = pos+1

    else:
        fig,ax = plt.subplots(figsize=(15,15))
        pos= 0
        key_array =[]
        w=0.25
        bar_heights = grouped[depth_column].mean()
        bar_error =(grouped[depth_column].std()*100.0).div(np.sqrt(grouped[depth_column].size()),axis=0)

        #positions for bars
        r1 = np.arange(len(grouped[depth_column[0]]))
        r2 = [x + w for x in r1]
        r3 = [x + w for x in r2]

        #names for groups
        for key,group in grouped:
            key_array.append(nlcd_dic[str(round(key))])

        #plotting
        plt.bar(r1,bar_heights[depth_column[0]]*100.0,yerr=bar_error[depth_column[0]],width=w,label='6 cm')
        plt.bar(r2,bar_heights[depth_column[1]]*100.0,yerr=bar_error[depth_column[1]],width=w,label='12 cm')
        plt.bar(r3,bar_heights[depth_column[2]]*100.0,yerr=bar_error[depth_column[2]],width=w,label='20 cm')


    #setting label names
    plt.xticks(r2,key_array,rotation=90,fontsize=16)

    #making the plot smaller so we can read text
    plt.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.48)

    #insert horizontal line
    plt.axhline(y=0.0,c='k',ls='--',lw=2.0)

    #setting up axes
    plt.ylabel(r'$\Delta$ '+'Volumetric Water Content (Oct-Aug) (%)',fontsize=18)
    plt.ylim(-2.0,80.0)

    #make legend
    plt.legend(loc=0,fontsize=20)

    #plot title
    plt.title('Comparing the (NLCD) with October 2017 Soil Moisture',fontsize=22)

    if len(depth_column) == 1:
        plt.savefig(path.figures+'Comparing_NLCD_GaplandFire/'+vegtype+'_'+depth_column+'_bar_all.png',dpi=500)
    else:
        plt.savefig(path.figures+'Comparing_NLCD_GaplandFire/'+vegtype+'_bar_all_oct_2017.png',dpi=500)
    plt.close()
    #plt.show()

def alt_sar_avg(above,in_situ,path,r_squared = False):
    #converting date for in_situ to datatime object
    in_situ['Date'] = pd.to_datetime(in_situ['Date'],errors='coerce')

    #drop nan values -9999
    in_situ =  in_situ[in_situ['Thaw_Depth'] != -9999]

    #groupby to make
    grouped_ab =above.groupby('SAR_Plot')
    grouped_insitu = in_situ.groupby('SAR_Plot')

    #get sar plot names
    sar_plots = above['SAR_Plot'].dropna().unique()

    if r_squared == False:
        fig,ax = plt.subplots(figsize=(20,10))
        pos = 0
        key_array = []
        sample_range = []
        #scan through to make boxplots of comparison
        for i in range(0,len(sar_plots)):
            ab_sar = grouped_ab.get_group(sar_plots[i])
            in_situ_sar = grouped_insitu.get_group(sar_plots[i])

            #plotting boxplots
            box_insitu= ax.boxplot([in_situ_sar['Thaw_Depth'],ab_sar['alt_x']],positions=[pos,pos+1],widths=0.5,patch_artist=True)

            colors = ['lightblue', 'salmon']
            #setting colors
            for patch, color in zip(box_insitu['boxes'], colors[:2]):
                patch.set_facecolor(color)
            pos=pos+2


        save_name = 'ALT/boxplot_in_situ_vs_above_ALT_2017.png'

        #legend
        ax.legend([box_insitu["boxes"][0],box_insitu["boxes"][1]], ['In-Situ TDR', 'ABoVE'], loc='upper right',fontsize=20)

        #plotting lines and sar plot text
        section = 2
        vert_lines =np.arange(section-0.5,pos+0.5,section)
        text_x = np.arange(section/2.0,pos+0.5,section)
        for j in range(0,len(vert_lines)):
            plt.axvline(x=vert_lines[j],linestyle='--')
            #plt.text(text_x[j]-0.5,105,sar_plots[j],ha='center',fontsize=16)

        #setting up the ticks and labels
        start = 0.5
        xticks = np.arange(start,pos,2)

        ax.set_xticks(xticks)
        ax.set_xticklabels(sar_plots,fontsize=16)
        #plt.xlabel('Depth (cm)',fontsize=20)
        plt.ylabel('Active Layer Depth (cm)',fontsize=20)
        plt.yticks(fontsize=16)
        plt.title('ABoVE vs In-Situ ALT',fontsize=24)
        plt.savefig(work_paths.figures+save_name,dpi=500)
        plt.close()
    else:
        fig,ax = plt.subplots(figsize=(15,10))

        save_name = 'ALT/sar_plot_r_squared_ALT_2017.png'
        x = grouped_insitu['Thaw_Depth'].mean().to_list()
        y = (grouped_ab['alt_x']).mean().to_list()
        y = [i * 100.0 for i in y]
        plt.scatter(x,y)

        r_square_append(x,y,np.array(sar_plots),0.0,ax)
        plt.xlim(20,80)
        plt.ylim(20,80)
        plt.xlabel('In-Situ ALT (cm)',fontsize=20)
        plt.ylabel('ABoVE ALT (cm)',fontsize=20)
        plt.title('In-Situ vs. ABoVE Active Layer Thickness (ALT)',fontsize=24)
        plt.savefig(work_paths.figures+save_name,dpi=500)
        plt.close()

def separate_depth_plot_avg(df_avg,path):
    grouped = df_avg.groupby(['Site','Depth'])
    fig,ax = plt.subplots(1,3,figsize=(15,5))
    ax=ax.flatten()
    print(grouped.describe())
    #plot Parameters
    symbol =['^','s','o']
    set = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999"]
    colors = {'Teller':set[2],'Kougarok':set[1],'Council':set[0]}
    size = [10,10,10]

    #append
    x6_array = []
    y6_array = []
    x6err_array =[]

    x12_array = []
    y12_array = []
    x12err_array =[]

    x20_array = []
    y20_array = []
    x20err_array =[]
    i = 0
    for key,group in grouped:
        if key[1] == 6:
            j = 0
        if key[1] == 12:
            j = 1
        if key[1] == 20:
            j = 2

        print(key)
        ax[j].errorbar(group['In_Situ'],group['Model'],xerr=group['In_Situ_std'],yerr=group['Model_std'],markersize=size[j],color=colors[key[0]],fmt=symbol[j])
        ax[j].set_ylim(0,1)
        ax[j].set_xlim(0,1)
        i=i+1
    grouped_depth = df_avg.groupby('Depth')
    l=0
    for key,group in grouped_depth:
        r2_other = r_square_append(group['In_Situ'],group['Model'],np.arange(0,len(group['In_Situ'])),0.0,ax[l])
        l=l+1

    #creating custom legends
    legend_elements = [ Line2D([0], [0], marker='^', color='w', label='6cm',
                          markerfacecolor='grey', markersize=10,markeredgecolor='k',markeredgewidth=2),
                      Line2D([0], [0], marker='s', color='w', label='12cm',
                            markerfacecolor='grey', markersize=15,markeredgecolor='k',markeredgewidth=2),
                        Line2D([0], [0], marker='o', color='w', label='20cm',
                            markerfacecolor='grey', markersize=25,markeredgecolor='k',markeredgewidth=2)]

    legend_elements2 = [ Line2D([0], [0], marker='o', color='w', label='Teller',
                          markerfacecolor=set[2], markersize=15,markeredgecolor='k'),
                      Line2D([0], [0], marker='o', color='w', label='Kougarok',
                            markerfacecolor=set[1], markersize=15,markeredgecolor='k'),
                        Line2D([0], [0], marker='o', color='w', label='Council',
                            markerfacecolor=set[0], markersize=15,markeredgecolor='k')]

    leg = fig.legend(handles=legend_elements, loc='upper right',prop={'size':16})
    fig.add_artist(leg)

    fig.legend(handles = legend_elements2, loc='upper left',prop={'size':16})
    plt.savefig(path.figures+'above_vs_hydro/Indi_depth_SAR_avg.png',dpi=500)
    plt.close()

def average_in_situ_plot_depth_year_comparison(in_situ_2017,in_situ_2019,path):
    """Function intended to make boxplots for each SAR_Plot at all depths.

    Parameters
    ----------
    in_situ : pandas DataFrame
        Description of parameter `in_situ`.
    path : type
        Description of parameter `path`.

    Returns
    -------
    plot
        Description of returned object.

    """
    groupby = ['SAR_Plot']

    #grouped_2017 = in_situ_2017.groupby(groupby)
    #grouped_2019 = in_situ_2019.groupby(groupby)

    #add year for groupby usage
    in_situ_2017['Year'] = 2017
    in_situ_2019['Year'] = 2019

    #subset each to esaily merge to one dataframe
    selected = ['Latitude','Longitude','SAR_Plot','VWC']
    in_situ_2017 = in_situ_2017[selected]
    in_situ_2019 = in_situ_2019[selected]

    print(in_situ_2017,in_situ_2019)
    #combine datasets into one
    df_merge =pd.concat(in_situ_2017,in_situ_2019)

    ax = in_situ.boxplot(by=['SAR_Plot'], column='VWC')

    set = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999"]
    colors = {'Teller':set[2],'Kougarok':set[1],'Council':set[0]}

    '''
    fig,ax = plt.subplots(figsize=(10,10))
    pos = 0
    for key,group in grouped:
        print(key)
        if len(groupby) > 1:
            grouped_depth = group.groupby(groupby[1])
            for key_2, group_2 in grouped_depth:
                print(key_2)
                print(pos)
                box = ax.boxplot(group_2['VWC'],positions=[pos],widths=0.5,patch_artist=True)

                #setting color
                colors=['palegreen']
                for patch, color in zip(box['boxes'], colors):
                    patch.set_facecolor(color)
                pos=pos+1
    '''
    plt.show()


def scatter_2017_vs_2019_comparison(df,path):
    grouped = df.groupby(['SAR_Plot'])

    x = []
    xerr = []
    y = []
    yerr =[]
    label = []
    for key,group in grouped:
        grouped_2 = group.groupby('Year')
        print(key)
        if len(grouped_2) >1:
            for key_2,group_2 in grouped_2:
                print(key_2)
                if group_2['Year'].unique() == 2017:
                    x.append(group_2['VWC'].mean())
                    xerr.append(group_2['VWC'].std())
                else:
                    y.append(group_2['VWC'].mean())
                    yerr.append(group_2['VWC'].std())
            label.append(key)
    print(x,y,len(x),len(y))
    plt.errorbar(x,y,marker='*')
    plt.xlabel('2017 Soil Moisture (%)')
    plt.ylabel('2019 Soil Moisture (%)')
    plt.title('Comparing 2017 and 2019 In-Situ Averages')
    plt.ylim(0,100)
    plt.xlim(0,100)
    plt.savefig(path.figures+'In_situ/2017_2019_comparison_scatterplot.png',dpi=500)
    plt.close()

def lidar_vs_roughness(data,path,depth='6cm'):
    '''
    Datasets used in this files
    filepath = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/combined_0.1m_may2017_3m_curv.csv'
    filepath_slope = 'Z:/AKSeward/Data/GIS/Teller/SAR/Lidar_In_situ_SM_analysis/STD_Slope/all_may_2017_slope_0.5mDTM_3mbuffer_Zstats_combined.csv'

    Purpose: Program meant to study the impact of surface roughness on soil moisture.

    '''
    print(list(data))
    SAR_plot_names = data['SAR_Plot'].unique()
    grouped = data.groupby('SAR_Plot')
    colors = {'TL_SAR_8':'b','TL_SAR_4':'r','TL_SAR_7':'g'}
    depths = [6,12,20]
    grouped = data.groupby(['SAR_Plot','VWC_Measurement_Depth'])
    for depth in depths:
        title='Surface Roughness vs. VWC at '+str(depth)+'cm'
        for key,group in grouped:
            #print(key,key[1])
            if key[1] == depth:
                print(key)
                plt.scatter(group['VWC'],group['STD'],color=colors[key[0]],label=key[0])
        print('Stop')
        plt.title(title)
        plt.xlabel('Volumetric Water Content (%)')
        plt.ylabel('Standard Deviation of Slope')
        plt.legend()
        plt.xlim(0,100)
        plt.savefig(path.figures+'Surface_Roughness/STD_0.5m_dem_slope_3m_'+str(depth)+'cm.png',dpi=500)
        plt.clf()
        plt.close()

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x'], point['y'], str(point['val']))

def flag_variance(df,sr_df):
    '''
    program to plot the average soil moisture against the standard deviation to discern how much variability is at each place
    SR_DF = file with surface roughness values
    '''
    #merge the two dataframes
    df = pd.merge(df,sr_df,on='SAR_Plot')
    #group the data for plotting
    grouped = df.groupby(['SAR_Plot','VWC_Measurement_Depth','Locale'])

    #average and std
    avg = grouped['VWC'].mean()
    std = grouped['VWC'].std()
    surf_rough = grouped['surf_rough'].mean()
    #merging dataframes
    merged = pd.merge(avg,std,on=['SAR_Plot','VWC_Measurement_Depth','Locale'])
    merged = pd.merge(merged,surf_rough,on=['SAR_Plot','VWC_Measurement_Depth','Locale'])
    #renaming dataframe
    merged = merged.rename(columns={'VWC_x':'mean','VWC_y':'std'})

    #reset index for manipulation
    merged.reset_index(inplace=True)

    #color dictionary
    colors={'Council':'green','Teller':'blue','Kougarok':'red'}
    color_map = [colors[site] for site in merged['Locale']]
    merged['colors'] = color_map

    #size dictionary
    sizes={'6':'^','12':'s','20':'o'}
    size_map = [sizes[str(size)] for size in merged['VWC_Measurement_Depth']]
    merged['mark'] = size_map

    #normalizing surface roughness
    merged['surf_rough_norm'] = merged['surf_rough']/merged['surf_rough'].max()

    #plotting routine
    fig,ax = plt.subplots(figsize=(10,10))
    grouped = merged.groupby(['VWC_Measurement_Depth'])
    for key,group in grouped:
        print(key)

        ax.scatter(group['mean'],group['std'],c=group['colors'],s=group['surf_rough_norm']*150.0,marker=group['mark'].iloc[0])

    #label points
    label_point(merged['mean'],merged['std'],merged['SAR_Plot'],ax)

    plt.xlabel('Volumetric Water Content (%)',fontsize=16)
    plt.ylabel('Standard Deviation of VWC',fontsize=16)




    #creating custom legends
    legend_elements = [ Line2D([0], [0], marker='^', color='w', label='6cm',
                          markerfacecolor='grey', markersize=20,markeredgecolor='k',markeredgewidth=2),
                      Line2D([0], [0], marker='s', color='w', label='12cm',
                            markerfacecolor='grey', markersize=20,markeredgecolor='k',markeredgewidth=2),
                        Line2D([0], [0], marker='o', color='w', label='20cm',
                            markerfacecolor='grey', markersize=20,markeredgecolor='k',markeredgewidth=2)]

    legend_elements2 = [ Line2D([0], [0], marker='o', color='w', label='Teller',
                          markerfacecolor='blue', markersize=15,markeredgecolor='k'),
                      Line2D([0], [0], marker='o', color='w', label='Kougarok',
                            markerfacecolor='red', markersize=15,markeredgecolor='k'),
                        Line2D([0], [0], marker='o', color='w', label='Council',
                            markerfacecolor='green', markersize=15,markeredgecolor='k')]

    leg = fig.legend(handles=legend_elements, loc='upper right',prop={'size':16})
    fig.add_artist(leg)

    fig.legend(handles = legend_elements2, loc='upper left',prop={'size':16})
    plt.savefig('Z:/JDann/Documents/Documents/Julian_Python/SAR_programs_20181003/Figures/In_situ/2020_02_18_variance_plots.png',dpi=500)
    plt.close()

def wavelet_analysis(raster):
    pass
