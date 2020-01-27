import pandas as pd
import numpy as np
import matplotlib.pyplot
import h5py



def main():
    filepath = 'Z:/AKSeward/EOW_Snow/2020_01_10_Baptiste_UAS_snow_cover/snow2019_assess.mat'
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


if __name__ == "__main__":
    main()
