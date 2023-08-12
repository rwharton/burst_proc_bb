import numpy as np
import time

def calc_mean_std(shape, N):
    """
    Generate random array of shape `shape` 
    then calculate the mean and standard dev
    across axis 0 N times
    """
    d = np.random.normal(size=shape)
    
    mt = np.zeros(N)
    st = np.zeros(N)

    for ii in range(N):
        print(ii)
        tm = time.time()
        mm = np.mean(d, axis=0)
        dtm = time.time() - tm
        print("mean : %.3f" %dtm)

        ts = time.time()
        ss = np.std(d, axis=0)
        dts = time.time() - ts
        print("std  : %.3f" %dts)
        print("\n")

        mt[ii] = dtm
        st[ii] = dts

    print("avg : %.2f +/- %.2f " %(np.mean(mt), np.std(mt)))
    print("std : %.2f +/- %.2f " %(np.mean(st), np.std(st)))

    return mt, st

        

def calc_mean_std_dat(d, N):
    """
    calculate the mean and standard dev
    across axis 0 N times for d
    """
    mt = np.zeros(N)
    st = np.zeros(N)

    for ii in range(N):
        print(ii)
        tm = time.time()
        mm = np.mean(d, axis=0)
        dtm = time.time() - tm
        print("mean : %.3f" %dtm)

        ts = time.time()
        ss = np.std(d, axis=0)
        dts = time.time() - ts
        print("std  : %.3f" %dts)
        print("\n")

        mt[ii] = dtm
        st[ii] = dts

    print("avg : %.2f +/- %.2f " %(np.mean(mt), np.std(mt)))
    print("std : %.2f +/- %.2f " %(np.mean(st), np.std(st)))

    return mt, st
