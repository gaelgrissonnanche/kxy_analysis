# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import numpy as np
import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def load_TS_file(folder, position, date, samplelabel = "", H = 0, measurement = "SYM", **trash):

    """
    Loads the TSweep files from my own VIs for the VTI, the data folder has to
    be put in 'data_raw'
    The argument 'bothH' if put to False, means that we are loading only a non-zero
    H file without wanting to load the opposite H file.
    If 'bothH' is True and H != 0, the function will try to load +H and -H
    data files, and will return a tuple (data_P, data_N)
    """

    if H != 0 and measurement[0:3] == "SYM":

        data_P = np.loadtxt("../data_raw/{0}/Data-TS-{1}T-{2}-{3}-{4}.dat".format(
                    folder, "{0:.1f}".format(abs(H)) , position,
                    samplelabel, date), dtype = "float", comments = "#")

        data_N = np.loadtxt("../data_raw/{0}/Data-TS--{1}T-{2}-{3}-{4}.dat".format(
                    folder, "{0:.1f}".format(abs(H)) , position,
                    samplelabel, date), dtype = "float", comments = "#")

    else:

        data = np.loadtxt("../data_raw/{0}/Data-TS-{1}T-{2}-{3}-{4}.dat".format(
                    folder, "{0:.1f}".format(H) , position,
                    samplelabel, date), dtype = "float", comments = "#")

        data_P = data
        data_N = data

    return (data_P, data_N)


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def load_pickle_geofactor(samplelabel):

    """
    Load the geofactor pickle for the corresponding samplelabel. \n
    It returns the geometric factor in a dictionnary of keys :\n
    "L" : the length between contacts on the sample, in meter \n
    "w" : the width of the sample, in meter \n
    "t" : the thicknes of the sample, in meter \n
    """

    with open("geofactor_" + samplelabel + ".pkl", "rb") as f: # 'r' means read

        unpickler_file = pickle.Unpickler(f)
        geofactor = unpickler_file.load() # keys of geofactor are "L", "w", "t"

    return geofactor


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def load_pickle_info(samplelabel):

    """
    Load the measurement pickle for the corresponding samplelabel. It contains
    the info regarding the measurements. The keys of info are (Field, Date), and
    return for each couple of (Field,Date) a particular dictionnary for each
    measurements. This dictionnary has for keys : \n

    "folder" : the name of the folder where data are stored in "data_raw" \n
    "position" : the position of the sample on the mount \n
    "comments" : comments regarding this particular set of points \n
    """

    with open("info_TS_" + samplelabel + ".pkl", "rb") as f: # 'w' means (over)write or create if doesn't exist

        unpickler_file = pickle.Unpickler(f)
        info = unpickler_file.load() # keys of measurements are (H,date)

    return info

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def load_pickle_FS(samplelabel):

    """
    Load the measurement pickle for the corresponding samplelabel. It contains
    a dictionnary file_FS of all isotherms, keys are (T0, date) and values
    [filepath, start_pause, end_pause]
    """

    with open("info_FS_" + samplelabel + ".pkl", "rb") as f: # 'w' means (over)write or create if doesn't exist

        unpickler_file = pickle.Unpickler(f)
        info = unpickler_file.load() # keys of measurements are (H,date)

    return info


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def save_TS_file(samplelabel, axe, TS_style, measurement, H, date, Data, Header):

  file_name = samplelabel + "_" + axe + "_" + TS_style + "_" + measurement + \
             "_" + str(H) + 'T_' + str(date) +".dat"

  np.savetxt('../data_analyzed/' + file_name,
             Data, fmt='%.7e',
             header = Header,
             comments = "#")

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def save_FS_files(samplelabel, axe, measurement, SYM_dict, data_type_list, Header):

    keys = sorted(SYM_dict["H"].keys())

    for Tav, date in keys:

        Data_tuple = ()
        for data_type in data_type_list:
            Data_tuple += (SYM_dict[data_type][Tav, date],)

        Data = np.vstack(Data_tuple)
        Data = Data.transpose()

        file_name =  samplelabel + "_" + axe + "_FS_" + measurement + "_" + \
                 "{0:g}".format(Tav) + 'K_' + str(date) + ".dat"

        np.savetxt('../Data_analyzed/'+ file_name,
                    Data, fmt='%.7e', header = Header, comments = "#")