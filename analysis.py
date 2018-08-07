# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
from .thermometry import Sther
from scipy import interpolate
import numpy as np
import itertools
import operator
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def TS_Kxx(Data_P, Data_N, Temperatures, geofactor, col_I,
                 R_heater = 5000, **trash):

    """
    Function returns Thermal Conductivity Kxx in a dictionnary : \n

    ThermalConductivity["I_P"] = I_P # in Amps \n
    ThermalConductivity["I_N"] = I_N # in Amps \n
    ThermalConductivity["I"] = I # in Amps \n

    ThermalConductivity["Kxx_P"] = Kxx_P # W / K m \n
    ThermalConductivity["Kxx_N"] = Kxx_N # W / K m \n
    ThermalConductivity["Kxx"] = Kxx # W / K m \n


    Data_P : data matrix at positive magnetic field \n
    Data_N : data matrix at negative magnetic field \n
    Temperatures : dictionnary out of temp_thermocouple/temp_thermometers function \n
    geofactor : dictionnary with keys "L", "w", "t", the geometric factor \n
    col_I : column of current in the heater \n
    col_T0 : column of T0 \n
    R_heater : resistance of the sample heater, usely 5 kOhm \n
    **trash : is put in order to absorbe any unecessary key of dictionnary input \n

    """

    ## First extract columns #####################
    ##############################################

    ## Current

    I_P = Data_P[:,col_I] #  current at positive field in Amps
    I_N = Data_N[:,col_I] # current at negative field in Amps

    I = ( I_P + I_N ) / 2  # average current (supposedly the same)

    ## dTx

    dTx_P = Temperatures["dTx_P"] # dTx at positive field
    dTx_N = Temperatures["dTx_N"] # dTx at negative field
    dTx = Temperatures["dTx"] # dTx symmetrized

    ## Geometric factor

    L = geofactor["L"]  # in meter
    w = geofactor["w"]  # in meter
    t = geofactor["t"] # in meter

    alpha=(t*w)/L # in meter


    ## Compute Kxx on positive field #############
    ##############################################

    Kxx_P = R_heater * I_P**2 / ( dTx_P * alpha ) # W / K m
    Kxx_N = R_heater * I_N**2 / ( dTx_N * alpha ) # W / K m

    Kxx = R_heater * I**2 / ( dTx * alpha ) # W / K m

    ## Create dictionnary of thermal conductivity
    ##############################################

    ThermalConductivity = {}

    ThermalConductivity["I_P"] = I_P
    ThermalConductivity["I_N"] = I_N
    ThermalConductivity["I"] = I

    ThermalConductivity["Kxx_P"] = Kxx_P
    ThermalConductivity["Kxx_N"] = Kxx_N
    ThermalConductivity["Kxx"] = Kxx


    return ThermalConductivity


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def TS_Kxy(Data_P, Data_N, ThermalConductivity, Temperatures, geofactor,
           col_dTy_0, col_dTy_Q, T_connected = "Tav", sign_dTy = 1,
                     Type = "E", Gain=1000, **trash):

    """
    Function returns the ThermalHallConductivity Kxy for a dTy the transverse
    gradient of heat current Jqx, with H applied along the z-axis, in the dictionnary : \n

    ThermalHallConductivity["V_dTy_P"] = V_dTy_P # in Volts \n
    ThermalHallConductivity["V_dTy_N"] = V_dTy_N # in Volts \n
    ThermalHallConductivity["V_dTy"] = V_dTy # in Volts \n

    ThermalHallConductivity["dTy_P"] = dTy_P # in Kelvin \n
    ThermalHallConductivity["dTy_N"] = dTy_N # in Kelvin \n
    ThermalHallConductivity["dTy"] = dTy # in Kelvin \n

    ThermalHallConductivity["Kxy_P"] = Kxy_P # W / K m \n
    ThermalHallConductivity["Kxy_N"] = Kxy_N # W / K m \n
    ThermalHallConductivity["Kxy"] = Kxy # W / K m \n


    !!!
    T_connected is the temperature corresponding to the location on the sample
    where the differential thermocouple measuring dTy is connected : \n
    - if connected on T+ : then T_connected = "T+" \n
    - if connected on T- : then T_connected = "T-" \n
    - if connected in diagonal : then T_connected = "Tav" \n
    !!!

    Data_P : data matrix at positive magnetic field \n
    Data_N : data matrix at negative magnetic field \n
    ThermalConductivity : dictionnary out of Kxx_analysis function \n
    Temperatures : dictionnary out of temp_thermocouple/temp_thermometers function \n
    geofactor : dictionnary with keys "L", "w", "t", the geometric factor \n
    Gain : the gain of the preamp, for homemade ones in Sherbrooke Gain is 1000 \n
    col_dTy_0 : column of the differential thermocouple voltage without heat \n
    col_dTy_Q : column of the differential thermocouple voltage with heat \n
    sign_dTy : 1 if sign of dTy ok, -1 if needed to be reversed \n
    Type : type of the thermocouple, by default type "E" \n
    **trash : is put in order to absorbe any unecessary key of dictionnary input \n


    """

    ## First extract columns #####################
    ##############################################

    V_dTy_P_0 = Data_P[:,col_dTy_0] / Gain  # in Volt
    V_dTy_N_0 = Data_N[:,col_dTy_0] / Gain  # in Volt

    V_dTy_P_Q = Data_P[:,col_dTy_Q] / Gain  # in Volt
    V_dTy_N_Q = Data_N[:,col_dTy_Q] / Gain  # in Volt

    ## Compute V_dTy

    V_dTy_P = sign_dTy *( V_dTy_P_Q - V_dTy_P_0 )
    V_dTy_N = sign_dTy *( V_dTy_N_Q - V_dTy_N_0 )

    V_dTy = (V_dTy_P - V_dTy_N) / 2

    ## Geometric factor

    L = geofactor["L"]  # in meter
    w = geofactor["w"]  # in meter

    ## Load Kxx from ThermalConductivity #########
    ##############################################

    Kxx_P = ThermalConductivity["Kxx_P"]
    Kxx_N = ThermalConductivity["Kxx_N"]
    Kxx = ThermalConductivity["Kxx"]

    ## Choose T_connected from Temperatures ######
    ##############################################

    Tav_P = Temperatures["Tav_P"]
    Tav_N = Temperatures["Tav_N"]
    Tav = Temperatures["Tav"]

    dTx_P = Temperatures["dTx_P"]
    dTx_N = Temperatures["dTx_N"]
    dTx = Temperatures["dTx"]

    Tp_P = Tav_P + dTx_P / 2
    Tp_N = Tav_N + dTx_N / 2
    Tp = Tav + dTx / 2

    Tm_P = Tav_P - dTx_P / 2
    Tm_N = Tav_N - dTx_N / 2
    Tm = Tav - dTx / 2


    ## Compute dTy

    if T_connected == "T+":

        dTy_P = V_dTy_P / Sther(Tp_P, Type)
        dTy_N = V_dTy_N / Sther(Tp_N, Type)

        dTy = V_dTy / Sther(Tp, Type)

    elif T_connected == "T-":

        dTy_P = V_dTy_P / Sther(Tm_P, Type)
        dTy_N = V_dTy_N / Sther(Tm_N, Type)

        dTy = V_dTy / Sther(Tm, Type)

    else :

        dTy_P = V_dTy_P / Sther(Tav_P, Type)
        dTy_N = V_dTy_N / Sther(Tav_N, Type)

        dTy = V_dTy / Sther(Tav, Type)


    ## Compute Kxy

    Kxy_P = Kxx_P * dTy_P / dTx_P * L / w
    Kxy_N = Kxx_N * dTy_N / dTx_N * L / w
    Kxy = Kxx * dTy / dTx * L / w


    ## Create dictionnary of thermal Hall conductivity
    ##############################################

    ThermalHallConductivity = {}

    ThermalHallConductivity["V_dTy_P"] = V_dTy_P
    ThermalHallConductivity["V_dTy_N"] = V_dTy_N
    ThermalHallConductivity["V_dTy"] = V_dTy

    ThermalHallConductivity["dTy_P"] = dTy_P
    ThermalHallConductivity["dTy_N"] = dTy_N
    ThermalHallConductivity["dTy"] = dTy

    ThermalHallConductivity["Kxy_P"] = Kxy_P
    ThermalHallConductivity["Kxy_N"] = Kxy_N
    ThermalHallConductivity["Kxy"] = Kxy


    return ThermalHallConductivity

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def FS_Kxy_discrete(files_FS, columns_FS, geofactor, sign_dTy = 1,
                    T_connected = "Tav", columns_calib = None,
                    R_heater = 5000, Gain=1000, **trash):

    """
    Function computes Kxy function of mangetic field for a dTy the transverse\n
    gradient of heat current Jqx, with H applied along the z-axis.

    Input:
    file_FS : dictionnary of all isotherms, keys are (T0, date) and values \n
            [filepath, start_pause, end_pause]
    if using thermometers for dTx, there is another element to the list:
    ..., file_calib] : path + filename towards the data_raw file which will be used
    for calibrating the Cernox.

    columns_FS : dictionnary that contains all columns for the FS files, \n
        col_H, col_Vy, col_T0, col_I, col_Rp, col_Rm, col_Vy, col_pause\n
    geofactor : dictionnary with keys "L", "w", "t", the geometric factor \n
    Gain : the gain of the preamp, for homemade ones in Sherbrooke Gain is 1000 \n
    sign_dTy : 1 if sign of dTy ok, -1 if needed to be reversed \n

    columns_calib : dictionnary that contains all columns for the calib file, \n
                    col_T0_calib, col_Rp_0_calib, col_Rm_0_calib \n
    R_heater : resistance of the sample heater, usely 5 kOhm \n

    !!!
    T_connected is the temperature corresponding to the location on the sample
    where the differential thermocouple measuring dTy is connected : \n
    - if connected on T+ : then T_connected = "T+" \n
    - if connected on T- : then T_connected = "T-" \n
    - if connected in diagonal : then T_connected = "Tav" \n
    !!!

    Output
    It returns the following dictionnaries \n
    that take for key the title of there content (e.g. "dTy") and have for values \n
    a dictionnary which takes (Tav, date) as a key (e.g. SYM_dict["dTy"] = dTy_SYM_dict \n
    which is a dictionnary of keys Tav, date) \n
    NSYM : means non-symmetrized in field data
    SYM : means symmetrized or antisymmetrized in field data

    NSYM_dict["Vy_stab"] = Vy_NSYM_pause_dict \n
    # in Volt, it is Vy paused @ field values function time \n
    NSYM_dict["H"] = H_NSYM_dict \n
    NSYM_dict["Vy"] = Vy_NSYM_dict \n
    NSYM_dict["T0"] = T0_NSYM_dict \n
    NSYM_dict["Tp"] = Tp_NSYM_dict \n
    NSYM_dict["Tm"] = Tm_NSYM_dict \n

    SYM_dict["H"] = H_SYM_dict \n
    SYM_dict["T0"] = T0_SYM_dict \n
    SYM_dict["Vy"] = Vy_SYM_dict \n
    SYM_dict["Tp"] = Tp_SYM_dict \n
    SYM_dict["Tm"] = Tm_SYM_dict \n
    SYM_dict["Tav"] = Tav_SYM_dict \n
    SYM_dict["dTx"] = dTx_SYM_dict \n
    SYM_dict["Kxx"] = Kxx_SYM_dict \n
    SYM_dict["dTy"] = dTy_SYM_dict \n
    SYM_dict["Kxy"] = Kxy_SYM_dict \n
    SYM_dict["I"] = I_SYM_dict \n
    """

    ## Longitudinal thermometry: thermometers or thermocouple? >>>>>>>>>>>>>>>>#
    if columns_calib == None:
        thermometry_xx = "thermocouples"
    else:
        thermometry_xx = "thermometers"

    ## Columns FS /////////////////////////////////////////////////////////////#
    col_H = columns_FS["col_H"]
    col_T0 = columns_FS["col_T0"]
    col_I = columns_FS["col_I"]
    col_Rp = columns_FS["col_Rp"]
    col_Rm = columns_FS["col_Rm"]
    col_Vy = columns_FS["col_Vy"]
    col_pause = columns_FS["col_pause"]

    ## Geometric factor ///////////////////////////////////////////////////////#
    L = geofactor["L"]
    w = geofactor["w"]
    t = geofactor["t"]
    alpha = w * t / L

    ## Loading FS Data ////////////////////////////////////////////////////////#
    keys_files_FS = sorted(files_FS.keys())

    ## Initialize variables ///////////////////////////////////////////////////#
    NSYM_dict = {}
    SYM_dict = {}
    Tav_list = []

    H_NSYM_dict = {}
    T0_NSYM_dict = {}
    Tp_NSYM_dict = {}
    Tm_NSYM_dict = {}
    Vy_NSYM_dict = {}
    Vy_diagnostic_dict = {}
    H_diagnostic_dict = {}

    H_SYM_dict = {}
    T0_SYM_dict = {}
    Tp_SYM_dict = {}
    Tm_SYM_dict = {}
    Tav_SYM_dict = {}
    dTx_SYM_dict = {}
    Kxx_SYM_dict = {}
    Vy_SYM_dict = {}
    dTy_SYM_dict ={}
    Kxy_SYM_dict = {}
    I_SYM_dict = {}

    dTy_p_SYM_dict = {}
    dTy_m_SYM_dict = {}
    Kxy_p_SYM_dict = {}
    Kxy_m_SYM_dict = {}

    ## FS loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
    for (T0, date) in keys_files_FS:

        # Load data for one isotherm
        list_elements = files_FS[T0, date]
        filename = list_elements[0]
        start_pause = list_elements[1]
        end_pause = list_elements[2]

        Data = np.loadtxt('../data_raw/' + filename,
                          dtype = "float", comments = "#")

        H_NSYM_raw = Data[:,col_H] # in Tesla
        T0_NSYM_raw = Data[:,col_T0] # in Kelvin
        Vy_NSYM_raw = sign_dTy * ( Data[:,col_Vy] ) / Gain # in Volt
        I = np.mean(Data[:, col_I]) # in Amps
        pause = Data[:,col_pause] # column giving 1 if H is paused, 0 otherwise

        # Thermometry along x axis
        if thermometry_xx == "thermometers":

            ## Thermometers calibration ///////////////////////////////////////#
            file_calib = list_elements[3]
            data_calib = np.loadtxt("../data_raw/" + file_calib, dtype = "float", comments = "#")
            col_T0_calib = columns_calib["col_T0_calib"]
            col_Rp_0_calib = columns_calib["col_Rp_0_calib"]
            col_Rm_0_calib = columns_calib["col_Rm_0_calib"]
            # T0
            T0_calib = data_calib[:,col_T0_calib]
            # T+
            Rp_0 = data_calib[:,col_Rp_0_calib]
            # T-
            Rm_0 = data_calib[:,col_Rm_0_calib]

            coeff_Rp = np.polyfit( np.log(Rp_0) - np.mean(np.log(Rp_0)), np.log(T0_calib), 8)
            coeff_Rm = np.polyfit( np.log(Rm_0) - np.mean(np.log(Rm_0)), np.log(T0_calib), 8)

            Tp_NSYM_raw = np.exp( np.polyval( coeff_Rp, np.log(Data[:,col_Rp]) - np.mean(np.log(Rp_0))) )
            Tm_NSYM_raw = np.exp( np.polyval( coeff_Rm, np.log(Data[:,col_Rm]) - np.mean(np.log(Rm_0))) )


        ## Seperate setps in field ////////////////////////////////////////////#
        H_NSYM_pause = [[float(H_NSYM_raw[i]) for i,value in it] for key,it in itertools.groupby(enumerate(pause), key=operator.itemgetter(1)) if key != 0]
        T0_NSYM_pause = [[float(T0_NSYM_raw[i]) for i,value in it] for key,it in itertools.groupby(enumerate(pause), key=operator.itemgetter(1)) if key != 0]
        Vy_NSYM_pause = [[float(Vy_NSYM_raw[i]) for i,value in it] for key,it in itertools.groupby(enumerate(pause), key=operator.itemgetter(1)) if key != 0]
        Tp_NSYM_pause = [[float(Tp_NSYM_raw[i]) for i,value in it] for key,it in itertools.groupby(enumerate(pause), key=operator.itemgetter(1)) if key != 0]
        Tm_NSYM_pause = [[float(Tm_NSYM_raw[i]) for i,value in it] for key,it in itertools.groupby(enumerate(pause), key=operator.itemgetter(1)) if key != 0]

        H_NSYM_mean = np.array([np.mean(ith_list[start_pause:end_pause]) for ith_list in H_NSYM_pause])
        T0_NSYM_mean = np.array([np.mean(ith_list[start_pause:end_pause]) for ith_list in T0_NSYM_pause])
        Vy_NSYM_mean = np.array([np.mean(ith_list[start_pause:end_pause]) for ith_list in Vy_NSYM_pause])
        Tp_NSYM_mean = np.array([np.mean(ith_list[start_pause:end_pause]) for ith_list in Tp_NSYM_pause])
        Tm_NSYM_mean = np.array([np.mean(ith_list[start_pause:end_pause]) for ith_list in Tm_NSYM_pause])

        ## Sort increasingly regarding H //////////////////////////////////////#
        index_ordered = np.argsort(H_NSYM_mean)
        H_NSYM = H_NSYM_mean[index_ordered]
        Vy_NSYM = Vy_NSYM_mean[index_ordered]
        T0_NSYM = T0_NSYM_mean[index_ordered]
        Tp_NSYM = Tp_NSYM_mean[index_ordered]
        Tm_NSYM = Tm_NSYM_mean[index_ordered]

        ## Makes sure there is no repetition in field values //////////////////#
        H_NSYM, index_unique = np.unique(H_NSYM, return_index = True)
        Vy_NSYM = Vy_NSYM[index_unique]
        T0_NSYM = T0_NSYM[index_unique]
        Tp_NSYM = Tp_NSYM[index_unique]
        Tm_NSYM = Tm_NSYM[index_unique]

        ## Power //////////////////////////////////////////////////////////////#
        Q = R_heater * I**2

        ## SYMMETRIZATION /////////////////////////////////////////////////////#
        ## Round field values to get same values on both +/-H because sometimes
        ## one can get 15.0000001 != 15
        H_NSYM = np.round(H_NSYM, 3)

        index_P = H_NSYM > 0
        index_N = H_NSYM < 0

        ## Seperate H > 0 and H < 0
        H_P  = H_NSYM[index_P]
        H_N  = H_NSYM[index_N]
        T0_P = T0_NSYM[index_P]
        T0_N = T0_NSYM[index_N]
        Vy_P = Vy_NSYM[index_P]
        Vy_N = Vy_NSYM[index_N]
        Tp_P = Tp_NSYM[index_P]
        Tp_N = Tp_NSYM[index_N]
        Tm_P = Tm_NSYM[index_P]
        Tm_N = Tm_NSYM[index_N]

        # Interpolate H < 0 measurements on H > 0
        Vy_N_interp = interpolate.interp1d(H_N, Vy_N)( - H_P )
        Tp_N_interp = interpolate.interp1d(H_N, Tp_N)( - H_P )
        Tm_N_interp = interpolate.interp1d(H_N, Tm_N)( - H_P )
        T0_N_interp = interpolate.interp1d(H_N, T0_N)( - H_P )

        # Symmetrize measurements
        Vy_SYM = ( Vy_P - Vy_N_interp) / 2.
        Tp_SYM = ( Tp_P + Tp_N_interp) / 2.
        Tm_SYM = ( Tm_P + Tm_N_interp) / 2.

        T0_SYM = ( T0_P + T0_N_interp) / 2.
        Tav_SYM = (Tp_SYM + Tm_SYM) / 2.
        dTx_SYM =  Tp_SYM - Tm_SYM

        ## Static value of Tav ////////////////////////////////////////////////#
        Tav = np.round( np.mean(Tav_SYM), 2)
        Tav_list.append(Tav)

        print("T0 = " + r"{0:.2f}".format(T0) + "K -- > Tav = " + \
              r"{0:.2f}".format(Tav) + " K @ date = " + str(date))

        ## Compute dTy & Kxy from thermocouples ///////////////////////////////#
        if T_connected == "T+":
            dTy_SYM = Vy_SYM / Sther(Tp_SYM)
        elif T_connected == "T-":
            dTy_SYM = Vy_SYM / Sther(Tm_SYM)
        elif T_connected == "Tav":
            dTy_SYM = Vy_SYM / Sther(Tav_SYM)

        Kxx_SYM = Q / dTx_SYM / alpha
        Kxy_SYM = Kxx_SYM * dTy_SYM / dTx_SYM * ( L / w )

        ## Compute dTy & Kxy from thermometers ////////////////////////////////#
        if thermometry_xx == "thermometers":
            dTy_p_SYM = Tp_P - Tp_N # dTy computed from T+ antisym signal
            dTy_m_SYM = Tm_P - Tm_N # dTy computed from T- antisym signal
            Kxy_p_SYM = Kxx_SYM * dTy_p_SYM / dTx_SYM * ( L / w )
            Kxy_m_SYM = Kxx_SYM * dTy_m_SYM / dTx_SYM * ( L / w )


        ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## Diagnostic: keep only paused H points at each steps ////////////////#
        ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

        # Invert order of the isotherms that start from H > 0 (meaning the ones
        # that go from H to -H), in order to plot all isotherms from -H to H.
        if H_NSYM_raw[0] > 0:
            # Invert the order of the plateaus previously from H to -H, now from -H to H
            Vy_NSYM_pause = Vy_NSYM_pause[::-1]
            H_NSYM_pause = H_NSYM_pause[::-1]

        # Flatten the list of list in one list and then a numpy array
        Vy_diagnostic_list = Vy_NSYM_pause
        Vy_diagnostic = np.array(list(itertools.chain.from_iterable(Vy_diagnostic_list)))


        # Create a H array fit for diagnostic, meaning each plateau is centered
        # on H_mean value and extend from H_mean - step_H / 2 to H_mean + step_H / 2
        step_H = H_NSYM[1] - H_NSYM[0]
        H_diagnostic_list = [np.linspace(np.mean(ith_list)-step_H/2, np.mean(ith_list)+step_H/2, len(ith_list)) for ith_list in H_NSYM_pause]
        H_diagnostic = np.array(list(itertools.chain.from_iterable(H_diagnostic_list)))

        ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

        ## Place in dictionary for each Tav ///////////////////////////////////#
        H_NSYM_dict[Tav, date] = H_NSYM
        T0_NSYM_dict[Tav, date] = T0_NSYM
        Tp_NSYM_dict[Tav, date] = Tp_NSYM
        Tm_NSYM_dict[Tav, date] = Tm_NSYM
        Vy_NSYM_dict[Tav, date] = Vy_NSYM
        Vy_diagnostic_dict[Tav, date] = Vy_diagnostic
        H_diagnostic_dict[Tav, date] = H_diagnostic

        H_SYM_dict[Tav, date] = H_P
        T0_SYM_dict[Tav, date] = T0_SYM
        Tp_SYM_dict[Tav, date] = Tp_SYM
        Tm_SYM_dict[Tav, date] = Tm_SYM
        Tav_SYM_dict[Tav, date] = Tav_SYM
        dTx_SYM_dict[Tav, date] = dTx_SYM
        Kxx_SYM_dict[Tav, date] = Kxx_SYM
        Vy_SYM_dict[Tav, date] = Vy_SYM
        dTy_SYM_dict[Tav, date] = dTy_SYM
        Kxy_SYM_dict[Tav, date] = Kxy_SYM
        I_SYM_dict[Tav, date] = I * np.ones(len(H_P))

        dTy_p_SYM_dict[Tav, date] = dTy_p_SYM
        dTy_m_SYM_dict[Tav, date] = dTy_m_SYM
        Kxy_p_SYM_dict[Tav, date] = Kxy_p_SYM
        Kxy_m_SYM_dict[Tav, date] = Kxy_m_SYM

    ## Place in output dictionnaries all measurements /////////////////////////#
    NSYM_dict["H_diagnostic"] = H_diagnostic_dict
    NSYM_dict["Vy_diagnostic"] = Vy_diagnostic_dict
    NSYM_dict["H"] = H_NSYM_dict
    NSYM_dict["Vy"] = Vy_NSYM_dict
    NSYM_dict["T0"] = T0_NSYM_dict
    NSYM_dict["Tp"] = Tp_NSYM_dict
    NSYM_dict["Tm"] = Tm_NSYM_dict

    SYM_dict["H"] = H_SYM_dict
    SYM_dict["T0"] = T0_SYM_dict
    SYM_dict["Vy"] = Vy_SYM_dict
    SYM_dict["Tp"] = Tp_SYM_dict
    SYM_dict["Tm"] = Tm_SYM_dict
    SYM_dict["Tav"] = Tav_SYM_dict
    SYM_dict["dTx"] = dTx_SYM_dict
    SYM_dict["Kxx"] = Kxx_SYM_dict
    SYM_dict["dTy"] = dTy_SYM_dict
    SYM_dict["Kxy"] = Kxy_SYM_dict
    SYM_dict["I"] = I_SYM_dict

    SYM_dict["dTy_p"] = dTy_p_SYM_dict
    SYM_dict["dTy_m"] = dTy_m_SYM_dict
    SYM_dict["Kxy_p"] = Kxy_p_SYM_dict
    SYM_dict["Kxy_m"] = Kxy_m_SYM_dict

    return NSYM_dict, SYM_dict

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def TS_Seebeck(Data_P, Data_N, Temperatures, col_S_0, col_S_Q, sign_S = 1,
                    Gain = 1000, **trash):

    """
    Function returns the Seebeck coefficient in a dictionnary : \n

    Seebeck["V_S_P"] = V_S_P # in Volts \n
    Seebeck["V_S_N"] = V_S_N # in Volts \n
    Seebeck["V_S"] = V_S # in Volts \n


    Seebeck["S_P"] = S_P # in microvolt / K
    Seebeck["S_N"] = S_N # in microvolt / K
    Seebeck["S"] = S # in microvolt / K


    Data_P : data matrix at positive magnetic field \n
    Data_N : data matrix at negative magnetic field \n
    Temperatures : dictionnary out of temp_thermocouple/temp_thermometers function \n
    Gain : the gain of the preamp, for homemade ones in Sherbrooke Gain is 1000 \n
    col_S_0 : column of the Seebeck voltage without heat \n
    col_S_Q : column of the Seebeck voltage with heat \n
    sign_S : 1 if sign of dTy ok, -1 if needed to be reversed \n
    **trash : is put in order to absorbe any unecessary key of dictionnary input \n


    """

    ## First extract columns #####################
    ##############################################

    V_S_P_0 = Data_P[:,col_S_0] / Gain  # in Volt
    V_S_N_0 = Data_N[:,col_S_0] / Gain  # in Volt

    V_S_P_Q = Data_P[:,col_S_Q] / Gain  # in Volt
    V_S_N_Q = Data_N[:,col_S_Q] / Gain  # in Volt

    ## Compute V_S

    V_S_P = sign_S *( V_S_P_Q - V_S_P_0 )
    V_S_N = sign_S *( V_S_N_Q - V_S_N_0 )

    V_S = (V_S_P + V_S_N) / 2

    ## Choose T_connected from Temperatures ######
    ##############################################

    dTx_P = Temperatures["dTx_P"] # in Kelvin
    dTx_N = Temperatures["dTx_N"] # in Kelvin
    dTx = Temperatures["dTx"] # in Kelvin

    ## Compute S #################################
    ##############################################

    S_P = V_S_P / dTx_P * 1e6 # in microvolt / K
    S_N = V_S_N / dTx_N * 1e6 # in microvolt / K
    S = V_S / dTx * 1e6 # in microvolt / K

    ## Create dictionnary of Seebeck #############
    ##############################################

    Seebeck = {}

    Seebeck["V_S_P"] = V_S_P
    Seebeck["V_S_N"] = V_S_N
    Seebeck["V_S"] = V_S

    Seebeck["S_P"] = S_P
    Seebeck["S_N"] = S_N
    Seebeck["S"] = S


    return Seebeck


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def TS_Nernst(Data_P, Data_N, Temperatures, geofactor, Field, col_N_0, col_N_Q,
              sign_N = 1, Gain=1000, **trash):

    """
    Function returns the Nernst coefficient in a dictionnary : \n

    Nernst["V_N_P"] = V_N_P  # volts
    Nernst["V_N_N"] = V_N_N # volts
    Nernst["V_N"] = V_N # volts

    Nernst["N_P"] = N_P # in microvolt / K
    Nernst["N_N"] = N_N # in microvolt / K
    Nernst["N"] = N # in microvolt / K

    Nernst["nu_P"] = nu_P # in microvolt / K T
    Nernst["nu_N"] = nu_N # in microvolt / K T
    Nernst["nu"] = nu # in microvolt / K T


    Data_P : data matrix at positive magnetic field \n
    Data_N : data matrix at negative magnetic field \n
    Temperatures : dictionnary out of temp_thermocouple/temp_thermometers function \n
    Gain : the gain of the preamp, for homemade ones in Sherbrooke Gain is 1000 \n
    col_N_0 : column of the Nernst voltage without heat \n
    col_N_Q : column of the Nernst voltage with heat \n
    sign_N : 1 if sign of dTy ok, -1 if needed to be reversed \n
    **trash : is put in order to absorbe any unecessary key of dictionnary input \n


    """

    ## First extract columns #####################
    ##############################################

    V_N_P_0 = Data_P[:,col_N_0] / Gain  # in Volt
    V_N_N_0 = Data_N[:,col_N_0] / Gain  # in Volt

    V_N_P_Q = Data_P[:,col_N_Q] / Gain  # in Volt
    V_N_N_Q = Data_N[:,col_N_Q] / Gain  # in Volt

    ## Compute V_N

    V_N_P = sign_N *( V_N_P_Q - V_N_P_0 )
    V_N_N = sign_N *( V_N_N_Q - V_N_N_0 )

    V_N = (V_N_P - V_N_N) / 2

    ## Geometric factor

    L = geofactor["L"]  # in meter
    w = geofactor["w"]  # in meter

    ## Choose T_connected from Temperatures ######
    ##############################################

    dTx_P = Temperatures["dTx_P"] # in Kelvin
    dTx_N = Temperatures["dTx_N"] # in Kelvin
    dTx = Temperatures["dTx"] # in Kelvin

    ## Compute N #################################
    ##############################################

    N_P = V_N_P / dTx_P * ( L / w ) *1e6 # in microvolt / K
    N_N = V_N_N / dTx_N * ( L / w ) *1e6 # in microvolt / K
    N = V_N / dTx * ( L / w ) *1e6 # in microvolt / K

    ## Compute nu ################################
    ##############################################

    nu_P = N_P / Field # in microvolt / K T
    nu_N = N_N / Field # in microvolt / K T
    nu = N / Field # in microvolt / K T

    ## Create dictionnary of Nernst #############
    ##############################################

    Nernst = {}

    Nernst["V_N_P"] = V_N_P  # volts
    Nernst["V_N_N"] = V_N_N # volts
    Nernst["V_N"] = V_N # volts

    Nernst["N_P"] = N_P # in microvolt / K
    Nernst["N_N"] = N_N # in microvolt / K
    Nernst["N"] = N # in microvolt / K

    Nernst["nu_P"] = nu_P # in microvolt / K T
    Nernst["nu_N"] = nu_N # in microvolt / K T
    Nernst["nu"] = nu # in microvolt / K T


    return Nernst


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def TS_CUT_discrete(SYM_dict, H_cut):

    """
    Function returns a cut of isotherms for a precise H_cut \n
    Variables in: \n
    - dictionnary that take for key the title of there content (e.g. "dTy") and have for values \n
    a dictionnary which takes (Tav, date) as a key (e.g. SYM_dict["dTy"] = dTy_SYM_dict \n
    which is a dictionnary of keys Tav, date) \n
    SYM : means symmetrized or antisymmetrized in field data
    - H_cut : the H value where to get the cut of the isotherm
    """

    TS_CUT_dict = {}

    # Initialize the dictionnaries values to list
    for label in SYM_dict.keys():
        TS_CUT_dict[label] = []
    TS_CUT_dict["date_list"] = []

    # Exactract the isotherm value at H_cut
    keys = sorted(SYM_dict["H"])
    for Tav, date in keys:

        H = SYM_dict["H"][Tav, date]
        index_cut = H == H_cut

        for label in SYM_dict.keys():
            TS_CUT_dict[label].append(float(SYM_dict[label][Tav, date][index_cut]))

        TS_CUT_dict["date_list"].append(date)

    # Convert list to numpy array
    for label in SYM_dict.keys():
        TS_CUT_dict[label] = np.array(TS_CUT_dict[label])
    TS_CUT_dict["date_list"] = np.array(TS_CUT_dict["date_list"])


    return TS_CUT_dict