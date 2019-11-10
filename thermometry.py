# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def Sther(T_Kelvin, Type = "E"):

    """
    This function returns the Seebeck coefficient of the thermocouple
    concerned (by default type "E") at a certain temperature. The input of the
    function is a temperature in Kelvin, but the coefficient below are for a
    polynomial function with T in Celsius. The output is S in [V / K]
    """


    if type(T_Kelvin) not in (np.ndarray, list) :

        T_Kelvin = np.array([T_Kelvin]) # convert in numpy array if T_Kevlin just a scalar

    elif  type(T_Kelvin) is list :

        T_Kelvin = np.array(T_Kelvin) # convert in numpy array if T_Kevlin just a scalar

    else :

        T_Kelvin = T_Kelvin



    ## Load coeff of thermocouples for E = f(T_Celsius) where E is in microVolts

    if Type == "E":

        coeff_E_below_270K = np.array([0,
                                      5.8665508708E1,
                                      4.5410977124E-2,
                                      -7.7998048686E-4,
                                      -2.5800160843E-5,
                                      -5.9452583057E-7,
                                      -9.3214058667E-9,
                                      -1.0287605534E-10,
                                      -8.0370123621E-13,
                                      -4.3979497391E-15,
                                      -1.6414776355E-17,
                                      -3.9673619516E-20,
                                      -5.5827328721E-23,
                                      -3.4657842013E-26])

        coeff_E_below_270K = coeff_E_below_270K[::-1] # ->reverse,
                                    # it is necessary to have cn, cn-1,...
                                    # for using poly1d and polyder


        coeff_E_above_270K = np.array([0,
                                      5.8665508710E1,
                                      4.5032275582E-2,
                                      2.8908407212E-5,
                                      -3.3056896652E-7,
                                      6.5024403270E-10,
                                      -1.9197495504E-13,
                                      -1.2536600497E-15,
                                      2.1489217569E-18,
                                      -1.4388041782E-21,
                                      3.5960899481E-25])

        coeff_E_above_270K = coeff_E_above_270K[::-1] # ->reverse,
                                    # it is necessary to have cn, cn-1,...
                                    # for using poly1d and polyder



    ## Convert T_Kelvin to Celsius

    T_Celsius = T_Kelvin - 273.15

    ## Selection of coefficient for temperature regime

    index_below = np.where(T_Celsius <= 0) # np.where returns the index of the condition
    index_above = np.where(T_Celsius > 0) # np.where returns the index of the condition

    S_values = np.zeros(np.size(T_Kelvin))


    E_below = np.poly1d(coeff_E_below_270K) # is a poly1d object in microVolt
    S_below = np.polyder(E_below) # is a poly1d object in microVolt / Celsius
    S_values[index_below] = S_below(T_Celsius[index_below])*1e-6 # is in Volt / K


    E_above = np.poly1d(coeff_E_above_270K) # is a poly1d object in microVolt
    S_above = np.polyder(E_above) # is a poly1d object in microVolt / Celsius
    S_values[index_above] = S_above(T_Celsius[index_above])*1e-6 # is in Volt / K


    return S_values


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def temp_thermocouple(Data_P, Data_N, col_T0, col_dTabs_0, col_dTabs_Q,
                      col_dTx_0, col_dTx_Q, Type = "E", Gain=1000, **trash):

    """
    Function returns the average temperature Tav of the sample,
    the longitudinal thermal gradient dTx, and dTy the transverse gradient if
    a magnetic field is applied. \n \n

    !!! The absolute thermocouple must be connected at T- !!! \n

    Data_P : data matrix at positive magnetic field \n
    Data_N : data matrix at negative magnetic field \n
    Gain : the gain of the preamp, for homemade ones in Sherbrooke Gain is 1000 \n
    col_T0 : column of temperature of the probe T0 \n
    col_dTabs_0 : column of the absolute thermocouple voltage without heat \n
    col_dTabs_Q : column of the absolute thermocouple voltage with heat \n
    col_dTx_0 : column of the differential thermocouple voltage without heat \n
    col_dTx_Q : column of the differential thermocouple voltage with heat \n
    Type : type of the thermocouple, by default type "E" \n
    **trash : is put in order to absorbe any unecessary key of dictionnary input \n

    The function returns the dictionnary "Temperatures" where : \n

    Temperatures["T0_P"] = T0_P \n
    Temperatures["T0_N"] = T0_N \n
    Temperatures["T0"] = T0 \n

    Temperatures["Tav_P"] = Tav_P \n
    Temperatures["Tav_N"] = Tav_N \n
    Temperatures["Tav"] = Tav \n

    Temperatures["dTx_P"] = dTx_P \n
    Temperatures["dTx_N"] = dTx_N \n
    Temperatures["dTx"] = dTx \n

    Temperatures["dTabs_P"] = dTabs_P \n
    Temperatures["dTabs_N"] = dTabs_N \n
    Temperatures["dTabs"] = dTabs \n

    Temperatures["dTy"] = dTy \n


    --> 'dTy' is here only the antisymmetrisation of dTabs, equivalent to the
    antisym of T-, it is not a differential value

    The function also returns the data_P, data_N matrix with the same number
    of lines, which is not always the case between H+ and H- measurements.
    """

    #!!! First make sure Data_P & Data_N have same number of T0 !!!#

    T0_P = Data_P[:,col_T0] # in Kelvin
    T0_N = Data_N[:,col_T0] # in Kelvin

    Data_P = Data_P[:min(len(T0_P), len(T0_N)), :]
    Data_N = Data_N[:min(len(T0_P), len(T0_N)), :]

    ## Extract columns ##

    T0_P = Data_P[:,col_T0] # in Kelvin
    T0_N = Data_N[:,col_T0] # in Kelvin

    T0 = (T0_P + T0_N) / 2 # the SYM_T0 (which should be the same)

    V_dTabs_P_0 = Data_P[:,col_dTabs_0] / Gain  # in Volt
    V_dTabs_N_0 = Data_N[:,col_dTabs_0] / Gain  # in Volt

    V_dTabs_P_Q = Data_P[:,col_dTabs_Q] / Gain  # in Volt
    V_dTabs_N_Q = Data_N[:,col_dTabs_Q] / Gain  # in Volt

    V_dTx_P_0 = Data_P[:,col_dTx_0] / Gain  # in Volt
    V_dTx_N_0 = Data_N[:,col_dTx_0] / Gain  # in Volt

    V_dTx_P_Q = Data_P[:,col_dTx_Q] / Gain  # in Volt
    V_dTx_N_Q = Data_N[:,col_dTx_Q] / Gain  # in Volt

    ## Compute dTabs

    dTabs_P = np.abs(V_dTabs_P_Q - V_dTabs_P_0) / Sther(T0_P, Type)
    dTabs_N = np.abs(V_dTabs_N_Q - V_dTabs_N_0) / Sther(T0_N, Type)

    error_P = np.ones_like(T0)
    error_N = np.ones_like(T0)
    while np.max(error_P) > 1e-8 and np.max(error_N) > 1e-8:

        dTabs_P_old = dTabs_P
        dTabs_N_old = dTabs_N

        dTabs_P = np.abs(V_dTabs_P_Q - V_dTabs_P_0) / Sther(T0_P + dTabs_P/2, Type)
        dTabs_N = np.abs(V_dTabs_N_Q - V_dTabs_N_0) / Sther(T0_N + dTabs_N/2, Type)

        error_P = np.abs(dTabs_P_old - dTabs_P)
        error_N = np.abs(dTabs_N_old - dTabs_N)



    ## Compute dTx

    dTx_P = abs(V_dTx_P_Q - V_dTx_P_0) / Sther(T0_P + dTabs_P, Type)
    dTx_N = abs(V_dTx_N_Q - V_dTx_N_0) / Sther(T0_N + dTabs_N, Type)

    error_P = np.ones_like(T0)
    error_N = np.ones_like(T0)
    while np.max(error_P) > 1e-8 and np.max(error_N) > 1e-8:

        dTx_P_old = dTx_P
        dTx_N_old = dTx_N

        dTx_P = np.abs(V_dTx_P_Q - V_dTx_P_0) / Sther(T0_P + dTabs_P + dTx_P/2, Type)
        dTx_N = np.abs(V_dTx_N_Q - V_dTx_N_0) / Sther(T0_N + dTabs_N + dTx_N/2, Type)

        error_P = np.abs(dTx_P_old - dTx_P)
        error_N = np.abs(dTx_N_old - dTx_N)


    ## Compute dTx, Tp, Tm, Tav, dTy

    Tav_P = T0_P + dTabs_P + dTx_P /2
    Tav_N = T0_N + dTabs_N + dTx_N/2

    Tp_P = Tav_P + dTx_P / 2
    Tp_N = Tav_N + dTx_N / 2
    Tp = ( Tp_P + Tp_N ) / 2

    Tm_P = Tav_P - dTx_P / 2
    Tm_N = Tav_N - dTx_N / 2
    Tm = ( Tm_P + Tm_N ) / 2

    dTabs = ( dTabs_P + dTabs_N ) / 2
    dTx = ( dTx_P + dTx_N ) / 2
    Tav = ( Tav_P + Tav_N ) / 2

    dTy = dTabs_P - dTabs_N

    ## Dictionnary output

    Temperatures = {}

    Temperatures["T0_P"] = T0_P
    Temperatures["T0_N"] = T0_N
    Temperatures["T0"] = T0

    Temperatures["Tav_P"] = Tav_P
    Temperatures["Tav_N"] = Tav_N
    Temperatures["Tav"] = Tav

    Temperatures["dTx_P"] = dTx_P
    Temperatures["dTx_N"] = dTx_N
    Temperatures["dTx"] = dTx

    Temperatures["dTabs_P"] = dTabs_P
    Temperatures["dTabs_N"] = dTabs_N
    Temperatures["dTabs"] = dTabs

    Temperatures["Tav_P"] = Tav_P
    Temperatures["Tav_N"] = Tav_N
    Temperatures["Tav"] = Tav

    Temperatures["dTx_P"] = dTx_P
    Temperatures["dTx_N"] = dTx_N
    Temperatures["dTx"] = dTx

    Temperatures["Tp_P"] = Tp_P
    Temperatures["Tp_N"] = Tp_N
    Temperatures["Tp"] = Tp

    Temperatures["Tm_P"] = Tm_P
    Temperatures["Tm_N"] = Tm_N
    Temperatures["Tm"] = Tm

    Temperatures["dTy"] = dTy

    return Temperatures, Data_P, Data_N


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def temp_thermometers(Data_P, Data_N, col_T0, col_Rp_0, col_Rp_Q,
        col_Rm_0, col_Rm_Q, deg_fit = 8, cut_off = True, graph = True, **trash):

    """
    Function returns the average temperature Tav of the sample,
    the longitudinal thermal gradient dTx, and dTy the transverse gradient if
    a magnetic field is applied. \n \n

    Data_P : data matrix at positive magnetic field \n
    Data_N : data matrix at negative magnetic field \n
    col_T0 : column of temperature of the probe T0 \n
    col_Rp_0 : column of the cernox resistance T+ without heat \n
    col_Rp_Q : column of the cernox resistance T+ with heat \n
    col_Rm_0 : column of the cernox resistance T- without heat \n
    col_Rm_Q : column of the cernox resistance T- with heat \n
    deg_fit : degree of fit for calibration of thermometers \n
    cut_off : if True, removes the points outside of the calibration range
    graph : True of display graphic calibration of thermometers \n
    **trash : is put in order to absorbe any unecessary key of dictionnary input \n

    The function returns the dictionnary "Temperatures" where : \n

    Temperatures["T0_P"] = T0_P \n
    Temperatures["T0_N"] = T0_N \n
    Temperatures["T0"] = T0 \n

    Temperatures["Tav_P"] = Tav_P \n
    Temperatures["Tav_N"] = Tav_N \n
    Temperatures["Tav"] = Tav \n

    Temperatures["dTx_P"] = dTx_P \n
    Temperatures["dTx_N"] = dTx_N \n
    Temperatures["dTx"] = dTx \n

    Temperatures["Tp_P"] = Tp_P \n
    Temperatures["Tp_N"] = Tp_N \n
    Temperatures["Tp"] = Tp \n

    Temperatures["Tm_P"] = Tp_P \n
    Temperatures["Tm_N"] = Tp_N \n
    Temperatures["Tm"] = Tp_N \n

    Temperatures["dTy_p"] = dTy_p \n
    Temperatures["dTy_m"] = dTy_m \n


    --> 'dTy' is here only the antisymmetrisation of dTabs, equivalent to the
    antisym of T-, it is not a differential value

    The function also returns the data_P, data_N matrix filtered from data
    points with Tp and Tm outside of the calibration range of the thermometers.
    Besides, data_P, data_N have now for certain the same number of lines, which
    is not always the case between H+ and H- measurements.
    """

    #!!! First make sure Data_P & Data_N have same number of T0 !!!#

    T0_P = Data_P[:,col_T0] # in Kelvin
    T0_N = Data_N[:,col_T0] # in Kelvin

    Data_P = Data_P[:min(len(T0_P), len(T0_N)), :]
    Data_N = Data_N[:min(len(T0_P), len(T0_N)), :]

    ## Extract columns

    T0_P = Data_P[:,col_T0] # in Kelvin
    t0_P = Data_P[:,col_T0] # in Kelvin for the plot of calibration
    T0_N = Data_N[:,col_T0] # in Kelvin
    t0_N = Data_N[:,col_T0] # in Kelvin for the plot of calibration

    # T+
    Rp_P_0 = Data_P[:,col_Rp_0]  # R+(H+,0)
    Rp_N_0 = Data_N[:,col_Rp_0]  # R+(H-,0)

    Rp_P_Q = Data_P[:,col_Rp_Q]  # R+(H+,Q)
    Rp_N_Q = Data_N[:,col_Rp_Q]  # R+(H-,Q)

    # T-
    Rm_P_0 = Data_P[:,col_Rm_0]  # R-(H+,0)
    Rm_N_0 = Data_N[:,col_Rm_0]  # R-(H-,0)

    Rm_P_Q = Data_P[:,col_Rm_Q]  # R-(H+,Q)
    Rm_N_Q = Data_N[:,col_Rm_Q]  # R-(H-,Q)

    #----------------------------------#
    ########### Positive H #############
    #----------------------------------#

    coeff_Rp_P = np.polyfit( np.log(Rp_P_0) - np.mean(np.log(Rp_P_0)), np.log(T0_P), deg_fit)
    coeff_Rm_P = np.polyfit( np.log(Rm_P_0) - np.mean(np.log(Rm_P_0)), np.log(T0_P), deg_fit)

    Tp_P = np.exp( np.polyval( coeff_Rp_P, np.log(Rp_P_Q) - np.mean(np.log(Rp_P_0))) ) # hot temperature T+
    Tm_P = np.exp( np.polyval( coeff_Rm_P, np.log(Rm_P_Q) - np.mean(np.log(Rm_P_0))) ) # cold temperature T-

    dTx_P = Tp_P - Tm_P # dTx(H+)
    Tav_P = (Tp_P + Tm_P) / 2. # Tav(H+)

    #----------------------------------#
    ########### Negative H #############
    #----------------------------------#

    coeff_Rp_N = np.polyfit( np.log(Rp_N_0) - np.mean(np.log(Rp_N_0)), np.log(T0_N), deg_fit)
    coeff_Rm_N = np.polyfit( np.log(Rm_N_0) - np.mean(np.log(Rm_N_0)), np.log(T0_N), deg_fit)

    Tp_N = np.exp( np.polyval( coeff_Rp_N, np.log(Rp_N_Q) - np.mean(np.log(Rp_N_0))) ) # hot temperature T+
    Tm_N = np.exp( np.polyval( coeff_Rm_N, np.log(Rm_N_Q) - np.mean(np.log(Rm_N_0))) ) # cold temperature T-

    dTx_N = Tp_N - Tm_N # dTx(H-)
    Tav_N = (Tp_N + Tm_N) / 2. # Tav(H-)

    if cut_off == True:

        #----------------------------------#
        ## To be in the calibration range ##
        #----------------------------------#

        # Keep Tp and Tm only in the range of calibration

        index_calib_P = ( Tp_P <= max(T0_P) ) * ( Tm_P <= max(T0_P) )
        index_calib_N = ( Tp_N <= max(T0_N) ) * ( Tm_N <= max(T0_N) )

        # To be sure _P and _N have the same size in the end, let's take:

        index_calib = index_calib_P * index_calib_N

        # Now cut the vectors to the right size

        T0_P = T0_P[index_calib]
        Tp_P = Tp_P[index_calib]
        Tm_P = Tm_P[index_calib]
        dTx_P = dTx_P[index_calib]
        Tav_P = Tav_P[index_calib]

        T0_N = T0_N[index_calib]
        Tp_N = Tp_N[index_calib]
        Tm_N = Tm_N[index_calib]
        dTx_N = dTx_N[index_calib]
        Tav_N = Tav_N[index_calib]

        Data_P = Data_P[index_calib,:]
        Data_N = Data_N[index_calib,:]

    else:

        #! If one does not want to remove the points outside of the calibration
        #! but make sure Data_P & Data_N have same number of T0

        Data_P = Data_P[:min(len(T0_P), len(T0_N)), :]
        Data_N = Data_N[:min(len(T0_P), len(T0_N)), :]

    #----------------------------------#
    ########### Symmetrize #############
    #----------------------------------#

    T0 = (T0_P + T0_N) / 2. # the SYM_T0 (which should be the same)

    Tp = ( Tp_P + Tp_N ) / 2.
    Tm = ( Tm_P + Tm_N ) / 2.

    dTx = ( dTx_P + dTx_N ) / 2.
    Tav = ( Tav_P + Tav_N ) / 2.

    dTy_p = ( Tp_P - Tp_N ) / 2. # dTy (T+)
    dTy_m = ( Tm_P - Tm_N ) / 2. # dTy (T-)

    ## Dictionnary output

    Temperatures = {}

    Temperatures["T0_P"] = T0_P
    Temperatures["T0_N"] = T0_N
    Temperatures["T0"] = T0

    Temperatures["Tav_P"] = Tav_P
    Temperatures["Tav_N"] = Tav_N
    Temperatures["Tav"] = Tav

    Temperatures["dTx_P"] = dTx_P
    Temperatures["dTx_N"] = dTx_N
    Temperatures["dTx"] = dTx

    Temperatures["Tp_P"] = Tp_P
    Temperatures["Tp_N"] = Tp_N
    Temperatures["Tp"] = Tp

    Temperatures["Tm_P"] = Tm_P
    Temperatures["Tm_N"] = Tm_N
    Temperatures["Tm"] = Tm

    Temperatures["dTy_p"] = dTy_p
    Temperatures["dTy_m"] = dTy_m

    #----------------------------------#
    ######## Plot calibration ##########

    if graph == True :

        mpl.rcdefaults()
        mpl.rcParams['font.size'] = 18. # change the size of the font in every figure
        mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
        mpl.rcParams['axes.labelsize'] = 18.
        mpl.rcParams['xtick.labelsize'] = 18
        mpl.rcParams['ytick.labelsize'] = 18
        mpl.rcParams['xtick.direction'] = "in"
        mpl.rcParams['ytick.direction'] = "in"
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True
        mpl.rcParams['xtick.major.width'] = 0.6
        mpl.rcParams['ytick.major.width'] = 0.6
        mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
        mpl.rcParams['pdf.fonttype'] = 42  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                            # editing the text in illustrator

        fig, axes = plt.subplots(2, 2, figsize = (11,8)) # figsize is w x h in inch of figure
        fig.subplots_adjust(left=0.02, right=0.96, bottom=0.05, top=0.9, wspace=0.3, hspace=0) # adjust the box of axes regarding the figure size


        ax1 = plt.subplot2grid((3,2), (0,0))
        ax2 = plt.subplot2grid((3,2), (1,0), rowspan=2)
        ax3 = plt.subplot2grid((3,2), (0,1))
        ax4 = plt.subplot2grid((3,2), (1,1), rowspan=2)

        bbox_props = dict(boxstyle="circle", fc="w", edgecolor="k", alpha=1, lw = 2)
        fig.text(0.44, 0.92, r"$H$ > 0", ha="center", va="center", size=15,
        bbox=bbox_props, color = "k")

        bbox_props = dict(boxstyle="circle", fc="w", edgecolor="k", alpha=1, lw = 2)
        fig.text(0.85, 0.92, r"$H$ < 0", ha="center", va="center", size=15,
        bbox=bbox_props, color = "k")

        #----------------------------------#
        ########### Positive H #############
        #----------------------------------#

        T0_P = t0_P

        ############## T0_fit ##############
        ax1.plot(np.log(Rp_P_0), np.log(T0_P), marker="o", ls="", ms=6, mew= 0, mfc = "r", mec = "r")
        ax1.plot(np.log(Rm_P_0), np.log(T0_P), marker="o", ls="", ms=6, mew= 0, mfc = "#440154", mec = "#440154")
        ax1.plot(np.log(Rp_P_0), np.polyval( coeff_Rp_P, np.log(Rp_P_0) - np.mean(np.log(Rp_P_0))), ls="-", lw=1.5,  c="k")
        ax1.plot(np.log(Rm_P_0), np.polyval( coeff_Rm_P, np.log(Rm_P_0) - np.mean(np.log(Rm_P_0))), ls="-", lw=1.5,  c="k")
        ax1.set_ylabel(r"log ( $T_0$ )", labelpad = 15) # alpha is the transparancy
        ax1.locator_params(nbins=10) #gives the number of ticks maximum we want, here 6
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])

        ############ T0 - T0_fit ###########
        ax2.axhline(y=0, ls ="--", c ="k")
        ax2.plot(T0_P,(T0_P-np.exp(np.polyval( coeff_Rp_P, np.log(Rp_P_0) - np.mean(np.log(Rp_P_0)))))*1e3,
                        marker="o", ls="--", lw=3, c="#DBDBDB", ms=8, mew = 0, mec="r", mfc = "r")
        ax2.plot(T0_P,(T0_P-np.exp(np.polyval( coeff_Rm_P, np.log(Rm_P_0) - np.mean(np.log(Rm_P_0)))))*1e3,
                        marker="o", ls="--", lw=3, c="#DBDBDB", ms=8, mew = 0, mec = "#440154", mfc = "#440154")
        ax2.set_xlabel(r"$T_0$ ( K )")
        ax2.set_ylabel(r"$\Delta T_0$ (fit) ( mK )")
        ax2.locator_params(nbins=6)


        #----------------------------------#
        ########### Negative H #############
        #----------------------------------#

        T0_N = t0_N

        ############## T0_fit ##############
        ax3.plot(np.log(Rp_N_0), np.log(T0_N), marker="o", ls="", ms=6, mfc = "r", mec = "r")
        ax3.plot(np.log(Rm_N_0), np.log(T0_N), marker="o", ls="", ms=6, mfc = "#440154", mec = "#440154")
        ax3.plot(np.log(Rp_N_0), np.polyval( coeff_Rp_N, np.log(Rp_N_0) - np.mean(np.log(Rp_N_0))), ls="-", lw=1.5,  c="k")
        ax3.plot(np.log(Rm_N_0), np.polyval( coeff_Rm_N, np.log(Rm_N_0) - np.mean(np.log(Rm_N_0))), ls="-", lw=1.5,  c="k")
        ax3.set_ylabel(r"log ( $T_0$ )", rotation = 270, labelpad = 35) # alpha is the transparancy
        ax3.yaxis.set_label_position("right")
        ax3.locator_params(nbins=10) #gives the number of ticks maximum we want, here 6
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])

        ############ T0 - T0_fit ###########
        ax4.axhline(y=0, ls ="--", c ="k")
        ax4.plot(T0_N,(T0_N-np.exp(np.polyval( coeff_Rp_N, np.log(Rp_N_0) - np.mean(np.log(Rp_N_0)))))*1e3,
                        marker="o", ls="--", lw=3, c="#DBDBDB", ms=8, mew = 0, mec = "r", mfc = "r")
        ax4.plot(T0_N,(T0_N-np.exp(np.polyval( coeff_Rm_N, np.log(Rm_N_0) - np.mean(np.log(Rm_N_0)))))*1e3,
                        marker="o", ls="--", lw=3, c="#DBDBDB", ms=8, mew = 0, mec = "#440154", mfc = "#440154")
        ax4.set_xlabel(r"$T_0$ ( K )")
        ax4.set_ylabel(r"$\Delta T_0$ (fit) ( mK )", rotation = 270, labelpad = 22)
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position("right")
        ax4.locator_params(nbins=6)


        fig.text(0.14,0.73, r"$T^+$", color = "r", fontstyle = "italic", fontweight = "bold")
        fig.text(0.18,0.73, r"$T^-$", color = "#440154", fontstyle = "italic", fontweight = "bold")
        fig.text(0.555,0.73, r"fit order = " + str(deg_fit), color = "k", fontstyle = "italic", fontweight = "bold")

        fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.98,wspace=0.1,hspace=0.05) # adjust the box of axes regarding the figure size

        f_manager = plt.get_current_fig_manager()
        f_manager.window.move(25, 75)
        plt.show()

    ###################################

    return Temperatures, Data_P, Data_N


