# -*- coding: utf-8 -*-

theta = 0

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import numpy as np
import matplotlib.pyplot as plt
from thermal_thermoelectric.analysis import *
from thermal_thermoelectric.thermometry import *
from thermal_thermoelectric.load_save_files import *
from thermal_thermoelectric.figures import *
from thermal_thermoelectric.decorators_theta import *
fig_print_TS = theta_TS_figure(theta)(fig_print_TS)
save_TS_file = theta_save_TS_file(theta)(save_TS_file)
save_TS_figure = theta_save_figure(theta)(save_TS_figure)
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

diagnostic_fig = False # shows the figures for diagnostic
print_fig = True # print in pdf the figures

## Which measurement? /////////////////////////////////////////////////////////#
H = 0 # is H_m, the field from the magnet
date = "2018-07-06"

## Info sample ////////////////////////////////////////////////////////////////#
axe = "ab"
samplename =  r"RuCl$_3$ 1804C"
samplelabel = r"RuCl3_1804C"
measurement = "TS_NSYM_Kxx"

## Data columns ///////////////////////////////////////////////////////////////#
columns = {
          "col_T0": 0, "col_I" : 2,
          "col_Rp_0": 3, "col_Rp_Q": 4,
          "col_Rm_0": 5, "col_Rm_Q": 6,
          "col_dTy_0": 7, "col_dTy_Q": 8
          }

## Load Pickles ///////////////////////////////////////////////////////////////#
geofactor = load_pickle_geofactor(samplelabel)
info = load_pickle_info(samplelabel)

## Load data //////////////////////////////////////////////////////////////////#
(data_P, data_N) = load_TS_file(samplelabel = samplelabel, H = H,
                measurement = measurement, date = date, **info[H, theta, date])

## Thermometers Analysis //////////////////////////////////////////////////////#
Temperatures, data_P, data_N = temp_thermometers(data_P, data_N,
                             graph = diagnostic_fig, cut_off = True, **columns)

T0 = Temperatures["T0"] # K
Tav = Temperatures["Tav"] # K
dTx = Temperatures["dTx"] # K
Tp = Temperatures["Tp"] # K
Tm = Temperatures["Tm"] # K

# Compute Kxx /////////////////////////////////////////////////////////////////#
ThermalConductivity = TS_Kxx(data_P, data_N, Temperatures, geofactor, **columns)

I = ThermalConductivity["I"] # Amps
Kxx = ThermalConductivity["Kxx"] # W / K m

## Save data //////////////////////////////////////////////////////////////////#
Data = np.vstack((Tav, dTx, Kxx, Tp, Tm, T0, I))
Data = Data.transpose()
Header = "Tav[K]\tdTx[K]\tKxx[W/Km]\tT+[K]\tT-[K]\tT0[K]\tI[A]\n"

save_TS_file(samplelabel = samplelabel, axe = axe, measurement = measurement,
             H = H, date = date, Data = Data, Header = Header)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
# Figures Diagnostic //////////////////////////////////////////////////////////#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

if diagnostic_fig is True:
    fig_list = []
    # Thermometry & Kxx
    fig = fig_diagno_thermometry(Temperatures, ThermalConductivity, H, measurement, samplename = samplename)
    #//////////////////////////////////////////////////////////////////////////#
    plt.show()


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
# Figures PRINT ///////////////////////////////////////////////////////////////#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

if print_fig is True:
    fig_list = []

    # dT/T
    fig = fig_print_TS(x = Tav, y = dTx / Tav * 100, H = H,
                ymin = 0, xlabel = r"$T$ ( K )", ylabel = r"dT$_x$ / T ( % )",
                samplename = samplename, samplename_loc = 3, color = "#515151")
    fig_list.append(fig)

    # dTx
    fig = fig_print_TS(x = Tav, y = dTx, H = H,
                xlabel = r"$T$ ( K )", ylabel = r"dT$_x$ ( K )",
                samplename = samplename, samplename_loc = 4, color = "#515151")
    fig_list.append(fig)

    # Kxx
    fig = fig_print_TS(x = Tav, y = Kxx, H = H,
                ymin = 0, ymax = 4, xlabel = r"$T$ ( K )", ylabel = r"$\kappa_{\rm xx}$ ( W / K m )",
                samplename = samplename, samplename_loc = 4, color = "r")
    fig_list.append(fig)

    #//////////////////////////////////////////////////////////////////////////#

    plt.show()

    ## Create PDF file for all figures >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
    save_TS_figure(samplelabel = samplelabel, axe = axe, measurement = measurement,
                   H = H, date = date, fig_list = fig_list)