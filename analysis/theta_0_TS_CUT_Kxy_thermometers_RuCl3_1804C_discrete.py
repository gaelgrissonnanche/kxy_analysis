# -*- coding: utf-8 -*-

theta = 0

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import numpy as np
import matplotlib.pyplot as plt
from thermal_thermoelectric.analysis import *
from thermal_thermoelectric.load_save_files import *
from thermal_thermoelectric.figures import *
from thermal_thermoelectric.decorators_theta import *
fig_print_TS = theta_TS_figure(theta)(fig_print_TS)
save_TS_file = theta_save_TS_file(theta)(save_TS_file)
save_TS_figure = theta_save_figure(theta)(save_TS_figure)
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

H_cut = 15

## Info sample >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

axe = "ab"
samplename =  r"RuCl$_3$ 1804C"
samplelabel = r"RuCl3_1804C"
measurement = r"TS_CUT_Kxx_Kxy"

sign_dTy = -1

## Columns ////////////////////////////////////////////////////////////////////#
columns_FS = {"col_H" : 0,
             "col_T0" : 1,
             "col_I" : 3,
             "col_Rp" : 4,
             "col_Rm" : 5,
             "col_Vy" : 6,
             "col_pause" : 7,
             }

columns_calib = {"col_T0_calib" : 0,
                "col_Rp_0_calib" : 3,
                "col_Rm_0_calib" : 5,
                }

## Geofactor //////////////////////////////////////////////////////////////////#
geofactor = load_pickle_geofactor(samplelabel)

## Analysis FS ////////////////////////////////////////////////////////////////#
files_FS = load_pickle_FS(samplelabel)
NSYM_dict, SYM_dict = FS_Kxy_discrete(files_FS, columns_FS, geofactor,
    sign_dTy = sign_dTy, columns_calib = columns_calib, T_connected = "T+")

## Analysis TS ////////////////////////////////////////////////////////////////#
TS_CUT = TS_CUT_discrete(SYM_dict, H_cut)

## Save data //////////////////////////////////////////////////////////////////#
Data = np.vstack((TS_CUT["Tav"], TS_CUT["dTx"], TS_CUT["Kxx"], TS_CUT["Kxy"],
                  TS_CUT["dTy"], TS_CUT["Tp"], TS_CUT["Tm"],
                  TS_CUT["T0"], TS_CUT["I"], TS_CUT["date_list"]))
Data = Data.transpose()
Header = "Tav[K]\tdTx[K]\tKxx[W/Km]\tKxy[W/Km]\tdTy[K]\tT+[K]\tT-[K]\tT0[K]\tI[A]\n"

# Get all dates of Fieldsweeps
dates = dates_function(SYM_dict)

save_TS_file(samplelabel = samplelabel, axe = axe,
measurement = measurement, H = H_cut, date = dates, Data = Data, Header = Header)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
# Figures PRINT ///////////////////////////////////////////////////////////////#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

fig_list = []

# dT/T
fig = fig_print_TS(x = TS_CUT["Tav"], y = TS_CUT["dTx"] / TS_CUT["Tav"] * 100, H = H_cut,
            ymin = 0, xlabel = r"T ( K )", ylabel = r"dT$_x$ / T$_{\rm av}$ ( % )",
            samplename = samplename, samplename_loc = 3, color = "#515151")
fig_list.append(fig)

# dTx
fig = fig_print_TS(x = TS_CUT["Tav"], y = TS_CUT["dTx"], H = H_cut,
            xlabel = r"T ( K )", ylabel = r"dT$_x$ ( K )",
            samplename = samplename, samplename_loc = 4, color = "#515151")
fig_list.append(fig)

# Kxy/T
fig = fig_print_TS(x = TS_CUT["Tav"], y = TS_CUT["Kxy"] / TS_CUT["Tav"] * 1e6, H = H_cut,
            ymin = 0, xlabel = r"T ( K )", ylabel = r"$\kappa_{\rm xy}$ / T ( $\mu$W / K$^2$ m )",
            samplename = samplename, samplename_loc = 2, color = "#440154")
fig_list.append(fig)

# Kxy
fig = fig_print_TS(x = TS_CUT["Tav"], y = TS_CUT["Kxy"] * 1e3, H = H_cut,
            ymin = 0, xlabel = r"T ( K )", ylabel = r"$\kappa_{\rm xy}$ ( mW / K m )",
            samplename = samplename, samplename_loc = 4, color = "#440154")
fig_list.append(fig)

# dTy
fig = fig_print_TS(x = TS_CUT["Tav"], y = TS_CUT["dTy"] * 1e3, H = H_cut,
            ymin = 0, xlabel = r"T ( K )", ylabel = r"dT$_y$ ( mK )",
            samplename = samplename, samplename_loc = 4, color = "#515151")
fig_list.append(fig)

# Kxx
fig = fig_print_TS(x = TS_CUT["Tav"], y = TS_CUT["Kxx"], H = H_cut,
            ymin = 0, xlabel = r"T ( K )", ylabel = r"$\kappa_{\rm xx}$ ( W / K m )",
            samplename = samplename, samplename_loc = 4, color = "r")
fig_list.append(fig)

#//////////////////////////////////////////////////////////////////////////#

plt.show()

## Create PDF file for all figures >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
save_TS_figure(samplelabel = samplelabel, axe = axe, measurement = measurement,
               H = H_cut, date = dates, fig_list = fig_list)

