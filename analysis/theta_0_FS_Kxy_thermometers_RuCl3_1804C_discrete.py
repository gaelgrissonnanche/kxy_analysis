# -*- coding: utf-8 -*-

theta = 0

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import numpy as np
import matplotlib.pyplot as plt
from thermal_thermoelectric.analysis import *
from thermal_thermoelectric.load_save_files import *
from thermal_thermoelectric.figures import *
from thermal_thermoelectric.decorators_theta import *
fig_print_FS = theta_FS_figure(theta)(fig_print_FS)
save_FS_files = theta_save_FS_files(theta)(save_FS_files)
save_FS_figure = theta_save_figure(theta)(save_FS_figure)
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

axe = "ab"
samplename =  r"RuCl$_3$ 1804C"
samplelabel = r"RuCl3_1804C"
measurement = r"FS_Kxx_Kxy"

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

## Save data //////////////////////////////////////////////////////////////////#
data_type_list = ["H", "dTy", "Kxy", "Kxx", "dTx", "Tav", "Vy", "I", "T0"]
Header = "H[T]\tdTy[K]\tKxy[W/Km]\tKxx[W/Km]\tdTx[K]\tTav[K]\tVy[V]\tI[A]\tT0[K]\n"

save_FS_files(samplelabel = samplelabel, axe = axe, measurement = measurement,
SYM_dict = SYM_dict, data_type_list = data_type_list, Header = Header)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
# Figures PRINT ///////////////////////////////////////////////////////////////#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

fig_list = []

# Vy paused for diagnostic
x = NSYM_dict["H_diagnostic"]
y = {t:v*1e6 for t, v in NSYM_dict["Vy_diagnostic"].items()}
fig = fig_print_FS(x, y, xmin = None, xmax = None, ymin = None, ymax = None,
    ylabel = r"$\Delta V_{\rm y}$ ( $\mu$V )", samplename = samplename,
    samplename_loc = 1,  colormap = "jet", vx = 0)
fig_list.append(fig)

# Vy SYM
x = SYM_dict["H"]
y = {t:v*1e9 for t, v in SYM_dict["Vy"].items()}
fig = fig_print_FS(x, y, xmin = 0, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\Delta V_{\rm y}$ ( nV )", samplename = samplename,
    samplename_loc = 1,  colormap = "jet")
fig_list.append(fig)

# Vy NONSYM
x = NSYM_dict["H"]
y = {t:(v - np.mean(v))*1e6 for t, v in NSYM_dict["Vy"].items()}
fig = fig_print_FS(x, y, xmin = -15, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\Delta V_{\rm y}$ ( $\mu$V )", samplename = samplename,
    samplename_loc = 2,  colormap = "jet", hl = 0, vx = 0)
fig_list.append(fig)

# T0 NONSYM
x = NSYM_dict["H"]
y = {t:(v - np.mean(v))*1e3 for t, v in NSYM_dict["T0"].items()}
fig = fig_print_FS(x, y, xmin = -15, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\Delta T_{\rm 0}$ ( mK )", samplename = samplename,
    samplename_loc = 3,  colormap = "jet", hl = 0, vx = 0)
fig_list.append(fig)

# Kxx SYM
x = SYM_dict["H"]
y = {t:(v - v[0]) for t, v in SYM_dict["Kxx"].items()}
fig = fig_print_FS(x, y, xmin = 0, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\Delta\kappa_{\rm xx}$ ( W / K m )", samplename = samplename,
    samplename_loc = 3,  colormap = "jet")
fig_list.append(fig)

# Kxy / H / T SYM
x = SYM_dict["H"]
y = {t:v*1e6/SYM_dict["H"][t]/SYM_dict["Tav"][t] for t, v in SYM_dict["Kxy"].items()}
fig = fig_print_FS(x, y, xmin = 0, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\kappa_{\rm xy}$ / $H$ $T$ ( $\mu$W / K$^2$ m T )", samplename = samplename,
    samplename_loc = 2,  colormap = "jet")
fig_list.append(fig)

# Kxy / H SYM
x = SYM_dict["H"]
y = {t:v*1e3/SYM_dict["H"][t] for t, v in SYM_dict["Kxy"].items()}
fig = fig_print_FS(x, y, xmin = 0, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\kappa_{\rm xy}$ / $H$ ( mW / K m T )", samplename = samplename,
    samplename_loc = 2,  colormap = "jet")
fig_list.append(fig)

# Kxy SYM
x = SYM_dict["H"]
y = {t:v*1e3 for t, v in SYM_dict["Kxy"].items()}
fig = fig_print_FS(x, y, xmin = 0, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\kappa_{\rm xy}$ ( mW / K m )", samplename = samplename,
    samplename_loc = 1,  colormap = "jet")
fig_list.append(fig)

# dTy SYM
x = SYM_dict["H"]
y = {t:v*1e3 for t, v in SYM_dict["dTy"].items()}
fig = fig_print_FS(x, y, xmin = 0, xmax = 15, ymin = None, ymax = None,
    ylabel = r"$\Delta T_{\rm y}$ ( mK )", samplename = samplename,
    samplename_loc = 1,  colormap = "jet")
fig_list.append(fig)


#//////////////////////////////////////////////////////////////////////////////#

plt.show()

## Create PDF file for all figures >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# Get all dates of Fieldsweeps
dates = dates_function(SYM_dict)
save_FS_figure(samplelabel = samplelabel, axe = axe, measurement = measurement,
                date = dates, fig_list = fig_list)





