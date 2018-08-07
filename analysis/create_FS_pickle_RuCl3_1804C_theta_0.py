# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample ////////////////////////////////////////////////////////////////#
axe = "a"
samplename =  r"RuCl$_3$ 1804C"
samplelabel = r"RuCl3_1804C"

## Calibration Dictionnary ////////////////////////////////////////////////////#
calib = {}
calib[20180706] = "2018-07-06-0deg-0T/Data-TS-0.0T-TOP-RuCl3_1804C-2018-07-06.dat"

## FieldSweep Dictionnary /////////////////////////////////////////////////////#

files_FS={} # keys are (T0, AAAAMMDD) the list ["filename", start_pause, end_pause_pause, calib_file]

## 2018-07-09
# files_FS[5, 20180709] = ["2018-07-09-0deg-FS/Data-FS-5.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[6, 20180709] = ["2018-07-09-0deg-FS/Data-FS-6.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[7, 20180709] = ["2018-07-09-0deg-FS/Data-FS-7.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[8, 20180709] = ["2018-07-09-0deg-FS/Data-FS-8.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[12, 20180709] = ["2018-07-09-0deg-FS/Data-FS-12.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[14, 20180709] = ["2018-07-09-0deg-FS/Data-FS-14.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[16, 20180709] = ["2018-07-09-0deg-FS/Data-FS-16.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[18, 20180709] = ["2018-07-09-0deg-FS/Data-FS-18.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[20, 20180709] = ["2018-07-09-0deg-FS/Data-FS-20.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[23, 20180709] = ["2018-07-09-0deg-FS/Data-FS-23.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[26, 20180709] = ["2018-07-09-0deg-FS/Data-FS-26.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[29, 20180709] = ["2018-07-09-0deg-FS/Data-FS-29.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[32, 20180709] = ["2018-07-09-0deg-FS/Data-FS-32.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[35, 20180709] = ["2018-07-09-0deg-FS/Data-FS-35.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[38, 20180709] = ["2018-07-09-0deg-FS/Data-FS-38.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[41, 20180709] = ["2018-07-09-0deg-FS/Data-FS-41.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[45, 20180709] = ["2018-07-09-0deg-FS/Data-FS-45.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[49, 20180709] = ["2018-07-09-0deg-FS/Data-FS-49.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[54, 20180709] = ["2018-07-09-0deg-FS/Data-FS-54.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[59, 20180709] = ["2018-07-09-0deg-FS/Data-FS-59.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[64, 20180709] = ["2018-07-09-0deg-FS/Data-FS-64.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[69, 20180709] = ["2018-07-09-0deg-FS/Data-FS-69.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[74, 20180709] = ["2018-07-09-0deg-FS/Data-FS-74.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]
# files_FS[79, 20180709] = ["2018-07-09-0deg-FS/Data-FS-79.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]]

## Tests
files_FS[15.96, 20180708] = ["2018-07-09-0deg-FS/Data-FS-15.96K-TOP-RuCl3_1804C-09-07-2018.dat", 20, None, calib[20180706]] # from 15T to -15T AUTOMATIC Needle Valve 0.3T/min, pause 400 seconde, wait for 120
files_FS[15.97, 20180709] = ["2018-07-09-0deg-FS/Data-FS-15.97K-TOP-RuCl3_1804C-09-07-2018.dat", 20, None, calib[20180706]] # from -15T to 15T AUTOMATIC Needle Valve 0.3T/min, pause 400 seconde, wait for 120
# files_FS[15.98, 20180701] = ["2018-07-09-0deg-FS/Data-FS-15.98K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # from 15T to -15T MANUAL Needle Valve 0.3T/min, pause 400 seconde, wait for 90
# files_FS[15.99, 20180709] = ["2018-07-09-0deg-FS/Data-FS-15.99K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # from -15T to 15T MANUAL Needle Valve 0.7T/min, pause 400 seconde, wait for 120
# files_FS[16.01, 20180709] = ["2018-07-09-0deg-FS/Data-FS-16.01K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # from 15T to -15T MANUAL Needle Valve 0.7T/min, pause 400 seconde, wait for 120
# files_FS[16.02, 20180710] = ["2018-07-09-0deg-FS/Data-FS-16.02K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from -15T to 15T MANUAL Needle Valve 0.7T/min, pause 400 seconde, wait for 120
# files_FS[16.03, 20180711] = ["2018-07-09-0deg-FS/Data-FS-16.03K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from 15T to -15T MANUAL Needle Valve 0.7T/min, pause 400 seconde, wait for 120
# files_FS[19.00, 20180712] = ["2018-07-09-0deg-FS/Data-FS-19.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from -15T to 15T MANUAL Needle Valve 0.3T/min, pause 400 seconde, wait for 120
# files_FS[19.01, 20180713] = ["2018-07-09-0deg-FS/Data-FS-19.01K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from 15T to -15T MANUAL Needle Valve 0.3T/min, pause 400 seconde, wait for 120
# files_FS[9.00, 20180714] = ["2018-07-09-0deg-FS/Data-FS-9.00K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from -15T to 15T MANUAL Needle Valve 0.7T/min, pause 400 seconde, wait for 120
# files_FS[9.01, 20180715] = ["2018-07-09-0deg-FS/Data-FS-9.01K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from 15T to -15T MANUAL Needle Valve 0.7T/min, pause 400 seconde, wait for 120
# files_FS[9.02, 20180716] = ["2018-07-09-0deg-FS/Data-FS-9.02K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from -15T to 15T MANUAL Needle Valve 0.3T/min, pause 400 seconde, wait for 120
# files_FS[9.03, 20180717] = ["2018-07-09-0deg-FS/Data-FS-9.03K-TOP-RuCl3_1804C-09-07-2018.dat", 0, None, calib[20180706]] # HEAT OFF from 15T to -15T MANUAL Needle Valve 0.3T/min, pause 400 seconde, wait for 120






## Save info_dict with Pickle /////////////////////////////////////////////////#
with open("info_FS_" + samplelabel + ".pkl", "wb") as f:
    pickler_file = pickle.Pickler(f)
    pickler_file.dump(files_FS)



