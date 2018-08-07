# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample ////////////////////////////////////////////////////////////////#
axe = "a"
samplename =  r"RuCl$_3$ 1804C"
samplelabel = r"RuCl3_1804C"

## Geometric factor ///////////////////////////////////////////////////////////#
## RuCl3 1804C
L = 956e-6 # in meter from L56
w = 3516e-6 # in meter
t = 80e-6 # in meter

## Save geofactor_dict with Pickle ////////////////////////////////////////////#
geofactor = { "L" : L, "w" : w, "t" : t }

with open("geofactor_" + samplelabel + ".pkl", "wb") as f:
    pickler_file = pickle.Pickler(f)
    pickler_file.dump(geofactor)


