# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import pickle
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

## Info sample ////////////////////////////////////////////////////////////////#
axe = "a"
samplename =  r"RuCl$_3$ 1804C"
samplelabel = r"RuCl3_1804C"

## Create dictionnaries of measurements info //////////////////////////////////#
info = {}     # keys are (H_ext, theta, date)

info[15, 0, "2018-04-15"] = {
                        "folder" : "2018-04-15-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, from 2K to 35K, dT/T 8% to 5%" ,
                        }

info[15, 0, "2018-04-17"] = {
                        "folder" : "2018-04-17-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 15% to 9%, but pump noise removed only at +15T",
                        }

info[15, 0, "2018-04-21"] = {
                        "folder" : "2018-04-21-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 15% to 9%, pump noise removed at both field",
                        }

info[0, 0, "2018-04-24"] = {
                        "folder" : "2018-04-24-0deg-0T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 4% to 2%",
                        }

info[15, 0, "2018-04-26"] = {
                        "folder" : "2018-04-26-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 10% to 8%",
                        }

info[15, 0, "2018-05-03"] = {
                        "folder" : "2018-05-03-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 5% to 4%",
                        }

info[15, 0, "2018-05-06"] = {
                        "folder" : "2018-05-06-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 5% to 4%, from 64K to 90K",
                        }

info[15, 0, "2018-05-12"] = {
                        "folder" : "2018-05-12-0deg-15T",
                        "position" : "BOT",
                        "comments" : "Good data, dT/T 5% to 4%, field changed @ 90K",
                        }

info[0, 0, "2018-07-06"] = {
                        "folder" : "2018-07-06-0deg-0T",
                        "position" : "TOP",
                        "comments" : "Good data, dT/T 2% to 1%",
                        }

## Save info_dict with Pickle /////////////////////////////////////////////////#

with open("info_TS_" + samplelabel + ".pkl", "wb") as f:

    pickler_file = pickle.Pickler(f)
    pickler_file.dump(info)

