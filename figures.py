# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
import numpy as np
from numpy import cos, sin, arange, round, array
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def fig_print_TS(x, y, H, xmin = 0, xmax = None, ymin = None, ymax = None,
    xlabel = "", ylabel = "", samplename = "", samplename_loc = 1, color = "r"):

    """
    This function will plot one curve on one plot in the dimensions asked for
    a better figure format that called "print figures". It takes for arguments:
    - x : the x-axis data
    - y : the y-axis data
    - H : field from the magnet H_m
    - xmin : x mininum
    - xmax : x maximum
    - ymin : y mininum
    - ymax : y maximum
    - xlabel : the label for the x-axis, should be a string
    - ylabel : the label for the y-axis, should be a string
    - samplename : the name of the sample
    - samplenameloc : the name of the sample location. It should be 0,1,2 or 3
                that correspond to upper left, upper right, low left, low right
    - color : the color of the points of the curve

    Returns a figure object as "fig"
    """

    #///// RC Parameters //////#
    mpl.rcdefaults()
    mpl.rcParams['font.size'] = 24. # change the size of the font in every figure
    mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
    mpl.rcParams['axes.labelsize'] = 24.
    mpl.rcParams['xtick.labelsize'] = 24
    mpl.rcParams['ytick.labelsize'] = 24
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.width'] = 0.6
    mpl.rcParams['ytick.major.width'] = 0.6
    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
    mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                        # editing the text in illustrator

    #///// Create Figure //////#
    fig , axes = plt.subplots(1,1, figsize=(9.2, 5.6)) # figsize is w x h in inch of figure
    fig.subplots_adjust(left = 0.17, right = 0.81, bottom = 0.18, top = 0.95) # adjust the box of axes regarding the figure size


    #///// Allow to shift the label ticks up or down with set_pad /////#
    for tick in axes.xaxis.get_major_ticks():
        tick.set_pad(7)
    for tick in axes.yaxis.get_major_ticks():
        tick.set_pad(8)

    #///// Legend //////#
    loc_dict = {}  # Position the sample name
    loc_dict[1] = [0.21,0.87, "left"] # coordinates
    loc_dict[2] = [0.79,0.87, "right"]
    loc_dict[3] = [0.21,0.23, "left"]
    loc_dict[4] = [0.79,0.23, "right"]

    loc = loc_dict[samplename_loc] # select the right location among coordinates
    fig.text(loc[0], loc[1], samplename, ha = loc[2])

    fig.text(0.82,0.87, r"$H$ = ", ha = 'left', fontsize = 20)
    fig.text(0.885,0.87, "{0:g}".format(H) + " T", ha="left", fontsize = 20)

    #///// Choose if horizontal line y = 0  is needed /////#
    if min(y) * max(y) < 0:
        axes.axhline(y = 0, ls = "--", c = "k", linewidth = 0.6)

    #///// Plot /////#
    line = axes.plot(x, y)
    plt.setp(line, ls = "--", c = '#DBDBDB', lw = 3, marker = "o", mfc = 'w', ms = 6.5, mec = color, mew = 2.5)
    line = axes.plot(x, y)
    plt.setp(line, ls = "", c = '#DBDBDB', lw = 3, marker = "o", mfc = color, ms = 6.5, mec = color, mew = 0, alpha = 0.7)

    #///// Choosing the limits and ticks /////#
    axes.set_xlim(xmin, xmax)   # limit for xaxis
    axes.set_ylim(ymin, ymax) # leave the ymax auto, but fix ymin
    axes.locator_params(axis = 'x', nbins = 6)
    axes.locator_params(axis = 'y', nbins = 6)
    axes.set_xlabel(xlabel, labelpad = 8)
    axes.set_ylabel(ylabel)

    return fig


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def fig_print_FS(x_dict, y_dict, xmin = 0, xmax = None, ymin = None, ymax = None,
    xlabel = r"$H$ ( T )", ylabel = "", samplename = "", samplename_loc = 1,  markevery = 1,
    colormap = "jet", hl = None, vx = None):

    """
    This function will plot all the isotherms of some type of measurements.
    It takes for arguments:
    - x_dict : the x-axis dictionnary (keys are temperature and date) for the x-axis data
    - y_dict : the y-axis dictionnary (keys are temperature and date) data
    - xmin : x mininum
    - xmax : x maximum
    - ymin : y mininum
    - ymax : y maximum
    - ylabel : the label for the y-axis, should be a string
    - samplename : the name of the sample
    - samplenameloc : the name of the sample location. It should be 0,1,2 or 3
                that correspond to upper left, upper right, low left, low right
    - colormap : the colormap to display the curves
    - hl : the value of y for the horizontal dashline
    - vx : the value of x for the vertical dashline

    Returns a figure object as "fig"
    """

    #///// RC Parameters //////#
    mpl.rcdefaults()
    mpl.rcParams['font.size'] = 24. # change the size of the font in every figure
    mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
    mpl.rcParams['axes.labelsize'] = 24.
    mpl.rcParams['xtick.labelsize'] = 24
    mpl.rcParams['ytick.labelsize'] = 24
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.width'] = 0.6
    mpl.rcParams['ytick.major.width'] = 0.6
    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
    mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                        # editing the text in illustrator

    #///// Create Figure //////#
    fig , axes = plt.subplots(1,1, figsize=(9.2, 5.6)) # figsize is w x h in inch of figure
    fig.subplots_adjust(left = 0.17, right = 0.81, bottom = 0.18, top = 0.95) # adjust the box of axes regarding the figure size


    #///// Allow to shift the label ticks up or down with set_pad /////#
    for tick in axes.xaxis.get_major_ticks():
        tick.set_pad(7)
    for tick in axes.yaxis.get_major_ticks():
        tick.set_pad(8)

    ## Field components ///////////////////////////////////////////////////////////#

    #///// Legend //////#
    loc_dict = {}  # Position the sample name
    loc_dict[1] = [0.21,0.87, "left"] # coordinates
    loc_dict[2] = [0.79,0.87, "right"]
    loc_dict[3] = [0.21,0.23, "left"]
    loc_dict[4] = [0.79,0.23, "right"]

    loc = loc_dict[samplename_loc] # select the right location among coordinates
    fig.text(loc[0], loc[1], samplename, ha = loc[2])

    fig.text(0.905,0.92, r"$T$ ( K ) -", fontsize = 12, ha = "right")
    fig.text(0.93,0.92, r"date", fontsize = 12, ha = "left")

    #///// Choose if horizontal line y = hl  is needed /////#
    if hl != None:
        axes.axhline(y = hl, ls = "--", c = "k", linewidth = 0.6)

    #///// Choose if vertical line x = vx  is needed /////#
    if vx != None:
        axes.axvline(x = vx, ls = "--", c = "k", linewidth = 0.6)

    ## Dictionnary keys
    keys = sorted(x_dict.keys()) # (Tav, date)

    #///// Colormap //////#
    cmap = mpl.cm.get_cmap(colormap, len(keys))
    colors = cmap(arange(len(keys)))

    if colormap == "viridis":
        colors[-1] = (1, 0, 0, 1)

    #///// Plot /////#
    k = len(keys)-1
    for key in keys[::-1]:

        T = key[0]
        date = key[1]

        line = axes.plot(x_dict[key],y_dict[key], markevery = markevery)
        plt.setp(line, ls ="-", c = colors[k], lw = 2.5, marker = "o", ms = 7, mew= 0, zorder = k)

        fig.text(0.99, 0.89 - k * 0.025, "{0}".format(round(T,2)) +
                 " - " "{0}".format(round(date,0)),
                      color = colors[k], fontsize = 12, ha ="right")

        k -= 1

    #///// Choosing the limits and ticks /////#
    axes.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    axes.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    axes.set_xlim(xmin, xmax)   # limit for xaxis
    axes.set_ylim(ymin, ymax) # leave the ymax auto, but fix ymin
    axes.locator_params(axis = 'x', nbins = 6)
    axes.locator_params(axis = 'y', nbins = 6)
    axes.set_xlabel(xlabel, labelpad = 8)
    axes.set_ylabel(ylabel, labelpad = 8)

    return fig


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def fig_diagno_thermometry(Temperatures, ThermalConductivity, H, measurement,
                           xmin = 0, xmax = None, samplename = ""):

    """
    This function will plot the diagnostic of the entire thermometry, meaning:
    dTx, Tav, dT/T, Kxx. It takes for arguments:
    - Temperatures : the dictionnary of temperatures
    - ThermalConductivity : the dictionnary of thermal conductivity
    - H : the magnetic field value
    - measurement : the string that give "NSYM_..." or "SYM" to indicate which
        type of experiment it is and if we are asking symmetrizing the data or not
    - xmin : x mininum
    - xmax : x maximum
    - samplename : the name of the sample

    Returns a figure object as "fig"
    """

    #///// RC Parameters //////#
    mpl.rcdefaults()
    mpl.rcParams['font.size'] = 20. # change the size of the font in every figure
    mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
    mpl.rcParams['axes.labelsize'] = 20.
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.width'] = 0.6
    mpl.rcParams['ytick.major.width'] = 0.6
    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
    mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                        # editing the text in illustrator

    #///// Create Figure //////#
    fig, axes = plt.subplots(2, 2, figsize = (11,8)) # figsize is w x h in inch of figure
    fig.subplots_adjust(left = 0.09, right = 0.96, bottom = 0.1, top = 0.88, wspace = 0.3, hspace = 0.2) # adjust the box of axes regarding the figure size


    if (H >= 0) and measurement[0:7] == "TS_NSYM": # only plot H >= 0 if NSYM needed
        axes[0,0].plot(Temperatures["T0_P"], Temperatures["dTx_P"], ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1,0].plot(Temperatures["T0_P"], Temperatures["Tav_P"] - Temperatures["T0_P"], ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[0,1].plot(Temperatures["T0_P"], Temperatures["dTx_P"]/Temperatures["Tav_P"]*100, ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1,1].plot(Temperatures["T0_P"], ThermalConductivity["Kxx_P"], ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[0,0].annotate("H +" + str(abs(H)) + " T", xy = (0,0), xytext = (0.12,0.83), ha="left", textcoords = "figure fraction", color = "r", fontsize = 20, fontweight = "bold")

    if (H < 0) and measurement[0:7] == "TS_NSYM": # only plot H < 0 if NSYM needed
        axes[0,0].plot(Temperatures["T0_N"], Temperatures["dTx_N"], ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1,0].plot(Temperatures["T0_N"], Temperatures["Tav_N"] - Temperatures["T0_N"], ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[0,1].plot(Temperatures["T0_N"], Temperatures["dTx_N"]/Temperatures["Tav_N"]*100, ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1,1].plot(Temperatures["T0_N"], ThermalConductivity["Kxx_N"], ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[0,0].annotate("H -" + str(abs(H)) + " T", xy = (0,0), xytext = (0.12,0.83), ha="left", textcoords = "figure fraction", color = "#440154", fontsize = 20, fontweight = "bold")

    if measurement[0:6] == "TS_SYM": # plots H>0, H<0 and SYM when needed
        axes[0,0].plot(Temperatures["T0_P"], Temperatures["dTx_P"], ls = "", lw = 1, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[0,0].plot(Temperatures["T0_N"], Temperatures["dTx_N"], ls = "", lw = 1, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[0,0].plot(Temperatures["T0"], Temperatures["dTx"], ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

        axes[1,0].plot(Temperatures["T0_P"], Temperatures["Tav_P"] - Temperatures["T0_P"], ls = "", lw = 1, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1,0].plot(Temperatures["T0_N"], Temperatures["Tav_N"] - Temperatures["T0_N"], ls = "", lw = 1, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1,0].plot(Temperatures["T0"], Temperatures["Tav"] - Temperatures["T0"], ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

        axes[0,1].plot(Temperatures["T0_P"], Temperatures["dTx_P"]/Temperatures["Tav_P"]*100, ls = "", lw = 1, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[0,1].plot(Temperatures["T0_N"], Temperatures["dTx_N"]/Temperatures["Tav_N"]*100, ls = "", lw = 1, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[0,1].plot(Temperatures["T0"], Temperatures["dTx"]/Temperatures["Tav"]*100, ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

        axes[1,1].plot(Temperatures["T0_P"], ThermalConductivity["Kxx_P"], ls = "", lw = 1, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1,1].plot(Temperatures["T0_N"], ThermalConductivity["Kxx_N"], ls = "", lw = 1, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1,1].plot(Temperatures["T0"], ThermalConductivity["Kxx"], ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

        axes[0,0].annotate("H +" + str(abs(H)) + " T", xy = (0,0), xytext = (0.12,0.83), ha="left", textcoords = "figure fraction", color = "r", fontsize = 20, fontweight = "bold")
        axes[0,0].annotate("H -" + str(abs(H)) + " T", xy = (0,0), xytext = (0.12,0.79), ha="left", textcoords = "figure fraction", color = "#440154", fontsize = 20, fontweight = "bold")
        axes[0,0].annotate("SYM", xy = (0,0), xytext = (0.12,0.75), ha="left", textcoords = "figure fraction", color = "#A7A7A7", fontsize = 20, fontweight = "bold")

    #///// dTx //////#
    axes[0,0].set_ylabel(r"$dT_x$ ( K )")
    axes[0,0].locator_params(nbins = 6) #gives the number of ticks maximum we want, here 6
    axes[0,0].set_xlim(xmin, xmax)
    axes[0,0].set_ylim(bottom=0)

    #///// Tav //////#
    axes[1,0].set_xlabel(r"$T_0$ ( K )")
    axes[1,0].set_ylabel(r"$T_{\mathrm{av}}$ - $T_0$ ( K )")
    axes[1,0].locator_params(nbins = 6)
    axes[1,0].set_xlim(xmin, xmax)
    axes[1,0].set_ylim(bottom=0)

    #///// dTx/T //////#
    axes[0,1].set_ylabel(r"$dT_x$ / $T_{\mathrm{av}}$ ( % )")
    axes[0,1].locator_params(nbins = 6)
    axes[0,1].set_xlim(xmin, xmax)

    #///// Kxx //////#
    axes[1,1].set_xlabel(r"$T_0$ ( K )")
    axes[1,1].set_ylabel(r"$\kappa_{\mathrm{xx}}$ ( W / K m )")
    axes[1,1].locator_params(nbins = 6)
    axes[1,1].set_xlim(xmin, xmax)
    axes[1,1].set_ylim(bottom=0)


    fig.text(0.15, 0.93, "Thermometry", fontsize = 34)
    fig.text(0.96, 0.93, samplename, fontsize = 34, ha = "right", color = "#A7A7A7")

    f_manager = plt.get_current_fig_manager()
    f_manager.window.move(25, 75)

    return fig


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def fig_diagno_ThermalHallConductivity(Temperatures, ThermalHallConductivity, H,
                           xmin = 0, xmax = None, samplename = ""):

    """
    This function will plot the diagnostic of the entire thermal Hall effect
    analysis, meaning:
    V_dTy, dTy, Kxy, Kxy/T. It takes for arguments:
    - Temperatures : the dictionnary of temperatures
    - ThermalHallConductivity : the dictionnary of thermal Hall conductivity
    - H : the magnetic field value
    - xmin : x mininum
    - xmax : x maximum
    - samplename : the name of the sample

    Returns a figure object as "fig"
    """

    #///// RC Parameters //////#
    mpl.rcdefaults()
    mpl.rcParams['font.size'] = 20. # change the size of the font in every figure
    mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
    mpl.rcParams['axes.labelsize'] = 20.
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.width'] = 0.6
    mpl.rcParams['ytick.major.width'] = 0.6
    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
    mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                        # editing the text in illustrator

    #///// Create Figure //////#
    fig, axes = plt.subplots(2, 2, figsize = (11,8)) # figsize is w x h in inch of figure
    fig.subplots_adjust(left = 0.09, right = 0.96, bottom = 0.1, top = 0.88, wspace = 0.3, hspace = 0.2)


    axes[0,0].plot(Temperatures["T0_P"], ThermalHallConductivity["V_dTy_P"]*1e6, ls = "", lw = 1, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
    axes[0,0].plot(Temperatures["T0_N"], ThermalHallConductivity["V_dTy_N"]*1e6, ls = "", lw = 1, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")

    axes[1,0].plot(Temperatures["T0"], ThermalHallConductivity["dTy"]*1e3, ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

    axes[0,1].plot(Temperatures["T0"], ThermalHallConductivity["Kxy"]*1e3, ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

    axes[1,1].plot(Temperatures["T0"], ThermalHallConductivity["Kxy"]/Temperatures["Tav"]*1e3, ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")


    if min(ThermalHallConductivity["dTy"]) * max(ThermalHallConductivity["dTy"]) < 0: # choose if horizontal line y = 0  is needed /////#
        axes[1,0].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)
        axes[0,1].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)
        axes[1,1].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)

    #///// V_y //////#
    axes[0,0].set_ylabel(r"$V_y$ ( $\mu$V )")
    axes[0,0].locator_params(nbins = 6) #gives the number of ticks maximum we want, here 6
    axes[0,0].set_xlim(xmin, xmax)

    #///// dTy //////#
    axes[1,0].set_xlabel(r"$T_0$ ( K )")
    axes[1,0].set_ylabel(r"$dT_y$ ( mK )")
    axes[1,0].locator_params(nbins = 6)
    axes[1,0].set_xlim(xmin, xmax)

    #///// Kxy //////#
    axes[0,1].set_ylabel(r"$\kappa_{\mathrm{xy}}$ ( mW / K m )")
    axes[0,1].locator_params(nbins = 6)
    axes[0,1].set_xlim(xmin, xmax)

    #///// Kxy/T //////#
    axes[1,1].set_xlabel(r"$T_0$ ( K )")
    axes[1,1].set_ylabel(r"$\kappa_{\mathrm{xy}}$ / T ( mW / K$^2$ m )")
    axes[1,1].locator_params(nbins = 6)
    axes[1,1].set_xlim(xmin, xmax)

    fig.text(0.15, 0.93, "Thermal Hall", fontsize = 34)
    fig.text(0.96, 0.93, samplename, fontsize = 34, ha = "right", color = "#A7A7A7")

    f_manager = plt.get_current_fig_manager()
    f_manager.window.move(25, 75)

    return fig


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def fig_diagno_Seebeck(Temperatures, Seebeck, H, measurement,
                           xmin = 0, xmax = None, samplename = ""):

    """
    This function will plot the diagnostic of the entire Seebeck analysis, meaning:
    V_S and S/T. It takes for arguments:
    - Temperatures : the dictionnary of temperatures
    - Seebeck : the dictionnary of Seebeck
    - H : the magnetic field value
    - measurement : the string that give "NSYM_..." or "SYM" to indicate which
        type of experiment it is and if we are asking symmetrizing the data or not
    - xmin : x mininum
    - xmax : x maximum
    - samplename : the name of the sample

    Returns a figure object as "fig"
    """

    #///// RC Parameters //////#
    mpl.rcdefaults()
    mpl.rcParams['font.size'] = 20. # change the size of the font in every figure
    mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
    mpl.rcParams['axes.labelsize'] = 20.
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.width'] = 0.6
    mpl.rcParams['ytick.major.width'] = 0.6
    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
    mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                        # editing the text in illustrator

    #///// Create Figure //////#
    fig, axes = plt.subplots(2, 1, figsize = (6, 8)) # figsize is w x h in inch of figure
    fig.subplots_adjust(left = 0.19, right = 0.95, bottom = 0.1, top = 0.9,wspace = 0.3,hspace = 0.2)

    if (H >= 0) and measurement[0:4] == "NSYM": # only plot H >= 0 if NSYM needed
        axes[0].plot(Temperatures["T0_P"], Seebeck["V_S_P"] * 1e6, ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1].plot(Temperatures["T0_P"], Seebeck["S_P"] / Temperatures["Tav_P"], ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")

    if (H < 0) and measurement[0:4] == "NSYM": # only plot H < 0 if NSYM needed
        axes[0].plot(Temperatures["T0_N"], Seebeck["V_S_N"] * 1e6, ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1].plot(Temperatures["T0_N"], Seebeck["S_N"] / Temperatures["Tav_N"], ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")

    if measurement[0:3] == "SYM": # plots H>0, H<0 and SYM when needed
        axes[0].plot(Temperatures["T0_P"], Seebeck["V_S_P"] * 1e6, ls = "", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[0].plot(Temperatures["T0_N"], Seebeck["V_S_N"] * 1e6, ls = "", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[0].plot(Temperatures["T0"], Seebeck["V_S"] * 1e6, ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")

        axes[1].plot(Temperatures["T0_P"], Seebeck["S_P"] / Temperatures["Tav_P"], ls = "", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1].plot(Temperatures["T0_N"], Seebeck["S_N"] / Temperatures["Tav_N"], ls = "", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1].plot(Temperatures["T0"], Seebeck["S"] / Temperatures["Tav"], ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")



    if min(Seebeck["S"]) * max(Seebeck["S"]) < 0: # choose if horizontal line y = 0  is needed /////#
        axes[0].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)
        axes[1].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)

    #///// V_S //////#
    axes[0].set_ylabel(r"$V_S$ ( $\mu$V )")
    axes[0].locator_params(nbins = 6) #gives the number of ticks maximum we want, here 6
    axes[0].set_xlim(xmin, xmax)

    #///// S/T //////#
    axes[1].set_xlabel(r"$T_0$ ( K )")
    axes[1].set_ylabel(r"$S$ / $T$ ( $\mu$V  / K$^2$ )")
    axes[1].locator_params(nbins = 6)
    axes[1].set_xlim(xmin, xmax)

    fig.text(0.19, 0.95, "Seebeck", fontsize = 20)
    fig.text(0.955, 0.94, samplename, fontsize = 24, ha = "right", color = "#A7A7A7")

    f_manager = plt.get_current_fig_manager()
    f_manager.window.move(950, 75)

    return fig


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#


def fig_diagno_Nernst(Temperatures, Nernst, H, measurement,
                           xmin = 0, xmax = None, samplename = ""):

    """
    This function will plot the diagnostic of the entire Seebeck analysis, meaning:
    V_N and nu/T. It takes for arguments:
    - Temperatures : the dictionnary of temperatures
    - Nernst : the dictionnary of Nernst
    - H : the magnetic field value
    - measurement : the string that give "NSYM_..." or "SYM" to indicate which
        type of experiment it is and if we are asking symmetrizing the data or not
    - xmin : x mininum
    - xmax : x maximum
    - samplename : the name of the sample

    Returns a figure object as "fig"
    """

    #///// RC Parameters //////#
    mpl.rcdefaults()
    mpl.rcParams['font.size'] = 20. # change the size of the font in every figure
    mpl.rcParams['font.family'] = 'Arial' # font Arial in every figure
    mpl.rcParams['axes.labelsize'] = 20.
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['xtick.direction'] = "in"
    mpl.rcParams['ytick.direction'] = "in"
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.major.width'] = 0.6
    mpl.rcParams['ytick.major.width'] = 0.6
    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
    mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                        # editing the text in illustrator

    #///// Create Figure //////#
    fig, axes = plt.subplots(2, 1, figsize = (6, 8)) # figsize is w x h in inch of figure
    fig.subplots_adjust(left = 0.19, right = 0.95, bottom = 0.1, top = 0.9,wspace = 0.3,hspace = 0.3)

    if (H >= 0) and measurement[0:4] == "NSYM": # only plot H >= 0 if NSYM needed
        axes[0].plot(Temperatures["T0_P"], Nernst["V_N_P"] * 1e6, ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")
        axes[1].plot(Temperatures["T0_P"], Nernst["nu_P"] * 1e3 / Temperatures["Tav_P"], ls = "--", lw = 3, marker = "o", c = "r", ms = 8, mfc = "r", mew = 0, mec = "r")

    if (H < 0) and measurement[0:4] == "NSYM": # only plot H < 0 if NSYM needed
        axes[0].plot(Temperatures["T0_N"], Nernst["V_N_N"] * 1e6, ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")
        axes[1].plot(Temperatures["T0_N"], Nernst["nu_N"] * 1e3 / Temperatures["Tav_N"], ls = "--", lw = 3, marker = "o", c = "#440154", ms = 8, mfc = "#440154", mew = 0, mec = "#440154")

    if measurement[0:3] == "SYM": # plots H>0, H<0 and SYM when needed
        axes[0].plot(Temperatures["T0"], Nernst["V_N"] * 1e6, ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")
        axes[1].plot(Temperatures["T0"], Nernst["nu"] * 1e3 / Temperatures["Tav"], ls = "--", lw = 3, marker = "o", c = "#A7A7A7", ms = 8, mfc  = "#A7A7A7", mew = 0, mec = "#A7A7A7")


    if min(Nernst["nu"]) * max(Nernst["nu"]) < 0: # choose if horizontal line y = 0  is needed /////#
        axes[0].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)
        axes[1].axhline(y = 0, ls = "--", c = "k", linewidth = 0.5, zorder = -1)

    #///// V_N //////#
    axes[0].set_xlabel(r"$T_0$ ( K )")
    axes[0].set_ylabel(r"$V_N$ ( $\mu$V )")
    axes[0].locator_params(nbins = 6) #gives the number of ticks maximum we want, here 6
    axes[0].set_xlim(xmin, xmax)

    #///// nu/T //////#
    axes[1].set_xlabel(r"$T_0$ ( K )")
    axes[1].set_ylabel(r"$\nu$ / $T$ ( nV / K$^2$ T )")
    axes[1].locator_params(nbins = 6)
    axes[1].set_xlim(xmin, xmax)

    fig.text(0.19, 0.95, "Seebeck", fontsize = 20)
    fig.text(0.955, 0.94, samplename, fontsize = 24, ha = "right", color = "#A7A7A7")

    f_manager = plt.get_current_fig_manager()
    f_manager.window.move(950, 75)

    return fig

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def save_TS_figure(samplelabel, axe, measurement, H, date, fig_list):

    figure_name = samplelabel + "_" + axe + "_" + measurement + \
             "_" + str(H) + 'T_' + str(date) +".pdf"

    file_figures = PdfPages("../figures/" + figure_name)

    for fig in fig_list[::-1]:
        file_figures.savefig(fig)

    file_figures.close()

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def save_FS_figure(samplelabel, axe, measurement, date, fig_list):

    figure_name = samplelabel + "_" + axe + "_" + measurement + \
                 "_" + str(date) +".pdf"

    file_figures = PdfPages("../figures/" + figure_name)

    for fig in fig_list[::-1]:
        file_figures.savefig(fig)

    file_figures.close()