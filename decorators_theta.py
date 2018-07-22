# -*- coding: Latin-1 -*-

## Modules <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
from numpy import cos, sin
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def theta_TS_figure(theta):

    def decorator(function_in):

        def function_modified(H, **parameters):

            ## Field components ///////////////////////////////////////////////#

            H_para = H * sin(theta) # H // to (ab) plane
            H_perp = H * cos(theta) # H perp to (ab) plane

            fig = function_in(H = H, **parameters)

            fig.text(0.82,0.73, r"$\theta$   = ", ha = 'left', fontsize = 20)
            fig.text(0.895,0.73, "{0:g}".format(theta), ha="left", fontsize = 20)

            fig.text(0.82,0.66, r"$H_{||}$ = ", ha = 'left', fontsize = 20)
            fig.text(0.895,0.66, "{0:g}".format(H_para) + " T", ha="left", fontsize = 20)

            fig.text(0.82,0.59, r"$H_{\perp}$= ", ha = 'left', fontsize = 20)
            fig.text(0.895,0.59, "{0:g}".format(H_perp) + " T", ha="left", fontsize = 20)

            return fig
        return function_modified
    return decorator


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

def theta_FS_figure(theta):

    def decorator(function_in):

        def function_modified(x, y, **parameters):

            ## Field components ///////////////////////////////////////////////#

            fig = function_in(x, y, **parameters)

            props = dict(boxstyle='square', facecolor='#EBEBEB', linewidth = 0)
            fig.text(0.81,0.04, r"$\theta$ = " + "{0:g}".format(theta), ha = 'right', fontsize = 24, bbox=props)

            return fig
        return function_modified
    return decorator