#!/usr/bin/env python
## License

# Copyright 2025 Juliane Mai - juliane.mai@uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.

# The code creating this figure is based on Lamisa Malik's initial version coded in R with each
# subpanel being a separate figure which needed Canvas (or similar) to arrange them all in a
# single figure (including adjustemnts in labels and legends). 


# This script generates figure 2 of the following publication:
#      "Consistent and unified phosphorus mass budget dataset across Canada from 1961 to 2021: gTREND-P-Canada"
#      Lamisa Malik, Juliane Mai, Danyka K. Byrnes, Kimberly J. Van Meter, Shuyu Chang, Meghan McLeod, Nandita B. Basu
#      Under review with Nature Scientific Data. DOI forthcoming.

# pyenv activate env-3.11.9
# python figure_2.py -p ../Figures/figure_2 -t png -y '1961,1991,2021' -x 'area' -u

# -------------------------------------------------------------------------
# General settings
#
dobw      = False  # True: black & white
docomp    = True   # True: Print classification on top of modules
dosig     = False  # True: add signature to plot
dolegend  = False  # True: add legend to each subplot
doabc     = True   # True: add subpanel numbering
dotitle   = True   # True: add catchment titles to subpanels

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse
import numpy as np


plotname    = ''
outtype     = ''
usetex      = False
serif       = False
years       = ['1961', '2001', '2021']
hist_type   = 'area'  # 'county'

parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Plots results of hydrologic calibration experiments w/ and w/o screening of parameters.''')
parser.add_argument('-p', '--plotname', action='store',
                    default=plotname, dest='plotname', metavar='plotname',
                    help='Name of plot output file for types pdf, html or d3, '
                    'and name basis for type png (default: '+__file__[0:__file__.rfind(".")]+').')
parser.add_argument('-s', '--serif', action='store_true', default=serif, dest="serif",
                help="Use serif font; default sans serif.")
parser.add_argument('-t', '--type', action='store',
                    default=outtype, dest='outtype', metavar='outtype',
                    help='Output type is pdf, png, html, or d3 (default: open screen windows).')
parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf, png and html.")
parser.add_argument('-y', '--years', action='store', default=years, dest='years',
                help="Years to plot data for. Comma-separated. E.g., '1961,2021'")
parser.add_argument('-x', '--hist_type', action='store', default=hist_type, dest='hist_type',
                help="Histogram type. Either county count ('county') or by area ('area').")

args      = parser.parse_args()
plotname  = args.plotname
outtype   = args.outtype
serif     = args.serif
usetex    = args.usetex
years     = [ yy for yy in args.years.split(',') ]
hist_type = args.hist_type



# Check input
outtype = outtype.lower()
outtypes = ['', 'pdf', 'png', 'html', 'd3']
if outtype not in outtypes:
    print('\nError: output type must be in ', outtypes)
    import sys
    sys.exit()

dowhite    = False  # True: black background, False: white background
title      = False   # True: title on plots, False: no plot titles
textbox    = False  # if true: additional information is set as text box within plot

if dowhite:
    fgcolor = 'white'
    bgcolor = 'black'
else:
    fgcolor = 'black'
    bgcolor = 'white'

del parser, args

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import color                            # in lib/
from position      import position      # in lib/
from abc2plot      import abc2plot      # in lib/
from brewer        import get_brewer    # in lib/
from autostring    import astr          # in lib/
from str2tex       import str2tex       # in lib/
from fsread        import fsread        # in lib/

# import fiona          # some shapefile coordinate stuff
import numpy as np
import xarray as xr
import pandas as pd
import csv
import fiona
import copy                       # deep copy objects, arrays etc
import time
import os
import shapefile
import glob
import geopandas as gpd
t1 = time.time()


def log_base(x,b):
    return np.log10(x)/np.log10(b)

catchfile_shp_ca = [ os.path.join('..','Data', 'county_shapefiles', str(year), str(year)) for year in years ]

# -------------------------
# (a) Outline of CA counties
# -------------------------

coord_catch_ca = {}
ID = {}
areas = {}

for iyear,year in enumerate(years):
    icoord_catch_ca = []
    iID = []
    iareas = {}

    with fiona.open(catchfile_shp_ca[iyear]+'.shp') as src:
        for ii in range(len(src)):
            icoord_catch_ca.append(src[ii]['geometry']['coordinates'])
            iID.append(str(src[ii]['properties']['ID']))

    data = gpd.read_file(catchfile_shp_ca[iyear]+'.shp')
    data_copy = data.copy()
    data_copy = data_copy.to_crs({'proj':'cea'})
    data['poly_area'] = data_copy['geometry'].area/ 10**6 # in km2
    for index, row in data.iterrows():
        iareas[str(row['ID'])] = row['poly_area']

    iID_uniq = np.unique(iID)
    if len(iID) != len(iID_uniq):
        print("There are duplicate IDs in {}!!! This is really risky.".format(catchfile_shp_ca[iyear]))

    coord_catch_ca[year] = icoord_catch_ca
    ID[year] = iID
    areas[year] = iareas

print("")

# -------------------------
# (c) P components of that year
# -------------------------

components        = [ 'fertilizer_P', 'crop_and_pasture_P_rem.','livestock_manure_P', 'domestic_waste_P', 'P_surplus'  ]
components_lamisa = [ 'FERT',         'CROP',                      'LVSK',               'WASTE',            'Surplus'    ]

val_components = {}
minmax_components = {}
for icomponent,component in enumerate(components):

    filename = os.path.join('..', 'Data', 'County_TREND-P-Canada', 'P_'+components_lamisa[icomponent]+'.csv')

    minval = 9999.
    maxval = -9999.

    ff = open(filename, "r")
    content = ff.readlines()

    ids = content[0].strip().split(',')[1:]
    
    iyears = []
    for line in content[1:]:
        iyears.append( line.strip().split(',')[0] )

    # initialize all IDs that exist in shapefile with NaN
    vals = { yy: { ii: np.nan for ii in ID[yy] } for yy in years }

    # now set values of IDs that exist in CSV with according values
    for line in content[1:]:
        ll = line.strip().split(',')
        year = ll[0]

        if year in years:
            vv   = ll[1:]

            for iii,ii in enumerate(ids):
                if vv[iii] == 'NAN' or vv[iii] == 'NA':
                    vals[year][ii] = np.nan
                else:
                    vals[year][ii] = float(vv[iii])
                    if vals[year][ii] < minval: minval = vals[year][ii]
                    if vals[year][ii] > maxval: maxval = vals[year][ii]

    val_components[component] = vals
    minmax_components[component] = {'min': minval, 'max': maxval}

# -------------------------------------------------------------------------
# Customize plots
#

# -------------------------------------------------------------------------
# Plotting of results
# -------------------------------------------------------------------------
# Main plot
ncol        = 3           # number columns
nrow        = 5           # number of rows
textsize    = 10          # standard text size
textsize2   = textsize-4  # standard text size
dxabc       = 0.95          # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.92          # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for signature
dysig       = -0.075      # % of (max-min) shift up from lower x-axis for signature
dxtit       = 0           # % of (max-min) shift to the right from left y-axis for title
dytit       = 1.2         # % of (max-min) shift up from lower x-axis for title
hspace      = 0.07        # x-space between subplots
vspace      = 0.15        # y-space between subplots

lwidth      = 0.2         # linewidth
elwidth     = 0.5         # errorbar line width
alwidth     = 1.0         # axis line width
glwidth     = 0.5         # grid line width
msize       = 3.0         # marker size
mwidth      = 0.0         # marker edge width
mcol1       = '0.7'       # primary marker colour
mcol2       = '0.0'       # secondary
mcol3       = '0.0'       # third
#mcols       = color.colours(['blue','green','yellow','orange','red','darkgray','darkblue','black','darkgreen','gray'])
#lcol0       = color.colours('black')    # line colour
#lcol1       = color.colours('blue')     # line colour
#lcol2       = color.colours('red')    # line colour
#lcol3       = color.colours('darkgreen')   # line colour
#lcols       = color.colours(['darkgray','gray','blue'])
markers     = ['o','v','s','^']

# Legend
llxbbox     = 0.5        # x-anchor legend bounding box
llybbox     = -0.6        # y-anchor legend bounding box
llrspace    = 0.          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
llhtextpad  = 0.4         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

dpi = 600
transparent = True
bbox_inches = 'tight'
pad_inches  = 0.035

import matplotlib as mpl
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import matplotlib.dates as dates
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt

plt.close('all')
mpl.use('TkAgg')

if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    mpl.rc('figure', figsize=(10.97,7.0)) # adjusted # figsize=(10.97,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
        if not serif:
            #   r'\usepackage{helvet}',                             # use Helvetica
            mpl.rcParams['text.latex.preamble'] = ( r'\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}' + # use MyriadPro font
                                                    r'\renewcommand{\familydefault}{\sfdefault}' +                   # normal text font is sans serif
                                                    r'\figureversion{lining,tabular}' +
                                                    r'\usepackage{wasysym}' )                                        # for permil symbol (load after MyriadPro)
        else:
            mpl.rcParams['text.latex.preamble'] = r'\usepackage{wasysym}'                     # for permil symbol
    else:
        # mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        # mpl.rc('font',**{'family':'serif','serif':['times']})
        mpl.rcParams['font.family'] = 'DejaVu Sans'
    mpl.rc('text.latex') #, unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(10.97,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
        if not serif:
            #   r'\usepackage{helvet}',                             # use Helvetica
            mpl.rcParams['text.latex.preamble'] = ( r'\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}' + # use MyriadPro font
                                                    r'\renewcommand{\familydefault}{\sfdefault}' +                   # normal text font is sans serif
                                                    r'\figureversion{lining,tabular}' +
                                                    r'\usepackage{wasysym}' )                                        # for permil symbol (load after MyriadPro)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex') #, unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
mpl.rc('text.latex') #, unicode=True)
mpl.rc('font', size=textsize)
mpl.rc('path', simplify=False) # do not remove
# print(mpl.rcParams)
mpl.rc('axes', linewidth=alwidth, edgecolor=fgcolor, facecolor=bgcolor, labelcolor=fgcolor)
mpl.rc('figure', edgecolor=bgcolor, facecolor='grey')
mpl.rc('grid', color=fgcolor)
mpl.rc('lines', linewidth=lwidth, color=fgcolor)
mpl.rc('patch', edgecolor=fgcolor)
mpl.rc('savefig', edgecolor=bgcolor, facecolor=bgcolor)
mpl.rc('patch', edgecolor=fgcolor)
mpl.rc('text', color=fgcolor)
mpl.rc('xtick', color=fgcolor)
mpl.rc('ytick', color=fgcolor)

if (outtype == 'pdf'):
    pdffile = plotname+'.pdf'
    print('Plot PDF ', pdffile)
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', plotname)
else:
    print('Plot X')

t1  = time.time()
ifig = 0

figsize = mpl.rcParams['figure.figsize']
mpl.rcParams['axes.linewidth'] = lwidth



from matplotlib.patches import Rectangle, Circle, Polygon
from mpl_toolkits.basemap import Basemap


# Read catchments as shape files
#  --> during plotting

if (catchfile_shp_ca != ''):
    # catchfile_shp_ca  = catchfile_shp_ca.split(',')
    # catchfile_shp_ca  = [ i.replace(' ','') for i in catchfile_shp_ca ]
    ncatchfile_shp_ca = len(catchfile_shp_ca)
else:
    ncatchfile_shp_ca = 0

# colors
if dobw:
    c = np.linspace(0.2, 0.85, nmod)
    c = np.ones(nmod)*0.7
    c = [ str(i) for i in c ]
    ocean_color = '0.1'
else:
    # c = [(165./255.,  0./255., 38./255.), # interception
    #      (215./255., 48./255., 39./255.), # snow
    #      (244./255.,109./255., 67./255.), # soil moisture
    #      (244./255.,109./255., 67./255.), # soil moisture
    #      (253./255.,174./255., 97./255.), # direct runoff
    #      (254./255.,224./255.,144./255.), # Evapotranspiration
    #      (171./255.,217./255.,233./255.), # interflow
    #      (116./255.,173./255.,209./255.), # percolation
    #      ( 69./255.,117./255.,180./255.), # routing
    #      ( 49./255., 54./255.,149./255.)] # geology
    c = get_brewer('rdylbu11', rgb=True)
    tmp = c.pop(5)   # rm yellow
    #c.insert(2,c[2]) # same colour for both soil moistures
    ocean_color = (151/256., 183/256., 224/256.)

    cc = color.get_brewer('dark_rainbow_256', rgb=True)
    cc = color.get_brewer('temp_19lev', rgb=True)[2:]
    #cc = cc[::-1] # reverse colors
    cmap = mpl.colors.ListedColormap(cc)

cc = {}

# green
cc['crop_and_pasture_P_rem.'] = [        (237/256.,248/256.,233/256.),
                        (199/256.,233/256.,192/256.),
                        (161/256.,217/256.,155/256.),
                        (116/256.,196/256.,118/256.),
                        (49/256.,163/256.,84/256.),
                        (0/256.,109/256.,44/256.)]
# purple
cc['fertilizer_P'] = [    (242/256.,240/256.,247/256.),
                        (218/256.,218/256.,235/256.),
                        (188/256.,189/256.,220/256.),
                        (158/256.,154/256.,200/256.),
                        (117/256.,107/256.,177/256.),
                        (84/256.,39/256.,143/256.)]
# orange
cc['livestock_manure_P'] = [     (254/256.,237/256.,222/256.),
                        (253/256.,208/256.,162/256.),
                        (253/256.,174/256.,107/256.),
                        (253/256.,141/256.,60/256.),
                        (230/256.,85/256.,13/256.),        
                        (166/256.,54/256.,3/256.)]
# gray
cc['domestic_waste_P'] = [   (247/256.,247/256.,247/256.),
                        (217/256.,217/256.,217/256.),
                        (189/256.,189/256.,189/256.),
                        (150/256.,150/256.,150/256.),
                        (99/256.,99/256.,99/256.),
                        (37/256.,37/256.,37/256.)]

# blue - yellow - red
cc['P_surplus'] = [     (116/256.,173/256.,209/256.),
                        (255/256.,255/256.,191/256.),   # light red
                        (254/256.,224/256.,144/256.),  
                        (253/256.,174/256.,97/256.),
                        (244/256.,109/256.,67/256.),
                        (215/256.,48/256.,39/256.),
                        (165/256.,0/256.,38/256.)]      # dark red

    
# -------------------------------------------------------------------------
# Plot
#

if (outtype == 'pdf'):
    print('Plot PDF ', plotname)
    pdf_pages = PdfPages(plotname)
elif (outtype == 'png'):
    print('Plot PNG ', plotname)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']

ifig = 0


# -------------------------------------------------------------------------
# Fig 1
# -------------------------------------------------------------------------
#
ifig += 1
iplot = 0
fig = plt.figure(ifig)
        

# -----------------------------------------
# Plot shapes
# -----------------------------------------


for icomponent,component in enumerate(components):

    for iyear,year in enumerate(years):

        if iyear == 0:          
            iplot = icomponent*ncol

        iplot += 1

        print('')
        #print('iplot=',iplot)

        # [left, bottom, width, height]
        pos = position(nrow,ncol,iplot,hspace=hspace,vspace=vspace) + [0.02*iyear,0.07*icomponent,0.15,0.11]
        #print('pos=',pos)
        sub = fig.add_axes(pos) #, axisbg='none')

        # Map: Canada
        llcrnrlon =  -126.
        urcrnrlon =   -30.   # maybe -30 (more space lower right)
        llcrnrlat =   34.
        urcrnrlat =   62.
        lat_1     =   45.   #(llcrnrlat+urcrnrlat)/2.  # first  "equator"
        lat_2     =   55.   #(llcrnrlat+urcrnrlat)/2.  # second "equator"
        lat_0     =   50.   #(llcrnrlat+urcrnrlat)/2.  # center of the map
        lon_0     =   -90.  #(llcrnrlon+urcrnrlon)/2.  # center of the map
        map4 = Basemap(projection='lcc',
                    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                    lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
                    resolution='i',
                    area_thresh=5000) # Lambert conformal
            
        map4.drawmapboundary(fill_color='0.9', linewidth=0.0)
        map4.drawmapboundary(color=None)


        # add county shape to plot
        nshapes = len(coord_catch_ca[year])

        base = 1.2
        thres = 0.001

        # color map bounds
        if 'surplus' in component:
            delta = (minmax_components[component]['max']-minmax_components[component]['min'])/len(cc[component])
            bounds = np.array([ minmax_components[component]['min'] + ii*delta for ii in range(len(cc[component])+1) ])

            delta = (np.log10(minmax_components[component]['max'])-np.log10(max(thres,minmax_components[component]['min'])))/int(len(cc[component])-1)
            bounds = np.array([ 10**(np.log10(max(thres,minmax_components[component]['min'])) + ii*delta) for ii in range(int(len(cc[component]))-1+1) ])

            delta = (log_base(minmax_components[component]['max'],base)-log_base(max(thres,minmax_components[component]['min']),base))/int(len(cc[component])-1)
            bounds = np.array([ base**(log_base(max(thres,minmax_components[component]['min']),base) + ii*delta) for ii in range(int(len(cc[component]))-1+1) ])

            bounds[0] = 0.0
            bounds = np.append( np.array([minmax_components[component]['min']]),bounds )

        else:
            
            # log
            delta = (np.log10(minmax_components[component]['max'])-np.log10(max(thres,minmax_components[component]['min'])))/len(cc[component])
            bounds = np.array([ 10**(np.log10(max(thres,minmax_components[component]['min'])) + ii*delta) for ii in range(len(cc[component])+1) ])

            delta = (log_base(minmax_components[component]['max'],base)-log_base(max(thres,minmax_components[component]['min']),base))/len(cc[component])
            bounds = np.array([ base**(log_base(max(thres,minmax_components[component]['min']),base) + ii*delta) for ii in range(len(cc[component])+1) ])

            bounds[0] = 0.0

        ncounty_with_data = 0
        histogram_data = np.array(bounds[0:-1])*0.0
        for iishape,ishape in enumerate(np.arange(nshapes)):
            xy = coord_catch_ca[year][ishape]

            nsubshapes = len(xy)
            set_color = False
            for inshape in np.arange(nsubshapes):

                ixy = xy[inshape]
                try:
                    nn = np.shape(ixy)
                except:
                    ixy = xy[inshape][0]

                if np.isnan(val_components[component][year][ID[year][iishape]]):
                    icolor = None
                else:
                    idx = np.sum(val_components[component][year][ID[year][iishape]] > bounds) - 1
                    idx = np.min([np.max([0,idx]),len(cc[component])-1])
                    icolor = cc[component][idx]
                    set_color = True

                if not(icolor is None):
                    if (len(np.shape(ixy)) == 2) and (np.shape(ixy)[1] == 2):
                        # shape is (nn,2)
                        sub.add_patch(Polygon([ map4(ii[0],ii[1]) for ii in ixy ], facecolor=icolor, edgecolor='0.2', linewidth=0.1, zorder=300))
                    elif (len(np.shape(ixy)) == 3) and (np.shape(ixy)[0] == 1):
                        # shape is (1,nn,2)
                        sub.add_patch(Polygon([ map4(ii[0],ii[1]) for ii in ixy[0] ], facecolor=icolor, edgecolor='0.2', linewidth=0.1, zorder=300))
                    else:
                        print('ishape  =',ishape)
                        print('inshape =',inshape)
                        print('len(np.shape(ixy)) =',len(np.shape(ixy)))
                        print('np.shape(ixy)      =',np.shape(ixy))
                        raise ValueError('CA: Dont know what to do! #1')

            if set_color:
                ncounty_with_data += 1
                if hist_type == 'county':
                    histogram_data[idx] += 1
                elif hist_type == 'area':
                    histogram_data[idx] += areas[year][ID[year][iishape]]
                else:
                    raise ValueError('Histogram type not implemented yet.')

        # number of counties
        if icomponent == len(components)-1:
            sub.text(0.2,0.2, str2tex('N$_\mathrm{county}$ = '+str(ncounty_with_data), usetex=usetex),
                     fontsize=textsize2, horizontalalignment="left", verticalalignment="bottom",transform=sub.transAxes,
                     rotation=-17)

        if icomponent == 0:
            sub.text(0.5,0.8, str2tex(year, usetex=usetex),
                     fontsize=textsize, horizontalalignment="center", verticalalignment="bottom",transform=sub.transAxes,)

        if iyear == 0:
            sub.text(0.0,0.5, str2tex(' '.join(component.split('_')).title().replace('And','and'), usetex=usetex),
                     fontsize=textsize, horizontalalignment="right", verticalalignment="center",transform=sub.transAxes,rotation=90)
            sub.text(0.05,0.5, str2tex('in [kg-P ha$^{-1}$ yr$^{-1}$]',usetex=usetex),
                     fontsize=textsize2, horizontalalignment="right", verticalalignment="center",transform=sub.transAxes,rotation=90)

        print("year = {} component = {} --> ncounty_with_data = {} of {}".format(year,component,ncounty_with_data,nshapes))


        # add abc
        sub.text(0.9,0.65,"{}.".format(chr(96+iplot)),
                         verticalalignment='bottom',horizontalalignment='left',rotation=0,
                         fontweight='bold',
                         fontsize=textsize,
                         transform=sub.transAxes)
            

        # [left, bottom, width, height]
        if True: #iyear == len(years)-1:
            pos_legend = pos+[0.326,0.04,-0.340,-0.085]
            #print('pos legend = ',pos_legend)
            sub  = fig.add_axes(pos_legend,frameon=False) #, left=0.26, right=0.86, top=0.422, bottom=0.382))
            sub2 = sub.twinx()

            # no frame
            sub2.spines['top'].set_visible(False)
            sub2.spines['right'].set_visible(False)
            sub2.spines['bottom'].set_visible(False)
            sub2.spines['left'].set_visible(False)

            xvals = np.linspace(bounds[0],bounds[-1],num=len(cc[component])+1)[0:-1]
            delta = (xvals[1]-xvals[0])/2.
            xvals += delta
            yvals = histogram_data / np.sum(histogram_data) * 100. # in %
            sub2.barh(xvals,yvals, align='center', height=1.8*delta, color=cc[component])
            #xvals = np.arange(len(bounds)+1)

            sub.set_xticks([])
            sub.set_yticks([])
            sub.tick_params(width = 0.0)

            sub2.set_xlim(sub2.get_xlim()[::-1])
            if iyear == len(years)-1:
                sub2.tick_params(labelsize=textsize2)
                sub2.set_yticks(np.linspace(bounds[0],bounds[-1],num=len(cc[component])+1))
                sub2.set_yticklabels( [ str2tex(astr(bb,prec=3),usetex=usetex) for bb in bounds ] )
                sub2.tick_params(width = lwidth, size = 1.5 )
            else:
                sub2.set_yticks([])
            

            #xvals = np.linspace(bounds[0],bounds[-1],num=len(cc[component])+1)
            for ival, val in enumerate(yvals):
                sub2.text(yvals[ival],xvals[ival], str2tex(astr(yvals[ival],prec=0)+'%',usetex=usetex),
                     fontsize=textsize2, horizontalalignment="right", verticalalignment="center",rotation=0)

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = plotname+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)

# -------------------------------------------------------------------------
# Finished
#

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()

t2    = time.time()
strin = '[m]: '+astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+astr(t2-t1,0)
print('Time ', strin)




    
