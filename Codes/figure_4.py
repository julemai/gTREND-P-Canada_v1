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


# This script generates figure 4 of the following publication:
#      "Long-term pan-Canadian phosphorus mass budget dataset (g)TREND-P-CanadaCANADA: Consistent,
#      unified county-scale and gridded data products from 1961 to 2021"
#      Lamisa Malik, Juliane Mai, Danyka K. Byrnes, Kimberly J. Van Meter, Shuyu Chang, Meghan McLeod, Nandita B. Basu
#      Under review with Nature Scientific Data. DOI forthcoming.

# pyenv activate env-3.11.9
# run figure_4.py -g ../Figures/figure_4

if __name__ == '__main__':

    # -----------------------
    # add subolder scripts/lib to search path
    # -----------------------
    import sys
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(dir_path+'/lib')
    
    import argparse
    import numpy as np
    import pandas as pd
    import time
    import scipy
    
    from position      import position      # in lib/

    pngbase   = ''
    pdffile   = ''
    usetex    = False

    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Plots figure 3 of publication.''')
    parser.add_argument('-g', '--pngbase', action='store',
                        default=pngbase, dest='pngbase', metavar='pngbase',
                        help='Name basis for png output files (default: open screen window).')
    parser.add_argument('-p', '--pdffile', action='store',
                        default=pdffile, dest='pdffile', metavar='pdffile',
                        help='Name of pdf output file (default: open screen window).')
    parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")

    args      = parser.parse_args()
    pngbase   = args.pngbase
    pdffile   = args.pdffile
    usetex    = args.usetex
    dobw      = False
    donolabel = True

    if pdffile != '':
        outtype   = 'pdf'
    elif pngbase != '':
        outtype   = 'png'
    else:
        outtype   = 'X'

    del parser, args
    # Comment|Uncomment - End

    import numpy as np

    data_folder = "../Data/Figure_4/"
    
    census_years = np.arange(1961, 2019, 5)  # Wang goes until 2018
    
    provinces = {'ON': 'Ontario',
                 'QC': 'Quebec',
                 'BC': 'British Columbia',
                 'MB': 'Manitoba',
                 'SK': 'Saskatchewan',
                 'AB': 'Alberta',
                 'AP': 'Atlantic Provinces' # = NS+NB+PE+NL
                 #'NS': 'Nova Scotia',              # only availabe as AP in Wang and IPNI
                 #'NB': 'New Brunswick',            # only availabe as AP in Wang and IPNI
                 #'PE': 'Prince Edward Island',     # only availabe as AP in Wang and IPNI
                 #'NL': 'Newfoundland & Labrador',  # only availabe as AP in Wang and IPNI
                 }

    components = { 'fert':    {'name': 'Fertilizer P',        'color': (203/256.,201/256.,226/256.)},
                   'man':     {'name': 'Livestock Manure P',  'color': (253/256.,190/256.,133/256.)},
                   'crop':    {'name': 'Crop P Removal',      'color': (186/256.,228/256.,179/256.)},
                   'waste':   {'name': 'Domestic Waste P',    'color': (204/256.,204/256.,204/256.)},
                   'surplus': {'name': 'P Surplus',           'color': 'k'},
                 }

    sources = { 'WANG': { 'filename': data_folder+'IPNI_data.xlsx', 'title': 'Wang et al. (2022)'},
                'IPNI': { 'filename': data_folder+'IPNI_data.xlsx', 'title': 'IPNI (2014)'},
              }

    unit = "[ 1,000 metric tons-P yr$^{-1}$ ]"

    # ----------------------------------------------------------------------
    # read data
    # ----------------------------------------------------------------------

    data = {'TREND': {}, 'WANG': {}, 'IPNI': {}}

    # ------------------------
    # (1) read data: provincial-scale TREND-P-Canada
    # ------------------------
    for province in provinces:

        if province != 'AP':
            filename = "{}/{}_Kton_no-pasture-removal.csv".format(data_folder,province)
            tmp = pd.read_csv(filename)
            tmp = tmp.set_index('year')
        else:
            xprovinces = ['NS', 'NB', 'PE', 'NL']
            for ixprovince,xprovince in enumerate(xprovinces):
                filename = "{}/{}_Kton_no-pasture-removal.csv".format(data_folder,xprovince)
                itmp = pd.read_csv(filename)
                itmp = itmp.set_index('year')
                if ixprovince == 0:
                    tmp = itmp
                else:
                    tmp += itmp
                
        # select only relevant years
        tmp = tmp.loc[ [ ii for ii in census_years ] ]

        data['TREND'][province] = tmp.copy(deep=True)

    # ------------------------
    # (2) read data: Wang
    # ------------------------
    for province in provinces:

        filename = "{}/Wang_Paper_data.xlsx".format(data_folder)
        tmp = pd.read_excel(open(filename, 'rb'), sheet_name=province)
        tmp = tmp.rename(columns={'Year': 'year'})
        tmp = tmp.set_index('year')

        # select only census years
        tmp = tmp.loc[ [ ii for ii in census_years ] ]

        # generate components
        for component in components:

            if component == 'crop':
                # "Crop" and "Residue"
                tmp[component] = pd.Series(tmp['Crop']+tmp['Residue'], index=tmp.index) * -1.0
            elif component == 'man':
                # "Manure" and "Manure.1"/1000
                tmp[component] = pd.Series(tmp['Manure']+tmp['Manure.1']/1000., index=tmp.index)
            elif component == 'fert':
                # "Fertilizer"
                tmp[component] = pd.Series(tmp['Fertilizer'], index=tmp.index)
            elif component == 'waste':
                # "Sludge" * 5 since 20% of sludge goes to cropland and 80% to freshwater --> multiply with 5 to get 100% from 20%
                tmp[component] = pd.Series(tmp['Sludge']*5.0, index=tmp.index)
            elif component == 'surplus':
                # "Budget" and "Budget.1"/1000
                tmp[component] = pd.Series(tmp['Budget']+tmp['Budget.1']/1000., index=tmp.index)
            else:
                raise ValueError('Component {} not known.'.format(component))

        # remove all unneccesary columns
        tmp = tmp[list(components)]

        data['WANG'][province] = tmp.copy(deep=True)

    # ------------------------
    # (3) read data: IPNI
    # ------------------------
    
    # initialize with what we had for Wang as IPNI has less years, less provinces, and less components
    tmp_master = tmp
    for col in tmp_master.columns:
        tmp_master[col].values[:] = np.nan

    # get all sheet names
    filename = "{}/IPNI_data.xlsx".format(data_folder)
    avail_components = pd.ExcelFile(filename).sheet_names
    avail_provinces = pd.read_excel(open(filename, 'rb'), sheet_name=avail_components[0]).rename(columns={'Year': 'year'}).set_index('year').columns

    for province in provinces:

        tmp_prov = tmp_master.copy(deep=True)
        if province in avail_provinces:
        
            for component in components:

                if component in avail_components:

                    filename = "{}/IPNI_data.xlsx".format(data_folder)
                    tmp = pd.read_excel(open(filename, 'rb'), sheet_name=component)
                    tmp = tmp.rename(columns={'Year': 'year'})
                    tmp = tmp.set_index('year')

                    if province in tmp.columns:
                        
                        # data for available years for current province and component
                        # 1. conversion from phosphate (P2O5) to phosphorus --> 0.437
                        # 2. conversion from tons to kilo-tons
                        tmp = (tmp.loc[ [ ii for ii in census_years if ii in tmp.index] ][province] * 0.437) / 1000.

                        # collect data in provincial dataframe
                        tmp_prov[component] = tmp

        data['IPNI'][province] = tmp_prov


    # ----------------------------------------------------------------------
    # plot data
    # ----------------------------------------------------------------------

    # Main plot
    nrow        = 7           # # of rows of subplots per figure
    ncol        = 5           # # of columns of subplots per figure
    hspace      = 0.04         # x-space between subplots
    vspace      = 0.04        # y-space between subplots
    right       = 0.9         # right space on page
    textsize    = 6           # standard text size
    dxabc       = 1.0         # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
    # dyabc       = -13       # % of (max-min) shift up from lower x-axis for a,b,c,... labels
    dyabc       = 0.0         # % of (max-min) shift up from lower x-axis for a,b,c,... labels
    dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for signature
    dysig       = -0.05       # % of (max-min) shift up from lower x-axis for signature
    dxtit       = 0           # % of (max-min) shift to the right from left y-axis for title
    dytit       = 1.3         # % of (max-min) shift up from lower x-axis for title

    lwidth      = 0.25         # linewidth
    elwidth     = 1.0         # errorbar line width
    alwidth     = 0.5         # axis line width
    glwidth     = 0.5         # grid line width
    msize       = 3.0         # marker size
    mwidth      = 1.0         # marker edge width
    mcol1       = 'blue'      # primary marker colour
    mcol2       = 'red'       # secondary
    mcol3       = 'red'       # third
    mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0']
    lcol1       = 'blue'   # primary line colour
    lcol2       = '0.0'
    lcol3       = '0.0'
    lcols       = ['None', 'None', 'None', 'None', 'None', '0.0']
    hatches     = [None, None, None, None, None, '//']

    # Legend
    llxbbox     = 0.0         # x-anchor legend bounding box
    llybbox     = 1.0         # y-anchor legend bounding box
    llrspace    = 0.          # spacing between rows in legend
    llcspace    = 1.0         # spacing between columns in legend
    llhtextpad  = 0.4         # the pad between the legend handle and text
    llhlength   = 1.5         # the length of the legend handles
    frameon     = False       # if True, draw a frame around the legend. If None, use rc

    # PNG
    dpi         = 600         # 150 for testing
    transparent = True
    bbox_inches = 'tight'
    pad_inches  = 0.035

    # Clock options
    ymax = 0.6
    ntextsize   = 'medium'       # normal textsize
    # modules
    bmod        = 0.5            # fraction of ymax from center to start module colours
    alphamod    = 0.7            # alpha channel for modules
    fwm         = 0.05           # module width to remove at sides
    ylabel1     = 1.15           # position of module names
    ylabel2     = 1.35           # position of class names
    mtextsize   = 'large'        # 1.3*textsize # textsize of module labels
    # bars
    bpar        = 0.4            # fraction of ymax from center to start with parameter bars
    fwb         = [0.7,0.4,0.3]  # width of bars
    plwidth     = 0.5
    # parameters in centre
    bplabel     = 0.1            # fractional distance of ymax of param numbers in centre from 0-line
    ptextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of param numbers in centre
    # yaxis
    space4yaxis = 2              # space for y-axis (integer)
    ytextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of y-axis
    sig         = 'J Mai' # sign the plot

    import matplotlib as mpl
    import matplotlib.patches as patches
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    
    plt.close('all')
    mpl.use('TkAgg')
    
    if (outtype == 'pdf'):
        mpl.use('PDF') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        # Customize: http://matplotlib.sourceforge.net/users/customizing.html
        mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
        mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
        if usetex:
            mpl.rc('text', usetex=True)
        #else:
        #    mpl.rc('font',**{'family':'serif','serif':['Palatino']})
            # mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
            # mpl.rc('font',**{'family':'serif','serif':['times']})
        mpl.rc('text.latex') #, unicode=True)
    elif (outtype == 'png'):
        mpl.use('Agg') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
        if usetex:
            mpl.rc('text', usetex=True)
        else:
            #mpl.rc('font',**{'family':'serif','serif':['Palatino']})
            #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Barlow Condensed']}) #['Barlow Condensed,Barlow Condensed Light']}) #Barlow Condensed']})
            #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
            print('None')
        mpl.rc('text.latex') #, unicode=True)
        mpl.rc('savefig', dpi=dpi, format='png')
    else:
        import matplotlib.pyplot as plt
        mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
    mpl.rc('font', size=textsize)
    mpl.rc('lines', linewidth=lwidth, color='black')
    mpl.rc('axes', linewidth=alwidth, labelcolor='black')
    mpl.rc('path', simplify=False) # do not remove

    from matplotlib.patches import Rectangle, Circle, Polygon
    from mpl_toolkits.basemap import Basemap

    if (outtype == 'pdf'):
        # pdffile = plotname+'.pdf'
        print('Plot PDF ', pdffile)
        pdf_pages = PdfPages(pdffile)
    elif (outtype == 'png'):
        print('Plot PNG ', pngbase)
    else:
        print('Plot X')

    t1  = time.time()
    ifig = 0

    figsize = mpl.rcParams['figure.figsize']
    mpl.rcParams['axes.linewidth'] = lwidth


    ifig = 0

    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
        
    fig = plt.figure(ifig)

    for isource,source in enumerate(sources):

        for component in components:

            iplot += 1

            #     [left, bottom, width, height]
            pos = position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)
        
            sub = fig.add_axes(pos,frameon=True) #, axisbg='none')

            # get data points and plot
            minval =  9999.
            maxval = -9999.
            if component == 'fert':
                minval = 0.0
                maxval = 200.0
            elif component == 'man':
                minval = 0.0
                maxval = 100.0
            elif component == 'crop':
                minval = 0.0
                maxval = 220.0
            elif component == 'waste':
                minval = 0.0
                maxval = 15.0
            elif component == 'surplus':
                minval = -52.0
                maxval = 100.0
            xvals_all = []
            yvals_all = []
            for province in provinces:

                xvals  = np.array(data['TREND'][province][component].values)
                yvals  = np.array(data[source][province][component].values)

                xvals_all += list(xvals)
                yvals_all += list(yvals)

                # iminval = min(np.min(xvals),np.min(yvals))
                # imaxval = max(np.max(xvals),np.max(yvals))
                # if minval > iminval:
                #     minval = iminval
                # if maxval < imaxval:
                #     maxval = imaxval

                sub.plot(xvals, yvals,
                             linestyle='None',
                             alpha=0.3,
                             marker='o', markeredgecolor='0.1', markerfacecolor='k', 
                             markersize=3.0, markeredgewidth=0.0)

            # correlation coefficient
            idx = (~np.isnan(xvals_all) & ~np.isnan(yvals_all) )
            r = np.ma.corrcoef(np.ma.masked_invalid(xvals_all),np.ma.masked_invalid(yvals_all))[0,1]
            if np.any(idx):
                sub.text(0.95,0.05,"{} = {:3.2f}".format("r",r),
                             verticalalignment='bottom',horizontalalignment='right',rotation=0,
                             fontsize=textsize,
                             transform=sub.transAxes)
            else:
                sub.text(0.5,0.5,"n/a",
                             verticalalignment='center',horizontalalignment='center',rotation=0,
                             fontsize=textsize,
                             transform=sub.transAxes)

            # plot 1-1 line
            if np.any(idx):
                sub.plot([minval,maxval], [minval,maxval], linestyle='-', linewidth=lwidth*3, color='k')
        
            # add label with component name
            sub.text(0.0,1.0,"{}.".format(chr(96+iplot)),
                             verticalalignment='bottom',horizontalalignment='left',rotation=0,
                             fontweight='bold',
                             fontsize=textsize+2,
                             transform=sub.transAxes)
            sub.text(0.05,0.95,"{}".format(components[component]['name'].replace('Crop P Removal','Crop P\nRemoval').replace('Livestock Manure','Livestock\nManure').replace('Domestic Waste','Domestic\nWaste')),
                             verticalalignment='top',horizontalalignment='left',rotation=0,
                             fontsize=textsize,
                             transform=sub.transAxes)

            # Add x-title
            if iplot == 8:
                sub.text(0.5,-0.28,"TREND-P-Canada {}".format(unit),
                             verticalalignment='top',horizontalalignment='center',rotation=0,
                             fontsize=textsize,
                             transform=sub.transAxes)
                
            # Add y-title
            if iplot == 1 or iplot == 6:
                sub.text(-0.5,0.5,"{}\n{}".format(sources[source]['title'],unit),  # first value position is moving right-left
                             rotation=90,
                             fontsize=textsize,
                             verticalalignment='center',horizontalalignment='center',
                             transform=sub.transAxes)

            # set ticks
            if np.any(idx):
                xticks = sub.get_yticks()
                sub.set_xticks(xticks)
                sub.set_yticks(xticks)
            else:
                sub.set_xticks([])
                sub.set_yticks([])

            sub.set_xlim([minval*0.98,maxval*1.02])
            sub.set_ylim([minval*0.98,maxval*1.02])


    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+".png"
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
    strin = '[m]: '+astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+str(t2-t1)
    print('Time ', strin)






    

