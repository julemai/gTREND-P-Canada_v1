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


# This script generates figure 3 of the following publication:
#      "Consistent and unified phosphorus mass budget dataset across Canada from 1961 to 2021: gTREND-P-Canada"
#      Lamisa Malik, Juliane Mai, Danyka K. Byrnes, Kimberly J. Van Meter, Shuyu Chang, Meghan McLeod, Nandita B. Basu
#      Under review with Nature Scientific Data. DOI forthcoming.

# pyenv activate env-3.11.9
# run figure_3.py -g ../Figures/figure_3

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

    data_folder = "../Data/Figure_3/"
    
    census_years = np.arange(1961, 2022, 5)
    
    provinces = {'ON': 'Ontario',
                 'QC': 'Quebec',
                 'BC': 'British Columbia',
                 'MB': 'Manitoba',
                 'SK': 'Saskatchewan',
                 'AB': 'Alberta',
                 'NS': 'Nova Scotia',
                 'NB': 'New Brunswick',
                 'PE': 'Prince Edward Island',
                 'NL': 'Newfoundland & Labrador',                 
                 }

    components = { 'crop':    {'name': 'Crop and Pasture P Removal', 'color': (186/256.,228/256.,179/256.)},
                    'man':     {'name': 'Livestock Manure P',         'color': (253/256.,190/256.,133/256.)},
                   'fert':    {'name': 'Fertilizer P',               'color': (203/256.,201/256.,226/256.)},
                   'waste':   {'name': 'Domestic Waste P',           'color': (204/256.,204/256.,204/256.)},
                   'surplus': {'name': 'P Surplus',                  'color': 'k'},
                 }

    unit = "[ 1,000 metric tons-P yr$^{-1}$ ]"

    # ----------------------------------------------------------------------
    # read data
    # ----------------------------------------------------------------------
    
    data = {'TREND': {}}
    for province in provinces:

        filename = "{}/{}_Kton.csv".format(data_folder,province)
        tmp = pd.read_csv(filename)
        tmp = tmp.set_index('year')

        data['TREND'][province] = tmp


    # ----------------------------------------------------------------------
    # plot data
    # ----------------------------------------------------------------------

    # Main plot
    nrow        = 5           # # of rows of subplots per figure
    ncol        = 4           # # of columns of subplots per figure
    hspace      = 0.06         # x-space between subplots
    vspace      = 0.04        # y-space between subplots
    right       = 0.9         # right space on page
    textsize    = 9           # standard text size
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
    iiplot = 0
    print('Plot - Fig ', ifig)
        
    fig = plt.figure(ifig)

    for iprovince,province in enumerate(provinces):

        iplot += 1
        if province == 'BC':
            iiplot += 3
        else:
            iiplot += 1

        #     [left, bottom, width, height]
        pos = position(nrow,ncol,iiplot,hspace=hspace,vspace=vspace)
        
        sub = fig.add_axes(pos,frameon=True) #, axisbg='none')

        width = 4.0
        bottom = np.zeros(len(census_years))

        for component in components:
            years = np.array(data['TREND'][province][component].index)
            vals  = np.array(data['TREND'][province][component].values)
            if component != 'surplus':
                if component != 'crop':
                    sub.bar(years, vals, width,
                                label=years,
                                bottom=bottom,
                                edgecolor='k', linewidth=lwidth,
                                color=components[component]['color'])
                    bottom += vals
                else:
                    sub.bar(years, -1.0*vals, width,
                                label=years,
                                edgecolor='k', linewidth=lwidth,
                                color=components[component]['color'])
            else:
                sub.plot(years, vals, linestyle='-', linewidth=1.5, color=components[component]['color'])
        

        # add label with province name
        sub.text(0.0,1.0,"{}.".format(chr(96+iplot)),
                         verticalalignment='bottom',horizontalalignment='left',rotation=0,
                         fontweight='bold',
                         fontsize=textsize+2,
                         transform=sub.transAxes)
        sub.text(0.05,0.95,"{}".format(provinces[province].replace('&','\n&').replace(' Island','\nIsland')),
                         verticalalignment='top',horizontalalignment='left',rotation=0,
                         fontsize=textsize,
                         transform=sub.transAxes)

        # Add x-title
        if iplot == 8:
            sub.text(1.1,-0.2,"Census years".format(unit),
                         verticalalignment='top',horizontalalignment='center',rotation=0,
                         fontsize=textsize,
                         transform=sub.transAxes)
            
        # Add y-title
        if iplot == 3:
            sub.text(-0.5,0.5,"Phosphorus flows {}".format(unit),  # first value position is moving right-left
                         rotation=90,
                         fontsize=textsize,
                         verticalalignment='center',horizontalalignment='left',
                         transform=sub.transAxes)

        # set limits
        sub.set_ylim([np.min(np.array(data['TREND'][province]['crop'].values)*-1.0)*1.1,np.max(bottom)*1.25])

        # set ticks
        sub.set_xticks(census_years, labels=None, minor=True)
        sub.set_xticks([1961, 1991, 2021], labels=['1961', '1991', '2021'], minor=False)


    # legend
    #     [left, bottom, width, height]
    iiplot = 3
    pos = position(nrow,ncol,iiplot,hspace=hspace,vspace=vspace)
        
    sub = fig.add_axes(pos,frameon=False) #, axisbg='none')
    sub.set_xticks([])
    sub.set_yticks([])

    xshift=0.3
    yshift=0.3

    # Create custom artists
    #      (left, bottom), width, height
    for icomponent,component in enumerate(components):
        if component != 'surplus':
            boxSTi_1 = patches.Rectangle( (0.00+xshift, 0.08+icomponent*0.12+yshift), 0.08, 0.06,
                                           alpha=0.6,
                                           linewidth=lwidth,
                                           edgecolor='k',
                                           facecolor=components[component]['color'],
                                           fill  = True,
                                           transform=sub.transAxes, clip_on=False )
            sub.add_patch(boxSTi_1)
            sub.text(0.11+xshift, 0.1+icomponent*0.12+yshift, components[component]['name'],
                         fontsize=textsize, horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
        else:
            line = patches.Rectangle( (0.00+xshift, 0.03+(-1)*0.1+0.05+yshift), 0.08, 0.00,
                                          color = components['surplus']['color'],
                                          alpha=1.0,
                                          fill  = True,
                                          transform=sub.transAxes, clip_on=False )
            sub.add_patch(line)
            sub.text(0.11+xshift, 0.1+(-1)*0.12+yshift, components[component]['name'],
                         fontsize=textsize, horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)


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






    

