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


# This script generates figure 5 of the following publication:
#      "Long-term pan-Canadian phosphorus mass budget dataset (g)TREND-P-CanadaCANADA: Consistent,
#      unified county-scale and gridded data products from 1961 to 2021"
#      Lamisa Malik, Juliane Mai, Danyka K. Byrnes, Kimberly J. Van Meter, Shuyu Chang, Meghan McLeod, Nandita B. Basu
#      Under review with Nature Scientific Data. DOI forthcoming.

# pyenv activate env-3.11.9
# run figure_5.py -g ../Figures/figure_5

if __name__ == '__main__':
    
    import argparse
    import numpy as np

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
    import pandas as pd
    import time
    import scipy

    # -----------------------------------------------------
    # define some variables
    # -----------------------------------------------------

    data_folder_TREND  = "../Data/County_TREND-P-Canada/"
    data_folder_gTREND = "../Data/Figure_5/"  # contains gridded data reaggregated to county scale
    
    census_years = np.arange(1961, 2022, 5)

    components = {
        'P_FERT': 'Fertilizer P',
        'P_LVSK': 'Livestock Manure P',
        'P_CROP': 'Crop and Pasture P Removal',
        'P_WASTE': 'Domestic Waste P',
        'P_Surplus': 'P Surplus',
        }

    unit = "kg-P ha$^{-1}$ yr$^{-1}$"

    data = {'TREND': {}, 'gTREND': {}}

    # -----------------------------------------------------
    # read data
    # -----------------------------------------------------
    
    for component in components:
        
        filename = "{}/{}.csv".format(data_folder_TREND,component)
        tmp = pd.read_csv(filename,na_values='NAN')
        tmp = tmp.set_index('Year')

        data['TREND'][component] = tmp

        filename = "{}/{}.csv".format(data_folder_gTREND,component)
        tmp = pd.read_csv(filename,na_values='NAN')
        tmp = tmp.set_index('Year')

        data['gTREND'][component] = tmp

    # -----------------------------------------------------
    # plot data
    # -----------------------------------------------------

    # Main plot
    nrow        = 4           # # of rows of subplots per figure
    ncol        = 3           # # of columns of subplots per figure
    hspace      = -0.1         # x-space between subplots
    vspace      = 0.01        # y-space between subplots
    right       = 0.9         # right space on page
    textsize    = 10           # standard text size
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
    import matplotlib.pyplot as plt

    plt.close('all')
    mpl.use('TkAgg')

    if (outtype == 'pdf'):
        mpl.use('PDF') # set directly after import matplotlib
        
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

    for icomponent,component in enumerate(components):

        iplot += 1

        #     [left, bottom, width, height]
        # pos = position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)
        pos = [0.35*(icomponent%3), 0.6-0.22*(icomponent//3), 0.25, 0.17]
        
        sub = fig.add_axes(pos,frameon=True) #, axisbg='none')

        # set limits
        maxval = 27
        sub.set_xlim([0,maxval])
        sub.set_ylim([0,maxval])

        sub.set_xticks([0,5,10,15,20,25])
        sub.set_yticks([0,5,10,15,20,25])

        # plot 1-1 line
        sub.plot([0,maxval], [0,maxval], linestyle='-', linewidth=lwidth*3, color='k')

        # plot datapoints
        for year in census_years:
            xy = pd.concat([data['TREND'][component].loc[year],data['gTREND'][component].loc[year]],axis=1).values

            sub.plot(xy[:,0], xy[:,1],
                             linestyle='None',
                             alpha=0.3,
                             marker='o', markeredgecolor='0.1', markerfacecolor='k', 
                             markersize=5.0, markeredgewidth=0.0)

        # correlation coefficient
        r2 = np.ma.corrcoef(np.ma.masked_invalid(xy[:,0]),np.ma.masked_invalid(xy[:,1]))[0,1]
        # coefficient of determination R2
        idx = ( ~np.isnan(xy[:,0]) & ~np.isnan(xy[:,1]) )
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xy[idx,0],xy[idx,1])
        R2 = r_value**2
        if R2 > 0.99:
            sub.text(0.95,0.05,"{} > {:3.2f}".format("R$^2$",0.99),
                         verticalalignment='bottom',horizontalalignment='right',rotation=0,
                         fontsize=textsize,
                         transform=sub.transAxes)
        else:
            sub.text(0.95,0.05,"{} = {:3.2f}".format("R$^2$",R2),
                         verticalalignment='bottom',horizontalalignment='right',rotation=0,
                         fontsize=textsize,
                         transform=sub.transAxes)
        

        # add label with component name
        sub.text(0.0,1.0,"{}.".format(chr(96+iplot)),
                         verticalalignment='bottom',horizontalalignment='left',rotation=0,
                         fontweight='bold',
                         fontsize=textsize+2,
                         transform=sub.transAxes)
        sub.text(0.05,0.95,"{}".format(components[component].replace('P Removal','P\nRemoval')),
                         verticalalignment='top',horizontalalignment='left',rotation=0,
                         fontsize=textsize,
                         transform=sub.transAxes)

        # Add x-title
        if iplot == 5:
            sub.text(0.5,-0.2,"TREND-P-Canada [{}]".format(unit),
                         verticalalignment='top',horizontalalignment='center',rotation=0,
                         fontsize=textsize,
                         transform=sub.transAxes)
            
        # Add y-title
        if iplot == 4:
            sub.text(-0.3,1.2,"gTREND-P-Canada [{}]".format(unit),
                         rotation=90,
                         fontsize=textsize,
                         verticalalignment='center',horizontalalignment='left',
                         transform=sub.transAxes)


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




