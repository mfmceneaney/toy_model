import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def convert_graph_to_csv(
    filename,
    x,
    y,
    xerr=None,
    yerr=None,
    xerr_syst=None,
    yerr_syst=None,
    delimiter=",",
    header=None,
    fmt=None,
    comments='',
    ):

    data = []
    if xerr_syst is None or len(xerr_syst)==0: xerr_syst = [0.0 for el in x]
    if yerr_syst is None or len(yerr_syst)==0: yerr_syst = [0.0 for el in x]
    for i, el in enumerate(x):
        data.append([i, x[i], y[i], xerr[i], yerr[i], xerr_syst[i], yerr_syst[i]])

    data = np.array(data)

    print("DEBUGGING: np.shape(data) = ",np.shape(data))
    print("DEBUGGING: np.shape(fmt) = ",np.shape(fmt))#DEBUGGING

    header = "REPLACEMENT_HEADER"+header

    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(filename, 'w') as file:
        file.write(filedata)

def get_data_from_tgrapherror(
    path = 'HB_costheta1_Q2_1.000_10.000_sgasym_0.00_bgasym_0.00.root',
    name = "Graph",
    ):

    # Get TGraphErrors from ROOT file
    try:
        f = ur.open(path)
        g = [np.array(f[name].member('fX')), f[name].member('fY'), f[name].member('fEX'), f[name].member('fEY')]

        return g
    except FileNotFoundError:
        print("DEBUGGING: FileNotFoundError: ",path)
        print("\t Returning empty list")
        return []

def get_arrs(out_file_list):

    # Initialize output list
    glist = []
    
    # Loop files
    for filename in out_file_list:
        g = get_data_from_tgrapherror(filename)
        if len(g)>0: glist.append(g)

    if len(glist)==0:
        print("ERROR: len(glist)==0")
        return {
            'x_mean':[],
            'y_mean':[],
            'xerr_mean':[],
            'yerr_mean':[],
            'y_min':[],
            'y_max':[],
            }

    # Convert to numpy
    glist  = np.array(glist)
    print("DEBUGGING: gshape = ",np.shape(glist))#DEBUGGING
    glist  = np.swapaxes(glist,0,1)
    print("DEBUGGING: gshape new = ",np.shape(glist))

    # Get arrays
    x_mean    = np.mean(glist[0],axis=0) #NOTE: Get mean across different graphs (axis=0) but not across bins (axis=1)
    y_mean    = np.mean(glist[1],axis=0)
    x_std     = np.std(glist[0],axis=0) #NOTE: Get std across different graphs (axis=0) but not across bins (axis=1)
    y_std     = np.std(glist[1],axis=0)
    xerr_mean = np.sqrt(np.mean(np.square(glist[2]),axis=0)) #NOTE: Get mean across different graphs (axis=0) but not across bins (axis=1)
    yerr_mean = np.sqrt(np.mean(np.square(glist[3]),axis=0))
    y_min     = np.min(glist[1],axis=0)
    y_max     = np.max(glist[1],axis=0)

    return {
            'x_mean':x_mean,
            'y_mean':y_mean,
            'x_std':x_std,
            'y_std':y_std,
            'xerr_mean':xerr_mean,
            'yerr_mean':yerr_mean,
            'y_min':y_min,
            'y_max':y_max
            }

def get_file_name(asym_num=0,idx=0):
    return "MLFit_x_y_z_0.000_1.000_A"+str(asym_num)+"__fitidx"+str(idx)+".root"

def main(nasyms,nreps):

    # Get file names
    file_lists = {asym_num:[get_file_name(asym_num=asym_num,idx=idx) for idx in range(nreps)] for asym_num in range(nasyms)}

    # Get graph data arrays
    graph_arrs = {asym_num:get_arrs(file_lists[asym_num]) for asym_num in range(nasyms)}

    # Plot graphs
    for asym_num in range(nasyms):
        plot_graph(**graph_arrs[asym_num],asym_num=asym_num)


def plot_graph(
    x_mean = [],
    y_mean = [],
    xerr_mean = [],
    yerr_mean = [],
    xerr_syst = [],
    yerr_syst = [],
    y_min  = [],
    y_max  = [],
    y_std  = [],
    x_std  = [],
    xlims = [0.0,1.0],
    ylims = [-0.05,0.2],
    title = 'Injection Results',
    xvar  = 'z',
    xtitle = '$z_{p\pi^{-}}$',
    ytitles = {0:'$\mathcal{A}_{LUT}^{\cos{\phi_{\Lambda}-\phi_{S_{\Lambda}}}}$',1:'$\mathcal{A}_{LUT}^{\cos{\phi_{S_{\Lambda}}}}$',2:'$\mathcal{A}_{LUT}^{\cos{2\phi_{\Lambda}-\phi_{S_{\Lambda}}}}$'},
    sgasyms = [0.00, 0.10, 0.00],
    bgasym = 0.00,
    color  = 'blue', #NOTE: COLOR OF DATA POINTS
    bcolor = 'gray', #NOTE:
    path = 'graph',
    verbose = True,
    yaml_args = {},
    asym_num = 0,
    ):

    # Set font sizes
    plt.rc('font', size=25) #controls default text size
    plt.rc('axes', titlesize=50) #fontsize of the title
    plt.rc('axes', labelsize=50) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
    plt.rc('legend', fontsize=20) #fontsize of the legend

    # Get some nicer plot settings
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    ecolor='black'
    elinewidth=2.0
    capsize=18
    capthick=2.0
    marker='o'
    linestyle=None
    linewidth=0.0
    markersize=20
    gridlinewidth=0.5
    axlinewidth=1

    # Plot 
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitles[asym_num],usetex=True)

    #TODO: DEBUGGING MESSAGE FOR BINS SEE IF SOMETHING GETS MESSED UP THERE AND MAKE SURE YOU ARE SETTING CORRECTLY...
    if len(y_std)>0: fb = plt.fill_between(x_mean, y_mean-y_std, y_mean+y_std, alpha=0.2, label='$\pm1\sigma$ Band', color=bcolor)
    elif len(yerr_syst)>0: g1 = plt.errorbar(x_mean,y_mean,xerr=None,yerr=yerr_syst,
                ecolor='gray', elinewidth=elinewidth*20, capsize=0,
                color=color, marker='o', linestyle=linestyle, alpha=0.5,
                linewidth=0, markersize=0,label='Systematic error of '+ytitles[asym_num])
    g2 = plt.errorbar(x_mean,y_mean,xerr=xerr_mean,yerr=yerr_mean,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=color, marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label=ytitles[asym_num])
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    if len(sgasyms)>asym_num: ax1.axhline(sgasyms[asym_num], color='red',linestyle='-',linewidth=axlinewidth,label='Injected '+ytitles[asym_num])
    plt.legend(loc='best')
    outpath = path+'__A'+str(asym_num)+'.pdf'
    print("DEBUGGING: plt.savefig(outpath) -> ",outpath)
    f1.savefig(outpath)

    # Save plot data to csv
    delimiter = ","
    header    = delimiter.join(["bin","x","y","xerr","yerr","xerrsyst","yerrsyst"]) #NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    fmt       = ["%d","%.3g","%.3g","%.3g","%.3g","%.3g","%.3g"]
    comments  = ""

    convert_graph_to_csv(
        outpath+'.csv',
        x_mean,
        y_mean,
        xerr=xerr_mean,
        yerr=yerr_mean,
        xerr_syst=xerr_syst,
        yerr_syst=yerr_syst,
        delimiter=delimiter,
        header=header,
        fmt=fmt,
        comments=comments
        )

if __name__=="__main__":

    nasyms = 3
    nreps = 50
    main(nasyms,nreps)
    print("DONE")
