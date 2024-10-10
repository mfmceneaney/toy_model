import numpy as np
import matplotlib.pyplot as plt


# x = np.array(
#     [10**i for i in range(2,9)]
# )

# a = [0.5,0.1,0.5,0.1,0.17,0.07,0.15,0.09]
# y = np.array(
#     [[el * 0.9**i for el in a] for i in range(2,9)]
# )

x = np.array(
    [10**i for i in range(4,7)]
)

y = np.array(
    [
        [ 0.551031 , 0.0869382 , 0.497453 , 0.100504 , 0.551057 , 0.0869317 , 0.498099 , 0.100376 , 0.517395 , 0.833843 , 0.562287 , 0.389695 , 0.723325 , 0.901871 , 0.584717 , 0.256758 ],
        [ 0.51241 , 0.0973439 , 0.496867 , 0.100615 , 0.51251 , 0.0973219 , 0.496882 , 0.100612 , 0.742586 , 0.804167 , 0.724596 , 0.28054 , 0.728171 , 0.103736 , 0.653657 , 0.366208 ],
        [ 0.498767 , 0.100245 , 0.496933 , 0.100602 , 0.498752 , 0.100248 , 0.49694 , 0.100601 , 0.747183 , 0.138772 , 0.741992 , 0.125246 , 0.748747 , 0.1566 , 0.748201 , 0.145662 ],
    ]
)

# x = np.array(
#     [10**i for i in range(4,9)]
# )

# y = np.array(
#     [
#         # [  0.265921 , 0.0987385 , 0.225441 , 0.104488 , 0.509133 , 0.41917 , 0.81338 , 0.366232 ],
#         # [ 0.596488 , 0.0668538 , 0.421957 , 0.107756 , 0.713353 , 0.407455 , 0.676266 , 0.278485 ],
#         [ 0.530279 , 0.0928033 , 0.457618 , 0.106235 , 0.764393 , 0.135009 , 0.722773 , 0.134216 ],
#         [ 0.510379 , 0.0977965 , 0.506421 , 0.0986673 , 0.234577 , 0.0463371 , 0.231212 , 0.0463928 ],
#         [ 0.501018 , 0.0997998 , 0.505245 , 0.098919 , 0.174756 , 0.0353737 , 0.0753589 , 0.0147953 ],
#         [ 0.499981 , 0.100013 , 0.499399 , 0.100122 , 0.0557952 , 0.0111784 , 0.0694246 , 0.0139523 ],
#         [ 0.5 , 0.1 , 0.500022 , 0.0999972 , 0.0177115 , 0.00354285 , 0.0221537 , 0.0044315 ]
#     ]
# )

# nevents = 100
# 0.265921 , 0.0987385 , 0.225441 , 0.104488 , 0.509133 , 0.41917 , 0.81338 , 0.366232 ]
# DONE
# ----------------------------------------------------------------------
# nevents = 1000
# 0.596488 , 0.0668538 , 0.421957 , 0.107756 , 0.713353 , 0.407455 , 0.676266 , 0.278485 ]
# DONE
# ----------------------------------------------------------------------
# nevents = 10000
# 0.530279 , 0.0928033 , 0.457618 , 0.106235 , 0.764393 , 0.135009 , 0.722773 , 0.134216 ]
# DONE
# ----------------------------------------------------------------------
# nevents = 100000
# 0.510379 , 0.0977965 , 0.506421 , 0.0986673 , 0.234577 , 0.0463371 , 0.231212 , 0.0463928 ]
# DONE
# ----------------------------------------------------------------------
# nevents = 1000000
# 0.501018 , 0.0997998 , 0.505245 , 0.098919 , 0.174756 , 0.0353737 , 0.0753589 , 0.0147953 ]
# DONE
# ----------------------------------------------------------------------
# nevents = 10000000
# 0.499981 , 0.100013 , 0.499399 , 0.100122 , 0.0557952 , 0.0111784 , 0.0694246 , 0.0139523 ]
# DONE
# ----------------------------------------------------------------------
# nevents = 100000000
# 0.5 , 0.1 , 0.500022 , 0.0999972 , 0.0177115 , 0.00354285 , 0.0221537 , 0.0044315 ]

keys = [
    "a0x_2d_val",
   "a1x_2d_val",
   "a0y_2d_val",
   "a1y_2d_val",
   "a0x_1d_val",
   "a1x_1d_val",
   "a0y_1d_val",
   "a1y_1d_val",
   "a0x_2d_err",
   "a1x_2d_err",
   "a0y_2d_err",
   "a1y_2d_err",
   "a0x_1d_err",
   "a1x_1d_err",
   "a0y_1d_err",
   "a1y_1d_err",
]

labels = [
   "$A_{0, x, 2D}$",
   "$A_{1, x, 2D}$",
   "$A_{0, y, 2D}$",
   "$A_{1, y, 2D}$",
   "$A_{0, x, 1D}$",
   "$A_{1, x, 1D}$",
   "$A_{0, y, 1D}$",
   "$A_{1, y, 1D}$",
   "$\delta A_{0, x, 2D}/A_{0, x, 2D}$",
   "$\delta A_{1, x, 2D}/A_{1, x, 2D}$",
   "$\delta A_{0, y, 2D}/A_{0, y, 2D}$",
   "$\delta A_{1, y, 2D}/A_{1, y, 2D}$",
   "$\delta A_{0, x, 1D}/A_{0, x, 1D}$",
   "$\delta A_{1, x, 1D}/A_{1, x, 1D}$",
   "$\delta A_{0, y, 1D}/A_{0, y, 1D}$",
   "$\delta A_{1, y, 1D}/A_{1, y, 1D}$",
]
print("len(labels) = ",len(labels))

# Set default figure size
figsize = (16,10)

# Set font sizes
plt.rc('font', size=25) #controls default text size
plt.rc('axes', titlesize=40) #fontsize of the title
plt.rc('axes', labelsize=40) #fontsize of the x and y labels
plt.rc('xtick', labelsize=20) #fontsize of the x tick labels
plt.rc('ytick', labelsize=20) #fontsize of the y tick labels
plt.rc('legend', fontsize=25) #fontsize of the legend

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True

# Error bar settings
ecolor='black'
elinewidth=2.0
capsize=18
capthick=2.0
### color='black'
marker='o'
linestyle=None
linewidth=0.0
markersize=20
gridlinewidth=0.5
axlinewidth=1

# Marker settings
markersize = 20
alpha      = 0.5

for i in range(0,len(keys)//4):

    idx_2d_val = i
    idx_1d_val = i+2
    idx_2d = i+len(keys)//2
    idx_1d = i+len(keys)//2+4
    print("DEBUGGING: idx_2d = ",idx_2d)
    print("DEBUGGING: idx_1d = ",idx_1d)

    # Plot comparison 2D vs. 1D results and errors
    f, ax = plt.subplots(figsize=figsize)
    # plt.ylim(0.0,2.0)
    ax.set_xscale('log')
    scaling_factor = 1.1
    h_2d = plt.errorbar(x,y[:,idx_2d_val],yerr=y[:,idx_2d],ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='b', marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,alpha=alpha,label=labels[idx_2d_val])
    h_1d = plt.errorbar(x*scaling_factor,y[:,idx_1d_val],yerr=y[:,idx_1d],ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='r', marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,alpha=alpha,label=labels[idx_1d_val])
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    plt.ylabel("Error")
    plt.xlabel("Dataset Size")
    filename = keys[idx_2d_val]+"__"+keys[idx_1d_val]+".pdf"
    plt.legend(loc='best')
    f.savefig(filename)

    # Plot error over value comparison 2D vs. 1D
    f = plt.figure(figsize=figsize)
    # plt.ylim(0.0,2.0)
    h_2d = plt.semilogx(x,np.divide(y[:,idx_2d],y[:,idx_2d_val]),"b^",markersize=markersize,alpha=alpha,label=labels[idx_2d])
    h_1d = plt.semilogx(x,np.divide(y[:,idx_1d],y[:,idx_1d_val]),"rv",markersize=markersize,alpha=alpha,label=labels[idx_1d])
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    plt.ylabel("Error")
    plt.xlabel("Dataset Size")
    filename = keys[idx_2d]+"__"+keys[idx_1d]+".pdf"
    plt.legend(loc='best')
    f.savefig(filename)

plt.show()
