from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from numpy import array

import sys

plot_dir = "plots/"

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)

"""
0   name
1   type
2   blazhko
3   period
4   log10(P_f/P_0)
5   metallicity
6   prior_mu
7   prior_mu_err
8   SF_ebv
9   SF_ebv_err
10  m_U_obs
11  m_U_obs_err
12  m_B_obs
13  m_B_obs_err
14  m_hipp_obs
15  m_hipp_obs_err
16  m_V_obs
17  m_V_obs_err
18  m_R_obs
19  m_R_obs_err
20  m_I_obs
21  m_I_obs_err
22  m_z_obs
23  m_z_obs_err
24  m_J_obs
25  m_J_obs_err
26  m_H_obs
27  m_H_obs_err
28  m_K_obs
29  m_K_obs_err
30  m_W1_obs
31  m_W1_obs_err
32  m_W2_obs
33  m_W2_obs_err
34  m_W3_obs
35  m_W3_obs_err
36  post_mu
37  post_mu_err
38  post_ebv
39  post_ebv_err
40  M_U_fit
41  M_U_fit_err
42  M_B_fit
43  M_B_fit_err
44  M_hipp_fit
45  M_hipp_fit_err
46  M_V_fit
47  M_V_fit_err
48  M_R_fit
49  M_R_fit_err
50  M_I_fit
51  M_I_fit_err
52  M_z_fit
53  M_z_fit_err
54  M_J_fit
55  M_J_fit_err
56  M_H_fit
57  M_H_fit_err
58  M_K_fit
59  M_K_fit_err
60  M_W1_fit
61  M_W1_fit_err
62  M_W2_fit
63  M_W2_fit_err
64  M_W3_fit
65  M_W3_fit_err
"""

mu_prior = []
mu_prior_err = []
mu_post = []
mu_post_err = []
type = []
blazhko = []

star_fit_data_file = file("full_rrl_fit_table.txt", "r")
first_line = star_fit_data_file.readline()
for line in star_fit_data_file:
    mu_prior.append(float(line.split()[6]))
    mu_prior_err.append(float(line.split()[7]))
    mu_post.append(float(line.split()[36]))
    mu_post_err.append(float(line.split()[37]))
    type.append(line.split()[1])
    blazhko.append(line.split()[2]=="True")
star_fit_data_file.close()

mu_prior = array(mu_prior)
mu_prior_err = array(mu_prior_err)
mu_post = array(mu_post)
mu_post_err = array(mu_post_err)
type = array(type)
blazhko = array(blazhko)



fig = plt.figure(figsize=(3.3, 4.5))
ax1 = fig.add_subplot(1,1,1)
# plot_title = (r"$\mu_{\rm Prior}$ vs $\mu_{\rm Posterior}$ with Residual")
# fig.suptitle(plot_title, fontsize=30)
divider = make_axes_locatable(ax1)
axHistx = divider.append_axes("top", 1.25, pad=0.12, sharex=ax1)

axHistx.plot([6.8, 12.1], [0.0, 0.0], color="black")
ax1.plot([6.8, 12.1], [6.8, 12.1], color="black")

residual_vals = []
residual_errs = []

for n in range(len(blazhko)):
    
    rrl = ["filler", mu_prior[n], mu_post[n], mu_post_err[n], mu_prior_err[n], type[n], blazhko[n]]


    ## RRab and Not Blazhko
    if rrl[5] == "RRab" and not rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="blue",
            markeredgecolor="blue",
            marker="s",
            markersize=3, alpha=0.5)
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="blue",
            markeredgecolor="blue",
            marker="s",
            markersize=3, alpha=0.5)
    
    ## RRab and Blazhko
    if rrl[5] == "RRab" and rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="blue",
            markeredgecolor="blue",
            marker="D",
            markersize=3, alpha=0.5)
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="blue",
            markeredgecolor="blue",
            marker="D",
            markersize=3, alpha=0.5)

for n in range(len(blazhko)):
    
    rrl = ["filler", mu_prior[n], mu_post[n], mu_post_err[n], mu_prior_err[n], type[n], blazhko[n]]

    ## RRc and Not Blazhko
    if rrl[5] == "RRc" and not rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="red",
            markeredgecolor="red",
            marker="s",
            markersize=3, alpha=0.5)
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="red",
            markeredgecolor="red",
            marker="s",
            markersize=3, alpha=0.5)
    
    ## RRc and Blazhko   
    if rrl[5] == "RRc" and rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="red",
            markeredgecolor="red",
            marker="D",
            markersize=3, alpha=0.5)
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="red",
            markeredgecolor="red",
            marker="D",
            markersize=3, alpha=0.5)
            
    residual_vals.append(rrl[2]-rrl[1])
    residual_errs.append(sqrt(rrl[3]**2 + rrl[4]**2))

for tl in axHistx.get_xticklabels():
    tl.set_visible(False)
    



ax1.set_xlabel(r"$\mu_{\rm Prior}$")
ax1.set_ylabel(r"$\mu_{\rm Post}$")
axHistx.set_ylabel(r"$\mu_{\rm Post} - \mu_{\rm Prior}$")

ax1.set_xlim(6.8, 12.1)
ax1.set_ylim(6.8, 12.1)
axHistx.set_ylim(-0.64, 0.64)


# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(1)
minorLocator_x = MultipleLocator(0.2)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(1)
minorLocator_y1 = MultipleLocator(0.2)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

majorLocator_y2 = MultipleLocator(0.2)
minorLocator_y2 = MultipleLocator(0.1)
axHistx.yaxis.set_major_locator(majorLocator_y2)
axHistx.yaxis.set_minor_locator(minorLocator_y2)


# pos = [left, bottom, width, height]
axHistx.set_position([0.18, 0.88, 0.78, 0.15])
ax1.set_position(    [0.18, 0.12, 0.78, 0.85])

# ax1.set_title(plot_title)
canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "mu_mu_plot.pdf", dpi=300)
close("all")


residual_vals = array(residual_vals)
residual_errs = array(residual_errs)

print "Number with +/- 1 sigma is", sum(abs(residual_vals/residual_errs) < 1), "out of", len(blazhko)
print "Fraction within +/- 1 sigma = ", sum(abs(residual_vals/residual_errs) < 1) / float(len(blazhko))
"""
Number with +/- 1 sigma is 112 out of 134
Fraction within +/- 1 sigma =  0.835820895522
"""