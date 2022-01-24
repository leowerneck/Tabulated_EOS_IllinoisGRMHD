import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams["text.usetex"] = "true"
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
    print("Correct usage is: python plots.py paper_data_nue_file paper_data_anue_file nrpy_data_nue_file nrpy_data_anue_file")
    sys.exit(1)

# Read from the data files
paper_data_nue  = np.loadtxt(sys.argv[1]).T
paper_data_anue = np.loadtxt(sys.argv[2]).T
nrpy_data_nue   = np.loadtxt(sys.argv[3]).T
nrpy_data_anue  = np.loadtxt(sys.argv[4]).T

# Begin plotting
fig = plt.figure(figsize=(3.54,3.54))

# Set grid
plt.grid(lw=0.5,ls=':')

# Set axes limits
plt.xlim(0,0.5)
plt.ylim(-0.2,0.2)

# Set axes labels
plt.xlabel(r"$t\ [\rm s]$")
plt.ylabel(r"$Y_{\rm e}(t)-Y_{\rm e}(0)$")

# Plot NRPy data
plt.plot(nrpy_data_nue[0] ,nrpy_data_nue[1] -5e-1,label=r"NRPyLeakage, only $\beta$-processes with $\nu_{\rm e}$",c='blue',lw=1)
plt.plot(nrpy_data_anue[0],nrpy_data_anue[1]-5e-3,label=r"NRPyLeakage, only $\beta$-processes with $\bar{\nu}_{\rm e}$",c='orange',lw=1)

# Plot paper data
plt.plot(paper_data_nue[0] ,paper_data_nue[1]      ,label=r"Ari's paper, only $\beta$-processes with $\nu_{\rm e}$",c='red',lw=0,marker='o',ms=4,mfc='none')
plt.plot(paper_data_anue[0],paper_data_anue[1]-5e-3,label=r"Ari's paper, only $\beta$-processes with $\bar{\nu}_{\rm e}$",c='green',lw=0,marker='o',ms=4,mfc='none')

# Save figure
plt.legend(fontsize=8)
plt.savefig("NRPyLeakage_semi_analytic_results.png",dpi=600,bbox_inches="tight",facecolor="white")

# Close fig
plt.close(fig)
