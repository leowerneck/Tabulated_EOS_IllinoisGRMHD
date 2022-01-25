import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams["text.usetex"] = "true"
import matplotlib.pyplot as plt

def generate_plot(outfile,harm_data_nue,harm_data_anue,nrpy_data_nue,nrpy_data_anue,
                  xl1,xl2,yl1,yl2,yt1,ytl1,l1="",l2=""):
    # Begin plotting
    fig,ax = plt.subplots(figsize=(3.54*2.25,3.54),ncols=2)

    # Set basic parameters
    ax[0].grid(lw=0.5,ls=':')
    ax[1].grid(lw=0.5,ls=':')

    # Set axes limits
    ax[0].set_xlim(xl1)
    ax[1].set_xlim(xl2)
    ax[0].set_ylim(yl1)
    ax[1].set_ylim(yl2)

    # Set axes labels
    ax[0].set_xlabel(r"$t\ [\rm s]$")
    ax[1].set_xlabel(r"$t\ [\rm s]$")
    ax[0].set_ylabel(r"$T(t)-T(0)$")
    ax[1].set_ylabel(r"$Y_{\rm e}(t)-Y_{\rm e}(0)$")
    ax[0].set_yticks(yt1)
    ax[0].set_yticklabels(ytl1)

    lw = 1.5

    # Plot harm data
    ax[0].plot(harm_data_nue[0] ,harm_data_nue[3] -1   ,label=r"\texttt{HARM}"+l1,c='red',lw=lw)
    ax[0].plot(harm_data_anue[0],harm_data_anue[3]-1   ,label=r"\texttt{HARM}"+l2,c='green',lw=lw)
    ax[1].plot(harm_data_nue[0] ,harm_data_nue[1] -5e-1,label=r"\texttt{HARM}"+l1,c='red',lw=lw)
    ax[1].plot(harm_data_anue[0],harm_data_anue[1]-5e-3,label=r"\texttt{HARM}"+l2,c='green',lw=lw)

    # Plot NRPy data
    ax[0].plot(nrpy_data_nue[0] ,nrpy_data_nue[3] -1   ,label=r"\texttt{NRPyLeakage}"+l1,c='blue',lw=lw,ls="--")
    ax[0].plot(nrpy_data_anue[0],nrpy_data_anue[3]-1   ,label=r"\texttt{NRPyLeakage}"+l2,c='orange',lw=lw,ls="--")
    ax[1].plot(nrpy_data_nue[0] ,nrpy_data_nue[1] -5e-1,label=r"\texttt{NRPyLeakage}"+l1,c='blue',lw=lw,ls="--")
    ax[1].plot(nrpy_data_anue[0],nrpy_data_anue[1]-5e-3,label=r"\texttt{NRPyLeakage}"+l2,c='orange',lw=lw,ls="--")

    # Save figure
    for i in range(2):
        ax[i].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile,dpi=600,bbox_inches="tight",facecolor="white")
    # Close fig
    plt.close(fig)

# Read from the data files
harm_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_only_nue.dat")).T
harm_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_only_anue.dat")).T
nrpy_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_only_nue.dat")).T
nrpy_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_only_anue.dat")).T
l1 = r", only $\mathcal{R}_{\rm ec}^{\nu_{\rm e}}$, $\mathcal{Q}_{\rm ec}^{\nu_{\rm e}}$"
l2 = r", only $\mathcal{R}_{\rm pc}^{\bar{\nu}_{\rm e}}$, $\mathcal{Q}_{\rm ec}^{\bar{\nu}_{\rm e}}$"
generate_plot("NRPyLeakage_semi_analytic_results_beta_only.png",harm_data_nue,harm_data_anue,nrpy_data_nue,nrpy_data_anue,
              [0,0.5],[0,0.5],[-4.2e-4,0.2e-4],[-0.25,0.15],[-4e-4,-3e-4,-2e-4,-1e-4,0],
              [r"$-4{\times}10^{-4}$",r"$-3{\times}10^{-4}$",r"$-2{\times}10^{-4}$",r"$-1{\times}10^{-4}$",r"$0$"],
              l1=l1,l2=l2)

harm_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_all_5e-1.dat")).T
harm_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_all_5e-3.dat")).T
nrpy_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_all_5e-1.dat")).T
nrpy_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_all_5e-3.dat")).T
generate_plot("NRPyLeakage_semi_analytic_results_all.png",harm_data_nue,harm_data_anue,nrpy_data_nue,nrpy_data_anue,
              [0,0.5],[0,0.5],[-5.25e-3,0.25e-3],[-0.15,0.15],[-5e-3,-4e-3,-3e-3,-2e-3,-1e-3,0],
              [r"$-5{\times}10^{-3}$",r"$-4{\times}10^{-3}$",r"$-3{\times}10^{-3}$",r"$-2{\times}10^{-3}$",r"$-1{\times}10^{-3}$",r"$0$"])
