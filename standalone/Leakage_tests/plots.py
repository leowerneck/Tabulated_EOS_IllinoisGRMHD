import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams["text.usetex"] = "true"
import matplotlib.pyplot as plt
import astropy.constants as ct

def generate_plot(outfile,
                  harm_data_nue  ,harm_data_anue,
                  nrpy_data_nue  ,nrpy_data_anue,
                  igm_Ye_data_nue,igm_Ye_data_anue,
                  igm_T_data_nue ,igm_T_data_anue,
                  xl1,xl2,yl1,yl2,yt1,ytl1,yt2,ytl2,l1="",l2=""):
    # Begin plotting
    fig,ax = plt.subplots(figsize=(3*2.25,3),ncols=2)

    # Set basic parameters
    ax[0].grid(lw=0.5,ls=':')
    ax[1].grid(lw=0.5,ls=':')

    # Set axes limits
    ax[0].set_xlim(xl2)
    ax[1].set_xlim(xl1)
    ax[0].set_ylim(yl2)
    ax[1].set_ylim(yl1)

    # Set axes labels
    ax[0].set_xlabel(r"$t\ [\rm s]$")
    ax[1].set_xlabel(r"$t\ [\rm s]$")
    ax[0].set_ylabel(r"$Y_{\rm e}(t)-Y_{\rm e}(0)$")
    ax[1].set_ylabel(r"$T(t)-T(0)$")
    ax[0].set_yticks(yt2)
    ax[0].set_yticklabels(ytl2)
    ax[1].set_yticks(yt1)
    ax[1].set_yticklabels(ytl1)

    lw = 1.5

    # Plot harm data
    # ax[0].plot(harm_data_nue[0] ,harm_data_nue[1] -5e-1,label=r"\texttt{HARM}, $Y_{\rm e}(0) = 0.5$"+l1,c='red',lw=lw)
    # ax[0].plot(harm_data_anue[0],harm_data_anue[1]-5e-3,label=r"\texttt{HARM}, $Y_{\rm e}(0) = 0.005$"+l2,c='green',lw=lw)
    # ax[1].plot(harm_data_nue[0] ,harm_data_nue[3] -1   ,label=r"\texttt{HARM}, $Y_{\rm e}(0) = 0.5$"+l1,c='red',lw=lw)
    # ax[1].plot(harm_data_anue[0],harm_data_anue[3]-1   ,label=r"\texttt{HARM}, $Y_{\rm e}(0) = 0.005$"+l2,c='green',lw=lw)
    ax[0].plot(harm_data_nue[0] ,harm_data_nue[1] -5e-1,label=r"\texttt{HARM+NUC}"+l1,c='red',lw=lw)
    ax[0].plot(harm_data_anue[0],harm_data_anue[1]-5e-3,label=r"\texttt{HARM+NUC}"+l2,c='green',lw=lw)
    ax[1].plot(harm_data_nue[0] ,harm_data_nue[3] -1   ,label=r"\texttt{HARM+NUC}"+l1,c='red',lw=lw)
    ax[1].plot(harm_data_anue[0],harm_data_anue[3]-1   ,label=r"\texttt{HARM+NUC}"+l2,c='green',lw=lw)

    # Plot NRPy data
    ax[0].plot(nrpy_data_nue[0] ,nrpy_data_nue[1] -5e-1,label=r"\texttt{NRPyLeakage}"+l1,c='blue',lw=lw,ls="--")
    ax[0].plot(nrpy_data_anue[0],nrpy_data_anue[1]-5e-3,label=r"\texttt{NRPyLeakage}"+l2,c='orange',lw=lw,ls="--")
    ax[1].plot(nrpy_data_nue[0] ,nrpy_data_nue[3] -1   ,label=r"\texttt{NRPyLeakage}"+l1,c='blue',lw=lw,ls="--")
    ax[1].plot(nrpy_data_anue[0],nrpy_data_anue[3]-1   ,label=r"\texttt{NRPyLeakage}"+l2,c='orange',lw=lw,ls="--")

    # Plot IGM data
    units_geom_to_cgs_T = (ct.G.cgs.value * ct.M_sun.cgs.value / ct.c.cgs.value**3 )
    ax[0].plot(units_geom_to_cgs_T*igm_Ye_data_nue[1] ,igm_Ye_data_nue[2] -5e-1,label=r"\texttt{IllinoisGRMHD}"+l1,c='black',lw=lw,ls=":")
    ax[0].plot(units_geom_to_cgs_T*igm_Ye_data_anue[1],igm_Ye_data_anue[2]-5e-3,label=r"\texttt{IllinoisGRMHD}"+l2,c='cyan',lw=lw,ls=":")
    ax[1].plot(units_geom_to_cgs_T*igm_T_data_nue[1]  ,igm_T_data_nue[2]  -1   ,label=r"\texttt{IllinoisGRMHD}"+l1,c='black',lw=lw,ls=":")
    ax[1].plot(units_geom_to_cgs_T*igm_T_data_anue[1] ,igm_T_data_anue[2] -1   ,label=r"\texttt{IllinoisGRMHD}"+l2,c='cyan',lw=lw,ls=":")

    # Save figure
    for i in range(2):
        ax[i].legend(fontsize=9,handlelength=1.2,handletextpad=0.3,borderaxespad=0.15,labelspacing=0.15)

    plt.tight_layout()
    plt.rcParams["figure.figsize"]
    plt.savefig(outfile,dpi=600,bbox_inches="tight",facecolor="white")
    # plt.savefig(outfile)
    # Close fig
    plt.close(fig)

# Read from the data files
# harm_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_only_nue.dat")).T
# harm_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_only_anue.dat")).T
# nrpy_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_only_nue.dat")).T
# nrpy_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_only_anue.dat")).T
l1 = r", only $\mathcal{R}_{\rm ec}^{\nu_{\rm e}}$, $\mathcal{Q}_{\rm ec}^{\nu_{\rm e}}$"
l2 = r", only $\mathcal{R}_{\rm pc}^{\bar{\nu}_{\rm e}}$, $\mathcal{Q}_{\rm ec}^{\bar{\nu}_{\rm e}}$"
# generate_plot("NRPyLeakage_semi_analytic_results_beta_only.png",harm_data_nue,harm_data_anue,nrpy_data_nue,nrpy_data_anue,
#               [0,0.5],[0,0.5],[-4.2e-4,0.2e-4],[-0.25,0.15],
#               [-4e-4,-3e-4,-2e-4,-1e-4,0],
#               [r"$-4{\times}10^{-4}$",r"$-3{\times}10^{-4}$",r"$-2{\times}10^{-4}$",r"$-1{\times}10^{-4}$",r"$0$"],
#               l1=l1,l2=l2)

harm_data_nue    = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_all_5e-1.dat")).T
harm_data_anue   = np.loadtxt(os.path.join("optically_thin","results_igm_paper","harm_all_5e-3.dat")).T
nrpy_data_nue    = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_all_5e-1.dat")).T
nrpy_data_anue   = np.loadtxt(os.path.join("optically_thin","results_igm_paper","nrpy_all_5e-3.dat")).T
igm_Ye_data_nue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","igm_Ye.maxmin.large_Ye.asc")).T
igm_Ye_data_anue = np.loadtxt(os.path.join("optically_thin","results_igm_paper","igm_Ye.maxmin.small_Ye.asc")).T
igm_T_data_nue   = np.loadtxt(os.path.join("optically_thin","results_igm_paper","igm_temperature.maxmin.large_Ye.asc")).T
igm_T_data_anue  = np.loadtxt(os.path.join("optically_thin","results_igm_paper","igm_temperature.maxmin.small_Ye.asc")).T

generate_plot("NRPyLeakage_semi_analytic_results_all.png",
              harm_data_nue,harm_data_anue,nrpy_data_nue,nrpy_data_anue,
              igm_Ye_data_nue,igm_Ye_data_anue,igm_T_data_nue,igm_T_data_anue,
              [0,0.5],[0,0.5],[-5.25e-3,0.25e-3],[-0.10-0.0125,0.10+0.0125],
              [-5e-3,-4e-3,-3e-3,-2e-3,-1e-3,0],
              [r"$-0.005$",r"$-0.004$",r"$-0.003$",r"$-0.002$",r"$-0.001$",r"$0.000$"],
              [-0.1,-0.05,0,0.05,0.1],
              [r"$-0.10$",r"$-0.05$",r"$0.00$",r"$+0.05$",r"$+0.10$"])

# harm_data_nue  = np.loadtxt("harm_large_Ye.dat").T
# harm_data_anue = np.loadtxt("harm_small_Ye.dat").T
# nrpy_data_nue  = np.loadtxt("nrpy_large_Ye.dat").T
# nrpy_data_anue = np.loadtxt("nrpy_small_Ye.dat").T
# generate_plot("NRPyLeakage_semi_analytic_results_all.png",harm_data_nue,harm_data_anue,nrpy_data_nue,nrpy_data_anue,
#               [0,0.5],[0,0.5],[-5.25e-3,0.25e-3],[-0.15,0.15],[-5e-3,-4e-3,-3e-3,-2e-3,-1e-3,0],
#               [r"$-5{\times}10^{-3}$",r"$-4{\times}10^{-3}$",r"$-3{\times}10^{-3}$",r"$-2{\times}10^{-3}$",r"$-1{\times}10^{-3}$",r"$0$"])
