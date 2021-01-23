# Tabulated Equation of State Support

In this repository we will be giving `IllinoisGRMHD` tabulated equation of state (TEOS) support. Our numerical experiments indicate that to achieve the best results with the conservative-to-primitive (con2prim) routine we are using, we should also evolve the entropy $S$ along with the other GRMHD quantities.

# Einstein Toolkit Thorns in this repository

## `IllinoisGRMHD`

High level modifications:

1. TEOS support
1. New primitive variables added: electron fraction, temperature, and entropy (optional)
1. Evolution of the electron fraction, $Y_{\rm e}$
1. Option to evolve the entropy, $S$
1. Reconstruction of the electron fraction, $S$, during PPM
1. Option to reconstruct the specific internal energy, $\epsilon$, during PPM
1. Option to reconstruct the entropy, $S$, during PPM
1. New conservative-to-primitive routine due to Palenzuela *et al.* (see [Palenzuela *et al.*](https://arxiv.org/pdf/1505.01607.pdf) and also the excellent review by [Siegel *et al.*](https://arxiv.org/pdf/1712.07538.pdf)).

## `Convert_to_HydroBase`

High level modifications:

1. Added conversion of the new primitive variables added to `IllinoisGRMHD`

## `ID_converter_ILGRMHD`

High level modifications:

1. Added conversion of the new primitive variables added to `IllinoisGRMHD`

## `ID_tabEOS_HydroBase_Quantities`

This is a new thorn, inteded to be used alongside neutron star initial data thorns, such as `NRPyPlusTOVID` and `Meudon_Bin_NS`. The idea is that these initial data thorns, which have been successfully used with `IllinoisGRMHD` in the past, set some of the `HydroBase` thorn quantities but not all of them. Therefore, this thorn provides initial data capabilities for $Y_{\rm e}$, $T$, and $S$, in a concise, modular, and flexible fashion.

When designing the thorn, care was taken so that we could add maximum flexibility at the parameter file level, allowing the thorn to work with many different types of initial data, as we briefly describe now.

When the functions in this thorn are invoked, we expect that the `HydroBase` variables $rho_{\rm HB}$ and $v^{i}_{\rm HB}$ to have already been set. The thorn will not modify the velocities. We then set $Y_{e}$, $T$, and $S$, allowing the following possiblities. The initialization options for these variables are as follows.

* **Electron fraction**
  1. $Y_{e}$ read from a file. The file should contain $Y_{e}\left(\rho\right)$. Typically this is such that we have neutrino-free beta-equilibrium.

* **Temperature**
  1. $T$ read from a file. The file should contain $T(\left(\rho\right)$. Typically this is normally used for constant entropy initial data.
  1. $T$ constant everywhere.

* **Entropy**
  1. $S$ constant everywhere.
  1. $S$ computed using the EOS table via interpolation, i.e. $S\left(\rho,Y_{e},T\right)$