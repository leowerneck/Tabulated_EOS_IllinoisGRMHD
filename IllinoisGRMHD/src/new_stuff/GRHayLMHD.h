#ifndef GRHAYLET_GRHAYLMHD_H
#define GRHAYLET_GRHAYLMHD_H

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <math.h>
#include <stdbool.h>

#include "GRHayLib.h"

enum GRHayLMHD_bc_dir {
    GRHayLMHD_xmin,
    GRHayLMHD_xmax,
    GRHayLMHD_ymin,
    GRHayLMHD_ymax,
    GRHayLMHD_zmin,
    GRHayLMHD_zmax
};

void GRHayLMHD_Prim2Con_SinglePoint(CCTK_ARGUMENTS,
                                    const int ijk,
                                    const int ijkx,
                                    const int ijky,
                                    const int ijkz);

#define LOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax)                                   \
    for(int k = _kmin; k < _kmax; k++) {                                                   \
        for(int j = _jmin; j < _jmax; j++) {                                               \
            for(int i = _imin; i < _imax; i++) {                                           \
                CCTK_ATTRIBUTE_UNUSED const int ijk  = CCTK_GFINDEX3D(cctkGH, i, j, k);    \
                CCTK_ATTRIBUTE_UNUSED const int ijkx = CCTK_GFINDEX4D(cctkGH, i, j, k, 0); \
                CCTK_ATTRIBUTE_UNUSED const int ijky = CCTK_GFINDEX4D(cctkGH, i, j, k, 1); \
                CCTK_ATTRIBUTE_UNUSED const int ijkz = CCTK_GFINDEX4D(cctkGH, i, j, k, 2);

#define OMPLOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax) \
    _Pragma("omp parallel for") LOOP3D(_imin, _imax, _jmin, _jmax, _kmin, _kmax)

#define ENDLOOP3D \
    }             \
    }             \
    }

// These are useful macros for packing (load) and unpacking (write) GRHayL's structs.
#define GRHAYLMHD_LOAD_PRIMS_AT_INDEX(struct_name, _ijk, _ijkx, _ijky, _ijkz) \
    ghl_initialize_primitives(rho[_ijk],                                      \
                              press[_ijk],                                    \
                              eps[_ijk],                                      \
                              vx[_ijk],                                       \
                              vy[_ijk],                                       \
                              vz[_ijk],                                       \
                              Bvec[_ijkx],                                    \
                              Bvec[_ijky],                                    \
                              Bvec[_ijkz],                                    \
                              entropy[_ijk],                                  \
                              Y_e[_ijk],                                      \
                              temperature[_ijk],                              \
                              &struct_name);

#define GRHAYLMHD_WRITE_PRIMS_AT_INDEX(struct_name, _ijk, _ijkx, _ijky, _ijkz) \
    ghl_return_primitives(&struct_name,                                        \
                          &rho[_ijk],                                          \
                          &press[_ijk],                                        \
                          &eps[_ijk],                                          \
                          &vx[_ijk],                                           \
                          &vy[_ijk],                                           \
                          &vz[_ijk],                                           \
                          &Bvec[_ijkx],                                        \
                          &Bvec[_ijky],                                        \
                          &Bvec[_ijkz],                                        \
                          &entropy[_ijk],                                      \
                          &Y_e[_ijk],                                          \
                          &temperature[_ijk]);

#define GRHAYLMHD_LOAD_CONS_AT_INDEX(struct_name, _ijk) \
    ghl_initialize_conservatives(rho_tilde[_ijk],       \
                                 tau_tilde[_ijk],       \
                                 S_x_tilde[_ijk],       \
                                 S_y_tilde[_ijk],       \
                                 S_z_tilde[_ijk],       \
                                 ent_tilde[_ijk],       \
                                 Y_e_tilde[_ijk],       \
                                 &struct_name);

#define GRHAYLMHD_WRITE_CONS_AT_INDEX(struct_name, _ijk) \
    ghl_return_conservatives(&struct_name,               \
                             &rho_tilde[_ijk],           \
                             &tau_tilde[_ijk],           \
                             &S_x_tilde[_ijk],           \
                             &S_y_tilde[_ijk],           \
                             &S_z_tilde[_ijk],           \
                             &ent_tilde[_ijk],           \
                             &Y_e_tilde[_ijk]);

#define GRHAYLMHD_LOAD_METRIC_AT_INDEX(struct_name, _ijk) \
    ghl_initialize_metric(alp[_ijk],                      \
                          betax[_ijk],                    \
                          betay[_ijk],                    \
                          betaz[_ijk],                    \
                          gxx[_ijk],                      \
                          gxy[_ijk],                      \
                          gxz[_ijk],                      \
                          gyy[_ijk],                      \
                          gyz[_ijk],                      \
                          gzz[_ijk],                      \
                          &struct_name);

#define GRHAYLMHD_LOAD_METRIC_ENFORCE_DETGAMMAEQ1(struct_name) \
    ghl_enforce_detgtij_and_initialize_ADM_metric(alp[ijk],    \
                                                  betax[ijk],  \
                                                  betay[ijk],  \
                                                  betaz[ijk],  \
                                                  gxx[ijk],    \
                                                  gxy[ijk],    \
                                                  gxz[ijk],    \
                                                  gyy[ijk],    \
                                                  gyz[ijk],    \
                                                  gzz[ijk],    \
                                                  &struct_name);

#define GRHAYLMHD_LOAD_EXTRINSIC_CURVATURE(struct_name) \
    ghl_initialize_extrinsic_curvature(kxx[ijk],        \
                                       kxy[ijk],        \
                                       kxz[ijk],        \
                                       kyy[ijk],        \
                                       kyz[ijk],        \
                                       kzz[ijk],        \
                                       &struct_name);

#define GRHAYLMHD_LOAD_METRIC(struct_name) GRHAYLMHD_LOAD_METRIC_AT_INDEX(struct_name, ijk)
#define GRHAYLMHD_LOAD_CONS(struct_name)   GRHAYLMHD_LOAD_CONS_AT_INDEX(struct_name, ijk)
#define GRHAYLMHD_WRITE_CONS(struct_name)  GRHAYLMHD_WRITE_CONS_AT_INDEX(struct_name, ijk)
#define GRHAYLMHD_LOAD_PRIMS(struct_name) \
    GRHAYLMHD_LOAD_PRIMS_AT_INDEX(struct_name, ijk, ijkx, ijky, ijkz)
#define GRHAYLMHD_WRITE_PRIMS(struct_name) \
    GRHAYLMHD_WRITE_PRIMS_AT_INDEX(struct_name, ijk, ijkx, ijky, ijkz)

#endif // GRHAYLET_GRHAYLMHD_H
