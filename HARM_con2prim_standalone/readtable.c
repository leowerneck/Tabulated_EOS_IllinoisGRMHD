/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

This file reads the EOS table and saves it to memory. For units and more information, see Schneider, Roberts, Ott, 2017, and
the user_guide for the tables from https://stellarcollapse.org/SROEOS 
The format is in the stellarcollapse format

******************************************************************************/

#include "decs.h"
#define H5_USE_16_API 1
#include "hdf5.h"
#include "nuc_eos.h"

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                              \
  do {                                                                   \
    int _error_code = fn_call;                                           \
    if (_error_code < 0) {	       				         \
      printf("HDF5 call '%s' returned error code %d",               \
                  #fn_call, _error_code);                                \
    }                                                                    \
  } while (0)

static int file_is_readable(const char* filename);

static int file_is_readable(const char* filename)
{
    FILE* fp = NULL;
    fp = fopen(filename, "r");
    if(fp != NULL)
    {
        fclose(fp);
        return 1;
    }
    return 0;
}


// define the variables
//namespace nuc_eos {
//  double temp0, temp1;
//  double energy_shift;
//
//  double eos_rhomax, eos_rhomin;
//  double eos_tempmin, eos_tempmax;
//  double eos_yemin, eos_yemax;
//  double eos_epsmin, eos_epsmax;
//  
//  double c2p_tempmin;
//  double c2p_tempmax;
//
//}
//namespace nuc_eos_private {
//  int nrho;
//  int ntemp;
//  int nye;
//
//  double * restrict alltables;
//  double * restrict epstable;
//  double * restrict logrho;
//  double * restrict logtemp;
//  double dlintemp, dlintempi;
//  double drholintempi;
//  double dlintempyei;
//  double drholintempyei;
//  double * restrict yes;
//  double dtemp, dtempi;
//  double drho, drhoi;
//  double dye, dyei;
//  double drhotempi;
//  double drhoyei;
//  double dtempyei;
//  double drhotempyei;
//}



// TODO: replace with version in ET EOS_Omni. NOTE: table arrangement changed.

// Cactus calls this function. It reads in the table and calls a fortran
// function to setup values for the fortran eos module
void nuc_eos_C_ReadTable(char* nuceos_table_name)
{
//  double temp0, temp1;
//  double energy_shift;
//
//  double eos_rhomax, eos_rhomin;
//  double eos_tempmin, eos_tempmax;
//  double eos_yemin, eos_yemax;
//  double eos_epsmin, eos_epsmax;
//  
//  double c2p_tempmin;
//  double c2p_tempmax;
//
//  int nrho;
//  int ntemp;
//  int nye;
//
//  double * restrict alltables;
//  double * restrict epstable;
//  double * restrict logrho;
//  double * restrict logtemp;
//  double dlintemp, dlintempi;
//  double drholintempi;
//  double dlintempyei;
//  double drholintempyei;
//  double * restrict yes;
//  double dtemp, dtempi;
//  double drho, drhoi;
//  double dye, dyei;
//  double drhotempi;
//  double drhoyei;
//  double dtempyei;
//  double drhotempyei;

  //printf("*******************************");
  //printf("Reading nuc_eos table file:");
  //printf("%s",nuceos_table_name);
  //printf("*******************************");

  hid_t file;
  if (!file_is_readable(nuceos_table_name)) {
    printf("Could not read nuceos_table_name %s \n",
               nuceos_table_name);
  }
  HDF5_ERROR(file = H5Fopen(nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_EOS_HDF5(NAME,VAR,TYPE,MEM)                                      \
  do {                                                                        \
    hid_t dataset;                                                            \
    HDF5_ERROR(dataset = H5Dopen(file, NAME));                                \
    HDF5_ERROR(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR));       \
    HDF5_ERROR(H5Dclose(dataset));                                            \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_EOSTABLE_HDF5(NAME,OFF)                                     \
  do {                                                                   \
    hsize_t offset[2]     = {OFF,0};                                     \
    H5Sselect_hyperslab(mem3, H5S_SELECT_SET, offset, NULL, var3, NULL); \
    READ_EOS_HDF5(NAME,alltables_temp,H5T_NATIVE_DOUBLE,mem3);           \
  } while (0)

  // Read size of tables
  READ_EOS_HDF5("pointsrho",  &nrho,  H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("pointstemp", &ntemp, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("pointsye",   &nye,   H5T_NATIVE_INT, H5S_ALL);


  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(nrho * ntemp * nye * NTABLES * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }
  if (!(logrho = (double*)malloc(nrho * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }
  if (!(logtemp = (double*)malloc(ntemp * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }
  if (!(yes = (double*)malloc(nye * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NTABLES, (hsize_t)nrho * ntemp * nye};
  hsize_t var3[2]       = { 1, (hsize_t)nrho * ntemp * nye};
  hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

  // Read alltables_temp
  READ_EOSTABLE_HDF5("logpress",  0);
  READ_EOSTABLE_HDF5("logenergy", 1);
  READ_EOSTABLE_HDF5("entropy",   2);
  READ_EOSTABLE_HDF5("munu",      3);
  READ_EOSTABLE_HDF5("cs2",       4);
  READ_EOSTABLE_HDF5("dedt",      5);
  READ_EOSTABLE_HDF5("dpdrhoe",   6);
  READ_EOSTABLE_HDF5("dpderho",   7);
  // chemical potentials
  READ_EOSTABLE_HDF5("muhat",     8);
  READ_EOSTABLE_HDF5("mu_e",      9);
  READ_EOSTABLE_HDF5("mu_p",     10);
  READ_EOSTABLE_HDF5("mu_n",     11);
  // compositions
  READ_EOSTABLE_HDF5("Xa",       12);
  READ_EOSTABLE_HDF5("Xh",       13);
  READ_EOSTABLE_HDF5("Xn",       14);
  READ_EOSTABLE_HDF5("Xp",       15);
  // average nucleus
  READ_EOSTABLE_HDF5("Abar",     16);
  READ_EOSTABLE_HDF5("Zbar",     17);
  // Gamma
  READ_EOSTABLE_HDF5("gamma",    18);

  // Read additional tables and variables
  READ_EOS_HDF5("logrho",       logrho,        H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("logtemp",      logtemp,       H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("ye",           yes,            H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("energy_shift", &energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL);

  HDF5_ERROR(H5Sclose(mem3));
  HDF5_ERROR(H5Fclose(file));

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(alltables = (double*)malloc(nrho * ntemp * nye * NTABLES 
				    * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }
  for(int iv = 0;iv<NTABLES;iv++) 
    for(int k = 0; k<nye;k++) 
      for(int j = 0; j<ntemp; j++) 
	for(int i = 0; i<nrho; i++) {
	  int indold = i + nrho*(j + ntemp*(k + nye*iv));
	  int indnew = iv + NTABLES*(i + nrho*(j + ntemp*k));
	  alltables[indnew] = alltables_temp[indold];
	}

  // free memory of temporary array
  free(alltables_temp);

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  energy_shift = energy_shift * EPSGF;
  for(int i=0;i<nrho;i++) {
    // rewrite:
    //logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
    // by using log(a^b*c) = b*log(a)+log(c)
    logrho[i] = logrho[i] * log(10.) + log(RHOGF);
  }

  for(int i=0;i<ntemp;i++) {
    //logtemp[i] = log(pow(10.0,logtemp[i]));
    logtemp[i] = logtemp[i]*log(10.0)+log(TEMPCONV);
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(epstable = (double*)malloc(nrho * ntemp * nye  
				    * sizeof(double)))) {
    printf("Cannot allocate memory for eps table\n");
  }

  if (!(presstable = (double*)malloc(nrho * ntemp * nye  
				    * sizeof(double)))) {
    printf("Cannot allocate memory for press table\n");
  }
  // convert units
  for(int i=0;i<nrho*ntemp*nye;i++) {

    { // pressure
      int idx = 0 + NTABLES*i;
      alltables[idx] = alltables[idx] * log(10.0) + log(PRESSGF);
      presstable[i] = exp(alltables[idx]);
    }

    { // eps
      int idx = 1 + NTABLES*i;
      alltables[idx] = alltables[idx] * log(10.0) + log(EPSGF);
      epstable[i] = exp(alltables[idx]);
    }

    { // cs2
      int idx = 4 + NTABLES*i;
      alltables[idx] *= LENGTHGF*LENGTHGF/TIMEGF/TIMEGF;
    }

    { // dedT
      int idx = 5 + NTABLES*i;
      alltables[idx] *= EPSGF;
    }

    { // dpdrhoe
      int idx = 6 + NTABLES*i;
      alltables[idx] *= PRESSGF/RHOGF;
    }

    { // dpderho
      int idx = 7 + NTABLES*i;
      alltables[idx] *= PRESSGF/EPSGF;
    }

  }

  temp0 = exp(logtemp[0]);
  temp1 = exp(logtemp[1]);

  // set up some vars
  dtemp  = (logtemp[ntemp-1] - logtemp[0]) / (1.0*(ntemp-1));
  dtempi = 1.0/dtemp;

  dlintemp = temp1-temp0;
  dlintempi = 1.0/dlintemp;

  drho  = (logrho[nrho-1] - logrho[0]) / (1.0*(nrho-1));
  drhoi = 1.0/drho;

  dye  = (yes[nye-1] - yes[0]) / (1.0*(nye-1));
  dyei = 1.0/dye;

  drhotempi   = drhoi * dtempi;
  drholintempi = drhoi * dlintempi;
  drhoyei     = drhoi * dyei;
  dtempyei    = dtempi * dyei;
  dlintempyei = dlintempi * dyei;
  drhotempyei = drhoi * dtempi * dyei;
  drholintempyei = drhoi * dlintempi * dyei;

  eos_rhomax = exp(logrho[nrho-1]);
  eos_rhomin = exp(logrho[0]);
  
  eos_tempmax = exp(logtemp[ntemp-1]);
  eos_tempmin = exp(logtemp[0]);

  eos_yemax = yes[nye-1];
  eos_yemin = yes[0];
  
  double epsmax = epstable[0];
  double epsmin = epstable[0];
  for(int i=1;i<nrho*ntemp*nye;i++){
    if ((epstable[i] > epsmax) && (epstable[i] < 1.0e150)){
      epsmax = epstable[i];
    }
    if (epstable[i] < epsmin){
      epsmin = epstable[i];
    }
  }

  double pressmax = presstable[0];
  double pressmin = presstable[0];
  for(int i=1;i<nrho*ntemp*nye;i++){
    if ((presstable[i] > pressmax) && (presstable[i] < 1.0e150)){
      pressmax = presstable[i];
    }
    if (presstable[i] < pressmin){
      pressmin = presstable[i];
    }
  }

  eos_pressmax = pressmax;
  eos_pressmin = pressmin;
  eos_epsmax = epsmax; //- energy_shift;
  eos_epsmin = epsmin; //- energy_shift;
  //fprintf(stderr, "eos_epsmax %e\n", eos_epsmax);
  //fprintf(stderr, "eos_epsmin %e\n", eos_epsmin);
  //fprintf(stderr, "eos_tempmax %e\n", eos_tempmax);
  //fprintf(stderr, "eos_tempmin %e\n", eos_tempmin);
  //fprintf(stderr, "eos_rhomax %e\n", eos_rhomax);
  //fprintf(stderr, "eos_rhomin %e\n", eos_rhomin);
  eos_epsmax = epsmax- energy_shift;
  eos_epsmin = epsmin- energy_shift;

}

void nuc_eos_c_get_energy_shift(double *energy_shift_fortran,
					    double *eos_tempmin_fortran,
					    double *eos_tempmax_fortran,
					    double *eos_yemin_fortran,
					    double *eos_yemax_fortran,
					    double *eos_rhomin_fortran,
					    double *eos_rhomax_fortran,
					    double *eos_epsmin_fortran,
					    double *eos_epsmax_fortran) {


  double eos_rhomax, eos_rhomin;
  double eos_tempmin, eos_tempmax;
  double eos_yemin, eos_yemax;
  double eos_epsmin, eos_epsmax;

  *energy_shift_fortran = energy_shift;
  *eos_tempmin_fortran = eos_tempmin;
  *eos_tempmax_fortran = eos_tempmax;
  *eos_yemin_fortran = eos_yemin;
  *eos_yemax_fortran = eos_yemax;
  *eos_rhomin_fortran = eos_rhomin;
  *eos_rhomax_fortran = eos_rhomax;
  *eos_epsmin_fortran = eos_epsmin;
  *eos_epsmax_fortran = eos_epsmax;

}


