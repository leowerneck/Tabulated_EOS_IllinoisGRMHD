/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

#include "harm_neutrinos.h"
#include "harm_units.h"

/***************************************************************************/
/***************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     N E U T R I N O    P H Y S I C S    R O U T I N E S

*****************************************************************************/

inline void check_electron_fraction(double *pf);


/******************************************************************************/
/******************************************************************************
 neutrino_source_func():
 ---------------------------
  -- right hand side of the equations of motion that are related to neutrinos; 
  -- this should include 
       1) the heating/cooling terms for the energy-momentum EOM; 
       2) the source term that enters into the electron fraction EOM;

******************************************************************************/
void neutrino_source_func(const NRPyEOS_params *restrict eos_params,
                          double *ph, double *optical_depth, struct of_state *q, struct of_geom *geom,
                          int ii, int jj, int kk, double *dU)
{

#if( EVOLVE_ELECTRON_FRACTION )

  
  double R_function;
  double Q_function;

  neutrino_absorption_heating_rate(eos_params, ph, optical_depth, &R_function, &Q_function) ;

  if (ph[YE]>0.9){
    fprintf(stderr," inside ye checkbound R_function %e \n", R_function); 
    fprintf(stderr," inside ye checkbound Q_function %e \n", Q_function);
  }

  /* Do not forget to add to the source since there may be other contributions : */
  dU[ UU ]  += Q_function * q->ucov[0] ;
  dU[ U1 ]  += Q_function * q->ucov[1] ;
  dU[ U2 ]  += Q_function * q->ucov[2] ;
  dU[ U3 ]  += Q_function * q->ucov[3] ;
  dU[ YE ]  += R_function;

  //fprintf(stderr, "Q_function %e\n", Q_function);
  //fprintf(stderr, "R_fucntion %e\n", R_fucntion);

#endif

  return;

}

/******************************************************************************/
/******************************************************************************
 neutrino_absorption_heating_rate():
 ---------------------------
  -- returns the scalar (local fluid frame) net neutrino heating (if
     positive) or cooling (if negative) rate per unit volume :

  -- usually used in the right hand side of the energy-momentum
     equations of motion that are related to neutrinos


  -- returns the scalar (local fluid frame) net neutrino absorption
     (if positive) or emission (if negative) rate per unit volume :

  -- usually used in the right hand side of the electron-fraction
     density equations of motion;

  -- TODO:   FINISH!

******************************************************************************/
void neutrino_absorption_heating_rate(const NRPyEOS_params *restrict eos_params,
                                      double *ph, double *optical_depth, double *R_code_units, double *Q_code_units)
{

#if( N_OPTICAL_DEPTHS )
  double Q_tot,  R_tot_nu_e, R_tot_anti_nu_e, R_tot; 

  R_tot_nu_e=0.;
  R_tot_anti_nu_e=0.;
  Q_tot=0.;
  R_tot=0.;

  if(ph[RHO]<eos_params->eos_rhomin){
    ph[RHO]=1.1*eos_params->eos_rhomin;
  }
  if(ph[RHO]>eos_params->eos_rhomax){
    ph[RHO]=0.99*eos_params->eos_rhomax;
  }
  if(ph[YE]<eos_params->eos_yemin){
    ph[YE]=1.01*eos_params->eos_yemin;
  }
  if(ph[YE]>eos_params->eos_yemax){
    ph[YE]=0.99*eos_params->eos_yemax;
  }
  if(ph[TEMP]>eos_params->eos_tempmax){
    ph[TEMP]=0.99*eos_params->eos_tempmax;
  }
  if(ph[TEMP]<eos_params->eos_tempmin){
    ph[TEMP]=1.01*eos_params->eos_tempmin;
  }





  double temperature, rho, energy_density, electron_fraction;
  electron_fraction=ph[YE];
  rho=ph[RHO]*INVRHOGF;
  energy_density=(ph[UU]/ph[RHO])*INVEPSGF;
  temperature=ph[TEMP]*1e6;// the factor is to convert to eV 

  double muhat_degeneracy,neutrino_degeneracy, anti_neutrino_degeneracy, neutron_degeneracy,proton_degeneracy,electron_degeneracy;

  muhat_degeneracy=0.;
  proton_degeneracy=0.;
  neutron_degeneracy=0.;
  electron_degeneracy=0.;

  double eta__nu_e=0.;
  double eta__nu_e_bar=0.;
  double eta__nu_mu=0.;
  double eta__nu_mu_bar=0.;
  double eta__nu_tau=0.;
  double eta__nu_tau_bar=0.;

   

  get_degeneracy_all(eos_params, rho, energy_density, temperature, electron_fraction, &electron_degeneracy, &proton_degeneracy, 
                     &neutron_degeneracy, &muhat_degeneracy, &eta__nu_e, &eta__nu_e_bar, &eta__nu_mu, &eta__nu_mu_bar, &eta__nu_tau, &eta__nu_tau_bar);


  for(int i_optdepth=0; i_optdepth < N_OPTICAL_DEPTHS/2; i_optdepth++){


    int __attribute__((unused)) neutrino_flavor, anti_neutrino_flavor;

    neutrino_flavor=i_optdepth;

    if (i_optdepth%2 == 0){ //neutrino
      anti_neutrino_flavor=i_optdepth+1;
    }
    else if (i_optdepth%2 == 1){ //neutrino_bar
      anti_neutrino_flavor=i_optdepth-1;
    }

    switch (neutrino_flavor) {
    case nu__e:

      neutrino_degeneracy = eta__nu_e; 
      anti_neutrino_degeneracy = eta__nu_e_bar;
      break;

    case nu__e_bar:

      neutrino_degeneracy = eta__nu_e_bar; 
      anti_neutrino_degeneracy = eta__nu_e;
      break;

    case nu__mu:
            
      neutrino_degeneracy = eta__nu_mu; 
      anti_neutrino_degeneracy = eta__nu_mu_bar;
      break;

    case nu__mu_bar:

      neutrino_degeneracy = eta__nu_mu_bar; 
      anti_neutrino_degeneracy = eta__nu_mu;
      break;

    case nu__tau:

      neutrino_degeneracy = eta__nu_tau; 
      anti_neutrino_degeneracy = eta__nu_tau_bar;
      break;

    case nu__tau_bar:

      neutrino_degeneracy = eta__nu_tau_bar; 
      anti_neutrino_degeneracy = eta__nu_tau;
      break;
    }


    double neutron_fraction,proton_fraction;
    neutron_fraction=get_particle_fraction(NEUTRON, electron_fraction);
    proton_fraction=get_particle_fraction(PROTON, electron_fraction);


    double __attribute__((unused)) Tau_Q, Tau_R;
    Tau_Q=optical_depth[N_OPTICAL_DEPTHS/2+i_optdepth];
    Tau_R=optical_depth[i_optdepth];


    Q_tot=Q_tot+effective_emission_rate(eos_params, neutrino_flavor, NEUTRINO_ENERGY, rho, energy_density, temperature,
                                        electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction, neutrino_degeneracy,
                                        anti_neutrino_degeneracy, neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_Q);


    if(neutrino_flavor==0){
      R_tot_nu_e=effective_emission_rate(eos_params, neutrino_flavor, NEUTRINO_NUMBER, rho, energy_density, temperature,
                                         electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction, neutrino_degeneracy, 
                                         anti_neutrino_degeneracy, neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_R);
    }
    else if (neutrino_flavor==1){
      R_tot_anti_nu_e=effective_emission_rate(eos_params, neutrino_flavor, NEUTRINO_NUMBER, rho, energy_density, temperature,
                                              electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction, neutrino_degeneracy, 
                                              anti_neutrino_degeneracy, neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_R);
    }

  }

  Q_tot=-Q_tot; //ruffert B25
  *Q_code_units=Q_tot*Q_CONV; //need to change units from eV/(cm3s) to code units.
  R_tot=(-R_tot_nu_e + R_tot_anti_nu_e)*C_amu; // Ruffert eq. 14 and 15 considering that we are evolving Ye*rho this is assuming we are using SRO tables that have the density in terms of the mass of the neutron
  *R_code_units=R_tot*R_CONV;

  // printf("(harm) ***************************************************\n");
  // printf("(harm) Before unit changes (i.e., everything in cgs):\n");
  // printf("(harm) rho      : %.15e\n",rho);
  // printf("(harm) R_tot    : %.15e\n",R_tot);
  // printf("(harm) Q_tot    : %.15e\n",Q_tot);
  // printf("(harm) R_tot/rho: %.15e\n",R_tot/rho);
  // printf("(harm) Q_tot/rho: %.15e\n",Q_tot/rho);
  // printf("(harm) ***************************************************\n");
  return;

#endif

  return;

}

void neutrino_electron_Q(const NRPyEOS_params *restrict eos_params,
                         double *ph, double *optical_depth, double *Q_code_units )
{

#if( N_OPTICAL_DEPTHS )
  int i_optdepth;
  double __attribute__((unused)) Q_tot,  R_tot_nu_e, R_tot_anti_nu_e, R_tot; 

  Q_tot=0.;
  i_optdepth=0;



  int __attribute__((unused)) neutrino_flavor, anti_neutrino_flavor;
  
  neutrino_flavor=i_optdepth;

  if (i_optdepth%2 == 0){ //neutrino
    anti_neutrino_flavor=i_optdepth+1;
  }
  else if (i_optdepth%2 == 1){ //neutrino_bar
    anti_neutrino_flavor=i_optdepth-1;
  }

  double temperature, rho, energy_density, electron_fraction;
  electron_fraction=ph[YE];
  rho=ph[RHO]*INVRHOGF;
  energy_density=(ph[UU]/ph[RHO])*INVEPSGF;
  temperature=ph[TEMP]*1e6;// the factor is to convert to eV 

  double muhat_degeneracy,neutrino_degeneracy, anti_neutrino_degeneracy, neutron_degeneracy,proton_degeneracy,electron_degeneracy;

  muhat_degeneracy=0.;
  proton_degeneracy=0.;
  neutron_degeneracy=0.;
  electron_degeneracy=0.;

  double eta__nu_e=0.;
  double eta__nu_e_bar=0.;
  double eta__nu_mu=0.;
  double eta__nu_mu_bar=0.;
  double eta__nu_tau=0.;
  double eta__nu_tau_bar=0.;

   

  get_degeneracy_all(eos_params, rho, energy_density, temperature, electron_fraction, &electron_degeneracy, &proton_degeneracy, 
                     &neutron_degeneracy, &muhat_degeneracy, &eta__nu_e, &eta__nu_e_bar, &eta__nu_mu, &eta__nu_mu_bar, &eta__nu_tau, &eta__nu_tau_bar);
  switch (neutrino_flavor) {
  case nu__e:

    neutrino_degeneracy = eta__nu_e; 
    anti_neutrino_degeneracy = eta__nu_e_bar;
    break;

  case nu__e_bar:

    neutrino_degeneracy = eta__nu_e_bar; 
    anti_neutrino_degeneracy = eta__nu_e;
    break;

  case nu__mu:
            
    neutrino_degeneracy = eta__nu_mu; 
    anti_neutrino_degeneracy = eta__nu_mu_bar;
    break;

  case nu__mu_bar:

    neutrino_degeneracy = eta__nu_mu_bar; 
    anti_neutrino_degeneracy = eta__nu_mu;
    break;

  case nu__tau:

    neutrino_degeneracy = eta__nu_tau; 
    anti_neutrino_degeneracy = eta__nu_tau_bar;
    break;

  case nu__tau_bar:

    neutrino_degeneracy = eta__nu_tau_bar; 
    anti_neutrino_degeneracy = eta__nu_tau;
    break;
  }


  double neutron_fraction,proton_fraction;
  neutron_fraction=get_particle_fraction(NEUTRON, electron_fraction);
  proton_fraction=get_particle_fraction(PROTON, electron_fraction);


  double __attribute__((unused)) Tau_Q, Tau_R;
  Tau_Q=optical_depth[N_OPTICAL_DEPTHS/2+i_optdepth];

  


  Q_tot=Q_tot+effective_emission_rate(eos_params, neutrino_flavor, NEUTRINO_ENERGY, rho, energy_density, temperature,
                                      electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction, neutrino_degeneracy, 
                                      anti_neutrino_degeneracy, neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_Q);


  Q_tot=Q_tot; //ruffert B25
  *Q_code_units=Q_tot; //need to change units from eV/(cm3s) to code units.
  return;

#endif

  return;

}



void antineutrino_electron_Q(const NRPyEOS_params *restrict eos_params,
                             double *ph, double *optical_depth, double *Q_code_units )
{

#if( N_OPTICAL_DEPTHS )
  int i_optdepth;
  double __attribute__((unused)) Q_tot,  R_tot_nu_e, R_tot_anti_nu_e, R_tot; 

  Q_tot=0.;
  i_optdepth=1;



  int __attribute__((unused)) neutrino_flavor, anti_neutrino_flavor;
  
  neutrino_flavor=i_optdepth;

  if (i_optdepth%2 == 0){ //neutrino
    anti_neutrino_flavor=i_optdepth+1;
  }
  else if (i_optdepth%2 == 1){ //neutrino_bar
    anti_neutrino_flavor=i_optdepth-1;
  }

  double temperature, rho, energy_density, electron_fraction;
  electron_fraction=ph[YE];
  rho=ph[RHO]*INVRHOGF;
  energy_density=(ph[UU]/ph[RHO])*INVEPSGF;
  temperature=ph[TEMP]*1e6;// the factor is to convert to eV 

  double muhat_degeneracy,neutrino_degeneracy, anti_neutrino_degeneracy, neutron_degeneracy,proton_degeneracy,electron_degeneracy;

  muhat_degeneracy=0.;
  proton_degeneracy=0.;
  neutron_degeneracy=0.;
  electron_degeneracy=0.;

  double eta__nu_e=0.;
  double eta__nu_e_bar=0.;
  double eta__nu_mu=0.;
  double eta__nu_mu_bar=0.;
  double eta__nu_tau=0.;
  double eta__nu_tau_bar=0.;

   

  get_degeneracy_all(eos_params, rho, energy_density, temperature, electron_fraction, &electron_degeneracy, &proton_degeneracy, 
                     &neutron_degeneracy, &muhat_degeneracy, &eta__nu_e, &eta__nu_e_bar, &eta__nu_mu, &eta__nu_mu_bar, &eta__nu_tau, &eta__nu_tau_bar);
  switch (neutrino_flavor) {
  case nu__e:

    neutrino_degeneracy = eta__nu_e; 
    anti_neutrino_degeneracy = eta__nu_e_bar;
    break;

  case nu__e_bar:

    neutrino_degeneracy = eta__nu_e_bar; 
    anti_neutrino_degeneracy = eta__nu_e;
    break;

  case nu__mu:
            
    neutrino_degeneracy = eta__nu_mu; 
    anti_neutrino_degeneracy = eta__nu_mu_bar;
    break;

  case nu__mu_bar:

    neutrino_degeneracy = eta__nu_mu_bar; 
    anti_neutrino_degeneracy = eta__nu_mu;
    break;

  case nu__tau:

    neutrino_degeneracy = eta__nu_tau; 
    anti_neutrino_degeneracy = eta__nu_tau_bar;
    break;

  case nu__tau_bar:

    neutrino_degeneracy = eta__nu_tau_bar; 
    anti_neutrino_degeneracy = eta__nu_tau;
    break;
  }


  double neutron_fraction,proton_fraction;
  neutron_fraction=get_particle_fraction(NEUTRON, electron_fraction);
  proton_fraction=get_particle_fraction(PROTON, electron_fraction);


  double __attribute__((unused)) Tau_Q, Tau_R;
  Tau_Q=optical_depth[N_OPTICAL_DEPTHS/2+i_optdepth];


  Q_tot=Q_tot+effective_emission_rate(eos_params, neutrino_flavor, NEUTRINO_ENERGY, rho, energy_density, temperature,
                                      electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction, neutrino_degeneracy, 
                                      anti_neutrino_degeneracy, neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_Q);


  //   }

  Q_tot=Q_tot; //ruffert B25
  *Q_code_units=Q_tot; //need to change units from eV/(cm3s) to code units.
  return;

#endif

  return;

}



void neutrino_x_Q(const NRPyEOS_params *restrict eos_params,
                  double *ph, double *optical_depth, double *Q_code_units )
{

#if( N_OPTICAL_DEPTHS )
  int i_optdepth;
  double __attribute__((unused)) Q_tot,  R_tot_nu_e, R_tot_anti_nu_e, R_tot; 

  Q_tot=0.;
  i_optdepth=5;



  int __attribute__((unused)) neutrino_flavor, anti_neutrino_flavor;
  
  neutrino_flavor=i_optdepth;

  if (i_optdepth%2 == 0){ //neutrino
    anti_neutrino_flavor=i_optdepth+1;
  }
  else if (i_optdepth%2 == 1){ //neutrino_bar
    anti_neutrino_flavor=i_optdepth-1;
  }

  double temperature, rho, energy_density, electron_fraction;
  electron_fraction=ph[YE];
  rho=ph[RHO]*INVRHOGF;
  energy_density=(ph[UU]/ph[RHO])*INVEPSGF;
  temperature=ph[TEMP]*1e6;// the factor is to convert to eV 

  double muhat_degeneracy,neutrino_degeneracy, anti_neutrino_degeneracy, neutron_degeneracy,proton_degeneracy,electron_degeneracy;

  muhat_degeneracy=0.;
  proton_degeneracy=0.;
  neutron_degeneracy=0.;
  electron_degeneracy=0.;

  double eta__nu_e=0.;
  double eta__nu_e_bar=0.;
  double eta__nu_mu=0.;
  double eta__nu_mu_bar=0.;
  double eta__nu_tau=0.;
  double eta__nu_tau_bar=0.;

   

  get_degeneracy_all(eos_params, rho, energy_density, temperature, electron_fraction, &electron_degeneracy, &proton_degeneracy, 
                     &neutron_degeneracy, &muhat_degeneracy, &eta__nu_e, &eta__nu_e_bar, &eta__nu_mu, &eta__nu_mu_bar, &eta__nu_tau, &eta__nu_tau_bar);
  switch (neutrino_flavor) {
  case nu__e:

    neutrino_degeneracy = eta__nu_e; 
    anti_neutrino_degeneracy = eta__nu_e_bar;
    break;

  case nu__e_bar:

    neutrino_degeneracy = eta__nu_e_bar; 
    anti_neutrino_degeneracy = eta__nu_e;
    break;

  case nu__mu:
            
    neutrino_degeneracy = eta__nu_mu; 
    anti_neutrino_degeneracy = eta__nu_mu_bar;
    break;

  case nu__mu_bar:

    neutrino_degeneracy = eta__nu_mu_bar; 
    anti_neutrino_degeneracy = eta__nu_mu;
    break;

  case nu__tau:

    neutrino_degeneracy = eta__nu_tau; 
    anti_neutrino_degeneracy = eta__nu_tau_bar;
    break;

  case nu__tau_bar:

    neutrino_degeneracy = eta__nu_tau_bar; 
    anti_neutrino_degeneracy = eta__nu_tau;
    break;
  }


  double neutron_fraction,proton_fraction;
  neutron_fraction=get_particle_fraction(NEUTRON, electron_fraction);
  proton_fraction=get_particle_fraction(PROTON, electron_fraction);


  double __attribute__((unused)) Tau_Q, Tau_R;
  Tau_Q=optical_depth[N_OPTICAL_DEPTHS/2+i_optdepth];


  Q_tot=Q_tot+effective_emission_rate(eos_params, neutrino_flavor, NEUTRINO_ENERGY, rho, energy_density, temperature,
                                      electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction, neutrino_degeneracy, 
                                      anti_neutrino_degeneracy, neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_Q);



  Q_tot=Q_tot; //ruffert B25
  *Q_code_units=Q_tot; //need to change units from eV/(cm3s) to code units.
  return;

#endif

  return;

}



/******************************************************************************/



/******************************************************************************/
/******************************************************************************
 check_electron_fraction(): 
 ---------------------------

  -- ensures that the resultant electron fraction is within physical
  -- bounds, i.e. [0,1]  ;

  -- TODO:  
        -- is it ok if Ye == 0 ?   should there be a non-zero floor on it? 

******************************************************************************/
void check_electron_fraction( double *pf)
{

#if( EVOLVE_ELECTRON_FRACTION )
  register double ye = pf[YE];
  
  if(pf[RHO]<eos_params->eos_rhomin){
    pf[RHO]=1.1*eos_params->eos_rhomin;
  }
  if(pf[RHO]>eos_params->eos_rhomax){
    pf[RHO]=0.99*eos_params->eos_rhomax;
  }
  if(pf[YE]<eos_params->eos_yemin){
    pf[YE]=1.01*eos_params->eos_yemin;
  }
  if(pf[YE]>eos_params->eos_yemax){
    pf[YE]=0.99*eos_params->eos_yemax;
  }
  if(pf[TEMP]>eos_params->eos_tempmax){
    pf[TEMP]=0.99*eos_params->eos_tempmax;
  }
  if(pf[TEMP]<eos_params->eos_tempmin){
    pf[TEMP]=1.01*eos_params->eos_tempmin;
  }


  if( isfinite(ye) ) {
    if     ( ye < 0. ) {  ye = YEMIN;  }
    else if( ye > 1. ) {  
      fprintf(stderr," inside ye checkbound pf[ye] %e \n", pf[YE]);
      fprintf(stderr," inside ye checkbound pf[rho] %e \n", pf[RHO]);
      fprintf(stderr," inside ye checkbound pf[temp] %e \n", pf[TEMP]);
      ye = eos_params->eos_yemax;  
    }
  }
  else {
    ye = YEMIN ;
  }

  pf[YE] = ye; 
#endif

  return;
  
}

/******************************************************************************************/
/******************************************************************************************
----------------------------------    REFERENCES    ---------------------------------------

[1] "Coalescing neutron stars - a step towards physical models"
    Ruffert et. al.
    Astron. Astrophys. 311, 532-566 (1996)


[2] "Implementation of a simplified approach to radiative transfer in general relativity"
     Galeazzi et. al.
     2 Oct. 2013


[3] "Beta Transition Rates in Hot and Dense Matter"
    K. Takahashi, M. F. El Eid, and W. Hillebrandt
    Astron. Astrophys. 67, 185-197 (1978)


HISTORY: 
-- Code written by Michael Kolacki (2017)  v1 


-----------------------------------------------------------------------

-- Notes (by Scott Noble): 
    -- Adding equation numbers from different papers [N] to identify
       where the equations below come from;  Comments look like 
             [1](B10) + 
       If there is a "+" then I have confirmed the equation. 

    -- look out for ?? for equations I'm not sure about 
    -- look out for !! for equations that I think are wrong; 

    -- !! need to make sure that arguments from one routine are passed
          consistently to the next routine, for all routines; I'm only
          checking the math expressions right now;

  --ari:
   Q is in units of eV/(cm^3s) and R in units 1/(cm^3s)

*******************************************************************************************/


   //#include "neutrino_cooling.h"



   /******************************************************************************************/
   /******************************************************************************************
    get_Temperature():
    -----------------
    Expected Parameters:
        -- energy_density   - in CGS units
        -- rho              - in CGS units

    Return Value:
        -- determines the Temperature (in eV) based on energy density and density
           values, with the ability to implement other equations of state.

        -- EOS should be set in neutrino_cooling.h

   *******************************************************************************************/
     double get_Temperature(double energy_density, double rho, double electron_fraction) {
     
     /*
      * T = ((Gamma - 1) * u) / (rho * C_k)
      *
      * C_k   :=  Boltzmann constant
      * Gamma :=  5/3
      * u,rho :=  Energy density and density, respectively
      *           Measured in CGS units
      * T     :=  Temperature
      *           Measured in eV
      */

     return ((2.0 / 3.0) * (energy_density / rho)) / C_k;



     // default:
     //    return 0.0;
     //    break;

   }


/******************************************************************************************/
/******************************************************************************************
    get_particle_fraction():
    -----------------------
    Expected Parameters:
        -- particle_type    - either ELECTRON, PROTON, or NEUTRON (in decs.h)
        -- energy_density   - in CGS units
        -- rho              - in CGS units
        -- temperature      - in eV - obtained from the get_Temperature() function
        -- ye -primitive value

    Return Value:
        -- ideally, determines the particle fraction for the particular fermion based on the
           assumptions of the matter observed. In the case of matter which is assumed to be
           completely disassociated, the case for NEUTRON and PROTON particle types should
           not need adjusting, but currently, the ELECTRON particle type only returns a
           constant value.

           This wll need adjusting.

           Parameters which should allow for this adjustment are included.

*******************************************************************************************/
double get_particle_fraction(int particle_type, double ye) {
  switch (MATTER_DISTRIBUTION_ASSUMPTION) {
  case COMPLETELY_DISSOCIATED:
    switch (particle_type) {
    case ELECTRON:
      return ye;  //ari: changed it to ye
      break;
	
    case PROTON:
      return get_particle_fraction(ELECTRON, ye);
      break;
	
    case NEUTRON:
      return 1.0 - get_particle_fraction(ELECTRON, ye);
      break;
	
    default:
      break;
	
    }
      
    break;


  default:
    break;
  }
    
  return( 1. );  /* HEREHERE  what should the default value be?? */
}


/******************************************************************************************/
/******************************************************************************************
    get_degeneracy():
    ----------------
    Expected Parameters:
        -- particle_type    - any of the described particle types (in decs.h)
        -- temperature      - in eV - from the get_Temperature() function
        -- energy_density   - in CGS units
        -- rho              - in CGS units

    Return Value:
        -- Determines the degeneracy parameter for each of the different particle types.
           For particle types of neutrinos, degeneracy (eta) is computed as:

           eta = eta_i,CEQ 

           where eta_i,CEQ is the degeneracy of the particular neutrino flavor at chemical
           equilibrium and Tau is the optical depth for the point of the observation.

           For particle types of fermions, degeneracy (eta) is computed as:

           eta = mu_i / T

           where mu_i is the chemical potential for particle i, and T is the temperature at
           the point of the observation. Currently, the function returns only constant values
           for ELECTRON, PROTON, and NEUTRON degeneracy, based on average values across runs
           seen in various papers.

           This will need adjusting.

           Return values of 0.0 for the 4 heavy neutrinos however, is based on the paper of
           Ruffert et. al. referenced at the top of the file (as well as several others, and
           may not be necessary to adjust.

*******************************************************************************************/
double get_degeneracy(const NRPyEOS_params *restrict eos_params,
                      int particle_type, double rho, double energy_density, double temperature, double electron_fraction) {
  /*
   * temperature   :=  Measured in eV
   *
   */
  double eta__nu_e__ceq, eta_e, eta_p, eta_n, eta_muhat;
  double __attribute__((unused)) Q;


  double xeps = energy_density*EPSGF;  //to convert from cgs to table units
  double xtemp = temperature*1e-6; // to convert from eV to table units (MeV)
  double xP = 0.;
  double rho0=rho*RHOGF; //to convert from cgs to table units
  double xmu_e = 0.0;
  double xmu_p = 0.0;
  double xmu_n = 0.0;
  double xmu_hat = 0.0;

  //fprintf(stderr, "xtemp table in get_degeneracy %e\n", xtemp);
  //fprintf(stderr, "rho0 table in get_degeneracy %e\n", rho0);

  // EOS_press_mu_with_T(rho0, &xeps, &xtemp, electron_fraction, &xP, &xmu_hat, &xmu_e, &xmu_p, &xmu_n);
  NRPyEOS_P_eps_mue_mup_mun_and_muhat_from_rho_Ye_T(eos_params,rho0,electron_fraction,xtemp,
                                                    &xP,&xeps,&xmu_e,&xmu_p,&xmu_n,&xmu_hat);

  //fprintf(stderr, "rho0 table %e\n", rho0);
  //fprintf(stderr, "xtemp table %e\n", xtemp);
  //fprintf(stderr, "mu_e table %e\n", xmu_e);

  eta_e= xmu_e/xtemp;
  eta_p= xmu_p/xtemp;
  eta_n= xmu_n/xtemp;
  eta_muhat=xmu_hat/xtemp;

  switch (particle_type) {
  case nu__e:
    /*
     * Electron, proton, and neutron degeneracies.
     * Used to determine an equilibrium degeneracy for electron neutrinos.
     */

    //eta_e = get_degeneracy(ELECTRON, temperature, optical_depth);
    //eta_p = get_degeneracy(PROTON, temperature, optical_depth);
    //eta_n = get_degeneracy(NEUTRON, temperature, optical_depth);


    /*
     * Q := Rest-mass-energy difference between a neutron and a proton
     *      Measured in eV
     */

    Q = 1.2935e6; //[1] in text of appendix A ari+

    /*
     * eta__nu_e__ceq := degeneracy parameter for electron neutrinos at chemical equilibrium
     */

    eta__nu_e__ceq = eta_e - eta_muhat; //- (Q / temperature); // [1] A5 ari+

    return( eta__nu_e__ceq ); //[1] A3, and eta0_nu_e=0 ari+ Considering Galeazzi
    break;
  case nu__e_bar:
    /*
     * Electron, proton, and neutron degeneracies.
     * Used to determine an equilibrium degeneracy for electron anti-neutrinos.
     */

    //eta_e = get_degeneracy(ELECTRON, temperature, optical_depth);
    //eta_p = get_degeneracy(PROTON, temperature, optical_depth);
    //eta_n = get_degeneracy(NEUTRON, temperature, optical_depth);


    /*
     * Q := Rest-mass-energy difference between a neutron and a proton
     *      Measured in eV
     */

    Q = 1.2935e6;

    /*
     * eta__nu_e__ceq := degeneracy parameter for electron neutrinos at chemical equilibrium
     */

    eta__nu_e__ceq = eta_e - eta_muhat;// - (Q / temperature); // [1] A5 ari+

    return( -eta__nu_e__ceq ); //[1] A4, and eta0_nu_e_bar=0 ari+ Considering Galeazzi
    break;
  case nu__mu:
    return 0.0; //[1] A2 ari+
    break;
  case nu__mu_bar:
    return 0.0;//[1] A2 ari+
    break;
  case nu__tau:
    return 0.0; //[1] A2 ari+
    break;
  case nu__tau_bar:
    return 0.0; //[1] A2 ari+
    break;
  case NEUTRON:
    return eta_n;//0.8; // 
    break;
  case PROTON:
    return eta_p;//1.0;
    break;
  case ELECTRON:
    return eta_e;//2.5;
    break;
  case MUHAT:
    return eta_muhat;
    break;
  default:
    return 0.0; // [1] A2
    break;

  }

}



/******************************************************************************************/
/******************************************************************************************
    get_degeneracy_all():
    ----------------
    Expected Parameters:
        -- particle_type    - any of the described particle types (in decs.h)
        -- temperature      - in eV - from the get_Temperature() function
        -- energy_density   - in CGS units
        -- rho              - in CGS units

    Return Value:
        -- Determines the degeneracy parameter for each of the different particle types.
           For particle types of neutrinos, degeneracy (eta) is computed as:

           eta = eta_i,CEQ 

           where eta_i,CEQ is the degeneracy of the particular neutrino flavor at chemical
           equilibrium and Tau is the optical depth for the point of the observation.

           For particle types of fermions, degeneracy (eta) is computed as:

           eta = mu_i / T



           Return values of 0.0 for the 4 heavy neutrinos however, is based on the paper of
           Ruffert et. al. referenced at the top of the file (as well as several others, and
           may not be necessary to adjust.

*******************************************************************************************/
void get_degeneracy_all(const NRPyEOS_params *restrict eos_params,
                        double rho, double energy_density, double temperature, double electron_fraction, 
                        double * eta_e, double * eta_p, double * eta_n, double * eta_muhat, double * eta__nu_e, double * eta__nu_e_bar,
                        double * eta__nu_mu, double * eta__nu_mu_bar, double * eta__nu_tau, double * eta__nu_tau_bar) {
  /*
   * temperature   :=  Measured in eV
   *
   */
  double __attribute__((unused)) Q;


  double xeps = energy_density*EPSGF;  //to convert from cgs to table units
  double xtemp = temperature*1e-6; // to convert from eV to table units (MeV)
  double xP = 0.;
  double rho0=rho*RHOGF; //to convert from cgs to table units
  double xmu_e = 0.0;
  double xmu_p = 0.0;
  double xmu_n = 0.0;
  double xmu_hat = 0.0;

  //fprintf(stderr, "xtemp table in get_degeneracy %e\n", xtemp);
  //fprintf(stderr, "rho0 table in get_degeneracy %e\n", rho0);

  // EOS_press_mu_with_T(rho0, &xeps, &xtemp, electron_fraction, &xP, &xmu_hat, &xmu_e, &xmu_p, &xmu_n);
  NRPyEOS_P_eps_mue_mup_mun_and_muhat_from_rho_Ye_T(eos_params,rho0,electron_fraction,xtemp,
                                                    &xP,&xeps,&xmu_e,&xmu_p,&xmu_n,&xmu_hat);

  //fprintf(stderr, "rho0 table %e\n", rho0);
  //fprintf(stderr, "xtemp table %e\n", xtemp);
  //fprintf(stderr, "mu_e table %e\n", xmu_e);

  *eta_e= xmu_e/xtemp;
  *eta_p= xmu_p/xtemp;
  *eta_n= xmu_n/xtemp;
  *eta_muhat=xmu_hat/xtemp;

  Q = 1.2935e6; //[1] in text of appendix A ari+

  /*
   * eta__nu_e__ceq := degeneracy parameter for electron neutrinos at chemical equilibrium
   */

  *eta__nu_e = xmu_e/xtemp - xmu_hat/xtemp; //- (Q / temperature); // [1] A5 ari+
  *eta__nu_e_bar = -(xmu_e/xtemp - xmu_hat/xtemp);// - (Q / temperature); // [1] A5 ari+
  *eta__nu_mu = 0.0; //[1] A2 ari+
  *eta__nu_mu_bar = 0.0;//[1] A2 ari+
  *eta__nu_tau = 0.0; //[1] A2 ari+
  *eta__nu_tau_bar = 0.0;

  // fprintf(stderr,"(harm) mu_e                = %.15e\n",xmu_e);
  // fprintf(stderr,"(harm) mu_p                = %.15e\n",xmu_p);
  // fprintf(stderr,"(harm) mu_n                = %.15e\n",xmu_n);
  // fprintf(stderr,"(harm) muhat               = %.15e\n",xmu_hat);
  // fprintf(stderr,"(harm) eta_e               = %.15e\n",*eta_e);
  // fprintf(stderr,"(harm) eta_p               = %.15e\n",*eta_p);
  // fprintf(stderr,"(harm) eta_n               = %.15e\n",*eta_n);
  // fprintf(stderr,"(harm) etahat              = %.15e\n",*eta_muhat);
  // fprintf(stderr,"(harm) eta_nue             = %.15e\n",*eta__nu_e);
  // fprintf(stderr,"(harm) eta_anue            = %.15e\n",*eta__nu_e_bar);
  
  return;

}


/******************************************************************************************/
/******************************************************************************************
    get_Fermi_integral():
    --------------------
    Expected Parameters:
        -- order    - the order of the integral being approximated
        -- eta      - the degeneracy/relativistic chemical potential of the particular
                      species being investigated.

    Return Value:
        -- approximates the Fermi integral, based on the methods described by

           K. Takahashi, M. F. El Eid, and W. Hillebrandt, Astron. Astrophys. 67, 185 (1978).

           Used in papers by Ruffert et. al. (referenced at the top of the file) and
           "Implementation of a simplified approach to radiative transfer in general
           relativity" by Galeazzi et. al.

*******************************************************************************************/
double get_Fermi_integral(int order, double eta) {

  if (eta > 0.001) {
    /*  HEREHERE -- can optimize this by not using calls to pow() for integer powers */ 
    switch (order) {
    case 0:
      return eta + log(1.0 + exp(-eta)); //Takahashi A4
      break;
    case 1:
      return ((eta*eta / 2.0) + 1.6449) / (1.0 + exp(-1.6855 * eta)); // Takahashi A2 ari+
      break;
    case 2:
      return ((eta*eta*eta / 3.0) + (3.2899 * eta)) / (1.0 - exp(-1.8246 * eta)); // Takahashi A2 ari+
      break;
    case 3:
      return ((pow(eta, 4.0) / 4.0) + (4.9348 * eta*eta) + 11.3644) / (1.0 + exp(-1.9039 * eta)); // Takahashi A2 ari+
      break;
    case 4:
      return ((pow(eta, 5.0) / 5.0) + (6.5797 * eta*eta*eta) + (45.4576 * eta)) / (1.0 - exp(-1.9484 * eta)); // Takahashi A2 ari+
      break;
    case 5:
      return ((pow(eta, 6.0) / 6.0) + (8.2247 * (eta*eta)*(eta*eta)) + (113.6439 * eta*eta) + 236.5323) / (1.0 + exp(-1.9727 * eta)); // Takahashi A2 ari+
      break;
    default:
      /* HEREHERE -- should there be a "default order" in this switch ? */
      fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of order "); 
      fprintf(stdout,"\t order = %d      eta = %26.16e \n", order, eta ); fflush(stdout); 
      // fail(FAIL_BASIC,0);
      exit(1);
      return(0.);
      break;
    }
  } 
  else {
    switch (order) {
    case 0:
      return eta + log(1.0 + exp(-eta)); //Takahashi A4
      break;
    case 1:
      return exp(eta) / (1.0 + 0.2159 * exp(0.8857 * eta)); // Takahashi A3 ari+
      break;
    case 2:
      return (2.0 * exp(eta)) / (1.0 + 0.1092 * exp(0.8908 * eta)); // Takahashi A3 ari+
      break;
    case 3:
      return (6.0 * exp(eta)) / (1.0 + 0.0559 * exp(0.9069 * eta)); // Takahashi A3 ari+
      break;
    case 4:
      return (24.0 * exp(eta)) / (1.0 + 0.0287 * exp(0.9257 * eta)); // Takahashi A3 ari+
      break;
    case 5:
      return (120.0 * exp(eta)) / (1.0 + 0.0147 * exp(0.9431 * eta)); // Takahashi A3 ari+
      break;
    default:
      /* HEREHERE -- should there be a "default order" in this switch ? */
      fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of order "); 
      fprintf(stdout,"\t order = %d      eta = %26.16e \n", order, eta ); fflush(stdout); 
      // fail(FAIL_BASIC,0);
      exit(1);
      return(1.);
      break;
      
    }
    
  }
  
}



/******************************************************************************************/
/******************************************************************************************
    get_neutrino_number_density():
    -----------------------------
    Expected Parameters:
        -- neutrino_flavor  - self explanatory. The flavor of neutrino being observed
        -- degeneracy       - the degeneracy of the particular neutrino being observed, as
                              determined by the get_degeneracy() function.
        -- temperature      - in eV - as determined by the get_Temperature() function.

    Return Value:
        -- returns the neutrino number density for the particular point of observation

*******************************************************************************************/
double get_neutrino_number_density(int neutrino_flavor, double degeneracy, double temperature) {
  /*
   * g__nu_i := Weight for neutrino type
   */

  double g__nu_i, fermi_int;

  switch (neutrino_flavor) {
  case nu__e:
  case nu__e_bar:
    g__nu_i = 1.0;
    break;
  default:
    g__nu_i = 4.0;
    break;
  } /* after [1](B19) + */

  fermi_int = get_Fermi_integral(2, degeneracy);  

  return( g__nu_i * 4.0 * C_pi * (temperature*temperature*temperature/C_hc3_ev) * fermi_int);  /* [1](B18) + */

}


/******************************************************************************************/
/******************************************************************************************
    get_neutrino_energy_density():
    -----------------------------
    Expected Parameters:
        -- neutrino_flavor  - self explanatory. The flavor of neutrino being observed
        -- degeneracy       - the degeneracy of the particular neutrino being observed, as
                              determined by the get_degeneracy() function.
        -- temperature      - in eV - as determined by the get_Temperature() function.

    Return Value:
        -- returns the neutrino energy density for the particular point of observation

*******************************************************************************************/
double get_neutrino_energy_density(int neutrino_flavor, double degeneracy, double temperature) {
  /*
   * g__nu_i := Weight for neutrino type
   */

  double g__nu_i, fermi_int;

  switch (neutrino_flavor) {
  case nu__e:
  case nu__e_bar:
    g__nu_i = 1.0;
    break;
  default:
    g__nu_i = 4.0;
    break;
  }  /* after [1](B19) + */

  fermi_int = get_Fermi_integral(3, degeneracy);

  return g__nu_i * 4.0 * C_pi * (pow(temperature, 4.0) / C_hc3_ev) * fermi_int;  /* [1](B19) + */

}


/******************************************************************************************/
/******************************************************************************************
    get_energy_moment():
    -------------------
    Expected Parameters:
        -- positron_or_electron - The type of particle under consideration.
        -- symbol               - Based on the symbols used in the paper by Ruffert et. al.,
                                  referenced at the top of the file, and defined in
                                  neutrino_cooling.h
        -- temperature          - in eV - as determined by the get_Temperature() function.
        -- electron_degeneracy  - degeneracy of the electrons at the point of observation.
                                  This IS also the intended parameter in the case of a
                                  POSITRON particle type.

    Return Value:
        -- returns the particular energy moment of either POSITRONs or ELECTRONs at the point
           of observation, with the specified temperature and degeneracy.
           These energy moments are used in the calculation of neutrino emission rates for
           each of three different processes.

*******************************************************************************************/
double get_energy_moment(int positron_or_electron, int symbol, double temperature, double electron_degeneracy) {
  double fermi_int;

  switch (positron_or_electron) {
  case ELECTRON:
    switch (symbol) {

    case ORDINARY:
      fermi_int = get_Fermi_integral(3, electron_degeneracy);
      return (8.0 * C_pi) * (pow(temperature, 4.0) / C_hc3_ev) * fermi_int; // [1] B5 ari+
      break;

    case TILDE:
      fermi_int = get_Fermi_integral(4, electron_degeneracy);
      return (8.0 * C_pi) * (pow(temperature, 5.0) / C_hc3_ev) * fermi_int; // [1] B6 ari+
      break;

    case STAR:
      fermi_int = get_Fermi_integral(5, electron_degeneracy);
      return (8.0 * C_pi) * (pow(temperature, 6.0) / C_hc3_ev) * fermi_int; // [1] B7 ari+
      break;

    default:
      break;
    }
      
    break;

  case POSITRON:
    switch (symbol) {

    case ORDINARY:
      fermi_int = get_Fermi_integral(3, -electron_degeneracy);
      return (8.0 * C_pi) * (pow(temperature, 4.0) / C_hc3_ev) * fermi_int; // [1] B5 ari+
      break;

    case TILDE:
      fermi_int = get_Fermi_integral(4, -electron_degeneracy);
      return (8.0 * C_pi) * (pow(temperature, 5.0) / C_hc3_ev) * fermi_int; // [1] B6 ari+
      break;

    case STAR:
      fermi_int = get_Fermi_integral(5, -electron_degeneracy);
      return (8.0 * C_pi) * (pow(temperature, 6.0) / C_hc3_ev) * fermi_int; // [1] B7 ari+
      break;

    default:
      break;
    }
      
    break;

  default:
    break;
  }

  /* HEREHERE -- should there be a "default value " in this switch ? */
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of positron_or_electron  or symbol "); 
  fprintf(stdout,"\t positron_or_electron = %d      symbol = %d \n", positron_or_electron, symbol); fflush(stdout); 
  // fail(FAIL_BASIC,0);
  exit(1);
  return(0.); 

}


/******************************************************************************************/
/******************************************************************************************
    get_scattering_coefficient():
    ----------------------------
    Expected Parameters:
        -- nucleon  - determines the coefficient to return.

    Return Value:
        -- returns the coefficient used in the calculation of the spectrally averaged
           scattering opacity for the respective nucleon

*******************************************************************************************/
double get_scattering_coefficient(int nucleon) {

  switch (nucleon) {

  case NEUTRON:
    /*
     * (1 + 5 * C_alpha_neutrino_sq) / 24
     *  where C_alpha_neutrino_sq  = 1.25*1.25
     */
    return 0.3671875; // [1] In the beggining of appendix A ari+
    break;

  case PROTON:
    /*
     * (4 * (C_V - 1)^2 + 5 * C_alpha_neutrino_sq) / 24
     *  where C_alpha_neutrino_sq  = 1.25*1.25
     *  and C_V = 1/2 + 2 * sin^2 (theta_w)
     *  with sin^2 (theta_w) = .23
     */
    return 7.8189 / 24.0; // [1] In the beggining of appendix A , ari: fixed number, original: 11.4989
    break;

  default:
    break;
  }
  
  /* HEREHERE -- should there be a "default value " in this switch ? */
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of nucleon "); 
  fprintf(stdout,"\t nucleon = %d \n", nucleon); fflush(stdout); 
  // fail(FAIL_BASIC,0);
  exit(1);
  return(0.); 
  
}


/******************************************************************************************/
/******************************************************************************************
    get_phase_space_blocking_value():
    --------------------------------
    Expected Parameters:
        -- interaction          - type of interaction for which the blocking factor is desired.
        -- particle_type        - necessary for blocking factors from scattering and Beta
                                  process decay.
                                  Irrelevant for other blocking factors.
        -- electron_degeneracy  - degeneracy of electrons at the point of observation.
                                  Irrelevant for Beta process decay.
        -- neutrino_degeneracy  - degeneracy of the flavor of neutrinos under observation at
                                  the point of observation

    Return Value:
        -- returns the approximated phase space blocking of leptons during scattering and for
           each electron-positron pair annihilation, Beta process decay, and plasmon decay

*******************************************************************************************/
double get_phase_space_blocking_value(int interaction, int particle_type, double electron_degeneracy, double neutrino_degeneracy) {
  double fermi_int_order_5;
  double fermi_int_order_4;
  double electron_fermi_int_order_4;
  double electron_fermi_int_order_3;
  double positron_fermi_int_order_4;
  double positron_fermi_int_order_3;
  double gamma_0;
  double gamma;

  switch (interaction) {

  case SCATTERING:

    switch (particle_type) {

    case ELECTRON:
      fermi_int_order_5 = get_Fermi_integral(5, neutrino_degeneracy);
      fermi_int_order_4 = get_Fermi_integral(4, neutrino_degeneracy);
      const double bf   = 1.0 / (1.0 + exp(-((fermi_int_order_5 / fermi_int_order_4) - electron_degeneracy)));
      return( bf ); //[1] A15 ari+
      break;
	
    case POSITRON:
      fermi_int_order_5 = get_Fermi_integral(5, neutrino_degeneracy); 
      fermi_int_order_4 = get_Fermi_integral(4, neutrino_degeneracy);
      return 1.0 / (1.0 + exp(-((fermi_int_order_5 / fermi_int_order_4) + electron_degeneracy))); //[1] A16
      break;
	
    default:
      break;
    }
    break;
      

  case BETA_PROCESS:
    switch (particle_type) {
    case NEUTRINO:
      fermi_int_order_5 = get_Fermi_integral(5, electron_degeneracy);    /* [1](B3)  ??? Should this not use electron_degeneracy (ari: yes, it should, fixed it) */
      fermi_int_order_4 = get_Fermi_integral(4, electron_degeneracy);    /* [1](B3)  ??? Should this not use electron_degeneracy */
      // fprintf(stderr,"(harm) eta_e               = %.15e\n",electron_degeneracy);
      // fprintf(stderr,"(harm) eta_nue             = %.15e\n",neutrino_degeneracy);
      // fprintf(stderr,"(harm) FD_5_nue            = %.15e\n",fermi_int_order_5);
      // fprintf(stderr,"(harm) FD_4_nue            = %.15e\n",fermi_int_order_4);
      // fprintf(stderr,"(harm) Ratio               = %.15e\n",fermi_int_order_5/fermi_int_order_4);
      // fprintf(stderr,"(harm) arg of exp          = %.15e\n",-((fermi_int_order_5 / fermi_int_order_4) - electron_degeneracy));
      return 1.0 / (1.0 + exp(-((fermi_int_order_5 / fermi_int_order_4) - neutrino_degeneracy)));  /* [1](B3)  */
      break;
	
    case ANTI_NEUTRINO:
      fermi_int_order_5 = get_Fermi_integral(5, -electron_degeneracy);
      fermi_int_order_4 = get_Fermi_integral(4, -electron_degeneracy);
      return 1.0 / (1.0 + exp(-((fermi_int_order_5 / fermi_int_order_4) - neutrino_degeneracy)));  /* [1](B3) */
      break;
	
    default:
      break;
    }
    break;

      
  case ELECTRON_POSITRON_PAIR_ANNIHILATION:
    electron_fermi_int_order_4 = get_Fermi_integral(4, electron_degeneracy);
    electron_fermi_int_order_3 = get_Fermi_integral(3, electron_degeneracy);
    positron_fermi_int_order_4 = get_Fermi_integral(4, -electron_degeneracy);
    positron_fermi_int_order_3 = get_Fermi_integral(3, -electron_degeneracy);
    return 1.0 / (1.0 + exp(-((0.5 * (electron_fermi_int_order_4 / electron_fermi_int_order_3))
                              + (0.5 * (positron_fermi_int_order_4 / positron_fermi_int_order_3)) - neutrino_degeneracy))); // [1] B9 ari+
    break;
      

  case PLASMON_DECAY:
    /*
     * alpha     := Fine structure constant
     *
     * gamma_0   := Plasma frequency
     *            = h_bar * Omega_0 / (m_e * c^2) = 2 * sqrt(alpha / (3 * Pi))
     *            = 5.565 * 10^-2
     */
    gamma_0 = 0.05565;
    gamma = gamma_0 * sqrt((1.0 / 3.0) * (pow(C_pi, 2.0) + 3.0 * pow(electron_degeneracy, 2.0)));
    return 1.0 / (1.0 + exp(-(1.0 + 0.5 * (pow(gamma, 2.0) / (1.0 + gamma)) - neutrino_degeneracy)));  //[1] B13 ari+
    break;
      

  default:
    break;

  }

  /* HEREHERE -- should there be a "default value " in this switch ? */
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of interaction or particle_type "); 
  fprintf(stdout,"\t interaction  = %d      particle_type = %d \n", interaction, particle_type); fflush(stdout); 
  // fail(FAIL_BASIC,0);
  exit(1);
  return(0.); 

}


/******************************************************************************************/
/******************************************************************************************
    get_spectrally_averaged_absorption_opacity():
    --------------------------------------------
    Expected Parameters:
        -- neutrino_flavor      - the particular flavor of neutrino currently under
                                  observation
        -- nucleon              - the particular nucleon for which scattering is being
                                  computed for
        -- type_of_transport    - used to differentiate between whether it is desired to
                                  calculate the neutrino number or neutrino energy emission rate
        -- rho                  - in CGS units
        -- temperature          - in eV - from get_Temperature() function
        -- electron_fraction    - computed from the get_particle_fraction() function
        -- neutron_degeneracy   - computed from the get_degeneracy() function
        -- proton_degeneracy    - computed from the get_degeneracy() function
        -- electron_degeneracy  - computed from the get_degeneracy() function
        -- neutrino_degeneracy  - computed from the get_degeneracy() function, utilizing the
                                  previously calculated value for optical depth

    Return Value:
        -- returns explicitly the spectrally averaged absorption opacity for the particular
           flavor of neutrino provided, ensuring that no improper pairing occurs (i.e.: nu__e
           neutrinos with a PROTON nucleon).

*******************************************************************************************/
double get_spectrally_averaged_absorption_opacity(int neutrino_flavor, int nucleon, int type_of_transport,
                                                  double rho, double temperature, double electron_fraction,
                                                  double muhat_degeneracy, 
                                                  double electron_degeneracy, double neutrino_degeneracy) {
  double Y_np;
  double Y_pn;
  double fermi_int1;
  double fermi_int2;
  double blocking_factor;

  switch (neutrino_flavor) {
  case nu__e:
    if (nucleon == NEUTRON) {
      if(electron_fraction<0.5){
        Y_np = (2.0 * electron_fraction - 1.0) / (exp(-muhat_degeneracy) - 1.0);  /* [1](A13) + */
      }
      else{
        Y_np=get_particle_fraction(NEUTRON, electron_fraction);
      }
                

      if (Y_np<0) {
        Y_np=get_particle_fraction(NEUTRON, electron_fraction);
      }

      fermi_int1 = get_Fermi_integral((4 + type_of_transport), neutrino_degeneracy);
      fermi_int2 = get_Fermi_integral((2 + type_of_transport), neutrino_degeneracy);

      blocking_factor = get_phase_space_blocking_value(SCATTERING, ELECTRON, electron_degeneracy, neutrino_degeneracy);
      return ((1.0 + 3.0 * C_alpha_neutrino_sq) / 4.0) * C_A_sigma_0 * rho * Y_np * 
        (temperature/ (C_me_c2_ev))*(temperature/ (C_me_c2_ev)) *
        (fermi_int1 / fermi_int2) * blocking_factor;    /* [1](A11) + */

    }
    else {
      fprintf(stderr, "Bad combination of values for  nucleon and neutrino_flavor\n");
    }

    break;

  case nu__e_bar:
    if (nucleon == PROTON) {
      if(electron_fraction>0.5){
        Y_pn = exp(-muhat_degeneracy) *
          ((2.0 * electron_fraction - 1.0) / (exp(-muhat_degeneracy) - 1.0)); /* [1](A14) + */
      }
      else{
        Y_pn=get_particle_fraction(PROTON, electron_fraction);
      }

      if (Y_pn<0){
        Y_pn=get_particle_fraction(PROTON, electron_fraction);
      }

      fermi_int1 = get_Fermi_integral((4 + type_of_transport), neutrino_degeneracy);
      fermi_int2 = get_Fermi_integral((2 + type_of_transport), neutrino_degeneracy);

      blocking_factor = get_phase_space_blocking_value(SCATTERING, POSITRON, electron_degeneracy, neutrino_degeneracy);

      return ((1.0 + 3.0 * C_alpha_neutrino_sq) *0.25) * C_A_sigma_0 * rho * Y_pn * 
        (temperature/ (C_me_c2_ev))*(temperature/ (C_me_c2_ev)) *
        (fermi_int1 / fermi_int2) * blocking_factor;  /* [1](A12) + */

    }
    else {
      fprintf(stderr, "Bad combination of values for  nucleon and neutrino_flavor \n");
    }

    break;

  default:
    fprintf(stderr, "Bad value of neutrino_flavor\n");
    break;
  }

  fprintf(stderr, "Bad value of absorption opacity \n");
  return(0.);

}


/******************************************************************************************/
/******************************************************************************************
    get_spectrally_averaged_scattering_opacity():
    --------------------------------------------
    Expected Parameters:
        -- neutrino_flavor      - the particular flavor of neutrino currently under
                                  observation
        -- nucleon              - the particular nucleon involved in the scattering
        -- type_of_transport    - used to differentiate between whether it is desired to
                                  calculate the neutrino number or neutrino energy emission rate
        -- rho                  - in CGS units
        -- temperature          - in eV - from get_Temperature() function
        -- nucleon_fraction     - computed from the get_particle_fraction() function
        -- nucleon_degeneracy   - computed from the get_degeneracy() function
        -- neutrino_degeneracy  - computed from the get_degeneracy() function, utilizing the
                                  previously calculated value for optical depth

    Return Value:
        -- returns explicitly the spectrally averaged scattering opacity for the particular
           flavor of neutrino provided

*******************************************************************************************/
double get_spectrally_averaged_scattering_opacity(int neutrino_flavor, int nucleon, int type_of_transport,
                                                  double rho, double temperature, double nucleon_fraction,
                                                  double nucleon_degeneracy, double neutrino_degeneracy) {

  double C_N, fermi_int1, fermi_int2, Y_NN, max;

  C_N = get_scattering_coefficient(nucleon);  /*  after [1](A1)  */
  fermi_int1 = get_Fermi_integral((4 + type_of_transport), neutrino_degeneracy);
  fermi_int2 = get_Fermi_integral((2 + type_of_transport), neutrino_degeneracy);
  max = (nucleon_degeneracy >= 0.0) ? nucleon_degeneracy : 0.0;   /* [1](A8) + */
  Y_NN = nucleon_fraction / (1.0 + ((2.0 / 3.0) * max));   /* [1](A8) + */

  return C_N * C_A_sigma_0 * rho * Y_NN * 
    (temperature / (C_me_c2_ev))*(temperature / (C_me_c2_ev)) * (fermi_int1 / fermi_int2);  /* [1](A6) + */

}


/******************************************************************************************/
/******************************************************************************************
    get_total_transport_opacity():
    -----------------------------
    Expected Parameters:
        -- neutrino_flavor      - the particular flavor of neutrino currently under
                                  observation
        -- type_of_transport    - used to differentiate between whether it is desired to
                                  calculate the neutrino number or neutrino energy emission rate
        -- neutrino_degeneracy  - computed from the get_degeneracy() function, utilizing the
                                  previously calculated value for optical depth
        -- rho                  - in CGS units
        -- temperature          - in eV - from get_Temperature() function
        -- electron_fraction    - computed from the get_particle_fraction() function
        -- neutron_fraction     - computed from the get_particle_fraction() function
        -- proton_fraction      - computed from the get_particle_fraction() function
        -- electron_degeneracy  - computed from the get_degeneracy() function


    Return Value:
        -- returns the total transport opacity for the particular flavor of neutrino provided,
           considering scattering opacity due to each proton and neutron scattering, and where
           appropriate, proton or neutron absorption.

*******************************************************************************************/
double get_total_transport_opacity(const NRPyEOS_params *restrict eos_params,
                                   int neutrino_flavor, int type_of_transport, 
                                   double rho, double energy_density, double electron_fraction, double temperature) {
  double kappa_s_n, kappa_s_p, kappa_a_n, kappa_a_p;

  //double neutrino_degeneracy,electron_degeneracy,muhat_degeneracy, proton_degeneracy, neutron_degeneracy;
  double neutron_fraction,proton_fraction;

  //fprintf(stderr, "temperature get_total_transport_opacity %e\n", temperature);
  //fprintf(stderr, "ye %e\n", electron_fraction);
  //fprintf(stderr, "rho %e\n", rho);

  double muhat_degeneracy,neutrino_degeneracy=0, neutron_degeneracy,proton_degeneracy,electron_degeneracy;
  double __attribute__((unused)) anti_neutrino_degeneracy;

  muhat_degeneracy=0.;
  proton_degeneracy=0.;
  neutron_degeneracy=0.;
  electron_degeneracy=0.;

  double eta__nu_e=0.;
  double eta__nu_e_bar=0.;
  double eta__nu_mu=0.;
  double eta__nu_mu_bar=0.;
  double eta__nu_tau=0.;
  double eta__nu_tau_bar=0.;

  double SMALL_kappa=1e-15;
  
  get_degeneracy_all(eos_params, rho, energy_density, temperature, electron_fraction, &electron_degeneracy, &proton_degeneracy, 
                     &neutron_degeneracy, &muhat_degeneracy,&eta__nu_e, &eta__nu_e_bar, &eta__nu_mu, &eta__nu_mu_bar, &eta__nu_tau, &eta__nu_tau_bar);
  switch (neutrino_flavor) {
  case nu__e:

    neutrino_degeneracy = eta__nu_e; 
    anti_neutrino_degeneracy = eta__nu_e_bar;
    break;
  case nu__e_bar:

    neutrino_degeneracy = eta__nu_e_bar; 
    anti_neutrino_degeneracy = eta__nu_e;
    break;
  case nu__mu:
            
    neutrino_degeneracy = eta__nu_mu; 
    anti_neutrino_degeneracy = eta__nu_mu_bar;
    break;
  case nu__mu_bar:

    neutrino_degeneracy = eta__nu_mu_bar; 
    anti_neutrino_degeneracy = eta__nu_mu;
    break;
  case nu__tau:

    neutrino_degeneracy = eta__nu_tau; 
    anti_neutrino_degeneracy = eta__nu_tau_bar;
    break;
  case nu__tau_bar:

    neutrino_degeneracy = eta__nu_tau_bar; 
    anti_neutrino_degeneracy = eta__nu_tau;
    break;
  }
    
  //neutrino_degeneracy=get_degeneracy(neutrino_flavor, rho,  energy_density, temperature,electron_fraction);
  //proton_degeneracy=get_degeneracy(PROTON, rho,  energy_density, temperature,electron_fraction);
  //neutron_degeneracy=get_degeneracy(NEUTRON, rho,  energy_density, temperature,electron_fraction);
  //muhat_degeneracy=get_degeneracy(MUHAT, rho,  energy_density, temperature,electron_fraction);
  //electron_degeneracy=get_degeneracy(ELECTRON, rho,  energy_density, temperature,electron_fraction);


  neutron_fraction=get_particle_fraction(NEUTRON, electron_fraction);
  proton_fraction=get_particle_fraction(PROTON, electron_fraction);



  switch (neutrino_flavor) {
  case nu__e:
    kappa_s_n = get_spectrally_averaged_scattering_opacity(neutrino_flavor, NEUTRON, type_of_transport, rho, temperature,
                                                           neutron_fraction, neutron_degeneracy, neutrino_degeneracy); /* [1](A6) */
    kappa_s_p = get_spectrally_averaged_scattering_opacity(neutrino_flavor, PROTON, type_of_transport, rho, temperature,
                                                           proton_fraction, proton_degeneracy, neutrino_degeneracy);  /* [1](A6) */
    kappa_a_n = get_spectrally_averaged_absorption_opacity(neutrino_flavor, NEUTRON, type_of_transport, rho, temperature,
                                                           electron_fraction, muhat_degeneracy, electron_degeneracy, neutrino_degeneracy);  /* [1](A11) */
    if(!isfinite(kappa_s_n)){
      //fprintf(stderr, "kappa_s_n electron %e\n", kappa_s_n);
      //fprintf(stderr, "electron_fraction %e\n", electron_fraction);
      //fprintf(stderr, "muhat_degeneracy %e\n", muhat_degeneracy);
      //fprintf(stderr, "rho %e\n", rho);
      //fprintf(stderr, "temperature %e\n", temperature);
      kappa_s_n=SMALL_kappa;
    }
    if(!isfinite(kappa_a_n)){
      //fprintf(stderr, "kappa_a_n electron %e\n", kappa_a_n);
      //fprintf(stderr, "electron_fraction %e\n", electron_fraction);
      //fprintf(stderr, "muhat_degeneracy %e\n", muhat_degeneracy);
      //fprintf(stderr, "rho %e\n", rho);
      //fprintf(stderr, "temperature %e\n", temperature);
      kappa_a_n=SMALL_kappa;
    }
    if(!isfinite(kappa_s_p)){
      //fprintf(stderr, "kappa_s_p electron %e\n", kappa_s_p);
      //fprintf(stderr, "electron_fraction %e\n", electron_fraction);
      //fprintf(stderr, "muhat_degeneracy %e\n", muhat_degeneracy);
      //fprintf(stderr, "rho %e\n", rho);
      //fprintf(stderr, "temperature %e\n", temperature);
      kappa_s_p=SMALL_kappa;
    }
    return kappa_s_n + kappa_s_p + kappa_a_n;  /* [1](A17) + */

    break;

  case nu__e_bar:
    kappa_s_n = get_spectrally_averaged_scattering_opacity(neutrino_flavor, NEUTRON, type_of_transport, rho, temperature,
                                                           neutron_fraction, neutron_degeneracy, neutrino_degeneracy);  /* [1](A6) */
    kappa_s_p = get_spectrally_averaged_scattering_opacity(neutrino_flavor, PROTON, type_of_transport, rho, temperature,
                                                           proton_fraction, proton_degeneracy, neutrino_degeneracy);    /* [1](A6) */
    kappa_a_p = get_spectrally_averaged_absorption_opacity(neutrino_flavor, PROTON, type_of_transport, rho, temperature,
                                                           electron_fraction, muhat_degeneracy, electron_degeneracy, neutrino_degeneracy);   /* [1](A12) */

    //fprintf(stderr, "kappa_s_n %e\n", kappa_s_n);
    //fprintf(stderr, "kappa_s_p %e\n", kappa_s_p);
    //fprintf(stderr, "kappa_a_p %e\n", kappa_a_p);
    //fprintf(stderr, "kappa_tot %e\n", kappa_s_n + kappa_s_p + kappa_a_p);
    if(!isfinite(kappa_s_n)){
      //fprintf(stderr, "kappa_s_n anti %e\n", kappa_s_n);
      //fprintf(stderr, "electron_fraction %e\n", electron_fraction);
      //fprintf(stderr, "muhat_degeneracy %e\n", muhat_degeneracy);
      //fprintf(stderr, "rho %e\n", rho);
      //fprintf(stderr, "temperature %e\n", temperature);
      kappa_s_n=SMALL_kappa;
    }
    if(!isfinite(kappa_a_p)){
      //fprintf(stderr, "kappa_a_p anti %e\n", kappa_a_p);
      //fprintf(stderr, "electron_fraction %e\n", electron_fraction);
      //fprintf(stderr, "muhat_degeneracy %e\n", muhat_degeneracy);
      //fprintf(stderr, "rho %e\n", rho);
      //fprintf(stderr, "temperature %e\n", temperature);
      kappa_a_p=SMALL_kappa;
    }
    if(!isfinite(kappa_s_p)){
      // fprintf(stderr, "kappa_s_p anti %e\n", kappa_s_p);
      // fprintf(stderr, "electron_fraction %e\n", electron_fraction);
      // fprintf(stderr, "muhat_degeneracy %e\n", muhat_degeneracy);
      // fprintf(stderr, "rho %e\n", rho);
      // fprintf(stderr, "temperature %e\n", temperature);
      kappa_s_p=SMALL_kappa;
    }

    return kappa_s_n + kappa_s_p + kappa_a_p;  /* [1](A18) + */

    break;

  default:
    kappa_s_n = get_spectrally_averaged_scattering_opacity(neutrino_flavor, NEUTRON, type_of_transport, rho, temperature,
                                                           neutron_fraction, neutron_degeneracy, neutrino_degeneracy);  /* [1](A6) */
    kappa_s_p = get_spectrally_averaged_scattering_opacity(neutrino_flavor, PROTON, type_of_transport, rho, temperature,
                                                           proton_fraction, proton_degeneracy, neutrino_degeneracy);    /* [1](A6) */

    //fprintf(stderr, "kappa_s_n %e\n", kappa_s_n);
    //fprintf(stderr, "kappa_s_p %e\n", kappa_s_p);

    if(!isfinite(kappa_s_n)){
      kappa_s_n=SMALL_kappa;
    }
    if(!isfinite(kappa_s_p)){
      kappa_s_p=SMALL_kappa;
    }
    return kappa_s_n + kappa_s_p;  /* [1](A19) + */

    break;

  }

  /* HEREHERE -- should there be a "default value " in this switch ? */
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Inside get_total_transport_opacity\n");
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of neutrino_flavor "); 
  fprintf(stdout,"\t neutrino_flavor  = %d \n", neutrino_flavor); fflush(stdout); 
  // fail(FAIL_BASIC,0);
  exit(1);
  return(0.); 

}


/******************************************************************************************/
/******************************************************************************************
    beta_process_emission_rate__neutrino_number():
    ---------------------------------------------
    Expected Parameters:
        -- neutrino_flavor      - the particular flavor of neutrino currently under
                                  observation
        -- rho                  - in CGS units
        -- temperature          - in eV - from get_Temperature() function
        -- electron_fraction    - computed from the get_particle_fraction() function
        -- neutrino_degeneracy  - computed from the get_degeneracy() function, utilizing the
                                  previously calculated value for optical depth
        -- neutron_degeneracy   - computed from the get_degeneracy() function
        -- proton_degeneracy    - computed from the get_degeneracy() function
        -- electron_degeneracy  - computed from the get_degeneracy() function

    Return Value:
        -- returns the total number of neutrinos emitted from a particular location due to
           Beta process decay,  R_\beta from  [1](B1-B4)

*******************************************************************************************/
double beta_process_emission_rate__neutrino_number(int neutrino_flavor, double rho, double temperature,
                                                   double electron_fraction, double neutrino_degeneracy,
                                                   double muhat_degeneracy, double electron_degeneracy) {

  double Y_pn, Y_np, energy_moment, blocking_factor;


  if (neutrino_flavor == nu__e) {
    if(electron_fraction>0.5){
      Y_pn = exp(-muhat_degeneracy) *
        ((2.0 * electron_fraction - 1.0) / (exp(-muhat_degeneracy) - 1.0)); /* [1](A14) + */
    }
    else{
      Y_pn=get_particle_fraction(PROTON, electron_fraction);
    }/* [1](A14) + */

    // printf("(harm) ***************** R_beta_nue *****************\n");
    energy_moment = get_energy_moment(ELECTRON, TILDE, temperature, electron_degeneracy);
    blocking_factor = get_phase_space_blocking_value(BETA_PROCESS, NEUTRINO, electron_degeneracy, neutrino_degeneracy); /* [1](B3) + */
    const double R_beta = ((1.0 + 3.0 * C_alpha_neutrino_sq) / 8.0) * (C_A_sigma_0 / (C_me2_c3_ev)) * rho * Y_pn * energy_moment * blocking_factor;
    // fprintf(stderr,"(harm) blocking_factor     = %.15e\n",blocking_factor);
    // fprintf(stderr,"(harm) energy_moment_star  = %.15e\n",energy_moment_star);
    // fprintf(stderr,"(harm) energy_moment_tilde = %.15e\n",energy_moment_tilde);
    // fprintf(stderr,"(harm) R_beta_nue          = %.15e\n",R_beta);
    // fprintf(stderr,"(harm) **********************************************\n");
    
    return R_beta;  /* [1](B1) + */

  }
  else if (neutrino_flavor == nu__e_bar) {
    if(electron_fraction<0.5){
      Y_np = (2.0 * electron_fraction - 1.0) / (exp(-muhat_degeneracy) - 1.0);  /* [1](A13) + */
    }
    else{
      Y_np=get_particle_fraction(NEUTRON, electron_fraction);
    }  /* [1](A13) + */
    energy_moment = get_energy_moment(POSITRON, TILDE, temperature, electron_degeneracy);
    blocking_factor = get_phase_space_blocking_value(BETA_PROCESS, ANTI_NEUTRINO, electron_degeneracy, neutrino_degeneracy); /* [1](B4) + */

    /* !! I added a factor of rho in the next line because I think it was missing:  */ 
    return ((1.0 + 3.0 * C_alpha_neutrino_sq) / 8.0) * (C_A_sigma_0 / (C_me2_c3_ev)) * rho * Y_np * energy_moment * blocking_factor;  /* [1](B2) +  */

  }
  else {
    /* HEREHERE -- should there be a "default value " in this switch ? */
    // fprintf(stdout,"(harm3d_neutrinos - ERROR) Inside beta_process_emission_rate__neutrino_number\n");
    // fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of neutrino_flavor "); 
    // fprintf(stdout,"\t neutrino_flavor  = %d \n", neutrino_flavor); fflush(stdout); 
    // fail(FAIL_BASIC,0);
    // exit(1);
    return(0.); 
  }

}


/******************************************************************************************/
/******************************************************************************************
    beta_process_emission_rate__neutrino_energy():
    ---------------------------------------------
    Expected Parameters:
        -- neutrino_flavor      - the particular flavor of neutrino currently under
                                  observation
        -- rho                  - in CGS units
        -- temperature          - in eV - from get_Temperature() function
        -- electron_fraction    - computed from the get_particle_fraction() function
        -- neutrino_degeneracy  - computed from the get_degeneracy() function, utilizing the
                                  previously calculated value for optical depth
        -- neutron_degeneracy   - computed from the get_degeneracy() function
        -- proton_degeneracy    - computed from the get_degeneracy() function
        -- electron_degeneracy  - computed from the get_degeneracy() function

    Return Value:
        -- returns the total energy of neutrinos emitted from a particular location due to
           Beta process decay

*******************************************************************************************/
double beta_process_emission_rate__neutrino_energy(int neutrino_flavor, double rho, double temperature,
                                                   double electron_fraction, double neutrino_degeneracy,
                                                   double muhat_degeneracy, double electron_degeneracy) {

  double energy_moment_star, energy_moment_tilde;
  double R_beta;
  double __attribute__((unused)) Y_pn, Y_np, blocking_factor;


  if (neutrino_flavor == nu__e) {


    //Y_pn = exp(proton_degeneracy - neutron_degeneracy) *
    //                         ((2.0 * electron_fraction - 1.0) / (exp(proton_degeneracy - neutron_degeneracy) - 1.0)); /* [1](A14) + */
    energy_moment_star = get_energy_moment(ELECTRON, STAR, temperature, electron_degeneracy);
    energy_moment_tilde = get_energy_moment(ELECTRON, TILDE, temperature, electron_degeneracy);
    //blocking_factor = get_phase_space_blocking_value(BETA_PROCESS, NEUTRINO, electron_degeneracy, neutrino_degeneracy);

    /* !! I added a factor of rho in the next line because I think it was missing:  */ 
    /* !! why don't we just call beta_process_emission_rate__neutrino_number() because it's the same function modulo the energy_moment ari: done  */
    R_beta=beta_process_emission_rate__neutrino_number(neutrino_flavor, rho, temperature, electron_fraction, neutrino_degeneracy, muhat_degeneracy, electron_degeneracy);
    //return ((1.0 + 3.0 * C_alpha_neutrino_sq) / 8.0) * (C_A_sigma_0 / (C_me2_c3)) * rho * Y_pn * energy_moment * blocking_factor;  /* [1](B14) + */
    
    return R_beta*energy_moment_star/energy_moment_tilde; /* [1](B14) */

  }
  else if (neutrino_flavor == nu__e_bar) {
    //Y_np = (2.0 * electron_fraction - 1.0) / (exp(proton_degeneracy - neutron_degeneracy) - 1.0);  /* [1](A13) + */
    energy_moment_star = get_energy_moment(POSITRON, STAR, temperature, electron_degeneracy);
    energy_moment_tilde = get_energy_moment(POSITRON, TILDE, temperature, electron_degeneracy);
    //blocking_factor = get_phase_space_blocking_value(BETA_PROCESS, ANTI_NEUTRINO, electron_degeneracy, neutrino_degeneracy);

    R_beta=beta_process_emission_rate__neutrino_number(neutrino_flavor, rho, temperature, electron_fraction, neutrino_degeneracy, muhat_degeneracy, electron_degeneracy);

    /* !! I added a factor of rho in the next line because I think it was missing:  */ 
    /* !! why don't we just call beta_process_emission_rate__neutrino_number() because it's the same function modulo the energy_moment   ari: done*/
    //return ((1.0 + 3.0 * C_alpha_neutrino_sq) / 8.0) * (C_A_sigma_0 / (C_me2_c3)) * rho * Y_np * energy_moment * blocking_factor;  /* [1](B15) +  */

    // printf("(harm) R_beta_anue = %.15e\n",R_beta*energy_moment_star/energy_moment_tilde);
    
    return R_beta*energy_moment_star/energy_moment_tilde; /* [1](B15)  */

  }
  else {
    /* HEREHERE -- should there be a "default value " in this switch ? */
    // fprintf(stdout,"(harm3d_neutrinos - ERROR) Inside beta_process_emission_rate__neutrino_energy\n");
    // fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of neutrino_flavor "); 
    // fprintf(stdout,"\t neutrino_flavor  = %d \n", neutrino_flavor); fflush(stdout); 
    // fail(FAIL_BASIC,0);
    // exit(1);
    return(0.); 
  }

}


/******************************************************************************************/
/******************************************************************************************
    electron_positron_pair_annihilation_emission_rate__neutrino_number():
    --------------------------------------------------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- temperature              - in eV - from get_Temperature() function
        -- electron_degeneracy      - computed from the get_degeneracy() function
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- anti_neutrino_degeneracy - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth

    Return Value:
        -- returns the total number of neutrinos emitted from a particular location due to
           electron-positron pair annihilation

*******************************************************************************************/
double electron_positron_pair_annihilation_emission_rate__neutrino_number(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                                          double neutrino_degeneracy, double anti_neutrino_degeneracy) {

  double electron_energy_moment, positron_energy_moment, C_values, neutrino_blocking_factor, anti_neutrino_blocking_factor;

  electron_energy_moment = get_energy_moment(ELECTRON, ORDINARY, temperature, electron_degeneracy);
  positron_energy_moment = get_energy_moment(POSITRON, ORDINARY, temperature, electron_degeneracy);


  if ((neutrino_flavor == nu__e) || (neutrino_flavor == nu__e_bar)) {
    /*
     * C_values  := (C_1 + C_2)_(nu_e, nu_e_bar) = (C_V - C_A)^2 + (C_V + C_A)^2
     *            = 1.92      ???  I think this is wrong  ari: yes, I got 2.3432
     *
     * C_V       := 1/2 + 2 * sin^2 (theta_W)  = 0.96
     *              with sin^2 (theta_W) = .23
     * C_A       := 1/2
     */
    C_values = 2.3432;
    neutrino_blocking_factor =
      get_phase_space_blocking_value(ELECTRON_POSITRON_PAIR_ANNIHILATION, NEUTRINO, electron_degeneracy, neutrino_degeneracy);
    anti_neutrino_blocking_factor =
      get_phase_space_blocking_value(ELECTRON_POSITRON_PAIR_ANNIHILATION, ANTI_NEUTRINO, electron_degeneracy, anti_neutrino_degeneracy);
    //fprintf(stderr, "electron_energy_moment %e\n", electron_energy_moment);
    //fprintf(stderr, "positron_energy_moment %e\n", positron_energy_moment);
    //fprintf(stderr, "neutrino_blocking_factor %e\n", neutrino_blocking_factor);
    //fprintf(stderr, "anti_neutrino_blocking_factor %e\n", anti_neutrino_blocking_factor);

    return( (C_values / 36.0) * (C_sigma_0_cl/ C_me2_c4) *
            electron_energy_moment * positron_energy_moment * neutrino_blocking_factor * anti_neutrino_blocking_factor ); /* [1](B8) */

  }else if ((neutrino_flavor == nu__mu) || (neutrino_flavor == nu__mu_bar) || (neutrino_flavor == nu__tau) || (neutrino_flavor == nu__tau_bar)) {
    /*
     * C_values  := (C_1 + C_2)_(nu_x, nu_x_bar) = (C_V - C_A)^2 + (C_V + C_A - 2)^2
     *            = -0.08     ???  how is this negative ?? ari:  It is not that, it is 0.5032
     *
     * C_V       := 1/2 + 2 * sin^2 (theta_W)
     *              with sin^2 (theta_W) = .23
     * C_A       := 1/2
     */
    C_values = 0.5032;
    neutrino_blocking_factor =
      get_phase_space_blocking_value(ELECTRON_POSITRON_PAIR_ANNIHILATION, NEUTRINO, electron_degeneracy, neutrino_degeneracy);

    return (C_values / 9.0) * (C_sigma_0_cl/ C_me2_c4) *
      electron_energy_moment * positron_energy_moment * neutrino_blocking_factor * neutrino_blocking_factor; // [1] B10

  }
  else { 
    /* HEREHERE -- should there be a "default value " in this switch ? */
    fprintf(stdout,"(harm3d_neutrinos - ERROR) Inside electron_positron_pair_annihilation_emission_rate__neutrino_number\n");
    fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of neutrino_flavor "); 
    fprintf(stdout,"\t neutrino_flavor  = %d \n", neutrino_flavor); fflush(stdout); 
    // fail(FAIL_BASIC,0);
    exit(1);
    return(0.); 
  }

}


/******************************************************************************************/
/******************************************************************************************
    electron_positron_pair_annihilation_emission_rate__neutrino_energy():
    --------------------------------------------------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- temperature              - in eV - from get_Temperature() function
        -- electron_degeneracy      - computed from the get_degeneracy() function
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- anti_neutrino_degeneracy - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth

    Return Value:
        -- returns the total energy of neutrinos emitted from a particular location due to
           electron-positron pair annihilation

*******************************************************************************************/
double electron_positron_pair_annihilation_emission_rate__neutrino_energy(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                                          double neutrino_degeneracy, double anti_neutrino_degeneracy) {

  double ord_electron_energy_moment, ord_positron_energy_moment, tilde_electron_energy_moment, tilde_positron_energy_moment, energy_moment_average;
  double R_ee;
  double __attribute__((unused)) C_values, neutrino_blocking_factor, anti_neutrino_blocking_factor;

  ord_electron_energy_moment = get_energy_moment(ELECTRON, ORDINARY, temperature, electron_degeneracy);
  ord_positron_energy_moment = get_energy_moment(POSITRON, ORDINARY, temperature, electron_degeneracy);
  tilde_electron_energy_moment = get_energy_moment(ELECTRON, TILDE, temperature, electron_degeneracy);
  tilde_positron_energy_moment = get_energy_moment(POSITRON, TILDE, temperature, electron_degeneracy);
  energy_moment_average = .5 * (tilde_electron_energy_moment * ord_positron_energy_moment + ord_electron_energy_moment * tilde_positron_energy_moment);

  R_ee=electron_positron_pair_annihilation_emission_rate__neutrino_number(neutrino_flavor, temperature, electron_degeneracy, neutrino_degeneracy, anti_neutrino_degeneracy);

  return R_ee * energy_moment_average/(ord_electron_energy_moment*ord_positron_energy_moment); //[1] B16 //

}


/******************************************************************************************/
/******************************************************************************************
    plasmon_decay_emission_rate__neutrino_number():
    ----------------------------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- temperature              - in eV - from get_Temperature() function
        -- electron_degeneracy      - computed from the get_degeneracy() function
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- anti_neutrino_degeneracy - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth

    Return Value:
        -- returns the total number of neutrinos emitted from a particular location due to
           plasmon decay

*******************************************************************************************/
double plasmon_decay_emission_rate__neutrino_number(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                    double neutrino_degeneracy, double anti_neutrino_degeneracy) {
  double gamma_0, gamma, C_value, neutrino_blocking_factor, anti_neutrino_blocking_factor;
  /*
   * gamma_0   := Plasma frequency
   *            = h_bar * Omega_0 / (m_e * c^2) = 2 * sqrt(C_alpha_neutrino / (3 * Pi))
   *            = 5.565 * 10^-2
   */

  gamma_0 = 0.05565;
  gamma = gamma_0 * sqrt((1.0 / 3.0) * (pow(C_pi, 2.0) + 3.0 * pow(electron_degeneracy, 2.0))); // [1]below B12

  if ((neutrino_flavor == nu__e) || (neutrino_flavor == nu__e_bar)) {
    /*
     * C_value   := (C_V)^2
     *            = 0.9216 ari: this one is correct
     *
     * C_V       := 1/2 + 2 * sin^2 (theta_W)
     *              with sin^2 (theta_W) = .23
     */

    C_value = 0.9216;
    neutrino_blocking_factor =
      get_phase_space_blocking_value(PLASMON_DECAY, NEUTRINO, electron_degeneracy, neutrino_degeneracy);
    anti_neutrino_blocking_factor =
      get_phase_space_blocking_value(PLASMON_DECAY, ANTI_NEUTRINO, electron_degeneracy, anti_neutrino_degeneracy);

    return ((C_pi * C_pi * C_pi) / (3.0 * C_fine_struct)) * C_value * (C_sigma_0_cl/ C_me2_c4) * (pow(temperature, 8.0) / pow(C_hc_ev, 6.0)) *
      pow(gamma, 6.0) * exp(-gamma) * (1.0 + gamma) * neutrino_blocking_factor * anti_neutrino_blocking_factor; //ari: changed to fine structure constant [1] B11

  }else if ((neutrino_flavor == nu__mu) || (neutrino_flavor == nu__mu_bar) || (neutrino_flavor == nu__tau) || (neutrino_flavor == nu__tau_bar)) {
    /*
     * C_value  := (C_V - 1)^2
     *            = 0.0016 ari: correct too
     *
     * C_V       := 1/2 + 2 * sin^2 (theta_W)
     *              with sin^2 (theta_W) = .23
     */

    C_value = 0.0016;
    neutrino_blocking_factor =
      get_phase_space_blocking_value(PLASMON_DECAY, NEUTRINO, electron_degeneracy, neutrino_degeneracy);

    return ((4.0 * C_pi * C_pi * C_pi) / (3.0 * C_fine_struct)) * C_value * (C_sigma_0_cl/ C_me2_c4) * (pow(temperature, 8.0) / pow(C_hc_ev, 6.0)) *
      pow(gamma, 6.0) * exp(-gamma) * (1.0 + gamma) * neutrino_blocking_factor * neutrino_blocking_factor; // [1] B12 ari:changed constant (the fine structure one)

  }
  else { 
    /* HEREHERE -- should there be a "default value " in this switch ? */
    fprintf(stdout,"(harm3d_neutrinos - ERROR) Inside plasmon_decay_emission_rate__neutrino_number\n");
    fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of neutrino_flavor  "); 
    fprintf(stdout,"\t neutrino_flavor  = %d \n", neutrino_flavor); fflush(stdout); 
    // fail(FAIL_BASIC,0);
    exit(1);
    return(0.); 
  }

}


/******************************************************************************************/
/******************************************************************************************
    plasmon_decay_emission_rate__neutrino_energy():
    ----------------------------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- temperature              - in eV - from get_Temperature() function
        -- electron_degeneracy      - computed from the get_degeneracy() function
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- anti_neutrino_degeneracy - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth

    Return Value:
        -- returns the total energy of neutrinos emitted from a particular location due to
           plasmon decay

*******************************************************************************************/
double plasmon_decay_emission_rate__neutrino_energy(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                    double neutrino_degeneracy, double anti_neutrino_degeneracy) {
  double gamma_0, gamma, R_gamma;
  /*
   * gamma_0   := Plasma frequency
   *            = h_bar * Omega_0 / (m_e * c^2) = 2 * sqrt(alpha / (3 * Pi))
   *            = 5.565 * 10^-2
   */

  gamma_0 = 0.05565;
  gamma = gamma_0 * sqrt((1.0 / 3.0) * (pow(C_pi, 2.0) + 3.0 * pow(electron_degeneracy, 2.0)));
  R_gamma = plasmon_decay_emission_rate__neutrino_number(neutrino_flavor, temperature, electron_degeneracy,
                                                         neutrino_degeneracy, anti_neutrino_degeneracy);

  return 0.5 * temperature * (2.0 + ((gamma * gamma) / (1.0 + gamma))) * R_gamma; //[1] B17

}


/******************************************************************************************/
/******************************************************************************************
    get_total_transport_opacity():
    -----------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- type_of_transport        - used to differentiate between whether it is desired to
                                      calculate the neutrino number or neutrino energy emission rate
        -- previous_optical_depth   - optical depth for the particular flavor of neutrino
                                      under observation up to the current location, and for the
                                      particular type of transport
        -- ds                       - current spacetime element of integration across the
                                      investigation region
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- rho                      - in CGS units
        -- temperature              - in eV - from get_Temperature() function
        -- electron_fraction        - computed from the get_particle_fraction() function
        -- neutron_fraction         - computed from the get_particle_fraction() function
        -- proton_fraction          - computed from the get_particle_fraction() function
        -- electron_degeneracy      - computed from the get_degeneracy() function
        -- neutron_degeneracy       - computed from the get_degeneracy() function
        -- proton_degeneracy        - computed from the get_degeneracy() function

    Return Value:
        -- returns the total integrated optical depth

*******************************************************************************************/
//double get_optical_depth(int neutrino_flavor, int type_of_transport, double previous_optical_depth, double ds,
//                              double neutrino_degeneracy, double rho, double energy_density, double electron_fraction,
//                              double neutron_fraction, double proton_fraction, double electron_degeneracy,
//                              double neutron_degeneracy, double proton_degeneracy) {

//    double kappa_t_j;

//    kappa_t_j = get_total_transport_opacity(neutrino_flavor, type_of_transport,
//					    rho, energy_density, electron_fraction);  /* [1](A17-A19) */ 

//    return previous_optical_depth + (kappa_t_j * ds);  /* [1](A20) ?? Do we want to use "previous_optical_depth" and calculate the total optical depth ?? */ 

//}


/******************************************************************************************/
/******************************************************************************************
    get_diffusion_timescale():
    -------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- type_of_transport        - used to differentiate between whether it is desired to
                                      calculate the neutrino number or neutrino energy emission rate
        -- previous_optical_depth   - optical depth for the particular flavor of neutrino
                                      under observation up to the current location, and for the
                                      particular type of transport
        -- ds                       - current spacetime element of integration across the
                                      investigation region
        -- s_total                  - total "distance" traveled along a geodesic across
                                      spacetime up to the current region of investigation
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- rho                      - in CGS units
        -- temperature              - in eV - from get_Temperature() function
        -- electron_fraction        - computed from the get_particle_fraction() function
        -- neutron_fraction         - computed from the get_particle_fraction() function
        -- proton_fraction          - computed from the get_particle_fraction() function
        -- electron_degeneracy      - computed from the get_degeneracy() function
        -- neutron_degeneracy       - computed from the get_degeneracy() function
        -- proton_degeneracy        - computed from the get_degeneracy() function

    Return Value:
        -- returns the expected timescale for a particular flavor of neutrino to diffuse
           along the provided path

*******************************************************************************************/
double get_diffusion_timescale(const NRPyEOS_params *restrict eos_params,
                               int neutrino_flavor, int type_of_transport, double rho, 
                               double energy_density, double electron_fraction, double temperature,double Tau) {

  double kappa_tot=get_total_transport_opacity(eos_params, neutrino_flavor, type_of_transport, rho, energy_density, 
                                               electron_fraction, temperature);
  //    if (Tau>2){
  //	fprintf(stderr, "diff timescale %e\n", (2.*3.0 * Tau*Tau / C_c) /kappa_tot);		
  //	}
  return (2.*3.0 * Tau*Tau / C_c) /kappa_tot;    /* [1](A21) combined with Rosswog and Liebend�orfer 2003 A32 and O'Connor's 2010 factor of 2*/

}

/******************************************************************************************/
/******************************************************************************************
    get_inverse_emission_timescale():
    --------------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- type_of_transport        - used to differentiate between whether it is desired to
                                      calculate the neutrino number or neutrino energy emission rate
        -- rho                      - in CGS units
        -- temperature              - in eV - from get_Temperature() function
        -- electron_fraction        - computed from the get_particle_fraction() function
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- anti_neutrino_degeneracy - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- neutron_degeneracy       - computed from the get_degeneracy() function
        -- proton_degeneracy        - computed from the get_degeneracy() function
        -- electron_degeneracy      - computed from the get_degeneracy() function

    Return Value:
        -- returns the inverse of the expected timescale for a particular flavor of neutrino
           to emit along the provided path

*******************************************************************************************/
double get_inverse_emission_timescale(int neutrino_flavor, int type_of_transport, double rho, double temperature,
                                      double electron_fraction, double neutrino_degeneracy, double anti_neutrino_degeneracy,
                                      double muhat_degeneracy,
                                      double neutron_degeneracy, double proton_degeneracy, double electron_degeneracy) {

  double number_density, R_Beta, R_ee, R_plasmon;
  double energy_density, Q_Beta, Q_ee, Q_plasmon;
  double SMALL_tot=1e-10;
  //fprintf(stderr, "Temperature get_inverse emission %e\n", temperature);

  switch (type_of_transport) {
  case NEUTRINO_NUMBER:
    number_density = get_neutrino_number_density(neutrino_flavor, neutrino_degeneracy, temperature);  /* [1](B18) + */

    R_Beta = beta_process_emission_rate__neutrino_number(neutrino_flavor, rho, temperature, electron_fraction, neutrino_degeneracy,
                                                         muhat_degeneracy, electron_degeneracy);  /* [1](B1-B4) + */

    R_ee = electron_positron_pair_annihilation_emission_rate__neutrino_number(neutrino_flavor, temperature, electron_degeneracy,
                                                                              neutrino_degeneracy, anti_neutrino_degeneracy); /* [1](B8-B10) + */
    R_plasmon = plasmon_decay_emission_rate__neutrino_number(neutrino_flavor, temperature, electron_degeneracy,
                                                             neutrino_degeneracy, anti_neutrino_degeneracy);  /* R_\gamma  [1](B11-B13) + */

    if(!isfinite(R_Beta)){
      R_Beta=SMALL_tot;
    }           
    if(!isfinite(R_plasmon)){
      R_plasmon=SMALL_tot;
    }
    if(!isfinite(R_ee)){
      R_ee=SMALL_tot;
    }
    if(number_density<1e-20){
      number_density=1e-20; 
    }
    //fprintf(stderr, "R_Beta %e\n", R_Beta);
    //fprintf(stderr, "R_ee %e\n", R_ee);
    //fprintf(stderr, "R_plasmon %e\n", R_plasmon);
    //fprintf(stderr, "number_density %e\n", number_density);

    // if( neutrino_flavor == nu__e ) {
    //   fprintf(stderr,"(harm) Electron neutrino:\n");
    //   fprintf(stderr,"(harm) R_beta_nue         = %.15e\n",R_Beta);
    //   fprintf(stderr,"(harm) R_pair_nue_anue    = %.15e\n",R_ee);
    //   fprintf(stderr,"(harm) R_plasmon_nue_anue = %.15e\n",R_plasmon);
    // }
    // else if( neutrino_flavor == nu__e_bar ) {
    //   fprintf(stderr,"(harm) Electron antineutrino:\n");
    //   fprintf(stderr,"(harm) R_beta_anue        = %.15e\n",R_Beta);
    //   fprintf(stderr,"(harm) R_pair_nue_anue    = %.15e\n",R_ee);
    //   fprintf(stderr,"(harm) R_plasmon_nue_anue = %.15e\n",R_plasmon);
    // }
    // else {
    //   fprintf(stderr,"(harm) Heavy lepton neutrinos:\n");
    //   fprintf(stderr,"(harm) R_beta_nux_anux    = %.15e\n",R_Beta);
    //   fprintf(stderr,"(harm) R_pair_nux_anux    = %.15e\n",R_ee);
    //   fprintf(stderr,"(harm) R_plasmon_nux_anux = %.15e\n",R_plasmon);
    // }

    return (R_Beta + R_ee + R_plasmon) / number_density; /* [1](B20) + */

    break;

  case NEUTRINO_ENERGY:
    energy_density = get_neutrino_energy_density(neutrino_flavor, neutrino_degeneracy, temperature);  /* [1](B19) + */

    Q_Beta = beta_process_emission_rate__neutrino_energy(neutrino_flavor, rho, temperature, electron_fraction, neutrino_degeneracy,
                                                         muhat_degeneracy, electron_degeneracy); /* [1](B14-B15) + */

    Q_ee = electron_positron_pair_annihilation_emission_rate__neutrino_energy(neutrino_flavor, temperature, electron_degeneracy,
                                                                              neutrino_degeneracy, anti_neutrino_degeneracy);  /* [1](B16) + */
    Q_plasmon = plasmon_decay_emission_rate__neutrino_energy(neutrino_flavor, temperature, electron_degeneracy,
                                                             neutrino_degeneracy, anti_neutrino_degeneracy);  /* Q_\gamma [1](B17) */

    if(!isfinite(Q_Beta)){
      Q_Beta=SMALL_tot;
    }
    if(!isfinite(Q_plasmon)){
      Q_plasmon=SMALL_tot;
    }
    if(!isfinite(Q_ee)){
      Q_ee=SMALL_tot;
    }

    if(energy_density<1e-20){
      energy_density=1e-20;
    }

    // if( neutrino_flavor == nu__e ) {
    //   fprintf(stderr,"(harm) Electron neutrino:\n");
    //   fprintf(stderr,"(harm) Q_beta_nue         = %.15e\n",Q_Beta);
    //   fprintf(stderr,"(harm) Q_pair_nue_anue    = %.15e\n",Q_ee);
    //   fprintf(stderr,"(harm) Q_plasmon_nue_anue = %.15e\n",Q_plasmon);
    // }
    // else if( neutrino_flavor == nu__e_bar ) {
    //   fprintf(stderr,"(harm) Electron antineutrino:\n");
    //   fprintf(stderr,"(harm) Q_beta_anue        = %.15e\n",Q_Beta);
    //   fprintf(stderr,"(harm) Q_pair_nue_anue    = %.15e\n",Q_ee);
    //   fprintf(stderr,"(harm) Q_plasmon_nue_anue = %.15e\n",Q_plasmon);
    // }
    // else {
    //   fprintf(stderr,"(harm) Heavy lepton neutrinos:\n");
    //   fprintf(stderr,"(harm) Q_beta_nux_anux    = %.15e\n",Q_Beta);
    //   fprintf(stderr,"(harm) Q_pair_nux_anux    = %.15e\n",Q_ee);
    //   fprintf(stderr,"(harm) Q_plasmon_nux_anux = %.15e\n",Q_plasmon);
    // }

    return (Q_Beta + Q_ee + Q_plasmon) / energy_density;  /* [1](B21) + */

    break;

  default:
    break;

  }

  /* HEREHERE -- should there be a "default value " in this switch ? */
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of type_of_transport  "); 
  fprintf(stdout,"\t type_of_transport   = %d \n", type_of_transport ); fflush(stdout); 
  // fail(FAIL_BASIC,0);
  exit(1);
  return(0.); 
    
}


/******************************************************************************************/
/******************************************************************************************
    effective_emission_rate():
    -------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- type_of_transport        - used to differentiate between whether it is desired to
                                      calculate the neutrino number or neutrino energy emission rate
        -- ds                       - current spacetime element of integration across the
                                      investigation region
        -- s_total                  - total "distance" traveled along a geodesic across
                                      spacetime up to the current region of investigation
        -- previous_optical_depth   - optical depth for the particular flavor of neutrino
                                      under observation up to the current location, and for the
                                      particular type of transport
        -- rho                      - in CGS units
        -- temperature              - in eV - from get_Temperature() function
        -- electron_fraction        - computed from the get_particle_fraction() function
        -- neutron_fraction         - computed from the get_particle_fraction() function
        -- proton_fraction          - computed from the get_particle_fraction() function
        -- neutrino_degeneracy      - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- anti_neutrino_degeneracy - computed from the get_degeneracy() function, utilizing the
                                      previously calculated value for optical depth
        -- neutron_degeneracy       - computed from the get_degeneracy() function
        -- proton_degeneracy        - computed from the get_degeneracy() function
        -- electron_degeneracy      - computed from the get_degeneracy() function

    Return Value:
        -- returns the effective number or energy of a particular flavor of neutrinos emitted
           from a particular location, as visible by some observer at the beginning of the
           integrated path.

*******************************************************************************************/
double effective_emission_rate(const NRPyEOS_params *restrict eos_params,
                               int neutrino_flavor, int type_of_transport, 
                               double rho, double energy_density, double temperature,
                               double electron_fraction, double muhat_degeneracy,double neutron_fraction, double proton_fraction,
                               double neutrino_degeneracy, double anti_neutrino_degeneracy,
                               double neutron_degeneracy, double proton_degeneracy, double electron_degeneracy, double Tau) {

  double number_density, neutrino_energy_density, t_diff, t_emiss__inv, R_value, Q_value;

  switch (type_of_transport) {
  case NEUTRINO_NUMBER:
    number_density = get_neutrino_number_density(neutrino_flavor, neutrino_degeneracy, temperature);  /* [1](B19) + */

    t_diff = get_diffusion_timescale(eos_params, neutrino_flavor, type_of_transport, rho, 
                                     energy_density, electron_fraction, temperature, Tau);   /* [1](A21) + */

    t_emiss__inv = get_inverse_emission_timescale(neutrino_flavor, type_of_transport, rho, temperature,
                                                  electron_fraction, neutrino_degeneracy, anti_neutrino_degeneracy, muhat_degeneracy,
                                                  neutron_degeneracy, proton_degeneracy, electron_degeneracy);  /* [1](B20) + */
	    
    R_value = t_emiss__inv * number_density;           /* [1](B20) + */

    //	    if (Tau>2.){
    //		fprintf(stderr, "tneutrino_degeneracy %e\n", neutrino_degeneracy);
    //		fprintf(stderr, "temp %e\n", temperature);
    //		fprintf(stderr, "t_diff %e\n", t_diff);
    //		fprintf(stderr, "t_emiss__inv %e\n", t_emiss__inv);
    //		fprintf(stderr, "number_density %e\n", number_density);		
    //		}

    return  R_value / (1.0 + t_diff * t_emiss__inv);   /* [1](B22) + */

    break;

  case NEUTRINO_ENERGY:
    neutrino_energy_density = get_neutrino_energy_density(neutrino_flavor, neutrino_degeneracy, temperature);  /* [1](B19) + */

    t_diff = get_diffusion_timescale(eos_params, neutrino_flavor, type_of_transport, rho, 
                                     energy_density, electron_fraction, temperature,Tau);    /* [1](A21) + */

    t_emiss__inv = get_inverse_emission_timescale(neutrino_flavor, type_of_transport, rho, temperature,
                                                  electron_fraction, neutrino_degeneracy, anti_neutrino_degeneracy, muhat_degeneracy,
                                                  neutron_degeneracy, proton_degeneracy, electron_degeneracy);  /* [1](B21) + */

    Q_value = t_emiss__inv * neutrino_energy_density;          /* [1](B21) + */

    return Q_value / (1.0 + t_diff * t_emiss__inv);   /* [1](B23) + */

    break;

  default:
    break;

  }

  /* HEREHERE -- should there be a "default value " in this switch ? */
  fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of type_of_transport "); 
  fprintf(stdout,"\t type_of_transport   = %d \n", type_of_transport ); fflush(stdout); 
  // fail(FAIL_BASIC,0);
  exit(1);
  return(0.); 

}


/******************************************************************************************/
/******************************************************************************************
    average_neutrino_energy_emitted():
    ---------------------------------
    Expected Parameters:
        -- neutrino_flavor          - the particular flavor of neutrino currently under
                                      observation
        -- ds                       - current spacetime element of integration across the
                                      investigation region
        -- s_total                  - total "distance" traveled along a geodesic across
                                      spacetime up to the current region of investigation
        -- previous_optical_depth__number - optical depth for the particular flavor of neutrino
                                      under observation up to the current location for
                                      specifically NEUTRINO_NUMBER
        -- previous_optical_depth__energy - optical depth for the particular flavor of neutrino
                                      under observation up to the current location for
                                      specifically NEUTRINO_ENERGY
        -- rho                      - in CGS units
        -- energy_density           - in CGS units
        -- temperature              - in eV - from get_Temperature() function
        -- electron_fraction        - computed from the get_particle_fraction() function
        -- neutron_fraction         - computed from the get_particle_fraction() function
        -- proton_fraction          - computed from the get_particle_fraction() function

    Return Value:
        -- returns the average expected energy of a particular flavor of neutrinos emitted
           from a particular location, as visible by some observer at the beginning of the
           integrated path, given the local environment.

*******************************************************************************************/
double average_neutrino_energy_emitted(const NRPyEOS_params *restrict eos_params,
                                       int neutrino_flavor, double rho, 
                                       double energy_density, double temperature,
                                       double electron_fraction, double neutron_fraction, double proton_fraction,
                                       double Tau_number, double Tau_energy) {
  int __attribute__((unused)) anti_neutrino_flavor;
  //double neutrino_degeneracy, anti_neutrino_degeneracy, proton_degeneracy, muhat_degeneracy;
  //double neutron_degeneracy, electron_degeneracy;
  double R_effective, Q_effective;

  switch (neutrino_flavor) {
  case nu__e:
    anti_neutrino_flavor = nu__e_bar;
    break;
  case nu__e_bar:
    anti_neutrino_flavor = nu__e;
    break;
  case nu__mu:
    anti_neutrino_flavor = nu__mu_bar;
    break;
  case nu__mu_bar:
    anti_neutrino_flavor = nu__mu;
    break;
  case nu__tau:
    anti_neutrino_flavor = nu__tau_bar;
    break;
  case nu__tau_bar:
    anti_neutrino_flavor = nu__tau;
    break;
  default:
    /* HEREHERE -- should there be a "default value " in this switch ? */
    fprintf(stdout,"(harm3d_neutrinos - ERROR) Inside average_neutrino_energy_emitted\n");
    fprintf(stdout,"(harm3d_neutrinos - ERROR) Invalid value of neutrino_flavor "); 
    fprintf(stdout,"\t neutrino_flavor  = %d \n", neutrino_flavor); fflush(stdout);
    exit(1);
    // fail(FAIL_BASIC,0);
    break;
  }


  double muhat_degeneracy,neutrino_degeneracy, anti_neutrino_degeneracy, neutron_degeneracy,proton_degeneracy,electron_degeneracy;

  muhat_degeneracy=0.;
  proton_degeneracy=0.;
  neutron_degeneracy=0.;
  electron_degeneracy=0.;

  double eta__nu_e=0.;
  double eta__nu_e_bar=0.;
  double eta__nu_mu=0.;
  double eta__nu_mu_bar=0.;
  double eta__nu_tau=0.;
  double eta__nu_tau_bar=0.;

   

  get_degeneracy_all(eos_params, rho, energy_density, temperature, electron_fraction, &electron_degeneracy, &proton_degeneracy, 
                     &neutron_degeneracy, &muhat_degeneracy, &eta__nu_e, &eta__nu_e_bar, &eta__nu_mu, &eta__nu_mu_bar, &eta__nu_tau, &eta__nu_tau_bar);
  switch (neutrino_flavor) {
  case nu__e:

    neutrino_degeneracy = eta__nu_e; 
    anti_neutrino_degeneracy = eta__nu_e_bar;
    break;
  case nu__e_bar:

    neutrino_degeneracy = eta__nu_e_bar; 
    anti_neutrino_degeneracy = eta__nu_e;
    break;
  case nu__mu:
            
    neutrino_degeneracy = eta__nu_mu; 
    anti_neutrino_degeneracy = eta__nu_mu_bar;
    break;
  case nu__mu_bar:

    neutrino_degeneracy = eta__nu_mu_bar; 
    anti_neutrino_degeneracy = eta__nu_mu;
    break;
  case nu__tau:

    neutrino_degeneracy = eta__nu_tau; 
    anti_neutrino_degeneracy = eta__nu_tau_bar;
    break;
  case nu__tau_bar:

    neutrino_degeneracy = eta__nu_tau_bar; 
    anti_neutrino_degeneracy = eta__nu_tau;
    break;
  }

  //neutrino_degeneracy = get_degeneracy(neutrino_flavor, rho, energy_density, temperature, electron_fraction);//get_degeneracy(neutrino_flavor, temperature, previous_optical_depth__number);
  //anti_neutrino_degeneracy = get_degeneracy(anti_neutrino_flavor, rho, energy_density, temperature, electron_fraction);//get_degeneracy(anti_neutrino_flavor, temperature, previous_optical_depth__number);
  //proton_degeneracy = get_degeneracy(PROTON, rho, energy_density, temperature, electron_fraction);//get_degeneracy(PROTON, temperature, previous_optical_depth__number);
  //neutron_degeneracy = get_degeneracy(NEUTRON, rho, energy_density, temperature, electron_fraction);//get_degeneracy(NEUTRON, temperature, previous_optical_depth__number);
  //electron_degeneracy = get_degeneracy(ELECTRON, rho, energy_density, temperature, electron_fraction);//get_degeneracy(ELECTRON, temperature, previous_optical_depth__number);
  //muhat_degeneracy=get_degeneracy(MUHAT, rho, energy_density, temperature, electron_fraction);

  /* [1](B22) */
  R_effective = effective_emission_rate(eos_params,
                                        neutrino_flavor, NEUTRINO_NUMBER, 
                                        rho, energy_density, temperature,
                                        electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction,
                                        neutrino_degeneracy, anti_neutrino_degeneracy,
                                        neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_number);

  /* [1](B23) */
  Q_effective = effective_emission_rate(eos_params,
                                        neutrino_flavor, NEUTRINO_ENERGY, 
                                        rho, energy_density, temperature,
                                        electron_fraction, muhat_degeneracy, neutron_fraction, proton_fraction,
                                        neutrino_degeneracy, anti_neutrino_degeneracy,
                                        neutron_degeneracy, proton_degeneracy, electron_degeneracy, Tau_energy);

  return Q_effective / R_effective;

}
