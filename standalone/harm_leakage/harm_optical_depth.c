/***********************************************************************************
 ***********************************************************************************/
#include "harm_units.h"
#include "harm_neutrinos.h"
#include "NRPyEOS.h"

/***************************************************************************/
/***************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 
     Optical depth for neutrinos as estimated in Nielsen 2014

*****************************************************************************/

void update_optical_depth_in_ijk(double ***optical_depth_new, double ***optical_depth_old, double ****opacity_nu, int neutrino_flavor, int type_of_transport, int i, int j, int k);
void optical_depth_func(double ****ph, int initialization);

void bounds_call_optical_depth( double ***optical_depth);



/******************************************************************************/
/******************************************************************************
 update_depth_func():
 ---------------------------
  -- Updates the optical depth via Neilsen 2014

******************************************************************************/

void update_optical_depth_in_ijk(double ***optical_depth_new, double ***optical_depth_old,double ****opacity_nu, int neutrino_flavor, int type_of_transport, int i, int j, int k){


  struct of_geom *geom_face1_ijk; //to get the spatial metric in coordinate 1
  struct of_geom *geom_face1_ip1jk;//to get the spatial metric in coordinate 1 p1 means plus 1

  struct of_geom *geom_face2_ijk; //to get the spatial metric in coordinate 2
  struct of_geom *geom_face2_ijp1k;//to get the spatial metric in coordinate 2 p1 means plus 1

  struct of_geom *geom_face3_ijk; //to get the spatial metric in coordinate 3
  struct of_geom *geom_face3_ijkp1;//to get the spatial metric in coordinate 3 p1 means plus 1


  double kappav_ijk, kappav_ip1jk, kappav_im1jk, kappav_ijp1k, kappav_ijm1k, kappav_ijkp1, kappav_ijkm1;

  double dx1_ijk_ip1jk,dx2_ijk_ip1jk,dx3_ijk_ip1jk,ds_ijk_ip1jk;
  double dx1_ijk_im1jk,dx2_ijk_im1jk,dx3_ijk_im1jk,ds_ijk_im1jk;
  double dx1_ijk_ijp1k,dx2_ijk_ijp1k,dx3_ijk_ijp1k,ds_ijk_ijp1k;
  double dx1_ijk_ijm1k,dx2_ijk_ijm1k,dx3_ijk_ijm1k,ds_ijk_ijm1k;
  double dx1_ijk_ijkp1,dx2_ijk_ijkp1,dx3_ijk_ijkp1,ds_ijk_ijkp1;
  double dx1_ijk_ijkm1,dx2_ijk_ijkm1,dx3_ijk_ijkm1,ds_ijk_ijkm1;

  double tau_tot, tau_ijk_ip1jk,tau_ijk_im1jk,tau_ijk_ijp1k,tau_ijk_ijm1k,tau_ijk_ijkp1,tau_ijk_ijkm1;

  //ALL_LOOP {

  //if((i!=0)&&(j!=0)&&(k!=0)&&(i!=N1TOT-1)&&(j!=N2TOT-1)&&(k!=N3TOT-1)){

  get_geometry(i,j,k,FACE1,ncurr,geom_face1_ijk);
  get_geometry(i,j,k,FACE2,ncurr,geom_face2_ijk);
  get_geometry(i,j,k,FACE3,ncurr,geom_face3_ijk);
    
  get_geometry(i+1,j,k,FACE1,ncurr,geom_face1_ip1jk);
  get_geometry(i,j+1,k,FACE2,ncurr,geom_face2_ijp1k);
  get_geometry(i,j,k+1,FACE3,ncurr,geom_face3_ijkp1);



  dx1_ijk_ip1jk=dx[1];
  //dx2_ijk_ip1jk=0.;
  //dx3_ijk_ip1jk=0.;

  //this is how ds is calculated, but with dx2_ijk_ip1jk=dx3_ijk_ip1jk=0
  //ds_ijk_ip1jk=sqrt(geom_face1_ip1jk->gcov[1][1]*dx1_ijk_ip1jk*dx1_ijk_ip1jk+2.*geom_face1_ip1jk->gcov[1][2]*dx1_ijk_ip1jk*dx2_ijk_ip1jk+ 
  //	2.*geom_face1_ip1jk->gcov[1][3]*dx1_ijk_ip1jk*dx3_ijk_ip1jk+geom_face1_ip1jk->gcov[2][2]*dx2_ijk_ip1jk*dx2_ijk_ip1jk+
  //	2.*geom_face1_ip1jk->gcov[2][3]*dx2_ijk_ip1jk*dx3_ijk_ip1jk+geom_face1_ip1jk->gcov[3][3]*dx3_ijk_ip1jk*dx3_ijk_ip1jk);

  ds_ijk_ip1jk=sqrt(geom_face1_ip1jk->gcov[1][1]*dx1_ijk_ip1jk*dx1_ijk_ip1jk);


  dx1_ijk_im1jk=dx[1];
  //dx2_ijk_im1jk=0.;
  //dx3_ijk_im1jk=0.;
  //ds_ijk_im1jk=sqrt(geom_face1_ijk->gcov[1][1]*dx1_ijk_im1jk*dx1_ijk_im1jk+2.*geom_face1_ijk->gcov[1][2]*dx1_ijk_im1jk*dx2_ijk_im1jk+ 
  //	2.*geom_face1_ijk->gcov[1][3]*dx1_ijk_im1jk*dx3_ijk_im1jk+geom_face1_ijk->gcov[2][2]*dx2_ijk_im1jk*dx2_ijk_im1jk+
  //	2.*geom_face1_ijk->gcov[2][3]*dx2_ijk_im1jk*dx3_ijk_im1jk+geom_face1_ijk->gcov[3][3]*dx3_ijk_im1jk*dx3_ijk_im1jk);


  ds_ijk_im1jk=sqrt(geom_face1_ijk->gcov[1][1]*dx1_ijk_im1jk*dx1_ijk_im1jk);

  //dx1_ijk_ijp1k=0.;
  dx2_ijk_ijp1k=dx[2];
  //dx3_ijk_ijp1k=0.;
  //ds_ijk_ijp1k=sqrt(geom_face2_ijp1k->gcov[1][1]*dx1_ijk_ijp1k*dx1_ijk_ijp1k+2.*geom_face2_ijp1k->gcov[1][2]*dx1_ijk_ijp1k*dx2_ijk_ijp1k+ 
  //	2.*geom_face2_ijp1k->gcov[1][3]*dx1_ijk_ijp1k*dx3_ijk_ijp1k+geom_face2_ijp1k->gcov[2][2]*dx2_ijk_ijp1k*dx2_ijk_ijp1k+
  //	2.*geom_face2_ijp1k->gcov[2][3]*dx2_ijk_ijp1k*dx3_ijk_ijp1k+geom_face2_ijp1k->gcov[3][3]*dx3_ijk_ijp1k*dx3_ijk_ijp1k);

  ds_ijk_ijp1k=sqrt(geom_face2_ijp1k->gcov[2][2]*dx2_ijk_ijp1k*dx2_ijk_ijp1k);

  //dx1_ijk_ijm1k=0.;
  dx2_ijk_ijm1k=dx[2];
  //dx3_ijk_ijm1k=0.;
  //ds_ijk_ijm1k=sqrt(geom_face2_ijk->gcov[1][1]*dx1_ijk_ijm1k*dx1_ijk_ijm1k+2.*geom_face2_ijk->gcov[1][2]*dx1_ijk_ijm1k*dx2_ijk_ijm1k+ 
  //	2.*geom_face2_ijk->gcov[1][3]*dx1_ijk_ijm1k*dx3_ijk_ijm1k+geom_face2_ijk->gcov[2][2]*dx2_ijk_ijm1k*dx2_ijk_ijm1k+
  //	2.*geom_face2_ijk->gcov[2][3]*dx2_ijk_ijm1k*dx3_ijk_ijm1k+geom_face2_ijk->gcov[3][3]*dx3_ijk_ijm1k*dx3_ijk_ijm1k);
      
  ds_ijk_ijm1k=sqrt(geom_face2_ijk->gcov[2][2]*dx2_ijk_ijm1k*dx2_ijk_ijm1k);

  //dx1_ijk_ijkp1=0.;
  //dx2_ijk_ijkp1=0.;
  dx3_ijk_ijkp1=dx[3];
  //ds_ijk_ijkp1=sqrt(geom_face3_ijkp1->gcov[1][1]*dx1_ijk_ijkp1*dx1_ijk_ijkp1+2.*geom_face3_ijkp1->gcov[1][2]*dx1_ijk_ijkp1*dx2_ijk_ijkp1+ 
  //  2.*geom_face3_ijkp1->gcov[1][3]*dx1_ijk_ijkp1*dx3_ijk_ijkp1+geom_face3_ijkp1->gcov[2][2]*dx2_ijk_ijkp1*dx2_ijk_ijkp1+
  //  2.*geom_face3_ijkp1->gcov[2][3]*dx2_ijk_ijkp1*dx3_ijk_ijkp1+geom_face3_ijkp1->gcov[3][3]*dx3_ijk_ijkp1*dx3_ijk_ijkp1);

  ds_ijk_ijkp1=sqrt(geom_face3_ijkp1->gcov[3][3]*dx3_ijk_ijkp1*dx3_ijk_ijkp1);


  //dx1_ijk_ijkm1=0.;
  //dx2_ijk_ijkm1=0.;
  dx3_ijk_ijkm1=dx[3];
  //ds_ijk_ijkm1=sqrt(geom_face3_ijk->gcov[1][1]*dx1_ijk_ijkm1*dx1_ijk_ijkm1+2.*geom_face3_ijk->gcov[1][2]*dx1_ijk_ijkm1*dx2_ijk_ijkm1+ 
  //  2.*geom_face3_ijk->gcov[1][3]*dx1_ijk_ijkm1*dx3_ijk_ijkm1+geom_face3_ijk->gcov[2][2]*dx2_ijk_ijkm1*dx2_ijk_ijkm1+
  //  2.*geom_face3_ijk->gcov[2][3]*dx2_ijk_ijkm1*dx3_ijk_ijkm1+geom_face3_ijk->gcov[3][3]*dx3_ijk_ijkm1*dx3_ijk_ijkm1);

  ds_ijk_ijkm1=sqrt(geom_face3_ijk->gcov[3][3]*dx3_ijk_ijkm1*dx3_ijk_ijkm1);


  int opt_depth_helper;

  if (type_of_transport==NEUTRINO_NUMBER){
    opt_depth_helper=neutrino_flavor;
  }
  else if (type_of_transport==NEUTRINO_ENERGY){
    opt_depth_helper=neutrino_flavor+N_OPTICAL_DEPTHS/2;
  }


  kappav_ijk=opacity_nu[i][j][k][opt_depth_helper];

  kappav_ip1jk=opacity_nu[i+1][j][k][opt_depth_helper];
  kappav_im1jk=opacity_nu[i-1][j][k][opt_depth_helper];

  kappav_ijp1k=opacity_nu[i][j+1][k][opt_depth_helper];
  kappav_ijm1k=opacity_nu[i][j-1][k][opt_depth_helper];

  kappav_ijkp1=opacity_nu[i][j][k+1][opt_depth_helper];
  kappav_ijkm1=opacity_nu[i][j][k-1][opt_depth_helper];


      
  double kappa_ds_ijk_ip1jk,kappa_ds_ijk_im1jk, kappa_ds_ijk_ijp1k,kappa_ds_ijk_ijm1k;
  double kappa_ds_ijk_ijkp1, kappa_ds_ijk_ijkm1;

  kappa_ds_ijk_ip1jk=ds_ijk_ip1jk*INVLENGTH*(kappav_ijk+kappav_ip1jk)*0.5;
  kappa_ds_ijk_im1jk=ds_ijk_im1jk*INVLENGTH*(kappav_ijk+kappav_im1jk)*0.5;
  kappa_ds_ijk_ijp1k=ds_ijk_ijp1k*INVLENGTH*(kappav_ijk+kappav_ijp1k)*0.5;
  kappa_ds_ijk_ijm1k=ds_ijk_ijm1k*INVLENGTH*(kappav_ijk+kappav_ijm1k)*0.5;
  kappa_ds_ijk_ijkp1=ds_ijk_ijkp1*INVLENGTH*(kappav_ijk+kappav_ijkp1)*0.5;
  kappa_ds_ijk_ijkm1=ds_ijk_ijkm1*INVLENGTH*(kappav_ijk+kappav_ijkm1)*0.5;

  tau_ijk_ip1jk=kappa_ds_ijk_ip1jk+optical_depth_old[i+1][j][k];
  tau_ijk_im1jk=kappa_ds_ijk_im1jk+optical_depth_old[i-1][j][k];
  tau_ijk_ijp1k=kappa_ds_ijk_ijp1k+optical_depth_old[i][j+1][k];
  tau_ijk_ijm1k=kappa_ds_ijk_ijm1k+optical_depth_old[i][j-1][k];
  tau_ijk_ijkp1=kappa_ds_ijk_ijkp1+optical_depth_old[i][j][k+1];
  tau_ijk_ijkm1=kappa_ds_ijk_ijkm1+optical_depth_old[i][j][k-1];




  tau_tot = 1e100;
  tau_tot = MIN(tau_tot, tau_ijk_ijp1k);
  tau_tot = MIN(tau_tot, tau_ijk_ijm1k);
  tau_tot = MIN(tau_tot, tau_ijk_ip1jk);
  tau_tot = MIN(tau_tot, tau_ijk_im1jk);



  if(N3!=1){ 
    tau_tot = MIN(tau_tot, tau_ijk_ijkm1);
    tau_tot = MIN(tau_tot, tau_ijk_ijkp1);
  }


  optical_depth_new[i][j][k]=tau_tot;

  return;

}

void optical_depth_func(double ****ph, int initialization){

  int i_optdepth;
  double ***optical_depth1;
  double ***optical_depth2;
  int i,j,k; 
 
  ALLOC_3D_ARRAY(optical_depth1,N1TOT,N2TOT,N3TOT);
  ALLOC_3D_ARRAY(optical_depth2,N1TOT,N2TOT,N3TOT);
  


  /* Determine the number of iterations and whether we are starting the optical depth from scratch or not: */ 

  int n_iter_optdepth = 500; /* Default number of iterations : */

  if( initialization ){
    int zero_depth_anywhere = 0; 
    ALL_LOOP {
      if( optical_depth[i][j][k][0] < 0. ) {
	zero_depth_anywhere = 1;
	break;
      }
    }
    if( zero_depth_anywhere ) {
      /* clear the optical depth and start from scratch, making sure to iterate enough : */
      n_iter_optdepth = MAX(N1_glob,MAX(N3_glob,N2_glob))*20.; /* need to iterate more if we are initializing it */
      ALL_LOOP   for (i_optdepth=0; i_optdepth<N_OPTICAL_DEPTHS;i_optdepth++) {
	optical_depth[i][j][k][i_optdepth] = 0.;
      }
    }
  }
      

  for (i_optdepth=0; i_optdepth<N_OPTICAL_DEPTHS;i_optdepth++){
    int type_of_transport;
    int neutrino_flavor;
    //    int i,j,k;
    
    
    if (i_optdepth<N_OPTICAL_DEPTHS/2){
      type_of_transport=0; //neutrino number
      neutrino_flavor=i_optdepth;
    }
    else{
      type_of_transport=1; //neutrino energy
      neutrino_flavor=i_optdepth-N_OPTICAL_DEPTHS/2;
    }

    ALL_LOOP {
      //for(int i=0;i<N1TOT-1;i++) for(int j=0;j<N2TOT;j++) for(int k=0;k<N3TOT;k++) {
      optical_depth1[i][j][k] = optical_depth[i][j][k][i_optdepth];
      optical_depth2[i][j][k] = optical_depth[i][j][k][i_optdepth];
      
      opacity_nu[i][j][k][i_optdepth]=get_total_transport_opacity(neutrino_flavor, type_of_transport,
								  ph[i][j][k][RHO]*INVRHOGF, (ph[i][j][k][UU]/ph[i][j][k][RHO])*INVEPSGF,
								  ph[i][j][k][YE], ph[i][j][k][TEMP]*1e6);
    }
    
    if( initialization ){
      for(i=1;i<N1TOT-1;i++) for(j=1;j<N2TOT-1;j++) for( k=1;k<N3TOT-1;k++) {
	    update_optical_depth_in_ijk(optical_depth1, optical_depth2, opacity_nu, neutrino_flavor, type_of_transport, i, j, k);
	  }
      bounds_call_optical_depth(optical_depth1);
    }


    int step_int;
    int step_int_new;
    double change_helper_new;
    double change_rate_helper;
    double change_helper_old;
    
    step_int = 0;
    step_int_new = 1;
    change_rate_helper = 1.;
    change_helper_old = 1.;
 
    double change_helper_old_eps;
    double change_rate_helper_eps;


    if(initialization){
      change_helper_old_eps=0.;
      change_rate_helper_eps=0.;
    }
    else{
      change_helper_old_eps=1.e-3;
      change_rate_helper_eps=1e-4;
    }

    while ((change_helper_old>change_helper_old_eps)){//((fabs(sum_helper1-sum_helper2)>1e-15*sum_helper1)){ ||(change_helper_old>1e-4)
      if(change_rate_helper<change_rate_helper_eps){
	break;
      }
      
      change_helper_new = 0.;
      
      step_int_new = step_int % 2;
      if (step_int_new==0){
	for(i=1;i<N1TOT-1;i++) for( j=1;j<N2TOT-1;j++) for(k=1;k<N3TOT-1;k++) {
              update_optical_depth_in_ijk(optical_depth2, optical_depth1, opacity_nu, neutrino_flavor, type_of_transport, i, j, k);

              if((optical_depth1[i][j][k] > 1e-6)&(optical_depth2[i][j][k] > 1e-6)){

                change_helper_new=MAX(change_helper_new, fabs(optical_depth1[i][j][k]-optical_depth2[i][j][k])/optical_depth1[i][j][k]);
              }
            }
	bounds_call_optical_depth(optical_depth2);
      }
      else {
	for( i=1;i<N1TOT-1;i++) for(j=1;j<N2TOT-1;j++) for(k=1;k<N3TOT-1;k++) {
              update_optical_depth_in_ijk(optical_depth1, optical_depth2, opacity_nu, neutrino_flavor, type_of_transport, i, j, k);
      
              if((optical_depth1[i][j][k] > 1e-6)&(optical_depth2[i][j][k] > 1e-6)) {
                change_helper_new=MAX(change_helper_new, fabs(optical_depth1[i][j][k]-optical_depth2[i][j][k])/optical_depth2[i][j][k]);
              }

            }
	bounds_call_optical_depth(optical_depth1);
      }
      if (step_int > n_iter_optdepth){
	break;
      }
      step_int++;
      
      mpi_global_max(&change_helper_new); /* Need to make sure that we find the maximum across all processors */
      
      if (step_int > 1){
	change_rate_helper=fabs(change_helper_old - change_helper_new)/change_helper_old;
      }

      change_helper_old=change_helper_new;


    } /* End of  "while ((change_helper_old>1e-3))"  loop */

    
      /* Copy latest value of optical depth to the global array */ 
    if (step_int_new==0){
      ALL_LOOP{
	optical_depth[i][j][k][i_optdepth] = optical_depth2[i][j][k];
      }
    }
    else if (step_int_new==1){
      ALL_LOOP{
	optical_depth[i][j][k][i_optdepth] = optical_depth1[i][j][k];
      }
    }
 
  }
    
  DEALLOC_3D_ARRAY(optical_depth1,N1TOT,N2TOT,N3TOT);
  DEALLOC_3D_ARRAY(optical_depth2,N1TOT,N2TOT,N3TOT);
  
  return;
  
}




void bounds_call_optical_depth( double ***optical_depth)
{
  unsigned int i,j,k,l,g,j2,k2;

#if( BC_TYPE_CHOICE  == BC_SPHERICAL_OUTFLOW5 )
  /* 1-lower face */
  if( bc_pid[1][BCDN] == BC_PHYS ) {
    GLOOP     N2_LOOP   N3_LOOP  {
      optical_depth[N1S-1-g][j][k] = optical_depth[N1S][j][k] ;
    }
  }

  /* 1-upper face */
  if( bc_pid[1][BCUP] == BC_PHYS ) {
    GLOOP     N2_LOOP   N3_LOOP  {
      optical_depth[N1E+1+g][j][k] = optical_depth[N1E][j][k] ;
    }
  }

  /* 2-lower face */
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    N1ALL_LOOP      GLOOP    N3_LOOP   {
      optical_depth[i][N2S-1-g][k] = optical_depth[i][N2S+g][k] ;
    }
  }

  /* 2-upper face */
  if( bc_pid[2][BCUP] == BC_PHYS ) {
    N1ALL_LOOP      GLOOP    N3_LOOP   {
      optical_depth[i][N2E+1+g][k] = optical_depth[i][N2E-g][k] ;
    }
  }

  /* 3-lower face */
  if( bc_pid[3][BCDN] == BC_PHYS ) {
    N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
      optical_depth[i][j][N3S-1-g] = optical_depth[i][j][N3E-g] ;
    }
  }

  /* 3-upper face */
  if( bc_pid[3][BCUP] == BC_PHYS ) {
    N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
      optical_depth[i][j][N3E+1+g] = optical_depth[i][j][N3S+g] ;
    }
  }

#endif

#if( BC_TYPE_CHOICE  == BC_SPHERICAL_FULLSPHERE4 )

  /* Note that if the subdomain covers an entire dimension, then its bc_pid is set to BC_PHYS 
     and it must perform the symmetry operation using its own data of course;   */

  /* 1-lower face :   */
  if( bc_pid[1][BCDN] == BC_PHYS ) {
    GLOOP     N2_LOOP  N3_LOOP  {
      optical_depth[N1S-1-g][j][k] = optical_depth[N1S][j][k] ;
    }
  }

  /* 1-upper face */
  if( bc_pid[1][BCUP] == BC_PHYS ) {
    GLOOP     N2_LOOP   N3_LOOP  {
      optical_depth[N1E+1+g][j][k] = 0.;//optical_depth[N1E][j][k] ;
    }
  }

  /* 2-lower face */
  if( bc_pid[2][BCDN] == BC_PHYS ) {
    /* This BC only works for an even number of azimuthal zones and domains: */
    if( ((totalsize[PH] > 1) && ((totalsize[PH] % 2) != 0)) ||
        ((    ncpux[PH] > 1) && ((    ncpux[PH] % 2) != 0))  ) {
      ERROR_MESSAGE("totalsize[PH] and ncpux[PH] must be even!! "); fail(FAIL_BASIC,0) ; 
    }

    N1ALL_LOOP      GLOOP    N3_LOOP   {
      k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
      optical_depth[i][N2S-1-g][k] = optical_depth[i][N2S+g][k2] ;
    }
  }

  /* 2-upper face */
  if( bc_pid[2][BCUP] == BC_PHYS ) {
    /* This BC only works for an even number of azimuthal zones and domains: */
    if( ((totalsize[PH] > 1) && ((totalsize[PH] % 2) != 0)) ||
        ((    ncpux[PH] > 1) && ((    ncpux[PH] % 2) != 0))  ) {
      ERROR_MESSAGE("totalsize[PH] and ncpux[PH] must be even!! "); fail(FAIL_BASIC,0) ; 
    }

    N1ALL_LOOP      GLOOP    N3_LOOP   {
      k2 = ( k - NG + totalsize[3]/2 ) % totalsize[3] + NG; 
      optical_depth[i][N2E+1+g][k] = optical_depth[i][N2E-g][k2] ;
    }
  }

  /* 3-lower face */
  if( bc_pid[3][BCDN] == BC_PHYS ) {
    N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
      optical_depth[i][j][N3S-1-g] = optical_depth[i][j][N3E-g] ;
    }
  }

  /* 3-upper face */
  if( bc_pid[3][BCUP] == BC_PHYS ) {
    N1ALL_LOOP      N2ALL_LOOP     GLOOP   {
      optical_depth[i][j][N3E+1+g] = optical_depth[i][j][N3S+g] ;
    }
  }

#endif

  return; 
} 
