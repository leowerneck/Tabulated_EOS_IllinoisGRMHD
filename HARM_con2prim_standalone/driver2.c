
#include "driver.h"
#include "nuc_eos.h"

//#include "EOS_functions.c"

double temp_guess=5.;
double gam=4./3.;
double gam_m1_o_gam=0.25;//(gam-1.)/gam;
int treat_floor_as_failure = 0;
double ye=0.5;
double Kpoly=1.e2;
int eoskey=2;

//Note, also W and logPmagP are parameters


void calc_b2(const double * con, const double * prim,const double g_con[4][4], 
      const double g_cov[4][4], const double alpha, double * bsquared)
{
   
  double B1_cov, B2_cov, B3_cov;  
  const double lorentz = prim[WLORENTZ];
      
  // Lower indices - covariant
  B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  const double u1_cov=g_cov[1][1]*prim[u1_con]+g_cov[1][2]*prim[u2_con]+g_cov[1][3]*prim[u3_con];
  const double u2_cov=g_cov[2][1]*prim[u1_con]+g_cov[2][2]*prim[u2_con]+g_cov[2][3]*prim[u3_con];
  const double u3_cov=g_cov[3][1]*prim[u1_con]+g_cov[3][2]*prim[u2_con]+g_cov[3][3]*prim[u3_con];

  
  const double b0_con=1.0/alpha * (prim[B1_con]*u1_cov+prim[B2_con]*u2_cov+prim[B3_con]*u3_cov);
  
  // B^2 = B^i * B_i
  const double B_squared = con[B1_con]*B1_cov + con[B2_con]*B2_cov + con[B3_con]*B3_cov;
  *bsquared=(B_squared+pow((alpha*b0_con),2))/pow(lorentz,2);
 
}



void collect_metric(struct metric *h, const double hxx, const double hxy, const double hxz,
        const double hyy, const double hyz, const double hzz){

  // get the metric and inverse metric components
  
  for (int i=0;i<4;i++) {
    h->lo[0][i] = 0.0;
    h->lo[i][0] = 0.0;
    h->up[0][i] = 0.0;
    h->up[i][0] = 0.0;
  }
  
  h->lo[1][1] = hxx;
  h->lo[1][2] = hxy;
  h->lo[1][3] = hxz;
  h->lo[2][2] = hyy;
  h->lo[2][3] = hyz;
  h->lo[3][3] = hzz;
  h->lo[0][0] = -1.;
        
  h->lo[2][1] = h->lo[1][2];
  h->lo[3][1] = h->lo[1][3];
  h->lo[3][2] = h->lo[2][3];
  
  h->lo_det = hxx*hyy*hzz + 2.0*hxy*hyz*hxz - hxz*hxz*hyy - hxx*hyz*hyz - hxy*hxy*hzz;
  h->lo_sqrt_det = sqrt(h->lo_det);
  
  h->up_det = 1.0 / h->lo_det;
  h->alpha = 1.0 / h->lo_det;
          
  h->up[1][1] = (-hyz*hyz + hyy*hzz) * h->up_det;
  h->up[1][2] = (hxz*hyz - hxy*hzz) * h->up_det;
  h->up[2][2] = (-hxz*hxz + hxx*hzz) * h->up_det;
  h->up[1][3] = (-hxz*hyy + hxy*hyz) * h->up_det;
  h->up[2][3] = (hxy*hxz - hxx*hyz) * h->up_det;
  h->up[3][3] = (-hxy*hxy + hxx*hyy) * h->up_det;
  h->up[0][0] = -1.;
  
  h->up[2][1] = h->up[1][2];
  h->up[3][1] = h->up[1][3];
  h->up[3][2] = h->up[2][3];

  

}


void collect_metric_noble(struct of_geom *h, const double hxx, const double hxy, const double hxz,
        const double hyy, const double hyz, const double hzz){

  // get the metric and inverse metric components
  
  for (int i=0;i<4;i++) {
    h->gcov[0][i] = 0.0;
    h->gcov[i][0] = 0.0;
    h->gcon[0][i] = 0.0;
    h->gcon[i][0] = 0.0;
  }
  
  h->gcov[1][1] = hxx;
  h->gcov[1][2] = hxy;
  h->gcov[1][3] = hxz;
  h->gcov[2][2] = hyy;
  h->gcov[2][3] = hyz;
  h->gcov[3][3] = hzz;
  h->gcov[0][0] = -1.;
        
  h->gcov[2][1] = h->gcov[1][2];
  h->gcov[3][1] = h->gcov[1][3];
  h->gcov[3][2] = h->gcov[2][3];

  

  h->g = sqrt(hxx*hyy*hzz + 2.0*hxy*hyz*hxz - hxz*hxz*hyy - hxx*hyz*hyz - hxy*hxy*hzz);

  h-> g_inv = 1./h->g;


  double up_det = h-> g_inv;

  h->gcon[1][1] = (-hyz*hyz + hyy*hzz) * up_det;
  h->gcon[1][2] = (hxz*hyz - hxy*hzz) * up_det;
  h->gcon[2][2] = (-hxz*hxz + hxx*hzz) * up_det;
  h->gcon[1][3] = (-hxz*hyy + hxy*hyz) * up_det;
  h->gcon[2][3] = (hxy*hxz - hxx*hyz) * up_det;
  h->gcon[3][3] = (-hxy*hxy + hxx*hyy) *up_det;
  h->gcon[0][0] = -1.;
  
  h->gcon[2][1] = h->gcon[1][2];
  h->gcon[3][1] = h->gcon[1][3];
  h->gcon[3][2] = h->gcon[2][3];


  double lo_det = hxx*hyy*hzz + 2.0*hxy*hyz*hxz - hxz*hxz*hyy - hxx*hyz*hyz - hxy*hxy*hzz;
  
  h->alpha = 1.0 / lo_det;
  h->ncon[0] = 1./h->alpha;//sqrt(h->gcon[0][0]) ;
  //h->alpha = 1.0 / h->ncon[0];
  h->ncon[1] = h->alpha * h->gcon[0][1];
  h->ncon[2] = h->alpha * h->gcon[0][2];
  h->ncon[3] = h->alpha * h->gcon[0][3];
  h->beta[0] =  h->gcon[0][1]/h->gcon[0][0];
  h->beta[1] =  h->gcon[0][2]/h->gcon[0][0];
  h->beta[2] =  h->gcon[0][3]/h->gcon[0][0];
          
 //fprintf(stderr, "alpha%15.6E\n",h->alpha);
 //fprintf(stderr, "g_inv%15.6E\n",h-> g_inv);
}

void prim2con_MHD(const double * prim, double * U, const double g_con[4][4],
                  const double g_cov[4][4], double alpha, double *bsq ){
  

  
  const double lorentz = prim[WLORENTZ];
  
  U[RHO] = prim[RHO]*lorentz/alpha;

  U[BCON1] = prim[B1_con];
  U[BCON2] = prim[B2_con];
  U[BCON3] = prim[B3_con];
  
  const double uvec0_con=lorentz/alpha;
  const double uvec1_con=prim[u1_con]-alpha*lorentz*g_con[0][1];
  const double uvec2_con=prim[u2_con]-alpha*lorentz*g_con[0][2];
  const double uvec3_con=prim[u3_con]-alpha*lorentz*g_con[0][3];

  const double B1_cov = g_cov[1][1]*prim[B1_con]+g_cov[1][2]*prim[B2_con]+g_cov[1][3]*prim[B3_con];
  const double B2_cov = g_cov[2][1]*prim[B1_con]+g_cov[2][2]*prim[B2_con]+g_cov[2][3]*prim[B3_con];
  const double B3_cov = g_cov[3][1]*prim[B1_con]+g_cov[3][2]*prim[B2_con]+g_cov[3][3]*prim[B3_con];

  const double uvec0_cov=g_cov[0][0]*uvec0_con+g_cov[0][1]*uvec1_con+g_cov[0][2]*uvec2_con+g_cov[0][3]*uvec3_con;
  const double uvec1_cov=g_cov[1][0]*uvec0_con+g_cov[1][1]*uvec1_con+g_cov[1][2]*uvec2_con+g_cov[1][3]*uvec3_con;
  const double uvec2_cov=g_cov[2][0]*uvec0_con+g_cov[2][1]*uvec1_con+g_cov[2][2]*uvec2_con+g_cov[2][3]*uvec3_con;
  const double uvec3_cov=g_cov[3][0]*uvec0_con+g_cov[3][1]*uvec1_con+g_cov[3][2]*uvec2_con+g_cov[3][3]*uvec3_con;

  const double u0_cov=g_cov[0][1]*prim[u1_con]+g_cov[0][2]*prim[u2_con]+g_cov[0][3]*prim[u3_con];
  const double u1_cov=g_cov[1][1]*prim[u1_con]+g_cov[1][2]*prim[u2_con]+g_cov[1][3]*prim[u3_con];
  const double u2_cov=g_cov[2][1]*prim[u1_con]+g_cov[2][2]*prim[u2_con]+g_cov[2][3]*prim[u3_con];
  const double u3_cov=g_cov[3][1]*prim[u1_con]+g_cov[3][2]*prim[u2_con]+g_cov[3][3]*prim[u3_con];



  const double b0_con=1.0/alpha * (prim[B1_con]*u1_cov+prim[B2_con]*u2_cov+prim[B3_con]*u3_cov);

  //const double b0_cov=1./lorentz * alpha*u0_cov*b0_con;
  const double b1_cov=B1_cov/lorentz+1./lorentz * alpha*u1_cov*b0_con;
  const double b2_cov=B2_cov/lorentz+1./lorentz * alpha*u2_cov*b0_con;
  const double b3_cov=B3_cov/lorentz+1./lorentz * alpha*u3_cov*b0_con;

  const double b1_con=prim[B1_con]/lorentz+b0_con*(alpha*prim[u1_con]/lorentz-pow(alpha,2)*g_con[0][1]);
  const double b2_con=prim[B2_con]/lorentz+b0_con*(alpha*prim[u2_con]/lorentz-pow(alpha,2)*g_con[0][2]);
  const double b3_con=prim[B3_con]/lorentz+b0_con*(alpha*prim[u3_con]/lorentz-pow(alpha,2)*g_con[0][3]);

  const double Bsquared=prim[B1_con]*B1_cov + prim[B2_con]*B2_cov + prim[B3_con]*B3_cov;
  const double bsquared=(Bsquared+pow((alpha*b0_con),2))/pow(lorentz,2);
  
  const double w_enthalpy=prim[RHO]+prim[UU]+prim[PRESS];

  const double b0_cov=-b0_con;//for minkowski//1./b0_con * (b1_cov*b1_con+b2_cov*b2_con+b3_cov*b3_con);

  if (bsquared==0.){
  U[UU]=lorentz*(w_enthalpy)*uvec0_cov+alpha*(prim[PRESS])+prim[RHO]*lorentz/alpha;
  }

  else {
  U[UU]=lorentz*(w_enthalpy+bsquared)*uvec0_cov+alpha*(prim[PRESS]+bsquared/2.)+(-alpha*b0_con)*b0_cov+prim[RHO]*lorentz/alpha;
  }

  U[QCOV1]=lorentz*(w_enthalpy+bsquared)*uvec1_cov+(-alpha*b0_con)*b1_cov;
  U[QCOV2]=lorentz*(w_enthalpy+bsquared)*uvec2_cov+(-alpha*b0_con)*b2_cov;
  U[QCOV3]=lorentz*(w_enthalpy+bsquared)*uvec3_cov+(-alpha*b0_con)*b3_cov;

  U[YE] = prim[YE]*U[RHO];

  *bsq = bsquared;
}


int main(int argc, char *argv[]){
int npoints= pow(2,12);
int numprims = 13; // rho, vel123, eps, U, B123, ye, temp, press, lorentz, ent, abar
int numcons = 9;  //

const int ALIGNED = 1;
double tol_x = 1e-9;
double tol_prim = 10.*tol_x;
const double perturbation=0.05;
const double perturbation_W=0.05; // only for Lorentz factor
int max_iterations   = 300;

nuc_eos_C_ReadTable("eostable.h5");


//fprintf(stderr, "rho_min%15.6E\n", eos_rhomin/rho_gf);
//fprintf(stderr, "rho_max%15.6E\n", eos_rhomax/rho_gf);

//fprintf(stderr, "temp_min%15.6E\n", eos_tempmin/temp_gf);
//fprintf(stderr, "temp_max%15.6E\n", eos_tempmax/temp_gf);


double* prim = malloc(sizeof(double) * numprims);
double* prim_all = malloc(sizeof(double) * npoints * numprims);
double* original_prim = malloc(sizeof(double) * npoints * numprims);
double* con = malloc(sizeof(double) * numcons);
double* con_all= malloc(sizeof(double) * npoints * numcons);


for(int i=0;i<(npoints*numprims);i++){
    prim_all[i]=0.0;
  }
  for(int i=0;i<(numprims);i++){
    prim[i]=0.0;
  }
  for(int i=0;i<(npoints*numcons);i++){
    con_all[i]=0.0;
  }
  for(int i=0;i<(numcons);i++){
    con[i]=0.0;
  }

  int n1,n2;
  double v,vx,vy,vz;
  double rhoexp;
  double delrho;
  double exp1;
  double del1;
  double tempexp;
  double deltemp;
  double tempexp_i;
  double rhoexp_i;



rhoexp_i = log10(1.e3*eos_rhomin);
rhoexp = rhoexp_i;

delrho = (log10(0.01*eos_rhomax)-rhoexp_i)/(sqrt(1.0*npoints)-1.0);
//fprintf(stderr, "rho_min%15.6E\n", eos_rhomin);
//fprintf(stderr, "delrho %15.6E\n", delrho);


double logWminus1 = 0.0;
double W = 0.0;
double logPmagP = 0.0;
double rho_cgs = 0.0;
double temp_mev = 0.0;

W = 2;
logPmagP = -5.;


  // define the metric
  struct metric geom;
  collect_metric(&geom, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);

for(n1=0;n1<sqrt(1.0*npoints);n1++){
    
    double exp2;
    double del2;

      tempexp_i = log10(10.*eos_tempmin);
      tempexp = tempexp_i;
      deltemp = (log10(0.01*eos_tempmax)-tempexp_i)/(sqrt(1.0*npoints)-1.0);

    for(n2=0;n2<sqrt(1.0*npoints);n2++){
      
      // Translate 2D indexing to 1D index
      //fprintf(stderr,"rho: %15.6E eos_rhomax: %15.6E eos_rhomin: %15.6E deltarho: %15.6E\n",pow(10, rhoexp),eos_rhomax,eos_rhomin, delrho);
      //fprintf(stderr,"T: %15.6E eos_Tmax: %15.6E eos_Tmin: %15.6E deltT: %15.6E\n",pow(10, tempexp),eos_tempmax,eos_tempmin, deltemp);

      int i = n1*sqrt(1.0*npoints)+n2;
  
      rho_cgs  = pow(10, rhoexp);
      temp_mev = pow(10,tempexp);

      // Velocity magnitude
      v = sqrt(1.0-1.0/(W*W));
      // Randomly choose vx, vy, vz
      vx = v*((double)rand())/((double)RAND_MAX);
      vy = sqrt(v*v - vx*vx)*((double)rand())/((double)RAND_MAX);
      vz = sqrt(v*v - vx*vx - vy*vy);
      // Verify the magnitude
      if(fabs(v*v - (vx*vx+vy*vy+vz*vz))>1.0e-10){
        printf("v problem!!!\n");
        exit(1);
      }
      
      // Set primitive variables
      prim_all[i*numprims+RHO] = rho_cgs;
      prim_all[i*numprims+u1_con] = W*vx; //Note, this is utilde
      prim_all[i*numprims+u2_con] = W*vy;
      prim_all[i*numprims+u3_con] = W*vz;
      prim_all[i*numprims+TEMP] = temp_mev;
      prim_all[i*numprims+YE] = ye;
      prim_all[i*numprims+WLORENTZ] = W;
    
      for(int k=0;k<numprims;k++){
        prim[k] = prim_all[i*numprims+k];
      }
       double xeps,xprs,xrho,xtemp,xye;
       double pressure, epsilon;
      // compute eps, press
      if (eoskey==2){

        xrho = prim[RHO];
        xtemp = prim[TEMP];
        xeps = 0.0;
        xye = prim[YE];
        xprs = 0.0;

        EOS_press_with_T(xrho, &xeps, &xtemp, xye, &xprs);
      }

      else if (eoskey==1){
        xprs=Kpoly*pow(prim[RHO],gam);
        xeps=xprs/((gam-1.)*prim[RHO]);
        //fprintf(stderr, "rho %e\n", rho_cgs);
        //fprintf(stderr, "eps %e\n", xeps);
      }


      epsilon = xeps;
      prim_all[i*numprims+EPS] = epsilon; 
      prim_all[i*numprims+UU] = epsilon*prim[RHO]; 
      prim[EPS] = epsilon;
      prim[UU] = epsilon*prim[RHO];
      pressure = xprs;
      prim_all[i*numprims+PRESS] = pressure;
      prim[PRESS] = pressure;
      prim_all[i*numprims+EPS] = epsilon; 
      prim[EPS] = epsilon;
      prim_all[i*numprims+PRESS] = pressure;
      prim[PRESS] = pressure;
        
      // check P > 0
      if(pressure<0.0){
        printf("negative initial pressure: press = %e\n",pressure);
        printf("Terminating...\n");
        exit(1);
      }
      
      // set up magnetic field
      double B;
      B = sqrt(2.0*pow(10.0,logPmagP)*pressure);
      double Bxhat;
      double Byhat;
      double Bzhat;
      if(ALIGNED==1){

        double rand_temp = ((double)rand())/((double)RAND_MAX);
    
        if(rand_temp<0.5){
          Bxhat = vx/v;
          Byhat = vy/v;
          Bzhat = vz/v;
        }
        else{
          Bxhat = -1.0*vx/v;
          Byhat = -1.0*vy/v;
          Bzhat = -1.0*vz/v;
        }
      }
      else{
        Bxhat = ((double)rand())/((double)RAND_MAX);
        Byhat = sqrt(1.0 - Bxhat*Bxhat)*((double)rand())/((double)RAND_MAX);
        Bzhat = (-Bxhat*vx-Byhat*vy)/vz;
        double normB = sqrt(Bxhat*Bxhat + Byhat*Byhat + Bzhat*Bzhat);
        Bxhat = Bxhat/normB;
        Byhat = Byhat/normB;
        Bzhat = Bzhat/normB;
    
        // Verify orthogonality
        if(fabs((Bxhat*vx+Byhat*vy+Bzhat*vz))>1.0e-10){
          printf("Orthogonality problem!\n");
          printf("Terminating...\n");
          exit(1);
        }
      }
      prim_all[i*numprims+B1_con] = -Bxhat*B;
      prim_all[i*numprims+B2_con] = -Byhat*B;
      prim_all[i*numprims+B3_con] = -Bzhat*B;
      // Verify the magnitude
      if(fabs(1.0 -(Bxhat*Bxhat+Byhat*Byhat+Bzhat*Bzhat))>1.0e-10){
        printf("B problem!!!\n");
        printf("%e\n",Bxhat*Bxhat+Byhat*Byhat+Bzhat*Bzhat);
        printf("Terminating...\n");
        exit(1);
      }
      
      // Save original primitive variables
      // and set primitive variables for present prim2con conversion
      for(int k=0;k<numprims;k++){
        original_prim[i*numprims+k] = prim_all[i*numprims+k];
        prim[k] = prim_all[i*numprims+k];
       // printf("Save prims:\n");
        //printf("%e\t%e\n",prim_all[i*numprims+k],original_prim[i*numprims+k]);
      }
      
      // Set up conservatives by converting prim to con
      prim2con_MHD(prim, con, geom.up, geom.lo, geom.alpha);

      // update con_all
      for(int k=0;k<numcons;k++){
        con_all[i*numcons+k] = con[k];
      }

      tempexp = tempexp + deltemp;

     }

  rhoexp = rhoexp + delrho;

  }


  for(int i=0;i<npoints;i++){
    
    for(int k=0;k<numprims;k++){
      int randnum = ((int)rand() % 2); // 0 or 1
      randnum = randnum*2; // 0 or 2
      randnum = randnum - 1; // -1 or 1
      if (k == WLORENTZ){
        // perturb W-1
        prim_all[i*numprims+k] = 1.0 + (prim_all[i*numprims+k] - 1.0) * (1.0 + ((double)randnum)*perturbation_W);
      } else {
        if (k!=YE){
        // perturb all other primitives
        prim_all[i*numprims+k] *= (1.0 + ((double)randnum)*perturbation);}
      }
    }


      // checks on initial guesses
    // ensure they are physical
    
    // density

    if (prim_all[i*numprims+RHO] < eos_rhomin){
      fprintf(stderr, "min_condition\n");
      prim_all[i*numprims+RHO] = (1.0+perturbation)*eos_rhomin;
    }
    if (prim_all[i*numprims+RHO] > eos_rhomax){
      prim_all[i*numprims+RHO] = (1.0-perturbation)*eos_rhomax;
    }
    
    // velocities

    vx = prim_all[i*numprims+u1_con];
    vy = prim_all[i*numprims+u2_con];
    vz = prim_all[i*numprims+u3_con];
    
    // specific internal energy
    if (prim_all[i*numprims+EPS] < eos_epsmin){
      prim_all[i*numprims+EPS] = fabs((1.0+perturbation)*eos_epsmin - eos_epsmin) + eos_epsmin;
      prim_all[i*numprims+UU] = prim_all[i*numprims+EPS]*prim_all[i*numprims+RHO];
    }
    if (prim_all[i*numprims+EPS] > eos_epsmax){
      prim_all[i*numprims+EPS] = (1.0-perturbation)*eos_epsmax;
      prim_all[i*numprims+UU] = prim_all[i*numprims+EPS]*prim_all[i*numprims+RHO];
    }
    if (prim_all[i*numprims+TEMP] > eos_tempmax){
      prim_all[i*numprims+TEMP] = (1.0-perturbation)*eos_tempmax;
    }

      if (prim_all[i*numprims+TEMP] < eos_tempmin){
      prim_all[i*numprims+TEMP] = fabs((1.0+perturbation)*eos_tempmin - eos_tempmin) + eos_tempmin;
    }

  }


  




  double error[npoints]; 



  for(int i=0;i<npoints;i++){
//fprintf(stderr, "con2[0]%15.6E\n", con[0]);
    
    // set primitives (initial guesses) and conservatives
    for(int k=0;k<numprims;k++){
      prim[k] = prim_all[i*numprims+k];
    }
    for(int k=0;k<numcons;k++){
      con[k] = con_all[i*numcons+k];
    }

     bool excise = false;
     bool grace = false;
     struct report report;

    // call con2prim process
    //con2prim_MHD_(prim,con,g.up,g.lo,excise,c2p_grace,&report);

    struct metric g;
    collect_metric(&g, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);

    double bsq;
    calc_b2(con, prim, g.up, g.lo, g.alpha, &bsq);
   
    struct of_geom geom_noble;
    collect_metric_noble(&geom_noble, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
    //double inv_gdet = geom_noble->g_inv; 
    //fprintf(stderr, "g_inv%15.6E\n",geom_noble.g_inv);

 
    double gamma;
    gamma=W;
    //fprintf(stderr, "con[RHO]/prim[RHO]%e\n", con[RHO]/prim[RHO]);
    //fprintf(stderr, "rhomax%15.6E\n", eos_rhomax);
    //fprintf(stderr, "rhomin%15.6E\n", eos_rhomin);
    //fprintf(stderr, "tmax%15.6E\n", eos_tempmax);
    //fprintf(stderr, "tmin%15.6E\n", eos_tempmin);
    //fprintf(stderr, "outside den %e\n", prim[RHO]);
    //fprintf(stderr, "prim[RHO]%e\n", prim[RHO]);
    //fprintf(stderr, "prim_original[T]%e\n", original_prim[TEMP]);
    //fprintf(stderr, "prim_original[EPS]%e\n", original_prim[EPS]);
    //fprintf(stderr, "prim[T]%e\n", prim[TEMP]);
    //fprintf(stderr, "prim[EPS]%e\n", prim[EPS]);
   // fprintf(stderr, "prim[PRESS]%e\n", prim[PRESS]);
    // fprintf(stderr, "con[UU]%e\n", con[UU]);
     //fprintf(stderr, "con[RHO]%e\n", con[RHO]);
    //fprintf(stderr, "point\n");
    Utoprim_2d_eos_backup(con, &geom_noble, prim, &gamma, &bsq);
    //fprintf(stderr, "prim[UU]%e\n", prim[UU]);
    //fprintf(stderr, "prim_all[UU]%e\n", prim_all[i*numprims+UU]);
    
    //if (report.failed) {

    //  printf("FAILED: con2prim reconstruction failed for the following reason:\n");
    //  printf("Report: %s\n", report.err_msg);
    //} else {
    //  printf("SUCCESS: con2prim converged in %i iterations!\n",count[i]);
    //}
    
    // update primitives
    for(int k=0;k<numprims;k++){
      //fprintf(stderr, "prim_all[i*numprims+k]%e\n", prim_all[i*numprims+k]);
      prim_all[i*numprims+k] = prim[k];
    }
  }

  for(int i=0;i<npoints;i++){

    // check for discrepancy  

    error[i] = 0.0;
    for(int k=0;k<numprims;k++){
      //if ((k!=WLORENTZ)&&(k!=EPS)&&(k!=PRESS)&&(k!=YE)&&(k!=TEMP)){
      if ((k!=WLORENTZ)&&(k!=EPS)&&(k!=PRESS)&&(k!=YE)){
        error[i] += fabs((original_prim[i*numprims+k]-prim_all[i*numprims+k]) / original_prim[i*numprims+k]);
        //fprintf(stderr, "k%i\n", k);
        //fprintf(stderr, "prim_all[i*numprims+k]%e\n", prim_all[i*numprims+k]);
        //fprintf(stderr, "prim_original[i*numprims+k]%e\n", original_prim[i*numprims+k]);
        //fprintf(stderr, "error %e\n", fabs((original_prim[i*numprims+k]-prim_all[i*numprims+k]) / original_prim[i*numprims+k]));
      }
      //fprintf(stderr, "prim_all[i*numprims+k]%e\n", prim_all[i*numprims+k]);
      //fprintf(stderr, "prim_all[i*numprims+k]%e\n", prim_perturbed[i*numprims+k]);
      //fprintf(stderr, "prim_original[i*numprims+k]%e\n", original_prim[i*numprims+k]);
      //fprintf(stderr, "error %e\n", error_perturbed[i]);
      }
      //error[i] = 0.25 * error[i];

      //fprintf(stderr, "error[i]%e\n", error[i]);
    

    // check if eps<0 and v>c
    int epsneg = 0;
    if(prim_all[i*numprims+EPS]<0.0)
      epsneg = 1;

    
    // for each point, store its parameters (log(rho), log(temp)) or (log(W-1), log(Pmag/P))
    // a recovery characteristics in output file
    int n22 = (i%((int)sqrt(1.0*npoints)));
    int n11 = (i-n22)/sqrt(1.0*npoints);



      tempexp = tempexp_i+n22*deltemp;
      rhoexp = rhoexp_i + n11*delrho;

     //FILE *fp1;
     //char fname1[256];

     // sprintf(fname1, "output_original_logWminus1-%.2g_logPmagP-%.2g_YE-%.1f.dat",logWminus1, logPmagP, ye_test);

     // fp1 = fopen(fname1, "w");
     fprintf(stderr,"%e %e %e\n", rhoexp, tempexp, error[i]);

     //FILE *fp2;
     //char fname2[256];

     // sprintf(fname2, "output_perturbed_logWminus1-%.2g_logPmagP-%.2g_YE-%.1f.dat",logWminus1, logPmagP, ye_test);

     // fp2 = fopen(fname2, "w");
      //fprintf(fp2,"%e %e %e %e\n",logWminus1, rhoexp, tempexp, error_perturbed[i]);

  } 
}
