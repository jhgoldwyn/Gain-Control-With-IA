#include <stdio.h>
#include <math.h>  /*per sin i cos*/
#include <stdlib.h> /*pel malloc i calloc*/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "memory.h"

// what data to save
int writeData; /* command line input */

// simulation run time
double tStop; /* command line input */
 
// cable
double d; /* command line input */ //10; diffusion
int synCable; /* command line input */ // 5; input position

//capacitance
double c=1.;

// leak
double glk=1;
double glkCable=.1;
double vlk=-70;

// Na channel
double gna=37;
double vna=55;

// Na activation [m]
double thetam=30;
double sm=15;

// Na inactivation [h]
double thetah=-39;
double sh=3.1;
double tauh0=1;
double tauh1=500;
double thh=57;
double sigmah=3;

// K channel
double gk=45;
double vk=-80;

// K activation
double thetan=-32.;
double sn=-8.;
double taun0=1;
double taun1=100.;
double thn=80.;
double sigman=26.;

// A channel
double ga; /* command line input */

// A activation [a]
double taua   = 2;
double thetaa = -50; 
double sigmaa = 20;

// A inactivation [b]
double taub=150;
double thetab=-70;
double sigmab=-6;

// kinetics rescaling
double phi=0.75;

// synapse 
double vsyne = 0.;
double vsyni = -85.;
double gsyne, gsyni; /* command line input */
double rsyne, rsyni; /* command line input */
double betae = .2;
double betai = .18;
double te, ti, teNext,tiNext; /* synapse event times */

double v,m,b,n,a,se,si;
double *vCable;
double minf,hinf,ninf,binf,ainf;
double tauh,taun,taua;
double ina,ik,ilk,ia,isyne,isyni;
double alphae;
 
//double twopi=8*atan(1.);
double twopi=3.1415926535;

int i;

int func (double t, const double x[], double dx[], void *params);

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);


int main(int argc, char *argv[]){

    double *x,*dx;
    double vlast,tspike;
    double si0;
    double rand;
    int ndim=15; 
    FILE *output;
    double tol=1.0e-5;
    double mu;  
 
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_rk8pd;  // Runge-Kutta Prince-Dormand 8-9.
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_rk2imp; // Runge-Kutta 2 implicit.
    const gsl_odeiv_step_type * T  = gsl_odeiv_step_rk4imp; // Runge-Kutta 4 implicit.
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_bsimp; // Burlirsch-Stoer implicit.. Need jacobian!
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_gear1; // Gear 1 implicit.
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_gear2; // Gear 2 implicit.

    gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, ndim);
    gsl_odeiv_control * c  = gsl_odeiv_control_y_new (tol, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (ndim);

    gsl_odeiv_system sys = {func, jac, ndim, &mu};
    
    /* Random generator stuff */
    const gsl_rng_type *U;
    gsl_rng *r;
    gsl_rng_env_setup();
    U = gsl_rng_default;
    r = gsl_rng_alloc (U);
    

    /* Parameters From Command Line */ 
    writeData = atof(argv[1]);   
    tStop     = atof(argv[2]);   
    ga        = atof(argv[3]);   
    gsyne     = atof(argv[4]); 
    gsyni     = atof(argv[5]); 
    rsyne     = atof(argv[6]); 
    rsyni     = atof(argv[7]); 
    d         = atof(argv[8]);
    synCable  = atoi(argv[9]); 

    double t, tf;
    double h = 1.e-3;
    int status;
   
    char name1[100];
    int number1,number2,number3;
    int ncount=0;  
    int nite;
    double sum;
    if (writeData==1){  /* save all data */
      sprintf(name1, "data10cpt_ga%g_ge%.1f_gi%g_re%.1f_ri%g_d%g_syn%d.txt",ga,gsyne, gsyni, rsyne,rsyni,d,synCable);}
    else if (writeData==2){ /* save spikes */
      sprintf(name1, "spike10cpt_ga%g_ge%.1f_gi%g_re%.1f_ri%g_d%g_syn%d.txt",ga,gsyne, gsyni, rsyne,rsyni,d,synCable);}
    output=fopen(name1, "w");

   /*open the data file*/
    output=fopen(name1, "w");
    if (output==NULL){
    printf("Error in opening data file\n");
    exit(1);
    }

/*memory allocation*/
x  = new_vector(0, ndim-1);
dx = new_vector(0, ndim-1); 
vCable = new_vector(0,8);

/*Initial point*/
double v0=-70;
x[0]=v0; /* v */
x[1]=1./(1.+exp((v0-thetan)/sn));     /* n */
x[2]=1/(1+exp(-(v0-thetaa)/sigmaa));  /* a */
x[3]=1/(1+exp(-(v0-thetab)/sigmab));  /* b */
x[4]=0; /* se */
x[5]=0; /* si */

for (i=6;i<15;i++){x[i] = v0;}   /* cable */

/*integrate*/

t=0.;
nite=0.;
vlast=v0;
te = 0;
ti = -999.; tiNext = 1000./rsyni; si0 = 0;
te = -999.; rand =gsl_ran_flat(r, 0,1); teNext = -1000.*log(1-rand)/(rsyne);
if (teNext<tiNext){tf=teNext;}
else {tf=tiNext;}

while (t<tStop){

  while (t<tf){

   status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tf, &h, x);
	
   if (status != GSL_SUCCESS)
    break;
  
  if (writeData==1){ /* save everything */
  fprintf(output," %.5f ",t);
  for (i=0;i<ndim;i++){ 
   fprintf(output," %.5f ",x[i]);
   }  
  fprintf(output, "\n");
  }

  /* spike detection */
  tspike =0.;
  if ((x[0]>=-10)&&(vlast<-10)){ 
     ncount++;	
     tspike=t;
     if (writeData==2){ /* save spikes */
        fprintf(output," %.5f %.5f",tspike, si0);  fprintf(output, "\n");
     }
  }

  vlast=x[0];

  nite++;
  }

  // generate next excitatory event
  if   (tf>=teNext){   
   te = teNext;
   rand=gsl_ran_flat(r, 0,1); teNext+=-1000.*log(1-rand)/(rsyne); 

   // for instantaneous rise se
   x[4]=1.; // ceiling at 1 

   si0 =  x[5]; // si value at the excitatory event onset
   if (writeData==2){fprintf(output," %.5f %.5f ",-tf, x[5]);  fprintf(output, "\n");} /* save excitatory events */

  }

  // generate next inhibitory event
  if   (tf>=tiNext){tiNext+=1000./rsyni; x[5]=1.;} // ceiling at 1 
  if (teNext<tiNext){tf=teNext;}
  else {tf=tiNext;}

  }

  if (writeData==2){fprintf(output," %.5f  %.5f ",tStop , tStop);  fprintf(output, "\n");}
  free_vector(x,0);
  free_vector(dx,0);

  fclose(output);

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  
}

/* Vector field we want to integrate*/
 
int func (double t, const double x[], double dx[], void *params)
{
    double mu = *(double *)params;
    
   double m3,n4,a3;

//  general variables

v = x[0];
m = 1./(1.+exp(-(v+thetam)/sm));   // m—>minf instantaneous
n = x[1];
a = x[2] ; // 1/(1+exp(-(v-thetaa)/sigmaa)); // a->ainf instantaneous
b = x[3];
se= x[4];
si= x[5];
for (i=0; i<9; i++){vCable[i] = x[6+i];}

//  currents
m3=m*m*m;
n4=n*n*n*n;
a3=a*a*a;
ina= gna*m3*(1-n)*(v-vna);
ik = gk*n4*(v-vk);
ia = ga*a3*b*(v-vk);
ilk= glk*(v-vlk);
isyni = gsyni*si*(v-vsyni);


/* update voltage at soma */
dx[0] = - (ilk + ina + ik + ia + isyni)/c;

/* update n gating */
ninf=1./(1.+exp((v-thetan)/sn));
taun=taun0+taun1/(1+exp((v+thn)/sigman));
dx[1]=phi*(ninf-n)/taun;

/* update a gating */
ainf = 1/(1+exp(-(v-thetaa)/sigmaa));
dx[2] = (ainf-a)/taua;

/* update b gating */
binf = 1./(1+exp(-(v-thetab)/sigmab));
dx[3]=(binf-b)/taub;

/* update excitation */
dx[4] = -betae*se; // instantaneous rise, exponential decay

/* update inhibition */
dx[5] = -betai*si;

/* CABLE */
for (i=0; i<9; i++){
dx[6+i] = - (glkCable*(vCable[i]-vlk)  )/c;
}

dx[0] += (d)*(vCable[0]-v)/c;
dx[6] += d*(v-2*vCable[0]+vCable[1])/c;
for (i=1;i<8;i++){ 
   dx[6+i] += (d/c)*(vCable[i-1]-2*vCable[i]+vCable[i+1]);
   } 
dx[14] += (d)*(vCable[7]-vCable[8]) /c  ;

if (synCable==0){ 
  dx[0] -= (gsyne*se*(v-vsyne)) / c; // soma
  }
else{  
  dx[5+synCable] -= (gsyne*se*(vCable[synCable-1]-vsyne)) / c;  // cable
  }

  return GSL_SUCCESS;
}







/* The jacobian
 * In case we wanna use implict Burslish-Stoer methods
 */
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;

    // Si usamos Burlish-Stoer tenemos de dar el Jacobiano
    /*gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));*/

    // Si usamos otro método podemos dejar el jacobiano
    // todo a cero.
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 0, 0.0);
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
  
    return GSL_SUCCESS;
}


