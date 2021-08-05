/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* 
 * Example: Non-Standard-Interactions and user-defined priors
 * Compile with ``make example6''
 *
 * This example is similar to Chapter 4 of hep-ph/0502147
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>    /* Standard Library of Complex Numbers */

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"            /* my input-output routines */
#include "lambdasMi_fix1.h"     
//#include "chiSBND.h" 

#define PI 3.14159265358979323846
//#define GLB_SIGMA_E 8
#define GLB_Q 6         /* Index of non-standard parameter Q */

/* If filename given, write to file; for empty filename write to screen */

char MYFILE20[]="./chi-LED/LED20i.dat";
char MYFILE30[]="./chi-LED/LED30i.dat";
char MYFILE40[]="./chi-LED/LED40i.dat";
char MYFILE50[]="./chi-LED/LED50i.dat";
char MYFILE60[]="./chi-LED/LED60i.dat";
char MYFILE70[]="./chi-LED/LED70i.dat";
char MYFILE80[]="./chi-LED/LED80i.dat";
char MYFILE90[]="./chi-LED/LED90i.dat";
char MYFILE100[]="./chi-LED/LED100i.dat";

/* Square of real number */
inline double square(double x)
{
    return x*x;
}

/*Valor absoluto de número real */
inline double absol(double value) {
  if (value < 0) {
    return -value;
  }
  else {
    return value;  
  }
}

int nexp = 3;
int nexpmu = 6;

/* Gauss likelihood (this is sufficient here due to the large event numbers
 * in a reactor experiment; for other setups, one should use Poisson statistics) */
inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}

/***************************************************************************
 * Function prior                                                      *
 ***************************************************************************
 * Calculate prior term of the form ((x - x_center)/error)^2.              *
 ***************************************************************************/
inline double prior(double x, double center, double sigma)
{
  double tmp = (x - center)/sigma;
  return tmp*tmp;
}

/***************************************************************************
 * Calculate chi^2 using Pedro's function
 ***************************************************************************/

double chiGabiNorm(int exp, int rule, int np, double *x, double *errors, void* user_data)
 {
  double chi2 = 0.0;
  double *true_rates ;
  double *signal_fit_rates;
  double *true_lar;
  double *fit_lar;
  double *bg_fit_rates;
  double fit_rate, signal_norm, bg_norm;
  double total_true_lar = 0.0;
  double total_fit_lar = 0.0;
  int ew_low, ew_high, i,j,k, aux;
 
  
  /*Roda em todos os experimentos no canal do eletron */
  if(exp==0)
  {
   signal_norm = 1.0 + x[0];
   bg_norm     = 1.0 + x[1];
   for (j=0; j<nexp; j++)
   {
    true_rates = glbGetRuleRatePtr(j, rule);
    signal_fit_rates = glbGetSignalFitRatePtr(j, rule);
    bg_fit_rates = glbGetBGFitRatePtr(j, rule);
    glbGetEnergyWindowBins(j, rule, &ew_low, &ew_high);
   
   for (i=ew_low; i <= ew_high; i++)
     {
    
    fit_rate  = signal_norm * signal_fit_rates[i] + bg_norm * bg_fit_rates[i];
    chi2 += likelihood(true_rates[i], fit_rate, true_rates[i]);
    
     }
 
    for(k=0; k < np; k++) chi2 += square(x[k]/errors[k]);
   }
  }

  /*Roda em todos os experimentos no canal do muon */
  if(exp==nexp)
    {
   signal_norm = 1.0 + x[0];
   bg_norm     = 1.0 + x[1];
   fit_lar = glbGetSignalFitRatePtr(nexp, rule);
   true_lar = glbGetRuleRatePtr(nexp, rule); 

   for (j=nexp + 1; j<nexpmu; j++)
   { 
    true_rates = glbGetRuleRatePtr(j, rule);
    signal_fit_rates = glbGetSignalFitRatePtr(j, rule);
    bg_fit_rates = glbGetBGFitRatePtr(j, rule);
    glbGetEnergyWindowBins(j, rule, &ew_low, &ew_high);
  
    for (aux=ew_low; aux <= ew_high; aux++)
     {
     total_true_lar += true_lar[aux];
     total_fit_lar += fit_lar[aux];
     }

    for (i=ew_low; i <= ew_high; i++)
     {
    
    fit_rate  = (total_true_lar/total_fit_lar)*(signal_norm * signal_fit_rates[i] +  bg_norm * bg_fit_rates[i]);
    chi2 += likelihood(true_rates[i], fit_rate, true_rates[i]);

     }
 
    for(k=0; k < np; k++) chi2 += square(x[k]/errors[k]); 
    }

    }

  return chi2;
 }
/**************************************************************
 *                   LER ARQUIVO DE m1,R                      *
 **************************************************************/                
    
/*    float m1[1266];
    float R[1266];
double leissaqui(void){
    const char lam[] = "mR.dat"; 
    FILE *laica = fopen(lam, "r");
    int size = 1266;
    int q=0;
    float a, b;    
    while (q < size) {
    fscanf (laica ,"%g %g", &a, &b);
    m1[q] = a;
    R[q] = b;
 //  printf("%g %g %g \n", a[p], b[q], c[p][q]); 
    q++;
  }
    fclose(laica);
   return 0;
}
   
*/
/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

double th12;
double th13;
double th23;
double deltacp;
double sdm;
double ldm;
int q;
//double r;
//double m;

/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{
  th12    = glbGetOscParams(p, GLB_THETA_12);
  th13    = glbGetOscParams(p, GLB_THETA_13);
  th23    = glbGetOscParams(p, GLB_THETA_23);
  deltacp = glbGetOscParams(p, GLB_DELTA_CP);
  sdm     = glbGetOscParams(p, GLB_DM_21);
  ldm     = glbGetOscParams(p, GLB_DM_31);
  q       = glbGetOscParams(p, GLB_Q);
//  r       = glbGetOscParams(p, GLB_R);
//  m       = glbGetOscParams(p, GLB_M);

  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
int my_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, th12, GLB_THETA_12);
  glbSetOscParams(p, th13, GLB_THETA_13);
  glbSetOscParams(p, th23, GLB_THETA_23);
  glbSetOscParams(p, deltacp, GLB_DELTA_CP);
  glbSetOscParams(p, sdm, GLB_DM_21);
  glbSetOscParams(p, ldm, GLB_DM_31);
  glbSetOscParams(p, q, GLB_Q);
//  glbSetOscParams(p, r, GLB_R);
//  glbSetOscParams(p, m, GLB_M);

  return 0;
}

/*******************************************************************************************************
 *                            RESOLVER EQUAÇÃO DE AUTOVALORES DE LAMBDA        (desisti)               *
 *******************************************************************************************************/
 /*   double lamb(float x, float y, int n){
       int step=50000;
       double f,g,h; 
       for(f=n; f<n+0.5+1e-5; f=f+0.5/step){
       g=1/tan(PI*(f))*PI*pow(10,2*(x+y));
       if(abso(g-f) < 1e-5){
       h=f;}}
return h; 
} */

/************************************************************************************************************
 *                                             OUTRAS FUNÇÕES                                               *
 ************************************************************************************************************/ 


/*double lam2(float x, float y){
   return sqrt(square(lamb(x,y,0)) + sdm*pow(10,2*y));
}

double lam3(float x, float y){
   return sqrt(square(lamb(x,y,0)) + ldm*pow(10,2*y));
}*/

/*double mi(float x, float y, double dm){
   double l;
   double w = sqrt(square(lambda0(x,y)) + dm*pow(10,2*y));
   if(w <= 0.5){
    l = log10(sqrt(w*tan(PI*w)/(PI*pow(10,2*y))));
    }
    else {l=-6;}
   return l; 
} */

/*double m3(float x, float y){
    double l3;
    if(lam3(x,y) <= 0.5){
     l3 = log10(sqrt(lam3(x,y)*tan(PI*lam3(x,y))/(PI*pow(10,2*y))));
    }
    return l3;
}*/

double L(double x, double y, double z){
   return (2*pow(10,2*y)*pow(10,2*x))/(pow(10,2*y)*pow(10,2*x) + square(PI)*pow(10,4*y)*pow(10,4*x) + square(z));
}

/************************************************************************************************************
 *                                           MATRIZ PMNS                                                    *
 ************************************************************************************************************/
  double complex u(int i, int j){
    int bf = 1344;
    double Lth13 = asin(sin(th13)/(sqrt(L(a3[bf],b[bf],e[0][bf]))));
    double Lth12 = asin(cos(th13)*sin(th12)/(sqrt(L(a2[bf],b[bf],d[0][bf]))*cos(Lth13)));
    double Lth23 = asin(cos(th13)*sin(th23)/(sqrt(L(a3[bf],b[bf],e[0][bf]))*cos(Lth13)));
    
    double complex U[3][3];
  U[0][0] = cos(Lth12)*cos(Lth13);
  U[0][1] = sin(Lth12)*cos(Lth13);
  U[0][2] = sin(Lth13) * cexp(-I * deltacp);

  U[1][0] = - sin(Lth12)*cos(Lth23) - cos(Lth12)*sin(Lth23)*sin(Lth13) * cexp(I*deltacp);
  U[1][1] =  cos(Lth12)*cos(Lth23) - sin(Lth12)*sin(Lth23)*sin(Lth13) * cexp(I*deltacp);
  U[1][2] =  sin(Lth23)*cos(Lth13);

  U[2][0] =  sin(Lth12)*sin(Lth23) - cos(Lth12)*cos(Lth23)*sin(Lth13) * cexp(I*deltacp);
  U[2][1] = - cos(Lth12)*sin(Lth23) - sin(Lth12)*cos(Lth23)*sin(Lth13) * cexp(I*deltacp);
  U[2][2] =  cos(Lth23)*cos(Lth13);
  return U[i][j];
}

/************************************************************************************************************
 *                                            SOMA DAS TORRES                                               *
 ************************************************************************************************************/

//#include "func_todas.h"


double complex sum11(int k, double sigma, double l, double ener){
 float sum=0;
 int i,j;
 float phi1, phi2;

 for(i=0; i<ntorres; i++)
 {
   phi1=2.54*(l/ener)*square(c[i][k])/pow(10,2*b[k]); // this calculates the phase phi_i
  for(j=0; j<ntorres; j++)
  {
    phi2=2.54*(l/ener)*square(c[j][k])/pow(10,2*b[k]); // this calculates the phase phi_j
    
    sum += L(a[k],b[k],c[i][k])*L(a[k],b[k],c[j][k])*cexp(-I*(phi1-phi2)-(pow((phi1-phi2)*sigma,2))/(2*pow(ener,2)));
  } 
 }
 return sum;
}
 
double complex sum12(int k, double sigma, double l, double ener){
 float sum=0;
 int i,j;
 float phi1, phi2;

 for(i=0; i<ntorres; i++)
 {
   phi1=2.54*(l/ener)*square(c[i][k])/pow(10,2*b[k]); // this calculates the phase phi_i
  for(j=0; j<ntorres; j++)
  {
    phi2=2.54*(l/ener)*square(d[j][k])/pow(10,2*b[k]); // this calculates the phase phi_j
    
    sum += L(a[k],b[k],c[i][k])*L(a2[k],b[k],d[j][k])*cexp(-I*(phi1-phi2)-(pow((phi1-phi2)*sigma,2))/(2*pow(ener,2)));
  } 
 }
 return sum;
}

double complex sum13(int k, double sigma, double l, double ener){
 float sum=0;
 int i,j;
 float phi1, phi2;

 for(i=0; i<ntorres; i++)
 {
   phi1=2.54*(l/ener)*square(c[i][k])/pow(10,2*b[k]); // this calculates the phase phi_i
  for(j=0; j<ntorres; j++)
  {
    phi2=2.54*(l/ener)*square(e[j][k])/pow(10,2*b[k]); // this calculates the phase phi_j
    
    sum += L(a[k],b[k],c[i][k])*L(a3[k],b[k],e[j][k])*cexp(-I*(phi1-phi2)-(pow((phi1-phi2)*sigma,2))/(2*pow(ener,2)));
  } 
 }
 return sum;
}

double complex sum22(int k, double sigma, double l, double ener){
 float sum=0;
 int i,j;
 float phi1, phi2;

 for(i=0; i<ntorres; i++)
 {
   phi1=2.54*(l/ener)*square(d[i][k])/pow(10,2*b[k]); // this calculates the phase phi_i
  for(j=0; j<ntorres; j++)
  {
    phi2=2.54*(l/ener)*square(d[j][k])/pow(10,2*b[k]); // this calculates the phase phi_j
    
    sum += L(a2[k],b[k],d[i][k])*L(a2[k],b[k],d[j][k])*cexp(-I*(phi1-phi2)-(pow((phi1-phi2)*sigma,2))/(2*pow(ener,2)));
  } 
 }
 return sum;
}

double complex sum23(int k, double sigma, double l, double ener){
 float sum=0;
 int i,j;
 float phi1, phi2;

 for(i=0; i<ntorres; i++)
 {
   phi1=2.54*(l/ener)*square(d[i][k])/pow(10,2*b[k]); // this calculates the phase phi_i
  for(j=0; j<ntorres; j++)
  {
    phi2=2.54*(l/ener)*square(e[j][k])/pow(10,2*b[k]); // this calculates the phase phi_j
    
    sum += L(a2[k],b[k],d[i][k])*L(a3[k],b[k],e[j][k])*cexp(-I*(phi1-phi2)-(pow((phi1-phi2)*sigma,2))/(2*pow(ener,2)));
  } 
 }
 return sum;
}

double complex sum33(int k, double sigma, double l, double ener){
 float sum=0;
 int i,j;
 float phi1, phi2;

 for(i=0; i<ntorres; i++)
 {
   phi1=2.54*(l/ener)*square(e[i][k])/pow(10,2*b[k]); // this calculates the phase phi_i
  for(j=0; j<ntorres; j++)
  {
    phi2=2.54*(l/ener)*square(e[j][k])/pow(10,2*b[k]); // this calculates the phase phi_j
    
    sum += L(a3[k],b[k],e[i][k])*L(a3[k],b[k],e[j][k])*cexp(-I*(phi1-phi2)-(pow((phi1-phi2)*sigma,2))/(2*pow(ener,2)));
  } 
 }
 return sum;
}

/***************************************************************************
 * Calculate oscillation probabilities.                                    *
 * Since for our setup, only P_ee is required, all other entries of P are  *
 * set to zero for simplicity. Furthermore, we neglect matter effects and  *
 * the filter feature (parameter filter_sigma).                            *
 * The formula for P_ee is Eq. (36) from hep-ph/0502147.                   *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:            The buffer where the probabilities are to be stored     *
 *   cp_sign:      +1 if probalities for neutrinos are requested, -1 for   *
 *                 anti-neutrinos.                                         *
 *   E:            The neutrino energy in GeV                              *
 *   psteps:       Number of constant density layers in the matter profile *
 *   length:       The lengths of these layers in km                       *
 *   density:      The individual densities of these layers in g/cm^3      *
 *   filter_sigma: Width of low-pass filter as given in the AEDL file      *
 ***************************************************************************/
int my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  int i, j;
  double L;
  
 /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  /* Calculate total baseline */
  L = 0.0;
  for (i=0; i < psteps; i++)
    L += length[i];
//  L = GLB_KM_TO_EV(L) * 1.0e9;      /* Convert to GeV^{-1} */
 
 /* LED */
  double complex s11, s12,  s13, s21 ,s22, s23, s31, s32, s33;

  s11 = sum11(q,filter_sigma,L,E);
  s12 = sum12(q,filter_sigma,L,E);
  s13 = sum13(q,filter_sigma,L,E);
  s21 = conj(s12);
  s22 = sum22(q,filter_sigma,L,E);
  s23 = sum23(q,filter_sigma,L,E);
  s31 = conj(s13);
  s32 = conj(s23);
  s33 = sum33(q,filter_sigma,L,E);
  
  
 /* Compute P_mm */
  P[1][1] = creal(conj(u(1,0))*u(1,0)*(u(1,0)*conj(u(1,0))*s11 + u(1,1)*conj(u(1,1))*s12 + u(1,2)*conj(u(1,2))*s13) + conj(u(1,1))*u(1,1)*(u(1,0)*conj(u(1,0))*s21 + u(1,1)*conj(u(1,1))*s22 + u(1,2)*conj(u(1,2))*s23) + conj(u(1,2))*u(1,2)*(u(1,0)*conj(u(1,0))*s31 + u(1,1)*conj(u(1,1))*s32 + u(1,2)*conj(u(1,2))*s33));

  /* Compute P_ee */
  P[0][0] = creal(conj(u(0,0))*u(0,0)*(u(0,0)*conj(u(0,0))*s11 + u(0,1)*conj(u(0,1))*s12 + u(0,2)*conj(u(0,2))*s13) + conj(u(0,1))*u(0,1)*(u(0,0)*conj(u(0,0))*s21 + u(0,1)*conj(u(0,1))*s22 + u(0,2)*conj(u(0,2))*s23) + conj(u(0,2))*u(0,2)*(u(0,0)*conj(u(0,0))*s31 + u(0,1)*conj(u(0,1))*s32 + u(0,2)*conj(u(0,2))*s33)); 

 if (cp_sign >= 0)
  {
  /* Compute P_me */
  P[1][0] = creal(conj(u(1,0))*u(0,0)*(u(1,0)*conj(u(0,0))*s11 + u(1,1)*conj(u(0,1))*s12 + u(1,2)*conj(u(0,2))*s13) + conj(u(1,1))*u(0,1)*(u(1,0)*conj(u(0,0))*s21 + u(1,1)*conj(u(0,1))*s22 + u(1,2)*conj(u(0,2))*s23) + conj(u(1,2))*u(0,2)*(u(1,0)*conj(u(0,0))*s31 + u(1,1)*conj(u(0,1))*s32 + u(1,2)*conj(u(0,2))*s33));
  
  }
else
{
 P[1][0] = creal(conj(u(0,0))*u(1,0)*(u(0,0)*conj(u(1,0))*s11 + u(0,1)*conj(u(1,1))*s12 + u(0,2)*conj(u(1,2))*s13) + conj(u(0,1))*u(1,1)*(u(0,0)*conj(u(1,0))*s21 + u(0,1)*conj(u(1,1))*s22 + u(0,2)*conj(u(1,2))*s23) + conj(u(0,2))*u(1,2)*(u(0,0)*conj(u(1,0))*s31 + u(0,1)*conj(u(1,1))*s32 + u(0,2)*conj(u(1,2))*s33));

}

 return 0;
}


/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{ 
 
  char* MYFILE="ev3_LED.dat"; 
  FILE* stream;
  if(strlen(MYFILE)>0) stream=fopen(MYFILE, "w");
   else stream = stdout; 

 // leissaqui(); 
 // le1();
 // le2();
 // le3();
  loadtowers1();
  
  /* Define standard oscillation parameters (cf. hep-ph/0405172v5) */
  double true_theta12 = asin(sqrt(0.307));
  double true_theta13 = asin(sqrt(0.02227));
  double true_theta23 = asin(sqrt(0.554));
  double true_deltacp = 1.54*PI;
  double true_q = 1344;
  double true_sdm = 7.4e-05;
  double true_ldm = -2.46e-03;
//  double true_m = -3.0;
//  double true_r = -2.0;
 

  /* Define one non-standard parameter sigma_E (wave packet energy spread
   * responsible for wave packet decoherence) */
//  double true_sigma_E = 0.0;

  /* Initialize libglobes */
  glbInit(argv[0]);
  glbDefineChiFunction(&chiGabiNorm,     2,        "chiGabiNorm",     NULL);
 
  /* Register non-standard probability engine. This has to be done
   * before any calls to glbAllocParams() or glbAllocProjections() */
  glbRegisterProbabilityEngine(7,      /* Number of parameters */
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);

  /* Initialize reactor experiment */
 glbInitExperiment("./detectors/SBN_lar_e.glb",&glb_experiment_list[0],&glb_num_of_exps); 
 glbInitExperiment("./detectors/SBN_mboone_e.glb",&glb_experiment_list[0],&glb_num_of_exps);
 glbInitExperiment("./detectors/SBN_icarus_e.glb",&glb_experiment_list[0],&glb_num_of_exps);
 
  glbInitExperiment("./detectors/SBN_lar_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
  glbInitExperiment("./detectors/SBN_mboone_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
  glbInitExperiment("./detectors/SBN_icarus_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
 

  /* Initialize parameter and projection vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection gabi_projection = glbAllocProjection();

  glbDefineParams(true_values,true_theta12,true_theta13,true_theta23, true_deltacp,true_sdm,true_ldm);
  glbSetOscParams(true_values,true_q, GLB_Q);   /* Non-standard parameter */
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  glbDefineParams(test_values,true_theta12,true_theta13,true_theta23, true_deltacp,true_sdm,true_ldm);
  glbSetOscParams(test_values,true_q, GLB_Q);   /* Non-standard parameter */
  glbSetDensityParams(test_values,1.0,GLB_ALL);
   
  /* Set starting values and input errors for all projections */  
  glbDefineParams(input_errors, 0.0,0.0,0.0,0.0,0.0,0.0);
  glbSetOscParams(input_errors,0.0, GLB_Q); 
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* The simulated data are computed */
  glbSetOscillationParameters(true_values);
  glbSetRates();
  
/*                                      THETA_12,  THETA_13,  THETA_23,  DELTA_CP,  DM_21,    DM_31 */
  glbDefineProjection(gabi_projection, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED);
  glbSetDensityProjectionFlag(gabi_projection, GLB_FIXED, GLB_ALL);
  glbSetProjectionFlag(gabi_projection, GLB_FIXED, GLB_Q);
  glbSetProjection(gabi_projection);
   
   /*  int w;
    int n_bins = glbGetNumberOfBins(0);
    double *true_rates_N = glbGetSignalRatePtr(0,0);
    double *bg_rates_N = glbGetBGRatePtr(0,0);
    double *center_bin_N = glbGetBinCentersListPtr(0);
    double *size_bin_N = glbGetBinSizeListPtr(0); 
    for(w=0;w<n_bins;w++) fprintf(stream,"%g %g %g %g\n", center_bin_N[w], size_bin_N[w] , true_rates_N[w], bg_rates_N[w]);
    */

//  for(w=0;w<n_bins;w++) fprintf(stream,"%g %g %g %g \n", center_bin_N[w], size_bin_N[w] ,chan1_rates_N[w] +chan2_rates_N[w], bg_rates_N[w]);
  

  /* Compute chi^2 without correlations */
   
     int x;
     double res;
   InitOutput(MYFILE50, ""); 
    int o=0;
    for(x=0; x < 2090; x++) {       
          
      glbSetOscParams(test_values,x,GLB_Q);
     
      res=glbChiNP(test_values,NULL,GLB_ALL);
  
    AddToOutput(a3[x],b[x],res);
    printf("%d \n", o);
    o++;
    }

   /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors); 
  glbFreeProjection(gabi_projection);
  
  if(strlen(MYFILE)>0) fclose(stream);

  exit(0);
}

