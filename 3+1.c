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
 * Example: GLoBES tour
 * Compile with ``make example-tour''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>    /* GNU Scientific library (required for root finding) */
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"            /* my input-output routines */


#define GLB_THETA_EE 6        /* Index of non-standard parameter thee*/
#define GLB_THETA_MM 7        /* Index of non-standard parameter thmm*/
#define GLB_THETA_ME 8        /* Index of non-standard parameter thme*/
#define GLB_DM_41 9          /* Index of non-standard parameter stdm*/


/* Square of real number */
inline double square(double x)
{
    return x*x;
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
   fit_lar = glbGetSignalFitRatePtr(0, rule);
   true_lar = glbGetRuleRatePtr(0, rule); 

   for (j=1; j<nexp; j++)
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

  /*Roda em todos os experimentos no canal do muon */
  if(exp==nexp)
  {
   signal_norm = 1.0 + x[0];
   bg_norm     = 1.0 + x[1];
   for (j=nexp; j<nexpmu; j++)
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

  return chi2;
}


/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

double th12;
double th13;
double th23;
double deltacp;
double sdm;
double ldm;
double thee;
double thmm;
double thme;
double stdm;

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
  sdm     = glbGetOscParams(p, GLB_DM_21) /* 1.0e-18 */;   /* Convert to GeV^2 */
  ldm     = glbGetOscParams(p, GLB_DM_31) /* 1.0e-18 */;   /* Convert to GeV^2 */
  thee    = glbGetOscParams(p, GLB_THETA_EE);
  thmm    = glbGetOscParams(p, GLB_THETA_MM);
  thme    = glbGetOscParams(p, GLB_THETA_ME);
  stdm    = glbGetOscParams(p, GLB_DM_41) /* 1.0e-18 */;   /* Convert to GeV^2 */

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
  glbSetOscParams(p, sdm /*1.0e18*/, GLB_DM_21);  /* Convert to eV^2 */
  glbSetOscParams(p, ldm /*1.0e18*/, GLB_DM_31);  /* Convert to eV^2 */ 
  glbSetOscParams(p, thee, GLB_THETA_EE);
  glbSetOscParams(p, thmm, GLB_THETA_MM);
  glbSetOscParams(p, thme, GLB_THETA_ME);
  glbSetOscParams(p, stdm /* 1.0e18 */, GLB_DM_41);  /* Convert to eV^2 */

  return 0;
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
  double Delta41;
  double see, smm, sme;
  double t;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  /* Calculate total baseline */
  L = 0.0;
  for (i=0; i < psteps; i++)
    L += length[i];
 /* L = GLB_KM_TO_EV(L) * 1.0e9;  */    /* Convert to GeV^{-1} */

  /* Compute P_ee */
  see = sin(2*thee);
/* if(L<0.2){t = 0;
  }else{t = (1.27 * L)/(E);}*/
 
  t = (1.27 * L)/(E);
  Delta41 = stdm * t;
  P[0][0] = 1 - (pow(see,2) * pow(sin(Delta41),2));

  /* Compute P_mm */
  smm = sin(2*thmm);
  P[1][1] = 1 - (pow(smm,2))* pow(sin(Delta41),2)*exp(-square(2*Delta41)*square(filter_sigma)/(2*pow(E,2)));

  /* Compute P_me */
  sme = sin(2*thme);
  P[1][0] = pow(sme,2) * pow(sin(Delta41),2);//*exp(-square(2*Delta41)*square(filter_sigma)/(2*pow(E,2)));;

  /* Compute P_em */
  P[0][1] = pow(sme,2) * pow(sin(Delta41),2);
  
  return 0;
}

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/
int main(int argc, char *argv[])
{
  /* char* MYFILE=""; */ 
  char* MYFILE="./sig-bg-SBND/lar-m-sig.dat"; /* if empty, write to screen, otherwise to file name given here */
  FILE* stream;
  if(strlen(MYFILE)>0) stream=fopen(MYFILE, "w+");
   else stream = stdout;
 
  char* MYFILE1="chiME.dat"; /* if empty, write to screen, otherwise to file name given here */
  FILE* gabi;
  if(strlen(MYFILE1)>0) gabi=fopen(MYFILE1, "w+");
   else gabi = stdout;

  /* Define my standard oscillation parameters */
  double true_theta12 = 0.0;//0.5843;
  double true_theta13 = 0.0;//0.148;
  double true_theta23 = 0.0;//0.738;
  double true_deltacp = 0.0; //0.0;
  double true_sdm = 0.0;//7.5e-05;
  double true_ldm = 0.0;//2.4e-03;
  double true_thetaee = asin(sqrt(0.0))/2;
  double true_thetamm = asin(sqrt(0.02))/2;
  double true_thetame = asin(sqrt(0.01))/2;
  double true_stdm = 1.0;

 
  glbInit(argv[0]); /* Initialize GLoBES and define chi^2 functions */
  glbDefineChiFunction(&chiGabiNorm,     2,        "chiGabiNorm",     NULL);

  /* Register non-standard probability engine. This has to be done
   * before any calls to glbAllocParams() or glbAllocProjections() */
  glbRegisterProbabilityEngine(10,      /* Number of parameters */
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);

  /* Initialize one experiment .glb */
 glbInitExperiment("./detectors/SBN_lar_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
 glbInitExperiment("./detectors/SBN_mboone_m.glb",&glb_experiment_list[0],&glb_num_of_exps); 
 glbInitExperiment("./detectors/SBN_icarus_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
 glbInitExperiment("./detectors/SBN_lar_e.glb",&glb_experiment_list[0],&glb_num_of_exps); 
 glbInitExperiment("./detectors/SBN_mboone_e.glb",&glb_experiment_list[0],&glb_num_of_exps);
 glbInitExperiment("./detectors/SBN_icarus_e.glb",&glb_experiment_list[0],&glb_num_of_exps);  


  
  /* Initialize a number of parameter vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
    
  /* Assign standard oscillation parameters */
  glbDefineParams(true_values,true_theta12,true_theta13,true_theta23, true_deltacp,true_sdm,true_ldm);  
  glbSetOscParams(true_values,true_thetaee, GLB_THETA_EE);   /* Non-standard parameter */
  glbSetOscParams(true_values,true_thetamm, GLB_THETA_MM);   /* Non-standard parameter */    
  glbSetOscParams(true_values,true_thetame, GLB_THETA_ME);   /* Non-standard parameter */
  glbSetOscParams(true_values,true_stdm, GLB_DM_41);         /* Non-standard parameter */    
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  glbDefineParams(test_values,true_theta12,true_theta13,true_theta23, true_deltacp,true_sdm,true_ldm);  
  glbSetOscParams(test_values,true_thetaee, GLB_THETA_EE);   /* Non-standard parameter */
  glbSetOscParams(test_values,true_thetamm, GLB_THETA_MM);   /* Non-standard parameter */    
  glbSetOscParams(test_values,true_thetame, GLB_THETA_ME);   /* Non-standard parameter */
  glbSetOscParams(test_values,true_stdm, GLB_DM_41);         /* Non-standard parameter */    
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  /* The simulated data are computed with "true_values" */
  glbSetOscillationParameters(true_values);
  glbSetRates(); 


  /* Return some low level information */
/* double ener=0;
  double i;
  for(i=0.025;i<4;i=i+0.005){
  ener=i;
  double resul=glbVacuumProbability(2,2,+1,ener,0.6);
  fprintf(stream,"%f %f \n",ener,resul);

} */

/*  fprintf(stream,"Simulated rates, experiment 0, rule 0: \n"); */
   int w;
    int n_bins = glbGetNumberOfBins(0);
    double *true_rates_N = glbGetSignalRatePtr(0,0);
    double *bg_rates_N = glbGetBGRatePtr(0,0);
  //  double *chan1_rates_N = glbGetChannelRatePtr(0,2,GLB_POST);
  //  double *chan2_rates_N = glbGetChannelRatePtr(0,3,GLB_POST);
    double *center_bin_N = glbGetBinCentersListPtr(0);
    double *size_bin_N = glbGetBinSizeListPtr(0); 
   for(w=0;w<n_bins;w++) fprintf(stream,"%g %g %g %g \n", center_bin_N[w], size_bin_N[w] ,true_rates_N[w], bg_rates_N[w]);

   /* int i;
  for(i=0;i<n_bins;i++){
    fprintf(gabi,"%g %g %g %g %g \n", center_bin_N[i], glbProfileProbability(0,2,2,+1,center_bin_N[i]), glbFlux(0,0,center_bin_N[i],0.11,2,+1)+glbFlux(0,0,center_bin_N[i],0.11,2,-1), glbXSection(0, 0, center_bin_N[i], 2, +1), size_bin_N[i]);
    }*/
   
// for(w=0;w<n_bins;w++) fprintf(stream,"%g %g %g %g \n", center_bin_N[w], size_bin_N[w] ,chan1_rates_N[w] +chan2_rates_N[w], bg_rates_N[w]);

//int rule_rate = glbShowRuleRates(stream, 0, 0, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG); 
 
   
 double z,y,chi2; 
   int strovet=1;   
    
   for(z=-4.0;z < 0.0 - 0.001; z=z+4.0/80)
   for(y=-2.0;y< 2.0+0.01;y=y+4.0/80) 
  {
    glbSetOscParams(test_values, asin(sqrt(pow(10,z)))/2, GLB_THETA_ME);
    glbSetOscParams(test_values,pow(10,y), GLB_DM_41); 

 
      chi2 =glbChiSys(test_values,GLB_ALL,0);
      fprintf(gabi,"%g %g %g \n",z,y, chi2); 

      printf("%i\n",strovet);
      strovet=strovet+1;
      } 
  

  /* Destroy parameter vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors); 

  if(strlen(MYFILE)>0) fclose(stream); 
  if(strlen(MYFILE1)>0) fclose(gabi);

  exit(0);
}
