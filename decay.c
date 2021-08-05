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
#include <gsl/gsl_math.h>    /* GNU Scientific library (required for root finding) */
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"            /* my input-output routines */
#include "syst_func_SBN.h"

#define PI 3.14159265358979323846

#define GLB_U 6               /* Index of |U_\mu4|^2 */
#define GLB_MG 7              /* Index of m4\bar{g} */

#include "decay_sbn.h"

/* If filename given, write to file; for empty filename write to screen */

char MYFILE1[]="./results/chi-dec-sbn.dat";

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{ 
 
  char* MYFILE="./results/ev-decay-sbn.dat"; 
  FILE* stream;
  if(strlen(MYFILE)>0) stream=fopen(MYFILE, "w+");
  else stream = stdout; 

  loadflux();
  
  /* Define standard oscillation parameters (cf. hep-ph/0405172v5) */
  double true_theta12 = 0.0;
  double true_theta13 = 0.0;
  double true_theta23 = 0.0;
  double true_deltacp = 0.0; 
  double true_sdm = 0.0;
  double true_ldm = 0.0;

  double true_u  = 0.0;
  double true_mg  = 0.0;
  
  /* Initialize libglobes */
  glbInit(argv[0]);
  glbDefineChiFunction(&chiSBN,     2,        "chiSBN",     NULL);
 
  /* Register non-standard probability engine. This has to be done
   * before any calls to glbAllocParams() or glbAllocProjections() */
  glbRegisterProbabilityEngine(8,      /* Number of parameters */
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);

  /* Initialize experiments */
  glbInitExperiment("./detectors/SBN_lar_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
  glbInitExperiment("./detectors/SBN_mboone_m.glb",&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment("./detectors/SBN_icarus_m.glb",&glb_experiment_list[0],&glb_num_of_exps);
  glbInitExperiment("./detectors/SBN_lar_e.glb",&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment("./detectors/SBN_mboone_e.glb",&glb_experiment_list[0],&glb_num_of_exps);
  glbInitExperiment("./detectors/SBN_icarus_e.glb",&glb_experiment_list[0],&glb_num_of_exps); 

  /* Initialize parameter and projection vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
 
  glbDefineParams(true_values,true_theta12,true_theta13,true_theta23,true_deltacp,true_sdm,true_ldm);
  glbSetOscParams(true_values, true_u, GLB_U);   /* Extra parameter */
  glbSetOscParams(true_values, true_mg, GLB_MG);   /* Extra parameter */
  
  glbSetDensityParams(true_values,1.0,GLB_ALL);

  glbDefineParams(test_values,true_theta12,true_theta13,true_theta23,true_deltacp,true_sdm,true_ldm);
  glbSetOscParams(test_values, true_u, GLB_U);   /* Extra parameter */
  glbSetOscParams(test_values, true_mg, GLB_MG);   /* Extra parameter */

  glbSetDensityParams(test_values,1.0,GLB_ALL);
   
  /* Set starting values and input errors for all projections */  
  glbDefineParams(input_errors, 0.0,0.0,0.0,0.0,0.0,0.0);
  glbSetOscParams(input_errors,0.0,  GLB_U); 
  glbSetOscParams(input_errors,0.0,  GLB_MG);
    
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
 
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* The simulated data are computed */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /*
    int w;
    int n_bins = glbGetNumberOfBins(0);
    double *true_rates_N = glbGetSignalRatePtr(0,0);
    double *bg_rates_N = glbGetBGRatePtr(0,0);
    double *center_bin_N = glbGetBinCentersListPtr(0);
    double *size_bin_N = glbGetBinSizeListPtr(0); 
    for(w=0;w<n_bins;w++) fprintf(stream,"%g %g %g %g\n", center_bin_N[w], size_bin_N[w] , true_rates_N[w], bg_rates_N[w]);*/
    
    printf("teu cu\n");
  /* Compute chi^2 without correlations */
    double a,b;
    double res;
   
    InitOutput(MYFILE1, ""); 
    int o=0;

    for(a=-3; a < 0 - 0.001; a=a+0.05)
    for(b=-1; b < 2 - 0.001; b=b+0.05)
     {       
          
      glbSetOscParams(test_values,pow(10,a),GLB_U);
      glbSetOscParams(test_values,pow(10,b),GLB_MG);
     
      res=glbChiSys(test_values,0,GLB_ALL);
  
    AddToOutput(a,b,res);
    printf("%d \n", o);
    o++;
      } 

   /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors); 
    
  if(strlen(MYFILE)>0) fclose(stream);

  exit(0);
}

