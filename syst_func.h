
/*

The following code fragment from the much more extensive example5.c calculates
χ
2
for a reactor setup with two detectors including five different systematical errors:
the flux normalization of the reactor (x[0]), the fiducial mass uncertainties of the
far (x[1]) and near (x[2]) detector, and energy calibration errors for the far (x[3])
and near (x[4]) detectors. The calculation follows Eq. (3.3), but includes the energy
calibration.


Function 3.3 double (*glb_chi_function)(int exp, int rule, int dim, double *params, double *errors, void* user_data).
In this function type, 

exp is the experiment number, 
rule is the rule number,
dim is the number of systematics parameters 
params is an array of the systematics nuisance parameters themselves,
errors is an array with the systematical errors. 


The parameter user_data is set as defined with glbDefineChiFunction. 
Note the the central values for params are 0, which means that 0 corresponds to the un-modified rates.
The array indices run from 0 to dim-1.

*/


double chiDCNorm(int exp, int rule, int np, double *x, double *errors, void* user_data)
 {
  int const EXP_DUNE = 0;
  double chi2 = 0.0, chi2auxS=0.0, chi2auxB=0.0, chi2auxSB=0.0, chiq=0.0;
  double *true_rates_DUNE ;
  double *signal_fit_rates_DUNE;
  double *bg_fit_rates_DUNE;
  int ew_low, ew_high, i,j;
  double ws=0.0, wb=0.0, Q,Qaux;
  
/*this part is the Beam normalization minimization */
  if(exp==0)
  {
   for (j=0; j<nexp; j++)
   {
    true_rates_DUNE = glbGetRuleRatePtr(j, rule);
    signal_fit_rates_DUNE=glbGetSignalFitRatePtr(j, rule);
    bg_fit_rates_DUNE=glbGetBGFitRatePtr(j, rule);
    glbGetEnergyWindowBins(j, rule, &ew_low, &ew_high);
    for (i=ew_low; i <= ew_high; i++)
     {
       Qaux=(true_rates_DUNE[i]+pow(errors[2]*signal_fit_rates_DUNE[i],2)+pow(errors[3]*bg_fit_rates_DUNE[i],2));
       ws+=pow(signal_fit_rates_DUNE[i],2)/Qaux;
       wb+=pow(bg_fit_rates_DUNE[i],2)/Qaux;
       chi2auxS+=((true_rates_DUNE[i]-(signal_fit_rates_DUNE[i]+bg_fit_rates_DUNE[i]))*signal_fit_rates_DUNE[i])/Qaux;
       chi2auxB+=((true_rates_DUNE[i]-(signal_fit_rates_DUNE[i]+bg_fit_rates_DUNE[i]))*bg_fit_rates_DUNE[i])/Qaux;
       chi2auxSB+=signal_fit_rates_DUNE[i]*bg_fit_rates_DUNE[i]/Qaux;
     }
    }
   ws+=1/pow(errors[0],2);
   wb+=1/pow(errors[1],2)-pow(chi2auxSB,2)/ws;
   chi2 +=-(1/ws)*pow(chi2auxS,2)*(1+pow(chi2auxSB,2)/(ws*wb))-pow(chi2auxB,2)/wb+2*chi2auxSB*chi2auxB*chi2auxS/(ws*wb);
   }


  if(exp==nexp)
  {
   for (j=nexp; j<nexpmu; j++)
   {
    true_rates_DUNE = glbGetRuleRatePtr(j, rule);
    signal_fit_rates_DUNE=glbGetSignalFitRatePtr(j, rule);
    bg_fit_rates_DUNE=glbGetBGFitRatePtr(j, rule);
    glbGetEnergyWindowBins(j, rule, &ew_low, &ew_high);
    for (i=ew_low; i <= ew_high; i++)
     {
       Qaux=(true_rates_DUNE[i]+pow(errors[2]*signal_fit_rates_DUNE[i],2)+pow(errors[3]*bg_fit_rates_DUNE[i],2));
       ws+=pow(signal_fit_rates_DUNE[i],2)/Qaux;
       wb+=pow(bg_fit_rates_DUNE[i],2)/Qaux;
       chi2auxS+=((true_rates_DUNE[i]-(signal_fit_rates_DUNE[i]+bg_fit_rates_DUNE[i]))*signal_fit_rates_DUNE[i])/Qaux;
       chi2auxB+=((true_rates_DUNE[i]-(signal_fit_rates_DUNE[i]+bg_fit_rates_DUNE[i]))*bg_fit_rates_DUNE[i])/Qaux;
       chi2auxSB+=signal_fit_rates_DUNE[i]*bg_fit_rates_DUNE[i]/Qaux;
     }
    }
   ws+=1/pow(errors[0],2);
   wb+=1/pow(errors[1],2)-pow(chi2auxSB,2)/ws;
   chi2 +=-(1/ws)*pow(chi2auxS,2)*(1+pow(chi2auxSB,2)/(ws*wb))-pow(chi2auxB,2)/wb+2*chi2auxSB*chi2auxB*chi2auxS/(ws*wb);
   }


    true_rates_DUNE = glbGetRuleRatePtr(exp, rule);
    signal_fit_rates_DUNE=glbGetSignalFitRatePtr(exp, rule);
    bg_fit_rates_DUNE=glbGetBGFitRatePtr(exp, rule);
    glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  for (i=ew_low; i <= ew_high; i++)
   {
     Qaux=(true_rates_DUNE[i]+pow(errors[2]*signal_fit_rates_DUNE[i],2)+pow(errors[3]*bg_fit_rates_DUNE[i],2));
     chi2 += pow(true_rates_DUNE[i]-(signal_fit_rates_DUNE[i]+bg_fit_rates_DUNE[i]),2)/Qaux;
   }
   return chi2;
}







/*




double chiDCNorm(int exp, int rule, int np, double *x, double *errors, void* user_data)
 {
  const EXP_DUNE = 0; const EXP_NEAR = 1;
  int n_bins = glbGetNumberOfBins(EXP_FAR);
  double *true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 0);
  double *true_rates_F = glbGetRuleRatePtr(EXP_FAR, 0);
  double signal_fit_rates_N[n_bins]; double signal_fit_rates_F[n_bins];
  double signal_norm_N, signal_norm_F;
  int ew_low, ew_high, i;
  double emin, emax, fit_rate; double chi2 = 0.0;

  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0),
  signal_fit_rates_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0),
  signal_fit_rates_N, n_bins, emin, emax);

  signal_norm_F = 1.0 + x[0] + x[1];
  signal_norm_N = 1.0 + x[0] + x[2];
  for (i=ew_low; i <= ew_high; i++)
   {

    fit_rate = signal_norm_F * signal_fit_rates_F[i];
    chi2 += likelihood(true_rates_F[i], fit_rate, true_rates_F[i]);

    fit_rate = signal_norm_N * signal_fit_rates_N[i];
    chi2 += likelihood(true_rates_N[i], fit_rate, true_rates_N[i]);
   }

   for (i=0; i < np; i++) chi2 += square(x[i]/errors[i]);
   return chi2;
}



*/
