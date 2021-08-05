/* RATES AND CHI-SQUARED FUNCTIONS FOR LSND */

inline double likelihood(double true_rate, double fit_rate)
{
  double res;
  res = fit_rate - true_rate;
  if (true_rate > 0)
  {
    if (fit_rate <= 0.0)
      res = 1e100;
    else
      res += true_rate * log(true_rate/fit_rate);
  }
  else
    res = fabs(res);

  return 2.0 * res;
}

/* Square of real number */
inline double square(double x)
{
    return x*x;
}

int nexp = 3;
int nexpmu = 6;


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

double chiSBN(int exp, int rule, int np, double *x, double *errors, void* user_data)
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
    chi2 += likelihood(true_rates[i], fit_rate/*, true_rates[i]*/);
    
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
    chi2 += likelihood(true_rates[i], fit_rate/*, true_rates[i]*/);

     }
 
    for(k=0; k < np; k++) chi2 += square(x[k]/errors[k]); 
    }

    }

  return chi2;
 }

