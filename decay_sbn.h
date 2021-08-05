// This function loads SBN fluxes
  int en = 501, flu=7;
  double icarus[501][7],mboone[501][7],larnd[501][7];

// This function loads several functions from three main files

  int loadflux(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "ICARUSplus.dat"
   q=0;  //counter of lines of the file
   file = fopen("./flux/ICARUSplus.dat", "r");
   while (q < en) 
   {
     // state=fscanf(file,"%lf %lf", &a[q], &b[q]); //Scan first two parameters of each column (variables)
    for(qi=0; qi<flu;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &icarus[q][qi]); // scan the column functions
     }
    q++;
   }
   fclose(file);
   // Read second file "MICROplus.dat"
   q=0; //counter of lines of the file
   file = fopen("./flux/MICROplus.dat", "r");
   while (q < en) 
   {
     //  state=fscanf (file,"%lf", &a2[q]); //Scan first parameter of each column (variable). Notice that it has only 1 not 2 as the first loop.
    for(qi=0; qi<flu;qi++) // counter of columns of the file ntorres gives the number of columns
    { 
     state=fscanf(file,"%lf", &mboone[q][qi]); // scan the column functions
    }
    q++;
   }
   fclose(file);
   // Read second file "LARplus.dat"
   q=0; //counter of lines of the file
    file = fopen("./flux/LARplus.dat", "r");
    while (q < en) 
    {
 //    state=fscanf (file,"%lf", &a3[q]); //Scan first parameter of each column (variable). Notice that it has only 1 not 2 as the first loop.
     for(qi=0; qi<flu;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &larnd[q][qi]); // scan the column functions
     }
     q++;
    }
   fclose(file);

  return 0;
  }


/***************************************************************************
 *                     Return the interpolated flux                        *
 ***************************************************************************/

double icarus_flux(double E, int f /* 1 a 6 */)
{
  int col = f;
  int n_steps = 500;
  double result;

  int k=0;
  while (k <= n_steps  &&  icarus[k][0] < E)
    k++;
  if (k <= 0 || k > n_steps)
    return 0.0;
  else
  {
    double E_lo    = icarus[k-1][0];
    double E_up    = icarus[k][0];
    double flux_lo = icarus[k-1][col];
    double flux_up = icarus[k][col];
    result  = flux_lo + (E - E_lo)*(flux_up - flux_lo)/(E_up - E_lo);
  }
  return result;
}

double mboone_flux(double E, int f /* 1 a 6 */)
{
  int col = f;
  int n_steps = 500;
  double result;

  int k=0;
  while (k <= n_steps  &&  mboone[k][0] < E)
    k++;
  if (k <= 0 || k > n_steps)
    return 0.0;
  else
  {
    double E_lo    = mboone[k-1][0];
    double E_up    = mboone[k][0];
    double flux_lo = mboone[k-1][col];
    double flux_up = mboone[k][col];
    result  = flux_lo + (E - E_lo)*(flux_up - flux_lo)/(E_up - E_lo);
  }
  return result;
}

double larnd_flux(double E, int f /* 1 a 6 */)
{
  int col = f;
  int n_steps = 500;
  double result;

  int k=0;
  while (k <= n_steps  &&  larnd[k][0] < E)
    k++;
  if (k <= 0 || k > n_steps)
    return 0.0;
  else
  {
    double E_lo    = larnd[k-1][0];
    double E_up    = larnd[k][0];
    double flux_lo = larnd[k-1][col];
    double flux_up = larnd[k][col];
    result  = flux_lo + (E - E_lo)*(flux_up - flux_lo)/(E_up - E_lo);
  }
  return result;
}

/***************************************************************************
 *                   Probability functions                                 *
 ***************************************************************************/

double complex decay_icarus_1(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 3.565904 - 3.556016;
     double complex dPdE1=0;
             
for (En = Erec + binwidth; En < icarus[500][0]; En= En+binwidth)
{
  dPdE1 = 0.5*2*Erec/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE1*icarus_flux(En,f)*binwidth;
 }
     return decay;
}

double complex decay_icarus_2(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 3.565904 - 3.556016;
     double complex dPdE2=0;
     
for (En = Erec + binwidth; En < icarus[500][0]; En= En+binwidth)
{
  dPdE2 = 0.5*2*(En-Erec)/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE2*icarus_flux(En,f)*binwidth;
 }
     return decay;
}

double complex decay_mboone_1(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 3.565904 - 3.556016;
     double complex dPdE1=0;
             
for (En = Erec + binwidth; En < mboone[500][0]; En= En+binwidth)
{
  dPdE1 = 0.5*2*Erec/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE1*mboone_flux(En,f)*binwidth;
 }
     return decay;
}

double complex decay_mboone_2(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 3.565904 - 3.556016;
     double complex dPdE2=0;
     
for (En = Erec + binwidth; En < mboone[500][0]; En= En+binwidth)
{
  dPdE2 = 0.5*2*(En-Erec)/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE2*mboone_flux(En,f)*binwidth;
 }
     return decay;
}

double complex decay_larnd_1(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 3.565904 - 3.556016;
     double complex dPdE1=0;
             
for (En = Erec + binwidth; En < larnd[500][0]; En= En+binwidth)
{
  dPdE1 = 0.5*2*Erec/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE1*larnd_flux(En,f)*binwidth;
 }
     return decay;
}

double complex decay_larnd_2(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 3.565904 - 3.556016;
     double complex dPdE2=0;
     
for (En = Erec + binwidth; En < larnd[500][0]; En= En+binwidth)
{
  dPdE2 = 0.5*2*(En-Erec)/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE2*larnd_flux(En,f)*binwidth;
 }
     return decay;
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

double u;
double mg;

/***************************************************************************
 * Store oscillation parameters in internal data structures.
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{

  th12    = glbGetOscParams(p, GLB_THETA_12);
  th13    = glbGetOscParams(p, GLB_THETA_13);
  th23    = glbGetOscParams(p, GLB_THETA_23);
  deltacp = glbGetOscParams(p, GLB_DELTA_CP);
  sdm     = glbGetOscParams(p, GLB_DM_21);   
  ldm     = glbGetOscParams(p, GLB_DM_31);

  u     = glbGetOscParams(p, GLB_U);
  mg    = glbGetOscParams(p, GLB_MG);

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

  glbSetOscParams(p, u, GLB_U);
  glbSetOscParams(p, mg, GLB_MG);

  return 0;
}

/***************************************************************************
 * Calculate oscillation probabilities.                                    *
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

  if (cp_sign >= 0)
   {
     if(fabs(L-0.11)<1.0e-8)
       {
     /*nu_e*/    
     P[1][0] = decay_larnd_1(u,mg,2,E,L) /*nu_mu*/ + decay_larnd_2(u,mg,5,E,L) /*anu_mu*/;
     /*nu_mu*/
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*larnd_flux(E,2); /*R_mu = 0*/
       }
     if(fabs(L-0.47)<1.0e-8)
       { 
     /*nu_e*/    
     P[1][0] = decay_larnd_1(u,mg,2,E,L) /*nu_mu*/ + decay_larnd_2(u,mg,5,E,L) /*anu_mu*/;
     /*nu_mu*/
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*larnd_flux(E,2); /*R_mu = 0*/
       }
      if(fabs(L-0.6)<1.0e-8)
      {
     /*nu_e*/    
     P[1][0] = decay_larnd_1(u,mg,2,E,L) /*nu_mu*/ + decay_larnd_2(u,mg,5,E,L) /*anu_mu*/;
     /*nu_mu*/
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*larnd_flux(E,2); /*R_mu = 0*/
      }
   }
 else
   {
    if(fabs(L-0.11)<1.0e-8)
       {
     /*anu_e*/    
     P[1][0] = decay_larnd_1(u,mg,5,E,L) /*anu_mu*/ + decay_larnd_2(u,mg,2,E,L) /*nu_mu*/;
     /*anu_mu*/
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*larnd_flux(E,5); /*R_mu = 0*/
       }
     if(fabs(L-0.47)<1.0e-8)
       {
      /*anu_e*/    
     P[1][0] = decay_larnd_1(u,mg,5,E,L) /*anu_mu*/ + decay_larnd_2(u,mg,2,E,L) /*nu_mu*/;
     /*anu_mu*/
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*larnd_flux(E,5); /*R_mu = 0*/
       }
      if(fabs(L-0.6)<1.0e-8)
      {
     /*anu_e*/    
     P[1][0] = decay_larnd_1(u,mg,5,E,L) /*anu_mu*/ + decay_larnd_2(u,mg,2,E,L) /*nu_mu*/;
     /*anu_mu*/
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*larnd_flux(E,5); /*R_mu = 0*/
      }
   }

 return 0;
}

