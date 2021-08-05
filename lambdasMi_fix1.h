
  int size = 2090, ntorres=41;
  double a[2090],a2[2090],a3[2090],b[2090],c[41][2090],d[41][2090],e[41][2090];

// This function loads several functions from three main files

  int loadtowers1(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "lm1.dat"
   q=0;  //counter of lines of the file
   file = fopen("./lambda_inv/L1.dat", "r");
   while (q < size) 
   {
    state=fscanf(file,"%lf", &a[q]); //Scan first two parameters of each column (variables)
    for(qi=0; qi<ntorres;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &c[qi][q]); // scan the column functions
     }
    q++;
   }
   fclose(file);
   // Read second file "lm2.dat"
   q=0; //counter of lines of the file
   file = fopen("./lambda_inv/L2.dat", "r");
   while (q < size) 
   {
    state=fscanf (file,"%lf", &a2[q]); //Scan first parameter of each column (variable). Notice that it has only 1 not 2 as the first loop.
    for(qi=0; qi<ntorres;qi++) // counter of columns of the file ntorres gives the number of columns
    { 
     state=fscanf(file,"%lf", &d[qi][q]); // scan the column functions
    }
    q++;
   }
   fclose(file);
   // Read second file "lm3.dat"
   q=0; //counter of lines of the file
    file = fopen("./lambda_inv/L3.dat", "r");
    while (q < size) 
    {
      state=fscanf (file,"%lf %lf", &a3[q], &b[q]); //Scan first parameter of each column (variable). Notice that it has only 1 not 2 as the first loop.
     for(qi=0; qi<ntorres;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &e[qi][q]); // scan the column functions
     }
     q++;
    }
   fclose(file);

   return 0;
  }
