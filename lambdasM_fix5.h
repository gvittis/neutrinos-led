
  int size = 2030, ntorres=41;
  double a[2030],a2[2030],a3[2030],b[2030],c[41][2030],d[41][2030],e[41][2030];

// This function loads several functions from three main files

  int loadtowers5(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "lm1.dat"
   q=0;  //counter of lines of the file
   file = fopen("./lambda/L1.dat", "r");
   while (q < size) 
   {
    state=fscanf(file,"%lf %lf", &a[q], &b[q]); //Scan first two parameters of each column (variables)
    for(qi=0; qi<ntorres;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &c[qi][q]); // scan the column functions
     }
    q++;
   }
   fclose(file);
   // Read second file "lm2.dat"
   q=0; //counter of lines of the file
   file = fopen("./lambda/L2.dat", "r");
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
    file = fopen("./lambda/L3.dat", "r");
    while (q < size) 
    {
     state=fscanf (file,"%lf", &a3[q]); //Scan first parameter of each column (variable). Notice that it has only 1 not 2 as the first loop.
     for(qi=0; qi<ntorres;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &e[qi][q]); // scan the column functions
     }
     q++;
    }
   fclose(file);

   return 0;
  }
