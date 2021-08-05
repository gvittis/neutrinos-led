
  int size = 2601, ntorres=41;
  double a[2601],a2[2601],a3[2601],b[2601],c[41][2601],d[41][2601],e[41][2601];

// This function loads several functions from three main files

  int loadtowers6(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "lm1.dat"
   q=0;  //counter of lines of the file
   file = fopen("./lambda/l1.dat", "r");
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
   file = fopen("./lambda/l2.dat", "r");
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
    file = fopen("./lambda/l3.dat", "r");
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
