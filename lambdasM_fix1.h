
  int size = 1255, ntorres=41;
  double a[1255],a2[1255],a3[1255],b[1255],c[41][1255],d[41][1255],e[41][1255];

// This function loads several functions from three main files

  int loadtowers(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "lambdam1.dat"
   q=0;  //counter of lines of the file
   file = fopen("./lambda/lambdam1.dat", "r");
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
   // Read second file "lambdam2.dat"
   q=0; //counter of lines of the file
   file = fopen("./lambda/lambdam2.dat", "r");
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
   // Read second file "lambdam2.dat"
   q=0; //counter of lines of the file
    file = fopen("./lambda/lambdam3.dat", "r");
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
