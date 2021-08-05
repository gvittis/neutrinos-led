
  int size = 1266, ntorres=3;
  double a[1266],a2[1266],a3[1266],b[1266],c[3][1266],d[3][1266],e[3][1266];

// This function loads several functions from three main files

  int loadtowers2(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "lambM1.dat"
   q=0;  //counter of lines of the file
   file = fopen("./lambda/lambM1.dat", "r");
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
   // Read second file "lambM2.dat"
   q=0; //counter of lines of the file
   file = fopen("./lambda/lambM2.dat", "r");
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
   // Read second file "lambM3.dat"
   q=0; //counter of lines of the file
    file = fopen("./lambda/lambM3.dat", "r");
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
