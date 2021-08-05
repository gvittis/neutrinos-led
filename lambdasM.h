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
 Some simple output functions for GLoBES examples
 */

/* Init standard output and write the string headline into it */
/* If filename given, write to file; for empty filename write to screen */

/**************************************************************
 *                FUNÇÃO LAMBDAS                              *
 **************************************************************/ 
   
    double a[1266];
    double b[1266];
    double c0[1266];    
    double c1[1266];   
    double c2[1266];
    double c3[1266];   
    double c4[1266];
    double c5[1266];   
    double c6[1266];
    double c7[1266];   
    double c8[1266];
    double c9[1266];   
    double c10[1266];
    double c11[1266];   
    double c12[1266];
    double c13[1266];   
    double c14[1266];
    double c15[1266];   
    double c16[1266];
    double c17[1266];   
    double c18[1266];
    double c19[1266];   
    double c20[1266];
    double c21[1266];   
    double c22[1266];
    double c23[1266];   
    double c24[1266];
    double c25[1266];   
    double c26[1266];
    double c27[1266];   
    double c28[1266];
    double c29[1266];   
    double c30[1266];
    double c31[1266];   
    double c32[1266];
    double c33[1266];   
    double c34[1266];
    double c35[1266];   
    double c36[1266];
    double c37[1266];   
    double c38[1266];
    double c39[1266];   
    double c40[1266];

double le1(void){
    const char lamb[] = "./lambda/lambdam1.dat"; 
    FILE *file = fopen(lamb, "r");
    int size = 1266;
    int q=0;
    double m, r, l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36, l37, l38, l39, l40;    
    while (q < size) {
    fscanf (file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &m, &r, &l0, &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8, &l9, &l10, &l11, &l12, &l13, &l14, &l15, &l16, &l17, &l18, &l19, &l20, &l21, &l22, &l23, &l24, &l25, &l26, &l27, &l28, &l29, &l30, &l31, &l32, &l33, &l34, &l35, &l36, &l37, &l38, &l39, &l40);
    a[q] = m;
    b[q] = r;
    c0[q] = l0;
    c1[q] = l1;
    c2[q] = l2;
    c3[q] = l3;
    c4[q] = l4;
    c5[q] = l5;
    c6[q] = l6;
    c7[q] = l7;
    c8[q] = l8;
    c9[q] = l9;
    c10[q] = l10;
    c11[q] = l11;
    c12[q] = l12;
    c13[q] = l13;
    c14[q] = l14;
    c15[q] = l15;
    c16[q] = l16;
    c17[q] = l17;
    c18[q] = l18;
    c19[q] = l19;
    c20[q] = l20;
    c21[q] = l21;
    c22[q] = l22;
    c23[q] = l23;
    c24[q] = l24;
    c25[q] = l25;
    c26[q] = l26;
    c27[q] = l27;
    c28[q] = l28;
    c29[q] = l29;
    c30[q] = l30;
    c31[q] = l31;
    c32[q] = l32;
    c33[q] = l33;
    c34[q] = l34;
    c35[q] = l35;
    c36[q] = l36;
    c37[q] = l37;
    c38[q] = l38;
    c39[q] = l39;
    c40[q] = l40;

 //  printf("%g %g %g \n", a[p], b[q], c[p][q]); 
    q++;
  }
    fclose(file);
   return 0;
}

    double a2[1266];
    double d0[1266];    
    double d1[1266];   
    double d2[1266];
    double d3[1266];   
    double d4[1266];
    double d5[1266];   
    double d6[1266];
    double d7[1266];   
    double d8[1266];
    double d9[1266];   
    double d10[1266];
    double d11[1266];   
    double d12[1266];
    double d13[1266];   
    double d14[1266];
    double d15[1266];   
    double d16[1266];
    double d17[1266];   
    double d18[1266];
    double d19[1266];   
    double d20[1266];
    double d21[1266];   
    double d22[1266];
    double d23[1266];   
    double d24[1266];
    double d25[1266];   
    double d26[1266];
    double d27[1266];   
    double d28[1266];
    double d29[1266];   
    double d30[1266];
    double d31[1266];   
    double d32[1266];
    double d33[1266];   
    double d34[1266];
    double d35[1266];   
    double d36[1266];
    double d37[1266];   
    double d38[1266];
    double d39[1266];   
    double d40[1266];

double le2(void){
    const char lamb[] = "./lambda/lambdam2.dat"; 
    FILE *file = fopen(lamb, "r");
    int size = 1266;
    int q=0;
    double m, l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36, l37, l38, l39, l40;    
    while (q < size) {
    fscanf (file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &m, &l0, &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8, &l9, &l10, &l11, &l12, &l13, &l14, &l15, &l16, &l17, &l18, &l19, &l20, &l21, &l22, &l23, &l24, &l25, &l26, &l27, &l28, &l29, &l30, &l31, &l32, &l33, &l34, &l35, &l36, &l37, &l38, &l39, &l40);
    a2[q] = m;
    d0[q] = l0;
    d1[q] = l1;
    d2[q] = l2;
    d3[q] = l3;
    d4[q] = l4;
    d5[q] = l5;
    d6[q] = l6;
    d7[q] = l7;
    d8[q] = l8;
    d9[q] = l9;
    d10[q] = l10;
    d11[q] = l11;
    d12[q] = l12;
    d13[q] = l13;
    d14[q] = l14;
    d15[q] = l15;
    d16[q] = l16;
    d17[q] = l17;
    d18[q] = l18;
    d19[q] = l19;
    d20[q] = l20;
    d21[q] = l21;
    d22[q] = l22;
    d23[q] = l23;
    d24[q] = l24;
    d25[q] = l25;
    d26[q] = l26;
    d27[q] = l27;
    d28[q] = l28;
    d29[q] = l29;
    d30[q] = l30;
    d31[q] = l31;
    d32[q] = l32;
    d33[q] = l33;
    d34[q] = l34;
    d35[q] = l35;
    d36[q] = l36;
    d37[q] = l37;
    d38[q] = l38;
    d39[q] = l39;
    d40[q] = l40;
 //  printf("%g %g %g \n", a[p], b[q], c[p][q]); 
    q++;
  }
    fclose(file);
   return 0;
}

   double a3[1266];
   double e0[1266];    
   double e1[1266];   
   double e2[1266];
   double e3[1266];   
   double e4[1266];
   double e5[1266];   
   double e6[1266];
   double e7[1266];   
   double e8[1266];
   double e9[1266];   
   double e10[1266];
   double e11[1266];   
   double e12[1266];
   double e13[1266];   
   double e14[1266];
   double e15[1266];   
   double e16[1266];
   double e17[1266];   
   double e18[1266];
   double e19[1266];   
   double e20[1266];
   double e21[1266];   
   double e22[1266];
   double e23[1266];   
   double e24[1266];
   double e25[1266];   
   double e26[1266];
   double e27[1266];   
   double e28[1266];
   double e29[1266];   
   double e30[1266];
   double e31[1266];   
   double e32[1266];
   double e33[1266];   
   double e34[1266];
   double e35[1266];   
   double e36[1266];
   double e37[1266];   
   double e38[1266];
   double e39[1266];   
   double e40[1266];

double le3(void){
    const char lamb[] = "./lambda/lambdam3.dat"; 
    FILE *file = fopen(lamb, "r");
    int size = 1266;
    int q=0;
    double m, l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36, l37, l38, l39, l40;    
    while (q < size) {
    fscanf (file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &m, &l0, &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8, &l9, &l10, &l11, &l12, &l13, &l14, &l15, &l16, &l17, &l18, &l19, &l20, &l21, &l22, &l23, &l24, &l25, &l26, &l27, &l28, &l29, &l30, &l31, &l32, &l33, &l34, &l35, &l36, &l37, &l38, &l39, &l40);
    a3[q] = m;
    e0[q] = l0;
    e1[q] = l1;
    e2[q] = l2;
    e3[q] = l3;
    e4[q] = l4;
    e5[q] = l5;
    e6[q] = l6;
    e7[q] = l7;
    e8[q] = l8;
    e9[q] = l9;
    e10[q] = l10;
    e11[q] = l11;
    e12[q] = l12;
    e13[q] = l13;
    e14[q] = l14;
    e15[q] = l15;
    e16[q] = l16;
    e17[q] = l17;
    e18[q] = l18;
    e19[q] = l19;
    e20[q] = l20;
    e21[q] = l21;
    e22[q] = l22;
    e23[q] = l23;
    e24[q] = l24;
    e25[q] = l25;
    e26[q] = l26;
    e27[q] = l27;
    e28[q] = l28;
    e29[q] = l29;
    e30[q] = l30;
    e31[q] = l31;
    e32[q] = l32;
    e33[q] = l33;
    e34[q] = l34;
    e35[q] = l35;
    e36[q] = l36;
    e37[q] = l37;
    e38[q] = l38;
    e39[q] = l39;
    e40[q] = l40;
 //  printf("%g %g %g \n", a[p], b[q], c[p][q]); 
    q++;
  }
    fclose(file);
   return 0;
}
