#include "udf.h"
 
/*Constants used in psat_h2o to calculate saturation pressure*/

#define PSAT_A 0.01
#define PSAT_TP 338.15
#define C_LOOP 8
#define H2O_PC 22.089E6
#define H2O_TC 647.286

double psat_h2o(double tsat)
/*                */
/* Computes saturation pressure of water vapor     */
/* as function of temperature            */
/* Equation is taken from THERMODYNAMIC PROPERTIES IN SI, */
/* by Reynolds, 1979                  */
/* Returns pressure in PASCALS, given temperature in KELVIN */
{
 int i;
 double var1,sum1,ans1,psat;
 double constants[8]={-7.4192420, 2.97221E-1, -1.155286E-1,
    8.68563E-3, 1.094098E-3, -4.39993E-3, 2.520658E-3, -5.218684E-4}; 

 /* var1 is an expression that is used in the summation loop */
 var1 = PSAT_A*(tsat-PSAT_TP);

 /* Compute summation loop */
 i = 0;
 sum1 = 0.0;
 while (i < C_LOOP){
     sum1+=constants[i]*pow(var1,i);
     ++i;
 }
 ans1 = sum1*(H2O_TC/tsat-1.0);

 /* compute exponential to determine result */
 /* psat has units of Pascals     */

 psat = H2O_PC*exp(ans1);
 return psat;
}