#include "udf.h"
///*Constants used in psat_h2o to calculate saturation pressure*/
//#define PSAT_A 0.01
//#define PSAT_TP 338.15
//#define C_LOOP 8
//#define H2O_PC 22.089E6
//#define H2O_TC 647.286
//
//double psat_h2o(double tsat)
///*                */
///* Computes saturation pressure of water vapor     */
///* as function of temperature            */
///* Equation is taken from THERMODYNAMIC PROPERTIES IN SI, */
///* by Reynolds, 1979                  */
///* Returns pressure in PASCALS, given temperature in KELVIN */
//{
// int i;
// double var1,sum1,ans1,psat;
// double constants[8]={-7.4192420, 2.97221E-1, -1.155286E-1,
//    8.68563E-3, 1.094098E-3, -4.39993E-3, 2.520658E-3, -5.218684E-4}; 
//
// /* var1 is an expression that is used in the summation loop */
// var1 = PSAT_A*(tsat-PSAT_TP);
//
// /* Compute summation loop */
// i = 0;
// sum1 = 0.0;
// while (i < C_LOOP){
//     sum1+=constants[i]*pow(var1,i);
//     ++i;
// }
// ans1 = sum1*(H2O_TC/tsat-1.0);
//
// /* compute exponential to determine result */
// /* psat has units of Pascals     */
//
// psat = H2O_PC*exp(ans1);
// return psat;
//}

/* Calculate the thermal conductivity of the membrane for its porosity higher than 60%
	 with the correlation of Garcia-Payo and lzquierdo-Gil [J Phys D 2004, 37(21): 3008-3016]
	 which is also commented by Hitsov [Sep Purification Tech 2015, 142: 48-64]
	 all the following reference numbers are in the review of Hitsov
*/
real ThermCond_Maxwell(real temp, real porosity, int opt) 
{
	real result;
	real kappa[2], beta;
	real const A[4] = {5.769, 5.769, 12.5, 4.167};
	real const B[4] = {0.9144, 8.914, -23.51, 1.452};
	int const solid = 0, gas = 1;
	int const PVDF = 0, PTFE = 1, PP = 2, PES = 3;
	int i = 0;
	kappa[solid] = 1.5e-3*fsqrt(temp); // the thermal conductivity of trapped air and steam by Jonsson [30]
	kappa[solid] = 2.72e-3+7.77e-5*temp; // correlated by Bahmanyar [36]
	switch(opt)
		case PVDF: i = PVDF; break;
		case PTFE: i = PTFE; break;
		case PP: i = PP; break;
		case PES: i = PES; break;
		default: Message("\n Function: ThermCond_Maxwell() has a wrong input argument of opt \n")
	kappa[gas] = A[i]*1.e-4*temp+B[i]*1.e-2;
	beta = (kappa[solid]-kappa[gas])/(kappa[solid]+2.*kappa[gas]);
	result = kappa[gas]*(1.+2.*beta*(1.-porosity))/(1.-beta*(1.-porosity));
	return result;
}