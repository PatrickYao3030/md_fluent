#include <math.h>
#include "udf.h"
#include "consts.h"
/*Constants used in psat_h2o to calculate saturation pressure*/
#define PSAT_A 0.01
#define PSAT_TP 338.15
#define C_LOOP 8
#define H2O_PC 22.089E6
#define H2O_TC 647.286

real psat_h2o(real tsat)
/* 
	 Computes saturation pressure of water vapor
	 as function of temperature            
	 Equation is taken from THERMODYNAMIC PROPERTIES IN SI, 
	 by Reynolds, 1979 
	 Returns pressure in PASCALS, given temperature in KELVIN
*/
{
 int i;
 real var1,sum1,ans1,psat;
 real constants[8]={-7.4192420, 2.97221E-1, -1.155286E-1,
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

real ThermCond_Maxwell(real temp, real porosity, int opt) 
/* Calculate the thermal conductivity of the membrane for its porosity higher than 60%
	 with the correlation of Garcia-Payo and lzquierdo-Gil [J Phys D 2004, 37(21): 3008-3016]
	 which is also commented by Hitsov [Sep Purification Tech 2015, 142: 48-64]
	 all the following reference numbers are in the review of Hitsov
*/
{
	real result;
	real kappa[2], beta;
	real const A[4] = {5.769, 5.769, 12.5, 4.167};
	real const B[4] = {0.9144, 8.914, -23.51, 1.452};
	int i = 0;
	kappa[GAS] = 1.5e-3*sqrt(temp); // the thermal conductivity of trapped air and steam by Jonsson [30]
	kappa[GAS] = 2.72e-3+7.77e-5*temp; // correlated by Bahmanyar [36]
	switch(opt){
		case 0: i = PVDF; break;
		case 1: i = PTFE; break;
		case 2: i = PP; break;
		case 3: i = PES; break;
		default: Message("\n Function: ThermCond_Maxwell() has a wrong input argument of opt \n");
	}
	kappa[SOLID] = A[i]*1.e-4*temp+B[i]*1.e-2;
	beta = (kappa[SOLID]-kappa[GAS])/(kappa[SOLID]+2.*kappa[GAS]);
	result = kappa[GAS]*(1.+2.*beta*(1.-porosity))/(1.-beta*(1.-porosity));
	return result;
}
real SatConc(real t) // saturated concentration for given temperature in term of the mass fraction of NaCl
{
	real result = 0.;
	result = .27;
	return result;
}

real LatentHeat(real t) // latent heat in water evaporation/condensation for given temperature (K) in 1 atm, Drioli, E.and M. Romano [IECR 40(5): 1277-1300]
{
	real result = 0.;
	result = 1.e+3*(1.7535*t+2024.3); // use SI unit (J/kg)
	return result;
}
real ThermCond_aq(real t,real c) // thermal conductivity for given temperature (K) and concentration (w%)
{
	real result=0.;
	result=(0.608+7.46e-4*(t-273.15))*(1.-0.98*(18.*c/(58.5-40.5*c)));
	return result;
}
/*
[Objectives] define the propertis in terms of temperature or mass fraction 
[methods]  1. obtain the temperature and mass fraction of the cell
		   2. calculate the properties with the given temperature/mass fraction and property function
[outputs]  the properties of materials
*/
//DEFINE_PROPERTY(ThermCond_aq0,c,t)//shuaitao
//{
//	real result;
//	real temp_tca = C_T(c,t);
//	real conc_tca = C_YI(c,t,1);
//	result = (0.608+(7.46e-4)*(temp_tca-273.15))*(1-0.98*(18*conc_tca/(58.5-40.5*conc_tca)));
//	return result;
//}

real ConvertX(int imat, int nmat, real MW[], real wi[])
/*
	[objs] convert the mass fraction of component imat into the molar fraction
	[meth] allocate a dynamic array for temporary storage
	       calculate the mole for each component
				 calculate the molar fraction for the specified index in the array
	[outs] molar fraction
*/
{
	real *r;
	real sum_N = 0., xi = 0.;
	int i;
	r = (real*)malloc(nmat*sizeof(real));
	for (i=0; i<nmat; i++) 
	{
		r[i] = wi[i]/MW[i];
		sum_N = sum_N+r[i];
	}
	xi = r[imat]/sum_N;
	free(r);
	return xi;
}

real ActivityCoefficient_h2o(x_nv)
/*
	[objs] correlate the activity coefficient in an aqueous sodium chloride solution
	[meth] use the correlation proposed by Lawson and Lloyd
	[outs] dimensionless activity coefficent
*/
{
	real alpha_h2o;
	alpha_h2o = 1.-0.5*x_nv-10.*pow(x_nv, 2.);
	return alpha_h2o;
}

real WaterVaporPressure_brine(real temperature, real mass_fraction_h2o)
/*
	[Objectives] calculate the vapor pressure for the specified component
	[methods] 1. convert the input mass fraction of water into the molar fraction of nonvolatile components
	          2. calculate the activity coefficient according to Lawson and Lloyd's correlation
						3. get the water vapor pressure by invoking psat_h2o
	[outputs] vapor pressure in SI (Pa)
*/
{
	real x_nv, alpha, vp, MW[2] = {18.01534, 58.4428}, wi[2];
	wi[1] = mass_fraction_h2o;
	wi[2] = 1.-mass_fraction_h2o;
	x_nv = 1.-ConvertX(0, 2, MW, wi);
	alpha = ActivityCoefficient_h2o(x_nv);
	vp = (1.-x_nv)*alpha*psat_h2o(temperature);
	return vp;
}

real Density_aqNaCl(real T, real w)
/*
	[objs] correlate the density of aqueous NaCl solution
	[meth] empirical correlation by B.S. Sparrow, Desalination 2003, 159(2): 161-170
	[outs] density in SI (kg/m3)
*/
{
	int i, j;
	real sum = 0.;
	real A[5][5] = {{1.001, 0.7666, -0.0149, 0.2663, 0.8845},
	                {-0.0214, -3.496, 10.02, -6.56, -31.37},
	                {-5.263, 39.87, 176.2, 363.5, -7.784},
	                {15.42,-167.0, 980.7, -2573.0, 876.6},
									{-0.0276, 0.2978, -2.017, 6.345, -3.914}};
	real B[5] = {0., 0., 0., 0., 0.};
	real C[5] = {1.e3, 1.e0, 1.e-3, 1.e-6, 1.e-6};
	if ((T<0.)||(T>300.) return sum;
	for (i=0; i<5; i++)
	{
		for (j=0; j<5; j++)
		{
			B[i] += A[i][j]*pow(w, j);
		}
		sum += B[i]*C[i]*pow(T, i);
	}
	return sum;
}

real Viscosity_aqNaCl(real T, real w_nv)
/*
	[objs] calculate the viscosity of aq NaCl solution for given T [C] and mass fraction of NaCl [?]
	[meth] ref. ??? from LST's codes
	[outs] viscosity [? SI]
*/
{
	real mu;
	mu = (17.02821-0.39206*T+0.188912*(conc_v*100)-0.00466*T*(conc_v*100)+0.003025*T*T+0.011738*(conc_v*100)*(conc_v*100))*0.001;
	return mu;
}