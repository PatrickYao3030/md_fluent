/********************************************************************
 This UDF is an implementation of the species sink/source near the
  wall boundary condition. 
 It's a preliminary use for membrane distillation.
********************************************************************/
/******************************************************************** 
 Allocate the appropriate number (1) of memory location(s) in the 
  User-Defined Memory dialog box in ANSYS FLUENT
********************************************************************/
/******************************************************************** 
 The codes are for running in serial only. 
********************************************************************/

#include "udf.h"
#include "mem.h"
#include "metric.h"

#include "consts.h"

// constants for defining the global array
#define MAXCELLNUM 9999
#define EPS 5.0e-4

FILE *fout0, *fout1, *fout2, *fout3, *fout4;

int UC_cell_index[MAXCELLNUM][2]; // index of corresponding thread of cells
real UC_cell_centroid[MAXCELLNUM][ND_ND][2]; // coordinate of centroid
real UC_cell_T[MAXCELLNUM][2]; // temperature
real UC_cell_WX[MAXCELLNUM][2]; // mass fraction of solute
real UC_cell_massflux[MAXCELLNUM]; // transmembrane mass flux toward the permeate
int gid = 0;

//struct CellInfo
//{
//	int seq;
//	int index;
//	real centroid[ND_ND];
//} WallCells[999], UCCell;

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

real SatConc(real t) // saturated concentration for given temperature in term of the mass fraction of NaCl
{
	real result = 0.;
	result = .27;
	return result;
}
real MassFlux(real TW0, real TW1, real WW0, real WW1)
{
	real result = 0.;
	real drv_force = 0., resistance = 0.;
	drv_force = psat_h2o(TW0)-psat_h2o(TW1);
	resistance = 1./3.0e-7;
	result = drv_force/resistance;
	return result;
}

DEFINE_INIT(idf_cells, domain)
// identify the cells, which are adjacent to both sides of the membrane
// use C_UDMI(0) to store the identifier, which marks as -1 or 1 for the adjacent cells and 0 (not set) for others 
{
	Domain *d_feed, *d_perm;
	cell_t i_cell, i_cell0, i_cell1;
	face_t i_face0, i_face1;
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	Thread *t_cell;
	real loc[ND_ND], loc0[ND_ND], loc1[ND_ND];
	int i = 0;
	real temp = 0.0;

	d_feed = Get_Domain(1);
	d_perm = Get_Domain(2);
	t_FeedFluid = Lookup_Thread(domain, 5);
	t_PermFluid = Lookup_Thread(domain, 6);
	t_FeedInterface = Lookup_Thread(domain, 13);
	t_PermInterface = Lookup_Thread(domain, 14);
	//fout0 = fopen("idf_cells0.out", "w");
	//fout1 = fopen("idf_cells1.out", "w");
	fout2 = fopen("idf_cell2.out", "w");
	fout3 = fopen("idf_cell3.out", "w");

	if(!Data_Valid_P()) 
	{
		Message("[idf_cells] Some accessing variables have not been allocated.\n");
		Message("[idf_cells] The wall cells have not been identified yet.\n");
		return;
	}

	begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	{
		i_cell0 = F_C0(i_face0, t_FeedInterface);
		C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
		C_UDMI(i_cell0, t_FeedFluid, 0) = +1; // mark the wall cells as +1, and others as 0 (no modification)
		UC_cell_index[gid][0] = i_cell0; // store the index of feed-side wall cell
		UC_cell_centroid[gid][0][0] = loc0[0];
		UC_cell_centroid[gid][1][0] = loc0[1];
		UC_cell_T[gid][0] = C_T(i_cell0, t_FeedFluid);
		UC_cell_WX[gid][0] = C_YI(i_cell0, t_FeedFluid, 0); // NOTE: this sentense is valid only if the mixture mode is used 
		begin_f_loop(i_face1, t_PermInterface) // search the symmetric cell (THE LOOP CAN ONLY RUN IN SERIAL MODE)
		{
			i_cell1 = F_C0(i_face1, t_PermInterface);
			C_CENTROID(loc1, i_cell1, t_PermFluid);
			if (fabs(loc0[0]-loc1[0])/loc0[0] < EPS) // In this special case, the pair of wall cells on both sides of membrane are symmetrical
			{
				fprintf(fout3, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
				C_UDMI(i_cell0, t_FeedFluid, 1) = i_cell1; // store the index of the found cell
				UC_cell_index[gid][1] = i_cell1; // store the index of permeate-side wall cell
				UC_cell_centroid[gid][0][1] = loc1[0];
				UC_cell_centroid[gid][1][1] = loc1[1];
				UC_cell_T[gid][1] = C_T(i_cell1, t_PermFluid);
				UC_cell_WX[gid][1] = C_YI(i_cell1, t_PermFluid, 0); // it'll lead to the error of ACCESS_VIOLATION if the domain is not a mixture
			}
		}
		end_f_loop(i_face1, t_PermInterface)
		gid++;
	}
	end_f_loop(i_face0, t_FeedInterface)

	//gid = 0; // reset the global index

	begin_f_loop(i_face1, t_PermInterface) // find the adjacent cells for the permeate-side membrane.
	{
		i_cell1 = F_C0(i_face1, t_PermInterface);
		C_CENTROID(loc1, i_cell1, t_PermFluid); // get the location of cell centroid
		begin_f_loop(i_face0, t_FeedInterface)
		{
			i_cell0 = F_C0(i_face0, t_FeedInterface);
			C_CENTROID(loc0, i_cell0, t_FeedFluid);
			if (fabs(loc1[0]-loc0[0])/loc1[0] < EPS)
			{
				fprintf(fout3, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
				C_UDMI(i_cell1, t_PermFluid, 1) = i_cell0;
			}
		}
		end_f_loop(i_face1, t_PermInterface)
		C_UDMI(i_cell1, t_PermFluid, 0) = +1; // mark the cell
	}
	end_f_loop(i_face1, t_PermInterface)
	
	for (i = 0; i<9999; i++)
		fprintf(fout2, "No.%d wall cell index %d located at %g %g with temperature of %g and mass fraction of %g, symmetric cell index %d at %g %g with temperature of %g and mass fraction of %g\n", i, UC_cell_index[i][0], UC_cell_centroid[i][0][0], UC_cell_centroid[i][1][0], UC_cell_T[i][0], UC_cell_WX[i][0], UC_cell_index[i][1], UC_cell_centroid[i][0][1], UC_cell_centroid[i][1][1], UC_cell_T[i][1], UC_cell_WX[i][1]);
	
	//temp = 353.15;
	//fprintf(fout0, "T = %g (K) saturated vapor pressure is %g (Pa)\n", temp, psat_h2o(temp));

	fclose(fout2);
	fclose(fout3);
	//fclose(fout0);
	//fclose(fout1);
}


DEFINE_ADJUST(calc_flux, domain)
// calculate the flux across the membrane
// the flux will also convert into the source of the adjacent cell, which will store in C_UDMI(1) for passing to the source term
{
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	face_t i_face0, i_face1;
	cell_t i_cell0, i_cell1;
	real loc0[ND_ND], loc1[ND_ND];
	real coeff = 3.e-7;
	real mass_flux, heat_flux; 
	int i = 0;

	//fout4 = fopen("idf_cell4.out", "w");

	t_FeedFluid = Lookup_Thread(domain, 5);
	t_PermFluid = Lookup_Thread(domain, 6);
	t_FeedInterface = Lookup_Thread(domain, 13);
	t_PermInterface = Lookup_Thread(domain, 14);

	for (i=0; i<9999; i++) // get the T and YI(0) of the wall cells
	{
		UC_cell_T[i][0] = C_T(UC_cell_index[i][0], t_FeedFluid);
		UC_cell_T[i][1] = C_T(UC_cell_index[i][1], t_PermFluid);
		//UC_cell_WX[i][0] = C_YI(UC_cell_index[i][0], t_FeedFluid, 0);
		//UC_cell_WX[i][1] = C_YI(UC_cell_index[i][1], t_PermFluid, 0);
		//fprintf(fout4, "Cell %d T = %g (K) psat(T) = %g (Pa)\n", UC_cell_index[i][0], UC_cell_T[i][0], psat_h2o(UC_cell_T[i][0]));
		//fprintf(fout4, "Feed-side wall cell %d T = %g and sat.P = %g, permeate-side wall cell %d T = %g and sat.P = %g\n", UC_cell_index[i][0], UC_cell_T[i][0], psat_h2o(UC_cell_T[i][0]), UC_cell_index[i][1], UC_cell_T[i][1], psat_h2o(UC_cell_T[i][1]));
		mass_flux = MassFlux(UC_cell_T[i][0], UC_cell_T[i][1], UC_cell_WX[i][0], UC_cell_WX[i][1]);
		UC_cell_massflux[i] = mass_flux;
		//fprintf(fout4, "No.%d membrane temperatures of feeding and permeating sides are %g and %g respectively, and its permeation flux is %g (kg/m2-s)\n", i, UC_cell_T[i][0], UC_cell_T[i][1],  mass_flux);
		C_UDMI(UC_cell_index[i][0], t_FeedFluid, 2) = -mass_flux; // store the permeation flux in the UDMI(2)
		C_UDMI(UC_cell_index[i][1], t_PermFluid, 2) = +mass_flux;
		if ((UC_cell_index[i][0] == 0) & (UC_cell_index[i][1] == 0)) return;
	}
	//fclose(fout4);
}

DEFINE_SOURCE(mass_source, i_cell, t_cell, dS, eqn)
{
	real conc, temp;
	real source; // returning result

	conc = 1.-C_YI(i_cell, t_cell, 0); // the mass fraction of NaCl
	temp = C_T(i_cell, t_cell);
	if (conc > SatConc(temp)) // calculation for the solution under saturation
	{
		Message("[mass_source] The solution has the saturated concentration.\n");
		source = 0.;
		return source;
	}
	else
	{
		source = C_UDMI(i_cell, t_cell, 0)*C_UDMI(i_cell, t_cell, 2); // mass source of the cell relates to the ratio of permeation flux and cell's height (0.5mm)
		source = C_UDMI(i_cell, t_cell, 0)*(-1.);
	}
  dS[eqn] = 0.;

  return source;
}
