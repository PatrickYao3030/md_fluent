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

FILE *fout0, *fout1, *fout2, *fout3, *fout4;

int UC_cell_index[MAXCELLNUM][2]; // index of corresponding thread of cells
real UC_cell_centroid[MAXCELLNUM][ND_ND][2]; // coordinate of centroid
real UC_cell_T[MAXCELLNUM][2]; // temperature
real UC_cell_WX[MAXCELLNUM][2]; // mass fraction of solute
real UC_cell_massflux[MAXCELLNUM]; // transmembrane mass flux toward the permeate
int gid = 0;
struct PorousMaterials membrane;
struct CellInfos WallCell[MAXCELLNUM][2];

void GetProp_Membrane(real init_temp) // Get the properties of the membrane for the given initial temperature
{
	membrane.thickness = 1.5e-6;
	membrane.porosity = 0.7;
	membrane.tortuosity = 1.2;
	membrane.conductivity = ThermCond_Maxwell(init_temp, membrane.porosity, PVDF);
	membrane.MDcoeff = 3.6e-7;
}

real MassFlux(real tw0, real tw1, real ww0, real ww1)
{
	extern real psat_h2o();
	real result = 0.;
	real drv_force = 0., resistance = 0.;
	real avg_temp = 0.;
	avg_temp = .5*(tw0+tw1);
	GetProp_Membrane(avg_temp);
	drv_force = psat_h2o(tw0)-psat_h2o(tw1);
	resistance = 1./membrane.MDcoeff;
	result = drv_force/resistance; // use SI unit (kg/m2-s)
	return result;
}

real HeatFlux(real tw0, real tw1, real mass_flux) // if tw0 > tw1, the mass_flux should be positive, then the output should also be positive, otherwise the negative result should be returned
{
	extern real LatentHeat();
	real avg_temp = 0., diff_temp = 0.;
	real latent_heat = 0.;
	real heat_flux_0 = 0., heat_flux_1 = 0.;
	real result = 0.;
	avg_temp = .5*(tw0+tw1);
	diff_temp = tw0-tw1;
	GetProp_Membrane(avg_temp);
	latent_heat = LatentHeat(avg_temp); // in the unit of (J/kg)
	heat_flux_0 = latent_heat*mass_flux; // latent heat flux
	heat_flux_1 = membrane.conductivity/membrane.thickness*diff_temp; // conductive heat flux
	result = heat_flux_0; // (W/m2) consider the latent heat flux only
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

	//d_feed = Get_Domain(32);
	//d_perm = Get_Domain(33);
	t_FeedFluid = Lookup_Thread(domain, 34);
	t_PermFluid = Lookup_Thread(domain, 35);
	t_FeedInterface = Lookup_Thread(domain, 30);
	t_PermInterface = Lookup_Thread(domain, 2);
	//fout0 = fopen("idf_cells0.out", "w");
	//fout1 = fopen("idf_cells1.out", "w");
	fout2 = fopen("idf_cell2.out", "w");
	fout3 = fopen("idf_cell3.out", "w");

	//if(!Data_Valid_P()) 
	//{
	//	Message("\n[idf_cells] Some accessing variables have not been allocated.\n");
	//	Message("[idf_cells] The wall cells have not been identified yet.\n");
	//	return;
	//}

	begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	{
		i_cell0 = F_C0(i_face0, t_FeedInterface);
		C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
		C_UDMI(i_cell0, t_FeedFluid, 0) = -1; // mark the wall cells as -1, and others as 0 (no modification)
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

DEFINE_ON_DEMAND(testGetDomain)
{
	Domain *d_feed = Get_Domain(32);
	Domain *d_perm = Get_Domain(33);
	Domain *domain = Get_Domain(1);
	Thread *t_FeedFluid = Lookup_Thread(domain, 34);
	Thread *t_PermFluid = Lookup_Thread(domain, 35);
	Thread *t_FeedInterface = Lookup_Thread(domain, 30);
	Thread *t_PermInterface = Lookup_Thread(domain, 2);
	//extern real ThermCond_Maxwell();
	//extern real psat_h2o();
	//real km = 0., tm = 333.15, porosty = .7;
	//km = ThermCond_Maxwell(tm, porosty, 1);
	//Message("\nThe membrane thermal conductivity is %g (W/m-K).\n", km);
	//Message("\nThe saturated vapor pressure is %g (Pa) for given temperature of %g (K).", psat_h2o(tm), tm);

}

DEFINE_ADJUST(calc_flux, domain)
// calculate the flux across the membrane
// the flux will also convert into the source of the adjacent cell, which will store in C_UDMI(1) for passing to the source term
{
	extern real SatConc();
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	face_t i_face0, i_face1;
	cell_t i_cell0, i_cell1;
	real loc0[ND_ND], loc1[ND_ND];
	//real coeff = 3.e-7;
	real mass_flux, heat_flux; 
	int i = 0;

	fout4 = fopen("idf_cell4.out", "w");

	t_FeedFluid = Lookup_Thread(domain, 34);
	t_PermFluid = Lookup_Thread(domain, 35);
	t_FeedInterface = Lookup_Thread(domain, 30);
	t_PermInterface = Lookup_Thread(domain, 2);

	for (i=0; i<9999; i++) // get the T and YI(0) of the wall cells
	{
		UC_cell_T[i][0] = C_T(UC_cell_index[i][0], t_FeedFluid);
		UC_cell_T[i][1] = C_T(UC_cell_index[i][1], t_PermFluid);
		UC_cell_WX[i][0] = C_YI(UC_cell_index[i][0], t_FeedFluid, 0);
		UC_cell_WX[i][1] = C_YI(UC_cell_index[i][1], t_PermFluid, 0);
		//fprintf(fout4, "Cell %d T = %g (K) psat(T) = %g (Pa)\n", UC_cell_index[i][0], UC_cell_T[i][0], psat_h2o(UC_cell_T[i][0]));
		//fprintf(fout4, "Feed-side wall cell %d T = %g and sat.P = %g, permeate-side wall cell %d T = %g and sat.P = %g\n", UC_cell_index[i][0], UC_cell_T[i][0], psat_h2o(UC_cell_T[i][0]), UC_cell_index[i][1], UC_cell_T[i][1], psat_h2o(UC_cell_T[i][1]));
		if (UC_cell_WX[i][0] > (1.-SatConc(UC_cell_T[i][0]))) // calculate the mass tranfer across the membrane only if the concentration is below the saturation
		{
			mass_flux = MassFlux(UC_cell_T[i][0], UC_cell_T[i][1], UC_cell_WX[i][0], UC_cell_WX[i][1]);
		}
		else
		{
			mass_flux = .0;
		}
		UC_cell_massflux[i] = mass_flux;
		heat_flux = HeatFlux(UC_cell_T[i][0], UC_cell_T[i][1], UC_cell_massflux[i]); // calculate the heat transfer across the membrane
		fprintf(fout4, "No.%d membrane temperatures of feeding and permeating sides are %g and %g respectively, and its permeation flux and heat flux are %g (kg/m2-s) and %g (J/m2-s) respectively.\n", i, UC_cell_T[i][0], UC_cell_T[i][1],  mass_flux);
		C_UDMI(UC_cell_index[i][0], t_FeedFluid, 1) = -heat_flux; // store the heat flux in the UDMI(1)
		C_UDMI(UC_cell_index[i][1], t_PermFluid, 1) = +heat_flux;
		C_UDMI(UC_cell_index[i][0], t_FeedFluid, 2) = -mass_flux; // store the permeation flux in the UDMI(2)
		C_UDMI(UC_cell_index[i][1], t_PermFluid, 2) = +mass_flux;
		if ((UC_cell_index[i][0] == 0) & (UC_cell_index[i][1] == 0)) return;
	}
	fclose(fout4);
}

DEFINE_SOURCE(mass_source, i_cell, t_cell, dS, eqn)
{
	real source; // returning result
	source = C_UDMI(i_cell, t_cell, 0)*C_UDMI(i_cell, t_cell, 2)/0.5e-3; // mass source of the cell relates to the ratio of permeation flux and cell's height (0.5mm)
  dS[eqn] = 0.;
  return source;
}

DEFINE_SOURCE(heat_source, i_cell, t_cell, dS, eqn)
{
	real source; // returning result
	source = fabs(C_UDMI(i_cell, t_cell, 0))*C_UDMI(i_cell, t_cell, 1)/0.5e-3; // heat source of the cell relates to the ratio of heat flux and cell's height (0.5mm)
	source = fabs(C_UDMI(i_cell, t_cell, 0))*C_UDMI(i_cell, t_cell, 1);
  dS[eqn] = 0.;
  return source;
}