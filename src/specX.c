/********************************************************************
 This UDF is an implementation of the species sink/source near the
  wall boundary condition. 
 It's a preliminary use for membrane distillation.
********************************************************************/
/******************************************************************** 
 Preliminary checklist before use:
 (1) Run the ANSYS FLUENT in 2-d serial mode
 (2) Allocate the three User-Defined Memories
 (3) Set the mixture for both cell zone conditions
********************************************************************/

#include "udf.h"
#include "mem.h"
#include "metric.h"

#include "consts.h"

FILE *fout0, *fout1, *fout2, *fout3, *fout4;

int gid = 0;
struct PorousMaterials membrane;
struct CellInfos WallCell[MAXCELLNUM][2];

void GetProp_Membrane(real temperature) // Get the properties of the membrane for the given temperature
{
	membrane.thickness = 1.5e-4;
	membrane.porosity = 0.7;
	membrane.tortuosity = 1.2;
	membrane.conductivity = ThermCond_Maxwell(temperature, membrane.porosity, PVDF);
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
	result = heat_flux_0; // (W/m2), consider the latent heat flux only
	return result;
}

real OverheatCheck(real q, real m, real cp, real t0, real t1, real tref) // if the overheat happens, it will return a revised heat flow (either being exothermal or endothermal) 
{
	real result = 0.;
	t1 = q/(m*cp)+t0;  
	if (q*(t1-tref)>0.) // with the absorbed heat (q>0), the calculated temperature should be lower than the refered one; with the released heat (q<0), t1 > tref
	{
		result = m*cp*(tref-t0);
	}
	else
	{
		result = m*cp*(t1-t0);
	}
}

DEFINE_INIT(idf_cells, domain)
/* 
   [objectives] identify the cell pairs, which are adjacent to both sides of the membrane
   [methods] 1. get a cell beside the feeding membrane boundary         
             2. find the corresponding cells with the same x-coordinate
   [outputs] 1. for all cells, the cell whose C_UDMI(0) = -1 means it belongs to the wall cell of feeding membrane
                                                          +1 means it belongs to the wall cell of permeating membrane
			 2. internal variables for recording the identified pairs of wall cells
*/
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

	t_FeedFluid = Lookup_Thread(domain, 32);
	t_PermFluid = Lookup_Thread(domain, 33);
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
		WallCell[gid][0].index = i_cell0;
		WallCell[gid][0].centroid[0] = loc0[0];
		WallCell[gid][0].centroid[1] = loc0[1];
		WallCell[gid][0].temperature = C_T(i_cell0, t_FeedFluid);
		WallCell[gid][0].massfraction.water = C_YI(i_cell0, t_FeedFluid, 0);
		begin_f_loop(i_face1, t_PermInterface) // search the symmetric cell (THE LOOP CAN ONLY RUN IN SERIAL MODE)
		{
			i_cell1 = F_C0(i_face1, t_PermInterface);
			C_CENTROID(loc1, i_cell1, t_PermFluid);
			if (fabs(loc0[0]-loc1[0])/loc0[0] < EPS) // In this special case, the pair of wall cells on both sides of membrane are symmetrical
			{
				fprintf(fout3, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
				//C_UDMI(i_cell0, t_FeedFluid, 1) = i_cell1; // store the index of the found cell
				WallCell[gid][1].index = i_cell1;
				WallCell[gid][1].centroid[0] = loc1[0];
				WallCell[gid][1].centroid[1] = loc1[1];
				WallCell[gid][1].temperature = C_T(i_cell1, t_PermFluid);
				WallCell[gid][1].massfraction.water = C_YI(i_cell1, t_PermFluid, 0);
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
				//C_UDMI(i_cell1, t_PermFluid, 1) = i_cell0;
			}
		}
		end_f_loop(i_face1, t_PermInterface)
		C_UDMI(i_cell1, t_PermFluid, 0) = +1; // mark the cell
	}
	end_f_loop(i_face1, t_PermInterface)
	
	for (i = 0; i<9999; i++)
		fprintf(fout2, "No.%d wall cell index %d located at %g %g with temperature of %g and mass fraction of %g, symmetric cell index %d at %g %g with temperature of %g and mass fraction of %g\n", i, WallCell[i][0].index, WallCell[i][0].centroid[0], WallCell[i][0].centroid[1], WallCell[i][0].temperature, WallCell[i][0].massfraction.water, WallCell[i][1].index, WallCell[i][1].centroid[0], WallCell[i][1].centroid[1], WallCell[i][1].temperature, WallCell[i][1].massfraction.water);
	
	fclose(fout2);
	fclose(fout3);
	//fclose(fout0);
	//fclose(fout1);
}

DEFINE_ON_DEMAND(testGetDomain)
{
	cell_t i_cell, i_cell0, i_cell1;
	face_t i_face0, i_face1;
	Thread *t_cell;
	real loc[ND_ND], loc0[ND_ND], loc1[ND_ND];
	int i = 0;
	real temp = 0.0;

	//Domain *d_feed = Get_Domain(32);
	//Domain *d_perm = Get_Domain(33);
	Domain *domain = Get_Domain(1);
	Thread *t_FeedFluid = Lookup_Thread(domain, 32);
	Thread *t_PermFluid = Lookup_Thread(domain, 33);
	Thread *t_FeedInterface = Lookup_Thread(domain, 30);
	Thread *t_PermInterface = Lookup_Thread(domain, 2);

	fout2 = fopen("idf_cell2.out", "w");
	fout3 = fopen("idf_cell3.out", "w");

	//begin_c_loop(i_cell, t_FeedFluid){
	//	Message("cell#%d\n", i_cell);
	//}
	//end_c_loop(i_cell, t_FeedFluid)

	begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	{
		i_cell0 = F_C0(i_face0, t_FeedInterface);
		C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
		Message("interface#%d's adjacent cell index is #%d, located at (%g,%g)\n", i_face0, i_cell0, loc0[0], loc0[1]);
		C_UDMI(i_cell0, t_FeedFluid, 0) = -1; // mark the wall cells as -1, and others as 0 (no modification)
		//UC_cell_index[gid][0] = i_cell0; // store the index of feed-side wall cell
		//UC_cell_centroid[gid][0][0] = loc0[0];
		//UC_cell_centroid[gid][1][0] = loc0[1];
		//UC_cell_T[gid][0] = C_T(i_cell0, t_FeedFluid);
		//UC_cell_WX[gid][0] = C_YI(i_cell0, t_FeedFluid, 0); // NOTE: this sentense is valid only if the mixture mode is used 
		//begin_f_loop(i_face1, t_PermInterface) // search the symmetric cell (THE LOOP CAN ONLY RUN IN SERIAL MODE)
		//{
		//	i_cell1 = F_C0(i_face1, t_PermInterface);
		//	C_CENTROID(loc1, i_cell1, t_PermFluid);
		//	if (fabs(loc0[0]-loc1[0])/loc0[0] < EPS) // In this special case, the pair of wall cells on both sides of membrane are symmetrical
		//	{
		//		fprintf(fout3, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
		//		C_UDMI(i_cell0, t_FeedFluid, 1) = i_cell1; // store the index of the found cell
		//		UC_cell_index[gid][1] = i_cell1; // store the index of permeate-side wall cell
		//		UC_cell_centroid[gid][0][1] = loc1[0];
		//		UC_cell_centroid[gid][1][1] = loc1[1];
		//		UC_cell_T[gid][1] = C_T(i_cell1, t_PermFluid);
		//		UC_cell_WX[gid][1] = C_YI(i_cell1, t_PermFluid, 0); // it'll lead to the error of ACCESS_VIOLATION if the domain is not a mixture
		//	}
		//}
		//end_f_loop(i_face1, t_PermInterface)
		gid++;
	}
	end_f_loop(i_face0, t_FeedInterface)

	//extern real ThermCond_Maxwell();
	//extern real psat_h2o();
	//real km = 0., tm = 333.15, porosty = .7;
	//km = ThermCond_Maxwell(tm, porosty, 1);
	//Message("\nThe membrane thermal conductivity is %g (W/m-K).\n", km);
	//Message("\nThe saturated vapor pressure is %g (Pa) for given temperature of %g (K).", psat_h2o(tm), tm);

}

DEFINE_ON_DEMAND(testGetVolume)
{
}

DEFINE_ADJUST(calc_flux, domain)
/*
	[objectives] calculate the flux across the membrane
	[methods] 1. get the temperatures and concentrations of the identified pair of wall cells
	          2. calculate the permeation flux according to the given temperature and concentration
			  3. calculate the permeative heat flux, here only latent heats are considered while the conjugated conductive heat transfer scheme is used
			  4. if the cell is overheated with the heat flux input, reset the permeation flux and go back to step 2
	[outputs] 1. C_UDMI(1) for mass flux
	          2. C_UDMI(2) for latent heat flux
*/
{
	extern real SatConc();
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	face_t i_face0, i_face1;
	cell_t i_cell0, i_cell1;
	real loc0[ND_ND], loc1[ND_ND];
	real mass_flux, heat_flux; 
	int i = 0;

	//fout4 = fopen("idf_cell4.out", "w");

	t_FeedFluid = Lookup_Thread(domain, 32);
	t_PermFluid = Lookup_Thread(domain, 33);
	t_FeedInterface = Lookup_Thread(domain, 30);
	t_PermInterface = Lookup_Thread(domain, 2);

	for (i=0; i<9999; i++) // get the T and YI(0) of the wall cells
	{
		WallCell[i][0].temperature = C_T(WallCell[i][0].index, t_FeedFluid);
		WallCell[i][1].temperature = C_T(WallCell[i][1].index, t_PermFluid);
		WallCell[i][0].massfraction.water = C_YI(WallCell[i][0].index, t_FeedFluid, 0);
		WallCell[i][1].massfraction.water = C_YI(WallCell[i][1].index, t_PermFluid, 0);

		if (WallCell[i][0].massfraction.water > (1.-SatConc(WallCell[i][0].temperature))) // calculate the mass tranfer across the membrane only if the concentration is below the saturation
		{
			mass_flux = MassFlux(WallCell[i][0].temperature, WallCell[i][1].temperature, WallCell[i][0].massfraction.water, WallCell[i][1].massfraction.water);
		}
		else
		{
			mass_flux = .0;
		}
		heat_flux = HeatFlux(WallCell[i][0].temperature, WallCell[i][1].temperature, mass_flux); // calculate the heat transfer across the membrane
		mass_flux = OverheatCheck(); 
		C_UDMI(WallCell[i][0].index, t_FeedFluid, 1) = -mass_flux; // store the permeation flux in the UDMI(1)
		C_UDMI(WallCell[i][1].index, t_PermFluid, 1) = +mass_flux;
		C_UDMI(WallCell[i][0].index, t_FeedFluid, 2) = -heat_flux; // store the heat flux in the UDMI(2)
		C_UDMI(WallCell[i][1].index, t_PermFluid, 2) = +heat_flux;
		if ((WallCell[i][1].index == 0) & (WallCell[i][1].index == 0)) return;
	}
	//fclose(fout4);
}

DEFINE_SOURCE(mass_source, i_cell, t_cell, dS, eqn)
/*
	[objectives] add the term of mass source for wall cells
	[methods] 1. convert the transmembrane mass transfer into the source term of the wall cell
	          2. set the mass source for each calculating cell
	[outputs] the mass source or sink
*/
{
	real source; // returning result
	source = fabs(C_UDMI(i_cell, t_cell, 0))*C_UDMI(i_cell, t_cell, 1)/0.5e-3; // mass source of the cell relates to the ratio of permeation flux and cell's height (0.5mm)
  dS[eqn] = 0.;
  return source;
}

DEFINE_SOURCE(heat_source, i_cell, t_cell, dS, eqn)
/*
	[objectives] add the term of latent heat source for wall cells
	[methods] 1. convert the latent heat flux into the source term of the wall cell
	          2. calculate the evaporation and condensation heat of the given cell
	[outputs] the latent heat source or sink
*/
{
	real source; // returning result
	source = fabs(C_UDMI(i_cell, t_cell, 0))*C_UDMI(i_cell, t_cell, 2)/0.5e-3; // heat source of the cell relates to the ratio of heat flux and cell's height (0.5mm)
  dS[eqn] = 0.;
  return source;
}