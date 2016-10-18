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

int gid = 0, rid = 0;
struct PorousMaterials membrane;
struct CellInfos WallCell[MAXCELLNUM][2];
struct MessageInfos CellPairInfo[MAXRECLINE];

extern real SatConc(), ThermCond_Maxwell(), WaterVaporPressure_brine(), LatentHeat(), Density_aqNaCl(), Viscosity_aqNaCl(), ThermCond_aq();

void GetProp_Membrane(real temperature) // Get the properties of the membrane for the given temperature
{
	membrane.thickness = 1.5e-4;
	membrane.porosity = 0.7;
	membrane.tortuosity = 1.2;
	membrane.conductivity = ThermCond_Maxwell(temperature, membrane.porosity, PVDF);
	membrane.MDcoeff = 9.5e-7;
}

void Monitor_CellPair(int opt, int rec_idx, int idx_cells)
/*
	[objectives] track the specified pair of cells adhered on the both sides of the membrane
	[methods] 1. put the cells' properties into the allocated storage space
	[outputs] 1. write the cell info into a data file
*/
{
	char str_line_buffer[79];
	//Message("Running Monitor_CellPair() for opt = %d, record index of %d\n", opt, rec_idx);
	if (opt == 0) return;
	if ((rec_idx >= 0) && (rec_idx <= MAXRECLINE))
	{
		CellPairInfo[rec_idx].flag = opt;
		sprintf(str_line_buffer, "T0 = %g, T1 = %g, JMLOC0 = %g, JH0 = %g", WallCell[idx_cells][0].temperature, WallCell[idx_cells][1].temperature, WallCell[idx_cells][0].flux.mass, WallCell[idx_cells][0].flux.heat);
		//Message("TF = %g, TP = %g", WallCell[idx_cells][0].temperature, WallCell[idx_cells][1].temperature);
		strcpy(CellPairInfo[rec_idx].content, str_line_buffer);
	}
	return;
}

real LocalMassFlux(real tw0, real tw1, real ww0, real ww1) // if tw0 > tw1, the mass_flux should be positive, or vice versa
{
	real JM = 0.;
	real drv_force = 0., resistance = 0.;
	real avg_temp = 0.;
	avg_temp = .5*(tw0+tw1);
	GetProp_Membrane(avg_temp);
	drv_force = WaterVaporPressure_brine(tw0, ww0)-WaterVaporPressure_brine(tw1, ww1);
	resistance = 1./membrane.MDcoeff;
	JM = drv_force/resistance; // use SI unit (kg/m2-s)
	return JM;
}

real LocalHeatFlux(int opt, real tw0, real tw1, real mass_flux) // if tw0 > tw1, the mass_flux should be positive, then the output should also be positive, otherwise the negative result should be returned
/*
	[objectives] calculate the heat transfer accompanying with the permeation
	[methods] 1. calculate the averaged membrane temperature
	          2. calculate the membrane properties and fill the workspace of membrane
						3. calculate the heat flux
	[outputs] 1. workspace of membrane with filled properties
	          2. latent heat flux, conductive heat flux and total heat flux (J/m2-s)
*/
{
	real avg_temp = 0., diff_temp = 0.;
	real latent_heat = 0.;
	real heat_flux_0 = 0., heat_flux_1 = 0.;
	real JH = 0.;
	if (((tw0>tw1) && (mass_flux>0)) || ((tw0<tw1) && (mass_flux<0))) // the sign of mass_flux should be same as tw0-tw1
	{
		avg_temp = .5*(tw0+tw1);
		diff_temp = tw0-tw1;
		GetProp_Membrane(avg_temp);
		latent_heat = LatentHeat(avg_temp); // in the unit of (J/kg)
		heat_flux_0 = latent_heat*mass_flux; // latent heat flux
		heat_flux_1 = membrane.conductivity/membrane.thickness*diff_temp; // conductive heat flux
		switch(opt)
		{
		case 0:
			JH = heat_flux_0; // (W/m2), consider the latent heat flux only
			break;
		case 1:
			JH = heat_flux_1; // (W/m2), consider the conductive heat flux only
			break;
		case 10:
			JH = heat_flux_0+heat_flux_1;
			break;
		default:
			JH = 0.;
			break;
		}
	}
	else if (tw0!=tw1)
	{
		Message("[ERROR] Mass flux has a wrong direction. ERR%g/%g/%g \n", mass_flux, tw0, tw1);
	}
	return JH;
}

void MembraneTransfer(int opt)
/*
	[objectives] calculate the md heat and mass transfer across the membrane
	[methods] 1. get the flow field adhered on the both sides of the membrane and store them into the constructed workspace "WallCell"
	          2. calculate the permeation flux according to the temperature differences
						3. calculate the heat flux according to the mass flux
						4. if the exerted heat causes the overcool/heat, revise the mass flux
						5. revise the latent heat flux according to the revised mass flux
	[outputs] 1. permeation flux in C_UDMI(1)
	          2. heat flux in C_UDMI(2)
*/
{
	int i = 0, iside = 0;
	real mass_flux, latent_heat_flux, conductive_heat_flux, total_heat_flux;
	cell_t i_cell[2];
	Thread *t_fluid[2];
	Domain *domain = Get_Domain(id_domain);
	t_fluid[0] = Lookup_Thread(domain, id_FeedFluid);
	t_fluid[1] = Lookup_Thread(domain, id_PermFluid);
	for (i=0; i<MAXCELLNUM; i++) 
	{
		if ((WallCell[i][0].index == 0) && (WallCell[i][1].index == 0)) // check the pre-existed workspace
		{
			if (i == 0) Message("Workspace WallCell is empty. Run the INITIATION first\n");
			break;
		}
		for (iside=0; iside<=1; iside++) // update the temperature and composition in workspace variable of WallCell
		{
			i_cell[iside] = WallCell[i][iside].index;
			WallCell[i][iside].temperature = C_T(i_cell[iside], t_fluid[iside]);
			WallCell[i][iside].massfraction.water = C_YI(i_cell[iside], t_fluid[iside], 0);
		}
		if (WallCell[i][0].massfraction.water > (1.-SatConc(WallCell[i][0].temperature))) // calculate the mass tranfer across the membrane only if the concentration is below the saturation
		{
			mass_flux = LocalMassFlux(WallCell[i][0].temperature, WallCell[i][1].temperature, WallCell[i][0].massfraction.water, WallCell[i][1].massfraction.water);
		}
		else
		{
			mass_flux = 0.;
		}
		latent_heat_flux = LocalHeatFlux(0, WallCell[i][0].temperature, WallCell[i][1].temperature, mass_flux); // calculate the heat transfer across the membrane
		conductive_heat_flux = LocalHeatFlux(1, WallCell[i][0].temperature, WallCell[i][1].temperature, mass_flux);
		total_heat_flux = LocalHeatFlux(10, WallCell[i][0].temperature, WallCell[i][1].temperature, mass_flux);
		if (id_message+opt > 2) Message("Cell pair of #%d and #%d: T = (%g, %g) K, with JM = %g, and JH_c = %g, JH_v = %g \n", i_cell[0], i_cell[1], WallCell[i][0].temperature, WallCell[i][1].temperature, mass_flux, conductive_heat_flux, latent_heat_flux);
		for (iside=0; iside<=1; iside++) // set the UDMI(1) and UDMI(2) as the mass and heat sources, respectively
		{
			WallCell[i][iside].flux.mass = mass_flux;
			WallCell[i][iside].flux.heat = total_heat_flux;
			C_UDMI(i_cell[iside], t_fluid[iside], 1) = mass_flux*WallCell[i][iside].area/WallCell[i][iside].volume; // SM = JM*A/V
			C_UDMI(i_cell[iside], t_fluid[iside], 2) = latent_heat_flux*WallCell[i][iside].area/WallCell[i][iside].volume; // SH = JH*A/V
		}
	}
	Monitor_CellPair(1, rid++, 17); // output the 17th cell pair for monitor
	return;
}

DEFINE_INIT(idf_cells_1007, domain)
/* 
   [objectives] 1. identify the cell pairs, which are adjacent to both sides of the membrane
                2. find the corresponding cells with the same x-coordinate
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
	real InterfaceArea[ND_ND];
	int i = 0;
	real temp = 0.0;

	t_FeedFluid = Lookup_Thread(domain, id_FeedFluid);
	t_PermFluid = Lookup_Thread(domain, id_PermFluid);
	t_FeedInterface = Lookup_Thread(domain, id_FeedInterface);
	t_PermInterface = Lookup_Thread(domain, id_PermInterface);
	fout0 = fopen("idf_cell0.out", "w");
	fout1 = fopen("idf_cell1.out", "w");

	//if(!Data_Valid_P()) 
	//{
	//	Message("\n[idf_cells] Some accessing variables have not been allocated.\n");
	//	Message("[idf_cells] The wall cells have not been identified yet.\n");
	//	return;
	//}

	begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	{
		i_cell0 = F_C0(i_face0, t_FeedInterface); // get the index for the cell adjacent the 
		C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
		C_UDMI(i_cell0, t_FeedFluid, 0) = -1; // mark the wall cells as -1, and others as 0 (no modification)
		WallCell[gid][0].index = i_cell0; // output the cell info to the workspace variable WallCell
		WallCell[gid][0].centroid[0] = loc0[0];
		WallCell[gid][0].centroid[1] = loc0[1];
		WallCell[gid][0].temperature = C_T(i_cell0, t_FeedFluid);
		WallCell[gid][0].volume = C_VOLUME(i_cell0, t_FeedFluid);
		F_AREA(InterfaceArea, i_face0, t_FeedInterface);
		WallCell[gid][0].area = NV_MAG(InterfaceArea);
		WallCell[gid][0].massfraction.water = C_YI(i_cell0, t_FeedFluid, 0);
		begin_f_loop(i_face1, t_PermInterface) // search the symmetric cell (THE LOOP CAN ONLY RUN IN SERIAL MODE)
		{
			i_cell1 = F_C0(i_face1, t_PermInterface);
			C_CENTROID(loc1, i_cell1, t_PermFluid);
			if (fabs(loc0[0]-loc1[0])/loc0[0] < EPS) // In this special case, the pair of wall cells on both sides of membrane are symmetrical
			{
				fprintf(fout0, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
				WallCell[gid][1].index = i_cell1;
				WallCell[gid][1].centroid[0] = loc1[0];
				WallCell[gid][1].centroid[1] = loc1[1];
				WallCell[gid][1].temperature = C_T(i_cell1, t_PermFluid);
				WallCell[gid][1].volume = C_VOLUME(i_cell1, t_PermFluid);
				F_AREA(InterfaceArea, i_face1, t_PermInterface);
				WallCell[gid][1].area = NV_MAG(InterfaceArea);
				WallCell[gid][1].massfraction.water = C_YI(i_cell1, t_PermFluid, 0);
			}
		}
		end_f_loop(i_face1, t_PermInterface)
		gid++;
	}
	end_f_loop(i_face0, t_FeedInterface)

	Message("The workspace WallCell[%d] has been created.\n", gid);
	Message("The identified wall cells, by using the feed-side enumeration, are summarized in idf_cell0.out.\n");

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
				fprintf(fout1, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
			}
		}
		end_f_loop(i_face1, t_PermInterface)
		C_UDMI(i_cell1, t_PermFluid, 0) = +1; // mark the cell
	}
	end_f_loop(i_face1, t_PermInterface)
	
	Message("The identified wall cells, by using the permeate-side enumeration, are summarized in redundant idf_cell1.out.\n");

	fclose(fout0);
	fclose(fout1);
}

DEFINE_ON_DEMAND(InterCellID_0923)
{
	cell_t i_cell, i_cell0, i_cell1;
	face_t i_face0, i_face1;
	Thread *t_cell;
	real loc[ND_ND], loc0[ND_ND], loc1[ND_ND];
	int i = 0;
	real temp = 0.0;
	Domain *domain = Get_Domain(id_domain);
	Thread *t_FeedFluid = Lookup_Thread(domain, id_FeedFluid);
	Thread *t_PermFluid = Lookup_Thread(domain, id_PermFluid);
	Thread *t_FeedInterface = Lookup_Thread(domain, id_FeedInterface);
	Thread *t_PermInterface = Lookup_Thread(domain, id_PermInterface);

	//fout2 = fopen("idf_cell2.out", "w");
	//fout3 = fopen("idf_cell3.out", "w");
	//begin_c_loop(i_cell, t_FeedFluid){
	//	Message("cell#%d\n", i_cell);
	//}
	//end_c_loop(i_cell, t_FeedFluid)

	begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	{
		i_cell0 = F_C0(i_face0, t_FeedInterface);
		C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
		Message("Feed-side interface#%d's adjacent cell index is #%d, located at (%g,%g)\n", i_face0, i_cell0, loc0[0], loc0[1]);
		C_UDMI(i_cell0, t_FeedFluid, 0) = -1; // mark the wall cells as -1, and others as 0 (no modification)
		gid++;
	}
	end_f_loop(i_face0, t_FeedInterface)
}

DEFINE_ON_DEMAND(WallCellProp_1018)
/*
	[objectives] check following properties of the wall cells: specific heat (cp)
	                                                           mass fraction (wx)
															                               density (rho)
															                               enthalpy (h)
															                               volume of the cell (vol)
																														 area vector (A[])
	[methods] get the properties by built-in macros 
	[outputs] FLUENT command-line output
*/
{
	int i = 0;
	cell_t i_cell;
	face_t i_face;
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	real cp[2], wx[2], rho[2], h[2], vol[2];
	real T[2] = {353.2, 303.2}, wi[2] = {0.965, 1.0};
	real A[ND_ND];
	Domain *domain = Get_Domain(id_domain);
	t_FeedFluid = Lookup_Thread(domain, id_FeedFluid);
	t_PermFluid = Lookup_Thread(domain, id_PermFluid);
	t_FeedInterface = Lookup_Thread(domain, id_FeedInterface);
	t_PermInterface = Lookup_Thread(domain, id_PermInterface);
	//for (i=0; i<MAXCELLNUM; i++)
	//{
	//	cp[0] = C_CP(WallCell[i][0].index, t_FeedFluid);
	//	cp[1] = C_CP(WallCell[i][1].index, t_PermFluid);
	//	wx[0] = C_YI(WallCell[i][0].index, t_FeedFluid, 0);
	//	wx[1] = C_YI(WallCell[i][1].index, t_PermFluid, 0);
	//	rho[0] = C_R(WallCell[i][0].index, t_FeedFluid);
	//	rho[1] = C_R(WallCell[i][1].index, t_PermFluid);
	//	h[0] = C_H(WallCell[i][0].index, t_FeedFluid);
	//	h[1] = C_H(WallCell[i][1].index, t_PermFluid);
	//	vol[0] = C_VOLUME(WallCell[i][0].index, t_FeedFluid);
	//	vol[1] = C_VOLUME(WallCell[i][1].index, t_PermFluid);
	//	Message("%d. Specific heat %g, mass fraction %g, density %g, enthalpy %g and cell volume %g\n", i, cp[0], wx[0], rho[0], h[0], vol[0]);
	//	if ((WallCell[i][1].index == 0) & (WallCell[i][1].index == 0)) return;
	//}
	//begin_f_loop(i_face, t_FeedInterface)
	//{
	//	i_cell = F_C0(i_face, t_FeedInterface);
	//	F_AREA(A, i_face, t_FeedInterface);
	//	Message("Cell#%d area vector is [%g, %g]\n", i_cell, A[0], A[1]);
	//}
	//end_f_loop(i_face, t_FeedInterface)
	Message("psat_h2o(%g) = %g\n", T[1], psat_h2o(T[1]));
	Message("WaterVaporPressure_brine(%g, %g) = %g\n", T[0], wi[0], WaterVaporPressure_brine(T[0], wi[0]));
}

DEFINE_ON_DEMAND(TProfile_0914)
/*
	[objectives] check the heat efficiency
	[Preliminary] run the initiation step in FLUENT to identify the wall cells
	[methods] 1. get the initial wall-temperature profiles 
	          2. calculate the permeation flux along the membrane (JM)
						3. calculate the heat flux related to the evaporation/condensation (JH_v)
						4. calculate the conductive heat flux (JH_c)
						5. eta = JH_v/(JH_v+JH_c)
	[outputs] FLUENT command-line output
*/
{
	MembraneTransfer(2); // opt = 2 indicates to show the debug messages
}

DEFINE_ON_DEMAND(OutputCells_0913)
/*
	[objectives] output the cells' info
	[Preliminary] run the case
	[methods] 1. check the flag of the records
	          2. output the messages of CellPairInfo
	[outputs] FLUENT command-line output
*/
{
	int idx = 0;
	Message("Output the info of the specified pair of cells ... \n");
	if (rid == 0)
	{
		Message("The CellPairInfo is empty! \n");
		return;
	}
	for (idx=0; idx<=rid; idx++)
	{
		Message("#%d: (flag %d) %s\n", idx, CellPairInfo[idx].flag, CellPairInfo[idx].content);
		if (CellPairInfo[idx].flag == 0) break;
	}
	return;
}

DEFINE_ADJUST(calc_flux_1007, domain)
/*
	[objectives] calculate the flux across the membrane
	[methods] revoke the function of MembraneTransfer with argument 0 (for normal run)
	                                                                1 (for debug run)
	[outputs] 1. C_UDMI(1) for mass flux
	          2. C_UDMI(2) for latent heat flux
*/
{
	MembraneTransfer(0);
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
	source = C_UDMI(i_cell, t_cell, 0)*C_UDMI(i_cell, t_cell, 1); // UDMI(0) = (+1/0/-1) and UDMI(1) stores the positive mass source
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
	source = C_UDMI(i_cell, t_cell, 0)*C_UDMI(i_cell, t_cell, 2); // UDMI(0) = (+1/0/-1) and UDMI(2) stores the positive heat source
  dS[eqn] = 0.;
  return source;
}

DEFINE_PROFILE(heat_flux_1008, t_face, SettingVariable)
/*
	[objectives] set the heat flux for either feed-side or permeate-side inteface between the membrane and feeding fluid
	[methods] 1. Get the index of adhered cell
	          2. Retrieve the cell number of workspace according to the cell index
						3. Get the temperature pair
						4. Calculate the mass flux and total heat flux
						5. Set the heat flux
	[outputs] the boundary condition of heat flux 
*/
{
	cell_t i_cell;
	face_t i_face;
	real Tw[2], wx[2];
	int iside, cell_num;
	real mass_flux, heat_flux;
	real dir = 0.0;
	if (THREAD_ID(t_face) == id_FeedInterface) 
	{
		dir = -1.0;
	}
	else if (THREAD_ID(t_face) == id_PermInterface)
	{
		dir = 1.0;
	}
	begin_f_loop(i_face, t_face)
	{
		i_cell = F_C0(i_face, t_face);
		cell_num = GetWID(i_cell);
		if (cell_num == -1) 
		{
			Message("The cell#%d hasn't matched with the workspace.\n", i_cell);
			return;
		}
		for (iside = 0; iside<=1; iside++)
		{
			Tw[iside] = WallCell[cell_num][iside].temperature;
			wx[iside] = WallCell[cell_num][iside].massfraction.water;
		}
		mass_flux = LocalMassFlux(Tw[0], Tw[1], wx[0], wx[1]);
		heat_flux = LocalHeatFlux(10, Tw[0], Tw[1], mass_flux);
		F_PROFILE(i_face, t_face, SettingVariable) = dir*heat_flux;
	}
	end_f_loop(i_face, t_face)
}

DEFINE_PROPERTY(density_aqNaCl_1017, i_cell, t_cell)
/*
	[objectives] set the fluid density
	[methods] 1. get the cell's temperature and mass fraction of NaCl
	          2. invoke the Density_aqNaCl()
						3. set the calculated density
	[outputs] density of aqueous NaCl solution
*/
{
	real T, w_nv, rho = 0.;
	T = C_T(i_cell, t_cell)-273.15; // celcius degree
	w_nv = C_YI(i_cell, t_cell, 1); // non-violatile component
	rho = Density_aqNaCl(T, w_nv);
	return rho;
}

DEFINE_PROPERTY(viscosity_aqNaCl_1017, i_cell, t_cell)
/*
	[objs] set the fluid viscosity
	[meth] 1. get the cell's temperature and mass fraction of NaCl
	       2. invoke the Viscosity_aqNaCl()
				 3. set the calculated viscosity with the empirical correlation
	[outs] density of aqueous NaCl solution
*/
{
	real T, w_nv, mu = 0.;
	T = C_T(i_cell, t_cell)-273.15; // celcius degree
	w_nv = C_YI(i_cell, t_cell, 1); // non-violatile component
	mu = Viscosity_aqNaCl(T, w_nv);
	return mu;
}

DEFINE_PROPERTY(thermal_conductivity_aqNaCl_1017, i_cell, t_cell)
/*
	[objs] set the fluid thermal conductivity
	[meth] 1. get the cell's temperature and mass fraction of NaCl
	       2. invoke the ThermCond_aq()
				 3. set the calculated thermal conductivity with the empirical correlation
	[outs] thermal conductivity of aqueous NaCl solution
*/
{
	real T, w_nv, lambda = 0.;
	T = C_T(i_cell, t_cell)-273.15; // celcius degree
	w_nv = C_YI(i_cell, t_cell, 1); // non-violatile component
	lambda = ThermCond_aq(T, w_nv);
	return lambda;
}

int GetWID(int searching_cell_index)
/*
	[objectives] Get the array index of workspace WallCell
	[methods] for each element of cell index in the workspace, check with the given index
	[outputs] if the index is equal to the given one, output the index; otherwise output the negative
*/
{
	int result=-1, i=0, iside=0;
	for (i=0; i<MAXCELLNUM; i++)
	{
		for (iside=0; iside<=1; iside++)
		{
			if (WallCell[i][iside].index == searching_cell_index)
			{
				result = i;
				return result;
			}
		}
	}
	if (result = -1) Message("[WARNING] None of cell index has been found in WallCell.\n");
	return result;
}