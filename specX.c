/********************************************************************
 This UDF is an implementation of the species sink/source near the
  wall boundary condition. 
 It's a preliminary use for membrane distillation.
********************************************************************/
/******************************************************************** 
 Allocate the appropriate number (1) of memory location(s) in the 
  User-Defined Memory dialog box in ANSYS FLUENT
********************************************************************/

#include "udf.h"
#include "mem.h"
#include "metric.h"

FILE *fout;

DEFINE_INIT(idf_cells, domain)
// identify the cells, which are adjacent to both sides of the membrane
// use C_UDMI(0) to store the identifier, which marks as -1 or 1 for the adjacent cells and 0 (not set) for others 
{
 	cell_t i_cell0, i_cell1; // global cell index
	cell_t i_cell;
	face_t i_face = -1; // global face index
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	Thread *t_cell;
	real loc0[ND_ND], loc1[ND_ND];
	int i_local;

	t_FeedFluid = Lookup_Thread(domain, 5);
	t_PermFluid = Lookup_Thread(domain, 6);
	t_FeedInterface = Lookup_Thread(domain, 13);
	t_PermInterface = Lookup_Thread(domain, 14);
	fout = fopen("idf_cells.out", "w");

	thread_loop_c(t_cell, domain)
	{
		begin_c_loop(i_cell, t_cell)
		{
			C_CENTROID(loc0, i_cell, t_cell); // get the location of current cell
			begin_f_loop(i_face, t_FeedInterface) // search the adjacent cells for the feed-side membrane. REMARKS: it has a very low efficiency.
			{
				i_cell0, i_cell1 = -1;
				i_cell0 = F_C0(i_face, t_FeedInterface);
				i_cell1 = F_C1(i_face, t_FeedInterface);
				if (i_cell1 == -1)
				{
					C_CENTROID(loc1, i_cell0, t_FeedFluid); // get the location of cell centroid
					if ((loc0[0] == loc1[0]) && (loc0[1] == loc1[1])) // do NOT use the sentence of ELSE because the cells will be evaluated repeatedly.
					{
						//fprintf(fout, "Cell#%d %g %g\n", i_cell0, loc1[0], loc1[1]);
						C_UDMI(i_cell, t_cell, 0) = -1;
					}
				}
			}
			end_f_loop(i_face, t_FeedInterface)
			begin_f_loop(i_face, t_PermInterface) // search the adjacent cells for the permeate-side membrane
			{
				i_cell0, i_cell1 = -1;
				i_cell0 = F_C0(i_face, t_PermInterface);
				i_cell1 = F_C1(i_face, t_PermInterface);
				if (i_cell1 == -1)
				{
					C_CENTROID(loc1, i_cell0, t_PermFluid);
					if ((loc0[0] == loc1[0]) && (loc0[1] == loc1[1]))
					{
						C_UDMI(i_cell, t_cell, 0) = 1;
					}
				}
			}
			end_f_loop(i_face, t_PermInterface)
		}
		end_c_loop(i_cell, t_cell)
	}
}

DEFINE_ADJUST(calc_flux, domain)
// calculate the flux across the membrane
// the flux will also convert into the source of the adjacent cell, which will store in C_UDMI(1) for passing to the source term
{
	Thread *t_cell;
	cell_t i_cell;
	real TW[2], XW[2]; // wall temperatures and mass fraction of component 0
	real loc0[ND_ND], loc1[ND_ND];
	real coeff = 3.e-7;
	real mass_flux, heat_flux; 

	thread_loop_c(t_cell, domain)
	{
		begin_c_loop(i_cell, t_cell)
		{
			if (C_UDMI(i_cell, t_cell, 0) == -1) // the cell is adjacent to the feed-side membrane wall
			{
				TW[0] = C_T(i_cell, t_cell); // the component of TW[0] indicates the feed-side temperature of the membrane surface
				XW[0] = C_YI(i_cell, t_cell, 0);
				C_CENTROID(loc0, i_cell, t_cell); // get the coordinate of the cell
			}
			if (C_UDMI(i_cell, t_cell, 0) == +1)
			{
				C_CENTROID(loc1, i_cell, t_cell);
				if (loc0[0] == loc1[0])
				{
					TW[1] = C_T(i_cell, t_cell); // the component of TW[1] means the permeate-side wall temperature
					XW[1] = C_YI(i_cell, t_cell, 0);
				}
			}
			mass_flux = coeff*(psat_h2o(TW[0])-psat_h2o(TW[1])); // use psat_h2o() in toolkits.c to calculate the saturated vapor pressure of h2o for the given temperature
			C_UDMI(i_cell, t_cell, 1) = mass_flux;
		}
		end_c_loop(i_cell, t_cell)
	}
}

DEFINE_SOURCE(evap_adj_membr, i_cell, t_cell, dS, eqn)
{
	real source; // returning result

	source = C_UDMI(i_cell, t_cell, 0)*.01;
  dS[eqn] = 0.;

  return source;
}
