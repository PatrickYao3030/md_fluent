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

FILE *fout0, *fout1, *fout2, *fout3;

struct CellInfo
{
	int seq;
	int index;
	real centroid[ND_ND];
} WallCells[999];

int gid = 0;

DEFINE_INIT(idf_cells, domain)
// identify the cells, which are adjacent to both sides of the membrane
// use C_UDMI(0) to store the identifier, which marks as -1 or 1 for the adjacent cells and 0 (not set) for others 
{
	cell_t i_cell0, i_cell1; // global cell index
	face_t i_face0, i_face1; // global face index
	Thread *t_FeedFluid, *t_PermFluid;
	Thread *t_FeedInterface, *t_PermInterface;
	real loc0[ND_ND], loc1[ND_ND];
	int i = 0;

	t_FeedFluid = Lookup_Thread(domain, 5);
	t_PermFluid = Lookup_Thread(domain, 6);
	t_FeedInterface = Lookup_Thread(domain, 13);
	t_PermInterface = Lookup_Thread(domain, 14);
	//fout0 = fopen("idf_cells0.out", "w");
	//fout1 = fopen("idf_cells1.out", "w");
	//fout2 = fopen("idf_cell2.out", "w");
	//fout3 = fopen("idf_cell3.out", "w");

	begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	{
		i_cell0, i_cell1 = -1;
		i_cell0 = F_C0(i_face0, t_FeedInterface);
		i_cell1 = F_C1(i_face0, t_FeedInterface);
		if (i_cell1 == -1)
		{
			C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
			C_UDMI(i_cell0, t_FeedFluid, 0) = -1; // mark the cell
		}
	}
	end_f_loop(i_face0, t_FeedInterface)

	begin_f_loop(i_face1, t_PermInterface) // find the adjacent cells for the permeate-side membrane.
	{
		i_cell0, i_cell1 = -1;
		i_cell0 = F_C0(i_face1, t_PermInterface);
		i_cell1 = F_C1(i_face1, t_PermInterface);
		if (i_cell1 == -1)
		{
			C_CENTROID(loc1, i_cell0, t_PermFluid); // get the location of cell centroid
			C_UDMI(i_cell0, t_PermFluid, 0) = +1; // mark the cell
		}
	}
	end_f_loop(i_face1, t_PermInterface)

	//for (i = 0; i<999; i++)
	//	fprintf(fout2, "%d Cell index %d and locates at %g %g\n", WallCells[i].seq, WallCells[i].index, WallCells[i].centroid[0], WallCells[i].centroid[1]);
	//fclose(fout2);
	//fclose(fout0);
	//fclose(fout1);
}

DEFINE_ADJUST(calc_flux, domain)
// calculate the flux across the membrane
// the flux will also convert into the source of the adjacent cell, which will store in C_UDMI(1) for passing to the source term
{
	Thread *t_FeedFluid, *t_PermFluid;
	cell_t i_cell0, i_cell1;
	real TW[2], XW[2]; // wall temperatures and mass fraction of component 0
	real loc0[ND_ND], loc1[ND_ND];
	real coeff = 3.e-7;
	real mass_flux, heat_flux; 

	t_FeedFluid = Lookup_Thread(domain, 5);
	t_PermFluid = Lookup_Thread(domain, 6);

	begin_c_loop(i_cell0, t_FeedFluid)
	{
		if (C_UDMI(i_cell0, t_FeedFluid, 0) == -1) // the cell is adjacent to the feed-side membrane wall
		{
			TW[0] = C_T(i_cell0, t_FeedFluid); // the component of TW[0] indicates the feed-side temperature of the membrane surface
			XW[0] = C_YI(i_cell0, t_FeedFluid, 0);
			C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the coordinate of the cell
			C_UDMI(i_cell0, t_FeedFluid, 1) = TW[0];
		}
	}
	end_c_loop(i_cell0, t_FeedFluid)
	begin_c_loop(i_cell1, t_PermFluid)
		if (C_UDMI(i_cell1, t_PermFluid, 0) == 1) 
		{
			TW[1] = C_T(i_cell1, t_PermFluid); // the component of TW[1] means the permeate-side wall temperature
			C_UDMI(i_cell1, t_PermFluid, 1) = TW[1];
		}
	end_c_loop(i_cell1, t_PermFluid)
}

DEFINE_SOURCE(evap_adj_membr, i_cell, t_cell, dS, eqn)
{
	real source; // returning result

	source = C_UDMI(i_cell, t_cell, 0)*.01;
  dS[eqn] = 0.;

  return source;
}
