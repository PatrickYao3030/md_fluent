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

#define MAXCELLNUM 9999

FILE *fout0, *fout1, *fout2, *fout3;

int UC_cell_index[MAXCELLNUM];
real UC_cell_centroid[MAXCELLNUM][ND_ND];
int gid = 0;

//struct CellInfo
//{
//	int seq;
//	int index;
//	real centroid[ND_ND];
//} WallCells[999], UCCell;


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
	extern int gid; // declare the global variables

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

	thread_loop_c(t_cell, domain)
	{
		begin_c_loop(i_cell, t_cell)
		{
			C_CENTROID(loc, i_cell, t_cell);
			fprintf(fout2, "%d cell#%d %g %g\n", gid, i_cell, loc[0], loc[1]);
			UC_cell_index[gid] = i_cell;
			UC_cell_centroid[gid][0] = loc[0];
			UC_cell_centroid[gid][1] = loc[1];
			fprintf(fout3, "Cell index %d and locates at %g %g\n", UC_cell_index[gid], UC_cell_centroid[gid][0], UC_cell_centroid[gid][1]);
			gid++;
		}
		end_c_loop(i_cell, t_cell)
	}

	//begin_f_loop(i_face0, t_FeedInterface) // find the adjacent cells for the feed-side membrane.
	//{
	//	i_cell0 = F_C0(i_face0, t_FeedInterface);
	//	C_CENTROID(loc0, i_cell0, t_FeedFluid); // get the location of cell centroid
	//	fprintf(fout2, "Feeding interface#%d, adj.cell#%d at %g %g\n", i_face0, i_cell0, loc0[0], loc0[1]);
	//	//begin_f_loop(i_face1, t_PermInterface)
	//	//{
	//	//	i_cell1 = F_C0(i_face1, t_PermInterface);
	//	//	C_CENTROID(loc1, i_cell1, t_PermFluid);
	//	//	fprintf(fout2, "i_cell0-%d, %g, %g, i_cell1-%d, %g, %g\n", i_cell0, loc0[0], loc0[1], i_cell1, loc1[0], loc1[1]);
	//	//}
	//	//end_f_loop(i_face1, t_PermInterface)
	//	C_UDMI(i_cell0, t_FeedFluid, 0) = -1; // mark the cell
	//}
	//end_f_loop(i_face0, t_FeedInterface)

	//begin_f_loop(i_face1, t_PermInterface) // find the adjacent cells for the permeate-side membrane.
	//{
	//	i_cell1 = F_C0(i_face1, t_PermInterface);
	//	C_CENTROID(loc1, i_cell1, t_PermFluid); // get the location of cell centroid
	//	fprintf(fout2, "Permeating interface#%d, adj.cell#%d at %g %g\n", i_face1, i_cell1, loc1[0], loc1[1]);
	//	C_UDMI(i_cell1, t_PermFluid, 0) = +1; // mark the cell
	//}
	//end_f_loop(i_face1, t_PermInterface)

	//for (i = 0; i<9999; i++)
	//	fprintf(fout3, "%d Cell index %d and locates at %g %g\n", Cells[i].seq, Cells[i].index, Cells[i].centroid[0], Cells[i].centroid[1]);
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
