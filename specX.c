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
			begin_f_loop(i_face, t_FeedInterface) // search the adjacent cells for the feed-side membrane
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

DEFINE_SOURCE(evap_adj_membr, i_cell, t_cell, dS, eqn)
{
	Domain *domain; 
	Thread *t_FeedInterface, *t_face; // pointer of the thread of faces
	face_t i_face = -1; // index of face
	cell_t i_cell0, i_cell1 = -1; // indexes of adjacent cells for boundary identification
	int i_local;
	real source; // returning result

	domain = Get_Domain(1); // explicit declaration of mixture
	t_FeedInterface = Lookup_Thread(domain, 13); // explicitly get the boundary, whose id (13) should be consisted with GUI display
	
	if (THREAD_ID(t_cell) == 5) // check the zone, whose id (5) should be shown at GUI display
	{
		c_face_loop(i_cell, t_cell, i_local) // loop all faces for the given cell of (i_cell, t_cell)
		{
			i_face = C_FACE(i_cell, t_cell, i_local); // global face index for the given cell
			t_face = C_FACE_THREAD(i_cell, t_cell, i_local); // return the thread of the face that is returned by above C_FACE 
			i_cell0 = F_C0(i_face, t_face); // return the adjacent cell index of the face, ref. "3.2.5 connectivity macro"
			i_cell1 = F_C1(i_face, t_face); // return none at the boundary
			if (i_cell1 == -1)
			{
				source = 0.01;
			}
			else
			{
				source = 0.;
			}
		}
	}
  dS[eqn] = 0.;
  return source;
}
