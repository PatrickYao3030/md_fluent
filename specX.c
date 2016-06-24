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

DEFINE_SOURCE(evap_adj_membr, i_cell, t_cell, dS, eqn)
{
	Domain *domain; 
	Thread *t_FeedInterface; // pointer of the thread of faces
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
			i_cell0 = F_C0(i_face, t_FeedInterface); // return the adjacent cell index of the face, ref. "3.2.5 connectivity macro"
			i_cell1 = F_C1(i_face, t_FeedInterface); // return none at the boundary
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
