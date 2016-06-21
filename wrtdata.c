/*****************************************************************
    Example of UDF for single phase that uses Get_Domain utility
******************************************************************/
// Study the features for identification of boundary by GGQ @ 2016/06/21

#include "udf.h"

FILE *fout;

void Print_Thread_Face_Centroids(Domain *domain, int id)
{
  real FC[ND_ND];
  face_t f;
  Thread *t = Lookup_Thread(domain, id);
  fprintf(fout,"thread id %d\n", id);
  begin_f_loop(f,t)
  {
    F_CENTROID(FC,f,t);
    if (ND_ND == 2)
		{
			fprintf(fout, "f%d %g %g\n", f, FC[0], FC[1]);
		}
		if (ND_ND == 3)
		{
			fprintf(fout, "f%d %g %g %g\n", f, FC[0], FC[1], FC[2]);
		}
  }
  end_f_loop(f,t)
  fprintf(fout, "\n");
}

DEFINE_ON_DEMAND(get_coords)
{
  Domain *domain;
  domain = Get_Domain(1);
  fout = fopen("faces.out", "w");
	fprintf(fout, "%d dimensional domain\n", ND_ND);
  Print_Thread_Face_Centroids(domain, 13);
	Print_Thread_Face_Centroids(domain, 14);
  fclose(fout);
}

DEFINE_ON_DEMAND(idf_bound)
// identify the adjacent cells for the specified boundary
// output the identified cells' info into data file
// the output process has a low efficiency
{
  Domain *domain;
	cell_t i_cell0, i_cell1 = -1; // global cell index
	face_t i_face = -1; // global face index
	Thread *t_FeedFluid;
	Thread *t_PermFluid;
	Thread *t_FeedInterface;
	Thread *t_PermInterface;
	real loc[ND_ND];

	domain = Get_Domain(1);
	t_FeedFluid = Lookup_Thread(domain, 5);
	t_PermFluid = Lookup_Thread(domain, 6);
	t_FeedInterface = Lookup_Thread(domain, 13);
	t_PermInterface = Lookup_Thread(domain, 14);
  fout = fopen("cells.out", "w");

	begin_f_loop(i_face, t_FeedFluid)
	{
		i_cell0 = F_C0(i_face, t_FeedInterface);
		i_cell1 = F_C1(i_face, t_FeedInterface);
		if (i_cell1 == -1)
		{
			C_CENTROID(loc, i_cell0, t_FeedFluid); // get the location of cell centroid
			fprintf(fout, "Cell#%d %g %g\n", i_cell0, loc[0], loc[1]);
		}
	}
	end_f_loop(i_face, t_FeedFluid)
	fprintf(fout, "Tranfer data are completed.\n");
}