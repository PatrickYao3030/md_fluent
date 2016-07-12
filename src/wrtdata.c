/*****************************************************************
    Example of UDF for single phase that uses Get_Domain utility
******************************************************************/
// Study the features for identification of boundary by GGQ @ 2016/06/21

#include "udf.h"

FILE *fout;

void Print_Thread_Face_Centroids(Thread *t_face)
{
  real FC[ND_ND];
  face_t i_face;
  
  fprintf(fout,"Face thread id %d\n", THREAD_ID(t_face));
  begin_f_loop(i_face,t_face)
  {
    F_CENTROID(FC,i_face,t_face);
    if (ND_ND == 2)
		{
			fprintf(fout, "f%d %g %g\n", i_face, FC[0], FC[1]);
		}
		if (ND_ND == 3)
		{
			fprintf(fout, "f%d %g %g %g\n", i_face, FC[0], FC[1], FC[2]);
		}
  }
  end_f_loop(i_face,t_face)
  fprintf(fout, "Completed face thread id %d\n");
	fprintf(fout, "\n");
}
