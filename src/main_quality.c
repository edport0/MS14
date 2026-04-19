#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri;
  double to, ti;
  const double clip_max = 10.0;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }

  //--- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 0);
  ti        = clock();

  if (!Msh)
    return 0;

  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));
  double* QalClip = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));
  if (!Qal || !QalClip) {
    fprintf(stderr, "## ERROR: failed to allocate quality fields\n");
    free(Qal);
    free(QalClip);
    return 1;
  }

  msh_quality(Msh, Qal, 1);

  double qmin = 1e30, qmax = 0.0, qmean = 0.0;
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
      if (Qal[iTri] < qmin) qmin = Qal[iTri];
      if (Qal[iTri] > qmax) qmax = Qal[iTri];
      qmean += Qal[iTri];
      QalClip[iTri] = fmin(Qal[iTri], clip_max);
  }
  qmean /= Msh->NbrTri;
  printf("Quality Q2 - min: %.4f  mean: %.4f  max: %.4f\n", qmin, qmean, qmax);

  if(!msh_write2dfield_Triangles("quality.solb", Msh->NbrTri, QalClip)){
    fprintf(stderr,"## ERROR: could not write clipped quality.solb\n");
    free(Qal);
    free(QalClip);
    return 1;
  }

  free(Qal);
  free(QalClip);

  return 0;
}
