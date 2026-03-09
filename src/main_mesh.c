#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri, iVer;
  double to, ti;

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

  //--- create neigbhors Q2 version
  to = clock();
  msh_neighborsQ2(Msh);
  ti = clock();
  printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors with hash table
  to = clock();
  msh_neighbors(Msh);
  ti = clock();
  printf("  time hash tab neigh.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));
  
  msh_quality(Msh, Qal, 0); // compute quality with Q1

  //---Printing stats
  double qmin = 1e30, qmax = 0.0, qmean = 0.0;
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
      if (Qal[iTri] < qmin) qmin = Qal[iTri];
      if (Qal[iTri] > qmax) qmax = Qal[iTri];
      qmean += Qal[iTri];
  }
  qmean /= Msh->NbrTri;
  printf("Quality Q1 - min: %.4f  mean: %.4f  max: %.4f\n", qmin, qmean, qmax);


  msh_write2dfield_Triangles("quality.solb", Msh->NbrTri, Qal);

  //--- TODO: compute metric field
  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }

  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);

  //--- Free memory
  if (Qal != NULL) {
    free(Qal);
    Qal = NULL;
  }
  if (Met != NULL) {
    free(Met);
    Met = NULL;
  }

  return 0;
}
