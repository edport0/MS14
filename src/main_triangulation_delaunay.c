#include "mesh.h"

static double rand01(void)
{
  return (double)rand()/(double)RAND_MAX;
}

int main(int argc, char* argv[])
{
   int i, iP;
   int npts = 20;
   double x, y;
   double to, ti;
   double tins0, tins1, tins;

   if(argc < 2) {
     printf(" usage : mesh file \n");
     return 0;
   }

   if(argc >= 3)
      npts = atoi(argv[2]);

    //--- read a mesh
    to = clock();
    Mesh* Msh = msh_read(argv[1], 0);
    ti = clock();

    if(!Msh)
        return 0;

    printf("Initial mesh \n");
    printf("  Vertices   %10d \n", Msh->NbrVer);
    printf("  Triangles  %10d \n", Msh->NbrTri);
    printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

    //---builds neighbor table for the 1st time
    if(!msh_neighbors(Msh)){
      fprintf(stderr, "## ERROR: msh_neighbors failed\n");
      return 1;
    }

    srand(7);
    tins0 = clock();

    for(i=0; i< npts; i++)
    {
      //--random point in the unit square
      x = rand01();
      y = rand01();

      iP = msh_insert_point_delaunay(Msh, x, y);
      if (iP == 0) {
          printf("Skipped point %3d at (%.6f, %.6f)\n", i + 1, x, y);
          continue;
      }

      printf("Inserted point %3d: vertex %3d at (%.6f, %.6f)\n", i + 1, iP, x, y);
    }
    tins1 = clock();
    tins = (tins1 - tins0) / CLOCKS_PER_SEC;


    printf("\nFinal Mesh\n");
    printf("  Vertices %10d\n", Msh->NbrVer);
    printf("  Triangles %10d\n", Msh->NbrTri);
    printf("  insertion time %lg (s)\n", tins);
    printf("TIMING delaunay %d %lg\n", npts, tins);

    if(!msh_write(Msh, "triangulation_delaunay_out.mesh")){
      fprintf(stderr,"## ERROR: could not write output mesh\n");
      return 1;
    }
   
    printf("Output written to triangulation_delaunay_out.mesh");

    return 0;
}
