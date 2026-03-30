#include "mesh.h"

static double rand01(void)
{
  return (double)rand()/(double)RAND_MAX;
}

int main(int argc, char* argv[])
{
   int i, iTri, iP, nzero, iEdg;
   double beta[3];
   int npts = 20;
   double x, y;
   double to, ti;

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

    for(i=0; i< npts; i++)
    {
      //--random point in the unit square
      x = rand01();
      y = rand01();


      //---find containing triangle
      if(!find_point_in_mesh(Msh, x, y, &iTri, beta)){
        printf("Point outside mesh\n");
      }

      nzero = 0;
      iEdg = -1;
      for(int k = 0; k < 3; k++){
        if(fabs(beta[k]) < 1e-12){
          nzero++;
          iEdg = k;
        }
      }

      if (nzero >= 2) {
            printf("Point %d = (%g,%g) already on an existing vertex\n", i + 1, x, y);
            continue;
      }

      iP = msh_add_vertex(Msh, x, y);
      if (iP == 0) {
          fprintf(stderr, "## ERROR: msh_add_vertex failed at insertion %d\n", i + 1);
          return 1;
      }

      if(nzero == 0){
          msh_split_triangle_naive(Msh, iTri, iP);
        }
      else if (nzero == 1){
        msh_split_edge_naive(Msh, iTri, iEdg, iP);
      }

      msh_neighbors(Msh);
      printf("Inserted point %3d: vertex %3d at (%.6f, %.6f) in triangle %d\n", i + 1, iP, x, y, iTri);
    }


    printf("\nFinal Mesh\n");
    printf("  Vertices %10d\n", Msh->NbrVer);
    printf("  Triangles %10d\n", Msh->NbrTri);

    if(!msh_write(Msh, "triangulation_naive_out.mesh")){
      fprintf(stderr,"## ERROR: could not write output mesh\n");
      return 1;
    }
   
    printf("Output written to triangulation_naive_out.mesh");

    return 0;
}
