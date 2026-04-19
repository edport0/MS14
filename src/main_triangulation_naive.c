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
   int use_walk = 1;
   double x, y;
   double to, ti;
   double tins0, tins1, tins;

   if(argc < 2) {
     printf(" usage : mesh file [npts] [walk|full]\n");
     return 0;
   }

   if(argc >= 3)
      npts = atoi(argv[2]);

   if(argc >= 4) {
      if(strcmp(argv[3], "walk") == 0)
        use_walk = 1;
      else if(strcmp(argv[3], "full") == 0)
        use_walk = 0;
      else {
        fprintf(stderr, "## ERROR: localization mode must be 'walk' or 'full'\n");
        return 1;
      }
   }

    to = clock();
    Mesh* Msh = msh_read(argv[1], 0);
    ti = clock();

    if(!Msh)
        return 0;

    printf("Initial mesh \n");
    printf("  Vertices   %10d \n", Msh->NbrVer);
    printf("  Triangles  %10d \n", Msh->NbrTri);
    printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);
    printf("  localization mode %s\n", use_walk ? "walk" : "full");

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

      if(use_walk){
        int iTriSeed = 1;
        if(!find_point_walk(Msh, x, y, iTriSeed, &iTri, beta)){
          if(!find_point_in_mesh(Msh, x, y, &iTri, beta)){
              printf("Point outside mesh\n");
              continue;
          }
        }
      }
      else{
        if(!find_point_in_mesh(Msh, x, y, &iTri, beta)){
            printf("Point outside mesh\n");
            continue;
        }
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
          if(!msh_split_triangle_naive(Msh, iTri, iP)) return 1;
      }
      else if (nzero == 1){
          if(!msh_split_edge_naive(Msh, iTri, iEdg, iP)) return 1;
      }

      if(!msh_neighbors(Msh)){
          fprintf(stderr, "## ERROR: msh_neighbors failed after split\n");
          return 1;
      }

      // Improve only the local star around the inserted vertex using a quality criterion.
      for(;;){
        int changed = 0;

        for(int jTri = 1; jTri <= Msh->NbrTri; ++jTri){
          int locP = -1;
          for(int k = 0; k < 3; ++k){
            if(Msh->Tri[jTri][k] == iP){
              locP = k;
              break;
            }
          }

          if(locP == -1)
            continue;

          if(msh_should_swap_edge_quality(Msh, jTri, locP, 1)){
            if(!msh_swap_edge(Msh, jTri, locP)){
              fprintf(stderr, "## ERROR: quality-based edge swap failed\n");
              return 1;
            }
            changed = 1;
            break;
          }
        }

        if(!changed)
          break;
      }

      if (iP == 0) {
          printf("Skipped point %3d at (%.6f, %.6f)\n", i + 1, x, y);
          continue;
      }

      printf("Inserted point %3d: vertex %3d at (%.6f, %.6f) in triangle %d\n", i + 1, iP, x, y, iTri);
    }
    tins1 = clock();
    tins = (tins1 - tins0) / CLOCKS_PER_SEC;


    printf("\nFinal Mesh\n");
    printf("  Vertices %10d\n", Msh->NbrVer);
    printf("  Triangles %10d\n", Msh->NbrTri);
    printf("  insertion time %lg (s)\n", tins);
    printf("TIMING naive %s %d %lg\n", use_walk ? "walk" : "full", npts, tins);

    if(!msh_write(Msh, "triangulation_naive_out.mesh")){
      fprintf(stderr,"## ERROR: could not write output mesh\n");
      return 1;
    }
   
    printf("Output written to triangulation_naive_out.mesh");

    return 0;
}
