#include "mesh.h"

static void print_triangle(Mesh* Msh, int iTri)
{
  printf("Triangle %d: (%d, %d, %d)\n",
         iTri,
         Msh->Tri[iTri][0],
         Msh->Tri[iTri][1],
         Msh->Tri[iTri][2]);
}

static void print_boundary_edge(BoundaryEdge* PilEdg, int k)
{
  printf("PilEdg[%d]: edge (%d, %d), from triangle %d edge %d\n",
         k,
         PilEdg[k].a,
         PilEdg[k].b,
         PilEdg[k].iTri,
         PilEdg[k].iEdg);
}

static void print_all_triangles(Mesh* Msh)
{
  for (int iTri = 1; iTri <= Msh->NbrTri; ++iTri) {
    if (Msh->Tri[iTri][0] <= 0)
      continue;
    print_triangle(Msh, iTri);
  }
}

static int run_cavity_case(Mesh* Msh, double x, double y)
{
  int iTriSeed, iP;
  int* cav;
  int nCav;
  BoundaryEdge* PilEdg;
  int nPilEdg;
  double beta[3];

  if (!msh_neighbors(Msh)) {
    fprintf(stderr, "## ERROR: failed to compute neighbors\n");
    return 0;
  }

  if (!find_point_in_mesh(Msh, x, y, &iTriSeed, beta)) {
    fprintf(stderr, "## ERROR: point (%.6f, %.6f) is outside the mesh\n", x, y);
    return 0;
  }

  iP = msh_add_vertex(Msh, x, y);
  if (iP == 0) {
    fprintf(stderr, "## ERROR: failed to add test point\n");
    return 0;
  }

  cav = malloc(Msh->NbrTri * sizeof(int));
  PilEdg = malloc((3 * Msh->NbrTri + 3) * sizeof(BoundaryEdge));
  if (!cav || !PilEdg) {
    fprintf(stderr, "## ERROR: failed to allocate cavity work arrays\n");
    free(cav);
    free(PilEdg);
    return 0;
  }

  if (!msh_build_cavity(Msh, iTriSeed, iP, cav, &nCav, PilEdg, &nPilEdg)) {
    fprintf(stderr, "## ERROR: cavity construction failed\n");
    free(cav);
    free(PilEdg);
    return 0;
  }

  printf("Point P = vertex %d at (%.6f, %.6f)\n", iP, x, y);
  printf("Seed triangle = %d\n", iTriSeed);
  printf("Cavity size = %d\n", nCav);
  for (int k = 0; k < nCav; ++k)
    print_triangle(Msh, cav[k]);

  printf("Boundary size = %d\n", nPilEdg);
  for (int k = 0; k < nPilEdg; ++k)
    print_boundary_edge(PilEdg, k);

  if (!msh_delete_cavity(Msh, cav, nCav)) {
    fprintf(stderr, "## ERROR: cavity deletion failed\n");
    free(cav);
    free(PilEdg);
    return 0;
  }

  if (!msh_star_cavity(Msh, iP, PilEdg, nPilEdg, cav, nCav)) {
    fprintf(stderr, "## ERROR: cavity starring failed\n");
    free(cav);
    free(PilEdg);
    return 0;
  }

  if (!msh_neighbors(Msh)) {
    fprintf(stderr, "## ERROR: failed to rebuild neighbors after starring\n");
    free(cav);
    free(PilEdg);
    return 0;
  }

  printf("After deletion + starring\n");
  print_all_triangles(Msh);

  free(cav);
  free(PilEdg);

  return 1;
}

int main(int argc, char* argv[])
{
  Mesh* Msh;
  double x = 0.9;
  double y = 0.2;

  if (argc < 2) {
    printf("usage : %s meshfile [x y]\n", argv[0]);
    return 0;
  }

  if (argc > 2) x = atof(argv[2]);
  if (argc > 3) y = atof(argv[3]);

  Msh = msh_read(argv[1], 0);
  if (!Msh) {
    fprintf(stderr, "## ERROR: failed to read mesh\n");
    return 1;
  }

  printf("== Case 1: expected one-triangle cavity ==\n");
  if (!run_cavity_case(Msh, x, y))
    return 1;

  Msh = msh_read(argv[1], 0);
  if (!Msh) {
    fprintf(stderr, "## ERROR: failed to reread mesh\n");
    return 1;
  }

  printf("\n== Case 2: expected two-triangle cavity ==\n");
  if (!run_cavity_case(Msh, 0.6, 0.3))
    return 1;

  return 0;
}
