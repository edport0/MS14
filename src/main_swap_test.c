#include "mesh.h"

static void print_triangle(Mesh* Msh, int iTri)
{
   printf("Triangle %d: (%d, %d, %d) \n", iTri, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2]);
}

static void print_neighbors(Mesh* Msh, int iTri)
{
   printf("TriVoi[%d]: (%d, %d, %d) \n", iTri, Msh->TriVoi[iTri][0], Msh->TriVoi[iTri][1], Msh->TriVoi[iTri][2]);
}

int main(int argc, char* argv[])
{
    int iTri = 1;
    int iEdg = 1;

    if(argc < 2) {
      printf(" usage :%s meshfile [iTri, iEdg] \n", argv[0]);
      return 0;
    }

    if(argc > 2) iTri = atoi(argv[2]);
    if(argc > 3) iEdg = atoi(argv[3]);

    Mesh* Msh = msh_read(argv[1], 0);
    if(!Msh){
        fprintf(stderr, "## ERROR: failed to read mesh\n");
        return 1;
    }

    if(!msh_neighbors(Msh)){
      fprintf(stderr, "## ERROR: failed to compute neighbors\n");
      return 1;
    }

    printf("Initial mesh\n");
    printf("  Vertices   %d \n", Msh->NbrVer);
    printf("  Triangles  %d \n", Msh->NbrTri);

    printf("Before edge swap \n");
    print_triangle(Msh, 1);
    if(Msh->NbrTri > 1) print_triangle(Msh, 2);
    print_neighbors(Msh, iTri);
    
    int jTri = Msh->TriVoi[iTri][iEdg];
    if(jTri == 0){
        printf("## ERROR: cannot swap a boundary edge\n");
        return 1;
    }

    printf("Swapping edge %d of triangle %d with neighbor triangle %d \n", iEdg, iTri, jTri);
    if(!msh_swap_edge(Msh, iTri, iEdg)){
        fprintf(stderr, "## ERROR: failed to swap edge\n");
        return 1;
    }

    printf("After edge swap \n");
    print_triangle(Msh, 1);
    if(Msh->NbrTri > 1) print_triangle(Msh, 2);

     if(!msh_write(Msh, "swap_test_out.mesh")){
      fprintf(stderr,"## ERROR: could not write output mesh\n");
      return 1;
    }
   
    printf("Output written to swap_test_out.mesh\n");

    return 0;
}
