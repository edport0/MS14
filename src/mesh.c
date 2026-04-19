#include "mesh.h"

int tri2edg[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

Mesh* msh_init()
{
  Mesh* Msh = malloc(sizeof(Mesh));
  if (!Msh) return NULL;

  Msh->Dim    = 0;
  Msh->NbrVer = 0;
  Msh->NbrTri = 0;
  Msh->NbrEfr = 0;
  Msh->NbrEdg = 0;

  Msh->NbrVerMax = 0;
  Msh->NbrTriMax = 0;
  Msh->NbrEfrMax = 0;
  Msh->NbrEdgMax = 0;

  Msh->Box[0] = 1.e30; // xmin
  Msh->Box[1] = -1.e30; // xmax
  Msh->Box[2] = 1.e30; // ymin
  Msh->Box[3] = -1.e30; // ymax

  //--- Data for the list of vertices
  Msh->Crd = NULL;

  //--- Data for the list of triangles
  Msh->Tri    = NULL;
  Msh->TriVoi = NULL;
  Msh->TriRef = NULL;
  Msh->TriMrk = NULL;

  //--- Data for the list of boundary edges
  Msh->Efr    = NULL;
  Msh->EfrVoi = NULL;
  Msh->EfrRef = NULL;

  //--- Data for the list of edges
  Msh->Edg = NULL;

  return Msh;
}

Mesh* msh_read(char* file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref;

  int fmsh = 0;

  if (!file) return NULL;

  Mesh* Msh = msh_init();

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".mesh")) {
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".meshb");
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".mesh");
      if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);

  Msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  Msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);

  //Msh->NbrVerMax = Msh->NbrVer;
  //Msh->NbrTriMax = Msh->NbrTri;
  Msh->NbrVerMax = Msh->NbrVer + 1000; //for insertion tests
  Msh->NbrTriMax = Msh->NbrTri + 2*1000; 
  

  //--- allocate arrays
  Msh->Crd    = calloc((Msh->NbrVerMax + 1), sizeof(double3d));
  Msh->Tri    = calloc((Msh->NbrTriMax + 1), sizeof(int3d));
  Msh->TriRef = calloc((Msh->NbrTriMax + 1), sizeof(int1d));
  Msh->TriMrk = calloc((Msh->NbrTriMax + 1), sizeof(int1d));


  //--- read vertices
  GmfGotoKwd(fmsh, GmfVertices);
  if (Msh->Dim == 2) {
    if (FilVer == GmfFloat) { // read 32 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        Msh->Crd[i][0] = (double)bufFlt[0];
        Msh->Crd[i][1] = (double)bufFlt[1];
      }
    }
    else { // read 64 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        Msh->Crd[i][0] = bufDbl[0];
        Msh->Crd[i][1] = bufDbl[1];
      }
    }
  }
  else {
    fprintf(stderr, "  ## ERROR: 3D is not implemented\n");
    exit(1);
  }

  //--- read triangles
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i = 1; i <= Msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    Msh->Tri[i][0] = bufTri[0];
    Msh->Tri[i][1] = bufTri[1];
    Msh->Tri[i][2] = bufTri[2];
    Msh->TriRef[i] = bufTri[3];
  }

  //--- read boundary edges
  if (readEfr == 1) {
    Msh->NbrEfr    = GmfStatKwd(fmsh, GmfEdges);
    Msh->NbrEfrMax = Msh->NbrEfr;

    Msh->Efr    = calloc((Msh->NbrEfrMax + 1), sizeof(int2d));
    Msh->EfrRef = calloc((Msh->NbrEfrMax + 1), sizeof(int1d));

    GmfGotoKwd(fmsh, GmfEdges);
    for (i = 1; i <= Msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      Msh->Efr[i][0] = bufEfr[0];
      Msh->Efr[i][1] = bufEfr[1];
      Msh->EfrRef[i] = bufEfr[2];
    }
  }

  GmfCloseMesh(fmsh);

  return Msh;
}

double* sol_read(char* file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[GmfMaxTyp];
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;

  int fsol = 0;

  if (!file) return NULL;

  double* sol = NULL;

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".sol")) {
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".solb");
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".sol");
      if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);

  SolTyp = GmfSolAtVertices; // read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);

  if (nbrSol == 0) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if (dim != mshDim) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if (nbrSol != mshNbrSol) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (NbrTyp != 1) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if (TypTab[0] != GmfSca) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }

  sol = (double*)calloc(nbrSol + 1, sizeof(double));

  GmfGotoKwd(fsol, SolTyp);

  for (i = 1; i <= nbrSol; ++i) {
    if (FilVer == GmfFloat) {
      GmfGetLin(fsol, SolTyp, &bufFlt);
      sol[i] = (double)bufFlt;
    }
    else {
      GmfGetLin(fsol, SolTyp, &bufDbl);
      sol[i] = bufDbl;
    }
  }

  if (!GmfCloseMesh(fsol)) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    // myexit(1);
  }

  return sol;
}

int msh_boundingbox(Mesh* Msh)
{
  int1d iVer;

  //--- compute bounding box
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    // TODO: Set Msh->Box
  }

  return 1;
}

int msh_write(Mesh* Msh, char* file)
{
  int iVer, iTri, iEfr;
  int FilVer = 2;

  if (!Msh) return 0;
  if (!file) return 0;

  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++)
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0);

  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++)
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);

  if (Msh->NbrEfr > 0) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr = 1; iEfr <= Msh->NbrEfr; iEfr++)
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }

  GmfCloseMesh(fmsh);

  return 1;
}

//--- Triangle quality functions, 1 for equilateral, +inf for degenerate
static double tri_area(double* A, double* B, double* C)
{
  return 0.5 * ((B[0] - A[0])*(C[1] - A[1]) - (C[0] - A[0])*(B[1] - A[1]));
}

static double incircle(double* A, double* B, double* C, double* D)
{
  double ax = A[0] - D[0];
  double ay = A[1] - D[1];
  double bx = B[0] - D[0];
  double by = B[1] - D[1];
  double cx = C[0] - D[0];
  double cy = C[1] - D[1];
  double det;

  det = (ax * ax + ay * ay) * (bx * cy - by * cx)
      - (bx * bx + by * by) * (ax * cy - ay * cx)
      + (cx * cx + cy * cy) * (ax * by - ay * bx);

  // Positive means D is inside circumcircle(ABC) when ABC is CCW.
  if (tri_area(A, B, C) < 0.0)
    det = -det;

  return det;
}

static int tri_contains_in_circumcircle(Mesh* Msh, int iTri, int iP)
{
  const double eps = 1.e-14;
  double* A;
  double* B;
  double* C;
  double* D;

  if (!Msh) return 0;
  if (iTri < 1 || iTri > Msh->NbrTri) return 0;
  if (iP < 1 || iP > Msh->NbrVer) return 0;

  A = Msh->Crd[Msh->Tri[iTri][0]];
  B = Msh->Crd[Msh->Tri[iTri][1]];
  C = Msh->Crd[Msh->Tri[iTri][2]];
  D = Msh->Crd[iP];

  return (incircle(A, B, C, D) > eps);
}

static int tri_find_edge_idx(Mesh* Msh, int iTri, int u, int v)
{
  int k;
  int a, b;

  if (!Msh) return -1;
  if (iTri < 1 || iTri > Msh->NbrTri) return -1;

  for (k = 0; k < 3; ++k) {
    a = Msh->Tri[iTri][tri2edg[k][0]];
    b = Msh->Tri[iTri][tri2edg[k][1]];

    if ((a == u && b == v) || (a == v && b == u))
      return k;
  }

  return -1;
}

static int boundary_add_or_cancel(BoundaryEdge* PilEdg, int* nPilEdg,
                                  int a, int b, int iTri, int iEdg)
{
  int k;

  if (!PilEdg || !nPilEdg) return 0;

  for (k = 0; k < *nPilEdg; ++k) {
    if ((PilEdg[k].a == a && PilEdg[k].b == b) ||
        (PilEdg[k].a == b && PilEdg[k].b == a)) {
      PilEdg[k] = PilEdg[*nPilEdg - 1];
      (*nPilEdg)--;
      return 1;
    }
  }

  PilEdg[*nPilEdg].a    = a;
  PilEdg[*nPilEdg].b    = b;
  PilEdg[*nPilEdg].iTri = iTri;
  PilEdg[*nPilEdg].iEdg = iEdg;
  (*nPilEdg)++;

  return 1;
}

static void tri_make_ccw(Mesh* Msh, int iTri)
{
    int a = Msh->Tri[iTri][0];
    int b = Msh->Tri[iTri][1];
    int c = Msh->Tri[iTri][2];

    if(tri_area(Msh->Crd[a], Msh->Crd[b], Msh->Crd[c]) < 0.0){
        int tmp = Msh->Tri[iTri][1];
        Msh->Tri[iTri][1] = Msh->Tri[iTri][2];
        Msh->Tri[iTri][2] = tmp;
    }
}

double msh_qualityQ1(Mesh* Msh, int iTri)
{
  const double alpha1 = sqrt(3.0)/12.0;

  double* A = Msh->Crd[Msh->Tri[iTri][0]];
  double* B = Msh->Crd[Msh->Tri[iTri][1]];
  double* C = Msh->Crd[Msh->Tri[iTri][2]];

  double l1sq = (B[0] - A[0])*(B[0] - A[0]) + (B[1] - A[1])*(B[1] - A[1]);
  double l2sq = (C[0] - B[0])*(C[0] - B[0]) + (C[1] - B[1])*(C[1] - B[1]);
  double l3sq = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1]);

  double area = fabs(tri_area(A, B, C));
  if(area < 1.e-30) 
    return 1.e30; // degenerate triangle

  return alpha1 * (l1sq + l2sq + l3sq) / area;
}

double msh_qualityQ2(Mesh* Msh, int iTri)
{
  const double alpha2 = 1.0 / (2.0 * sqrt(3.0));

  double* A = Msh->Crd[Msh->Tri[iTri][0]];
  double* B = Msh->Crd[Msh->Tri[iTri][1]];
  double* C = Msh->Crd[Msh->Tri[iTri][2]];

  double l1 = sqrt((B[0] - A[0])*(B[0] - A[0]) + (B[1] - A[1])*(B[1] - A[1]));
  double l2 = sqrt((C[0] - B[0])*(C[0] - B[0]) + (C[1] - B[1])*(C[1] - B[1]));
  double l3 = sqrt((A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1]));

  double hmax = fmax(l1, fmax(l2, l3));
  double perim = l1 + l2 + l3;
  double area = fabs(tri_area(A, B, C));

  if(area < 1.e-30) 
    return 1.e30; // degenerate triangle

  double rho = 2.0 * area / perim; // inradius

  return alpha2*hmax / rho;
}

int msh_quality(Mesh* Msh, double* field, int useQ2)
{
  int iTri;
  if(!Msh || !field) return 0;

  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    if (useQ2 == 0)
      field[iTri] = msh_qualityQ1(Msh, iTri);
    else
      field[iTri] = msh_qualityQ2(Msh, iTri);
  }

  return 1;
}


static double tri_quality_pts(double* A, double* B, double* C, int useQ2)
{
  double area = fabs(tri_area(A, B, C));
  if (area < 1.e-30)
    return 1.e30;

  if (useQ2 == 0) {
    const double alpha1 = sqrt(3.0) / 12.0;
    double l1sq = (B[0] - A[0])*(B[0] - A[0]) + (B[1] - A[1])*(B[1] - A[1]);
    double l2sq = (C[0] - B[0])*(C[0] - B[0]) + (C[1] - B[1])*(C[1] - B[1]);
    double l3sq = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1]);
    return alpha1 * (l1sq + l2sq + l3sq) / area;
  }
  else {
    const double alpha2 = 1.0 / (2.0 * sqrt(3.0));
    double l1 = sqrt((B[0] - A[0])*(B[0] - A[0]) + (B[1] - A[1])*(B[1] - A[1]));
    double l2 = sqrt((C[0] - B[0])*(C[0] - B[0]) + (C[1] - B[1])*(C[1] - B[1]));
    double l3 = sqrt((A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1]));
    double hmax = fmax(l1, fmax(l2, l3));
    double perim = l1 + l2 + l3;
    double rho = 2.0 * area / perim;
    return alpha2 * hmax / rho;
  }
}


int msh_neighborsQ2(Mesh* Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTriMax), sizeof(int3d));

  memset(Msh->TriVoi, 0, (Msh->NbrTriMax) * sizeof(int3d));

  //--- Compute the neighbors using a quadratic-complexity algorithm
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {

      if(Msh->TriVoi[iTri][iEdg] != 0) 
        continue; // already set, skip

      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices
      for (jTri = 1; jTri <= Msh->NbrTri; jTri++) {
        if (iTri == jTri)
          continue;

        for (jEdg = 0; jEdg < 3; jEdg++) {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];

          if(((iVer1 == jVer1) && (iVer2 == jVer2)) || ((iVer1 == jVer2) && (iVer2 == jVer1))){
            Msh->TriVoi[iTri][iEdg] = jTri; // set neighbor of iTri at edge iEdg
            Msh->TriVoi[jTri][jEdg] = iTri; // set neighbor of jTri at edge jEdg
            break;
          }
        }
      }
    }
  }

  return 1;
}

int msh_neighbors(Mesh* Msh)
{
  int iTri, iEdg, iVer1, iVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTriMax + 1), sizeof(int3d));

  memset(Msh->TriVoi, 0, sizeof(int3d) * (Msh->NbrTriMax + 1)); // reset!

  //--- initialize HashTable and set the hash table
  int SizHead = 3*Msh->NbrTri;
  int NbrMaxObj = 3*Msh->NbrTri/2 + Msh->NbrTri;
  HashTable* hsh = hash_init(SizHead, NbrMaxObj);

  //--- Compute the neighbors using the hash table
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      int iObj = hash_find(hsh, iVer1, iVer2);

      if(iObj == 0){
        hash_add(hsh, iVer1, iVer2, iTri);
      }
      else{
        int jTri = hsh->LstObj[iObj][2]; //1st triangle having iVer1-iVer2 as edge
        hsh ->LstObj[iObj][3] = iTri; //2nd triangle

        int jEdg;
        for(jEdg=0; jEdg<3; jEdg++){
          int jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          int jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];

          if((jVer1 == iVer1 && jVer2 == iVer2) || (jVer1 == iVer2 && jVer2 == iVer1))break;
        }

        Msh->TriVoi[iTri][iEdg] = jTri;
        Msh->TriVoi[jTri][jEdg] = iTri;
      }
    }
  }

  int NbrBnd = 0;
  for(int iObj = 1; iObj <= hsh->NbrObj; iObj++){
    if(hsh->LstObj[iObj][3] == 0) NbrBnd++;
  }
  printf("  Boundary edges:           %d\n", NbrBnd);

    // --- Edge count
  printf("  Number of edges: %d\n", hsh->NbrObj);

  // --- Collision diagnostics
  int    maxChain  = 0;
  int    emptySlots = 0;
  double avgChain  = 0.0;
  int    usedSlots = 0;

  for (int i = 0; i < hsh->SizHead; i++) {
    int len  = 0;
    int iObj = hsh->Head[i];
    if (iObj == 0) { emptySlots++; continue; }
    while (iObj != 0) {
      len++;
      iObj = hsh->LstObj[iObj][4];
    }
    if (len > maxChain) maxChain = len;
    avgChain += len;
    usedSlots++;
  }
  avgChain /= usedSlots;

  printf("  Hash key:          (iVer1 + iVer2) %% SizHead\n");
  printf("  Max chain length:  %d\n", maxChain);
  printf("  Avg chain length:  %.2f\n", avgChain);
  printf("  Empty slots:       %d / %d\n", emptySlots, hsh->SizHead);

  free(hsh->Head);
  free(hsh->LstObj);
  free(hsh);

  return 1;
}

HashTable* hash_init(int SizHead, int NbrMaxObj)
{
  HashTable* hsh = malloc(sizeof(HashTable));
  if (!hsh) return NULL;

  hsh->SizHead   = SizHead;
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj    = 0;

  // calloc zeroes memory: Head[i]=0 means no obj at key i
  hsh->Head   = calloc(SizHead,    sizeof(int));
  hsh->LstObj = calloc(NbrMaxObj + 1, sizeof(int5d));

  return hsh;
}

IntStack* initStack(int cap)
{
  IntStack* stk;

  if (cap <= 0)
    cap = 16;

  stk = malloc(sizeof(IntStack));
  if (!stk) return NULL;

  stk->data = malloc(cap * sizeof(int));
  if (!stk->data) {
    free(stk);
    return NULL;
  }

  stk->top = 0;
  stk->cap = cap;

  return stk;
}

int pushStack(IntStack* stk, int value)
{
  int* newData;
  int  newCap;

  if (!stk) return 0;

  if (stk->top >= stk->cap) {
    newCap = 2 * stk->cap;
    newData = realloc(stk->data, newCap * sizeof(int));
    if (!newData) return 0;

    stk->data = newData;
    stk->cap  = newCap;
  }

  stk->data[stk->top] = value;
  stk->top++;

  return 1;
}

int popStack(IntStack* stk, int* value)
{
  if (!stk || !value) return 0;
  if (stk->top <= 0) return 0;

  stk->top--;
  *value = stk->data[stk->top];

  return 1;
}

int sizeStack(IntStack* stk)
{
  if (!stk) return 0;

  return stk->top;
}

void freeStack(IntStack* stk)
{
  if (!stk) return;

  free(stk->data);
  free(stk);
}

int hash_find(HashTable* hsh, int iVer1, int iVer2)
{
  // compute key
  int key = (iVer1 + iVer2) % hsh->SizHead;

  // start at the first object with this key
  int iObj = hsh->Head[key];

  // follow the chain until it finds a match or reach the end
  while (iObj != 0) {
    if ( (hsh->LstObj[iObj][0] == iVer1 && hsh->LstObj[iObj][1] == iVer2) ||
         (hsh->LstObj[iObj][0] == iVer2 && hsh->LstObj[iObj][1] == iVer1) ) {
      return iObj; // found! return its index in LstObj
    }
    iObj = hsh->LstObj[iObj][4]; // follow chain to next object with same key
  }

  return 0; // not found
}

int hash_add(HashTable* hsh, int iVer1, int iVer2, int iTri)
{
 
  if (hsh->NbrObj >= hsh->NbrMaxObj) {
    printf("  ## ERROR: Hash table is full ! \n");
    return 0;
  }

  int key = (iVer1 + iVer2) % hsh->SizHead;

  //allocate a new slot
  int iObj = ++hsh->NbrObj;

  //store edge data
  hsh->LstObj[iObj][0] = iVer1;
  hsh->LstObj[iObj][1] = iVer2;
  hsh->LstObj[iObj][2] = iTri; // first triangle
  hsh->LstObj[iObj][3] = 0;    // second triangle, not known yet

  //add to the head of the chain
  hsh->LstObj[iObj][4] = hsh->Head[key]; // point to previous 1st obj
  hsh->Head[key] = iObj; // point to new 1st obj

  return iObj;
}

int hash_suppr(HashTable* hsh, int iVer1, int iVer2, int iTri)
{
 
  // to be implemented

  // ===> suppress this entry in the hash tab

  return 0;
}

int msh_write2dfield_Vertices(char* file, int nfield, double* field)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);

  for (iVer = 1; iVer <= nfield; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]);

  GmfCloseMesh(fmsh);

  return 1;
}

int msh_write2dfield_Triangles(char* file, int nfield, double* field)
{
  int iTri;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);

  for (iTri = 1; iTri <= nfield; iTri++)
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]);

  GmfCloseMesh(fmsh);

  return 1;
}

int msh_write2dmetric(char* file, int nmetric, double3d* metric)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSymMat;

  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);

  for (iVer = 1; iVer <= nmetric; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]);

  GmfCloseMesh(fmsh);

  return 1;
}


int msh_add_vertex(Mesh *Msh, double x, double y)
{
    if (!Msh) return 0;
    if (Msh->NbrVer + 1 > Msh->NbrVerMax) {
        fprintf(stderr, "Not enough space in vertex arrays\n");
        return 0;
    }

    Msh->NbrVer++;
    Msh->Crd[Msh->NbrVer][0] = x;
    Msh->Crd[Msh->NbrVer][1] = y;

    return Msh->NbrVer;
}

int compute_barycentrics(Mesh* Msh, double x, double y, int iTri, double beta[3])
{
    int iVer1, iVer2, iVer3;
    double den;

    if(!Msh || !beta) return 0;
    if(iTri < 1 || iTri > Msh->NbrTri) return 0;
  
    iVer1 = Msh->Tri[iTri][0];
    iVer2 = Msh->Tri[iTri][1];
    iVer3 = Msh->Tri[iTri][2];

    double* A = Msh->Crd[iVer1];
    double* B = Msh->Crd[iVer2];
    double* C = Msh->Crd[iVer3];

    den = tri_area(A,B,C);

    //--- barycentric coordinates of (x,y) in triangle A-B-C
    beta[0] = tri_area((double[]){x,y}, B, C) / den;
    beta[1] = tri_area(A, (double[]){x,y}, C) / den;
    beta[2] = tri_area(A, B, (double[]){x,y}) / den;

    return 1;
}



int find_point_in_mesh(Mesh* Msh, double x, double y, int *iTriFound, double beta[3])
{
  int iTri;
  const double eps = 1e-12;

  if(!Msh || !iTriFound || !beta)
    return 0;

  *iTriFound = 0;
  beta[0] = beta[1] = beta[2] = 0.0;

  for(iTri=1; iTri<=Msh->NbrTri; iTri++){
    
    compute_barycentrics(Msh, x, y, iTri, beta);

    if(beta[0] >= -eps && beta[1] >= -eps && beta[2] >= -eps) {
      if(fabs(beta[0])< eps) beta[0] = 0.0;
      if(fabs(beta[1])< eps) beta[1] = 0.0;
      if(fabs(beta[2])< eps) beta[2] = 0.0;

      *iTriFound = iTri;
      return 1;
    }
  }
  return 0;
}

int find_point_walk(Mesh* Msh, double x, double y, int iTriStart, int *iTriFound, double beta[3])
{
  int iTri, step, kmin, jTri;
  const double eps = 1e-12;

  if(!Msh || !iTriFound || !beta || !Msh->TriVoi)
    return 0;

  if(Msh->NbrTri < 1) return 0;

  if(iTriStart < 1 || iTriStart > Msh->NbrTri) iTriStart = 1;

  iTri = iTriStart;

  for(step = 0; step < 2*Msh->NbrTri + 10; ++step){

      if(!compute_barycentrics(Msh, x, y, iTri, beta))
        return 0;

      if(beta[0] >= -eps && beta[1] >= -eps && beta[2] >= -eps) {
        if(fabs(beta[0])< eps) beta[0] = 0.0;
        if(fabs(beta[1])< eps) beta[1] = 0.0;
        if(fabs(beta[2])< eps) beta[2] = 0.0;
        *iTriFound = iTri;
        return 1;
      }

      kmin = 0;
      if(beta[1] < beta[kmin]) kmin = 1;
      if(beta[2] < beta[kmin]) kmin = 2;

      jTri = Msh->TriVoi[iTri][kmin];
      if(jTri == 0) return 0;

      iTri = jTri; 
  }
  return 0;
}

int msh_split_triangle_naive(Mesh* Msh, int iTri, int iP)
{
    int A, B, C;
    int iTri2, iTri3;

    if(!Msh) return 0;
    
    if(Msh->NbrTri + 2 > Msh->NbrTriMax)
    {
      fprintf(stderr, "## ERROR: not enough space in triangle arrays\n");
      return 0;
    }

    A = Msh->Tri[iTri][0];
    B = Msh->Tri[iTri][1];
    C = Msh->Tri[iTri][2];


    //---The new triangle slots at the end
    iTri2 = Msh->NbrTri + 1;
    iTri3 = Msh->NbrTri + 2;

    //--- Replace the old (ABC) with (ABP)
    Msh->Tri[iTri][0] = A;
    Msh->Tri[iTri][1] = B;
    Msh->Tri[iTri][2] = iP;

    //--- Add (BCP)
    Msh->Tri[iTri2][0] = B;
    Msh->Tri[iTri2][1] = C;
    Msh->Tri[iTri2][2] = iP;

    //--- Add (CAP)
    Msh->Tri[iTri3][0] = C;
    Msh->Tri[iTri3][1] = A;
    Msh->Tri[iTri3][2] = iP;

    Msh->NbrTri += 2;

    return 1;
}

int msh_split_edge_naive(Mesh* Msh, int iTri, int iEdg, int iP)
{
  int A, B, C, D;
  int jTri, jEdg;
  int iTri2, iTri3;

  if(!Msh) return 0;

  A = Msh->Tri[iTri][tri2edg[iEdg][0]];
  B = Msh->Tri[iTri][tri2edg[iEdg][1]];
  C = Msh->Tri[iTri][iEdg];

  jTri = 0;
  if(Msh->TriVoi)
    jTri = Msh->TriVoi[iTri][iEdg];


  //--- Boundary edge case
  if(jTri == 0){
    iTri2 = Msh->NbrTri + 1;

    if (Msh->NbrTri + 1 > Msh->NbrTriMax) {
      fprintf(stderr, "## ERROR: not enough space in triangle arrays\n");
      return 0;
    }

    //---Replace the old triangle (ABC) by (APC)
    Msh->Tri[iTri][0] = A;
    Msh->Tri[iTri][1] = iP;
    Msh->Tri[iTri][2] = C;

    //---Add (PBC)
    Msh->Tri[iTri2][0] = iP;
    Msh->Tri[iTri2][1] = B;
    Msh->Tri[iTri2][2] = C;

    Msh->NbrTri += 1;
    return 1;
  }

  //--- Interior edge case
  jEdg = -1;
  for(int k = 0; k < 3; k++){
    if(Msh->TriVoi[jTri][k] == iTri){
      jEdg = k;
      break;
    }
 }

  if(jEdg == -1){
    fprintf(stderr,"## ERROR: inconsistent TriVoi\n");
    return 0;
  }

  D = Msh->Tri[jTri][jEdg]; //vertex of jTri opposite of the common edge

  if (Msh->NbrTri + 2 > Msh->NbrTriMax) {
        fprintf(stderr, "## ERROR: not enough space in triangle arrays\n");
        return 0;
  }

  iTri2 = Msh->NbrTri + 1;
  iTri3 = Msh->NbrTri + 2;

  //---Overwrite iTri
  Msh->Tri[iTri][0] = A;
  Msh->Tri[iTri][1] = iP;
  Msh->Tri[iTri][2] = C;

  //---Overwrite jTri
  Msh->Tri[jTri][0] = iP;
  Msh->Tri[jTri][1] = B;
  Msh->Tri[jTri][2] = C;

  //---New appended triangles
  Msh->Tri[iTri2][0] = B;
  Msh->Tri[iTri2][1] = iP;
  Msh->Tri[iTri2][2] = D;

  Msh->Tri[iTri3][0] = iP;
  Msh->Tri[iTri3][1] = A;
  Msh->Tri[iTri3][2] = D;

  Msh->NbrTri += 2;

  return 1;
}


int msh_should_swap_edge_quality(Mesh* Msh, int iTri, int iEdg, int useQ2)
{
    const double eps = 1.e-12;
    int A, B, C, D;
    int jTri, jEdg;
    double *pA, *pB, *pC, *pD;
    double sABC, sABD, sCDA, sCDB;
    double before1, before2, after1, after2;
    double beforeWorst, afterWorst, beforeSum, afterSum;

    if(!Msh) return 0;
    if(!Msh->TriVoi || !Msh->Tri) return 0;
    if(iTri < 1 || iTri > Msh->NbrTri) return 0;
    if(iEdg < 0 || iEdg > 2) return 0;

    jTri = Msh->TriVoi[iTri][iEdg];
    if(jTri == 0) return 0;

    A = Msh->Tri[iTri][tri2edg[iEdg][0]];
    B = Msh->Tri[iTri][tri2edg[iEdg][1]];
    C = Msh->Tri[iTri][iEdg];

    jEdg = -1;
    for(int k = 0; k < 3; ++k){
      if(Msh->TriVoi[jTri][k] == iTri){
        jEdg = k;
        break;
      }
    }
    if(jEdg == -1) return 0;

    D = Msh->Tri[jTri][jEdg];

    pA = Msh->Crd[A];
    pB = Msh->Crd[B];
    pC = Msh->Crd[C];
    pD = Msh->Crd[D];

    // Swap is valid only if the quadrilateral is strictly convex.
    sABC = tri_area(pA, pB, pC);
    sABD = tri_area(pA, pB, pD);
    sCDA = tri_area(pC, pD, pA);
    sCDB = tri_area(pC, pD, pB);
    if(sABC * sABD >= -eps) return 0;
    if(sCDA * sCDB >= -eps) return 0;

    before1 = tri_quality_pts(pA, pB, pC, useQ2);
    before2 = tri_quality_pts(pB, pA, pD, useQ2);
    after1  = tri_quality_pts(pC, pD, pB, useQ2);
    after2  = tri_quality_pts(pD, pC, pA, useQ2);

    beforeWorst = fmax(before1, before2);
    afterWorst  = fmax(after1, after2);
    beforeSum   = before1 + before2;
    afterSum    = after1 + after2;

    if(afterWorst < beforeWorst - eps) return 1;
    if(fabs(afterWorst - beforeWorst) <= eps && afterSum < beforeSum - eps) return 1;

    return 0;
}


int msh_swap_edge(Mesh* Msh, int iTri, int iEdg)
{
    int A, B, C, D;
    int jTri, jEdg;

    if(!Msh) return 0;
    if(!Msh->TriVoi || !Msh->Tri) return 0;

    //--- Neighbor triangle across edge iEdg of iTri
    jTri = Msh->TriVoi[iTri][iEdg];
    if (jTri == 0) {
        fprintf(stderr, "## ERROR: cannot swap a boundary edge\n");
        return 0;
    }


    //--- iTri = (ABC) with edge iEdg = (AB)
    A = Msh->Tri[iTri][tri2edg[iEdg][0]];
    B = Msh->Tri[iTri][tri2edg[iEdg][1]];
    C = Msh->Tri[iTri][iEdg];

    jEdg = -1;
    for(int k = 0; k < 3; k++){
      if(Msh->TriVoi[jTri][k] == iTri){
        jEdg = k;
        break;
      }
    }

    if(jEdg == -1){
      fprintf(stderr,"## ERROR: inconsistent TriVoi\n");
      return 0;
    }

    //--- jTri = (BAD) with edge jEdg = (BA)
    D = Msh->Tri[jTri][jEdg]; //vertex of jTri opposite of the common edge

    //--- iTri becomes (CDB) with edge (CD)
    Msh->Tri[iTri][0] = C;
    Msh->Tri[iTri][1] = D;
    Msh->Tri[iTri][2] = B;

    //--- jTri becomes (DCA) with edge (DC)
    Msh->Tri[jTri][0] = D;
    Msh->Tri[jTri][1] = C;
    Msh->Tri[jTri][2] = A;

    tri_make_ccw(Msh, iTri);
    tri_make_ccw(Msh, jTri);

    if(!msh_neighbors(Msh)){
      fprintf(stderr,"## ERROR: failed to update neighbors after edge swap\n");
      return 0;
    }

    return 1;
}


int msh_build_cavity(Mesh* Msh, int iTriSeed, int iP,
                     int* cav, int* nCav,
                     BoundaryEdge* PilEdg, int* nPilEdg)
{
  IntStack* PilElt;
  int iTri, jTri;
  int iEdg;
  int a, b;

  if (!Msh || !cav || !nCav || !PilEdg || !nPilEdg) return 0;
  if (!Msh->Tri || !Msh->TriVoi || !Msh->TriMrk) return 0;
  if (iTriSeed < 1 || iTriSeed > Msh->NbrTri) return 0;
  if (iP < 1 || iP > Msh->NbrVer) return 0;

  *nCav = 0;
  *nPilEdg = 0;

  // TriMrk convention for cavity construction:
  // 0 = not visited, 1 = pushed in stack, 2 = accepted in cavity, 3 = visited but outside cavity.
  memset(Msh->TriMrk, 0, (Msh->NbrTriMax + 1) * sizeof(int1d));

  PilElt = initStack(16);
  if (!PilElt) return 0;

  if (!pushStack(PilElt, iTriSeed)) {
    freeStack(PilElt);
    return 0;
  }
  Msh->TriMrk[iTriSeed] = 1;

  for (iEdg = 0; iEdg < 3; ++iEdg) {
    a = Msh->Tri[iTriSeed][tri2edg[iEdg][0]];
    b = Msh->Tri[iTriSeed][tri2edg[iEdg][1]];
    if (!boundary_add_or_cancel(PilEdg, nPilEdg, a, b, iTriSeed, iEdg)) {
      freeStack(PilElt);
      return 0;
    }
  }

  while (sizeStack(PilElt) > 0) {
    if (!popStack(PilElt, &iTri)) {
      freeStack(PilElt);
      return 0;
    }

    if (Msh->TriMrk[iTri] == 2 || Msh->TriMrk[iTri] == 3)
      continue;

    if (!tri_contains_in_circumcircle(Msh, iTri, iP)) {
      Msh->TriMrk[iTri] = 3;
      continue;
    }

    Msh->TriMrk[iTri] = 2;
    cav[*nCav] = iTri;
    (*nCav)++;

    for (iEdg = 0; iEdg < 3; ++iEdg) {
      jTri = Msh->TriVoi[iTri][iEdg];
      if (jTri == 0) continue;

      if (Msh->TriMrk[jTri] == 0) {
        if (tri_contains_in_circumcircle(Msh, jTri, iP)) {
          if (!pushStack(PilElt, jTri)) {
            freeStack(PilElt);
            return 0;
          }
          Msh->TriMrk[jTri] = 1;

          if (!boundary_add_or_cancel(PilEdg, nPilEdg,
                                      Msh->Tri[iTri][tri2edg[iEdg][0]],
                                      Msh->Tri[iTri][tri2edg[iEdg][1]],
                                      iTri, iEdg)) {
            freeStack(PilElt);
            return 0;
          }

          for (int jEdg = 0; jEdg < 3; ++jEdg) {
            if (Msh->TriVoi[jTri][jEdg] == iTri)
              continue;

            a = Msh->Tri[jTri][tri2edg[jEdg][0]];
            b = Msh->Tri[jTri][tri2edg[jEdg][1]];
            if (!boundary_add_or_cancel(PilEdg, nPilEdg, a, b, jTri, jEdg)) {
              freeStack(PilElt);
              return 0;
            }
          }
        }
        else {
          Msh->TriMrk[jTri] = 3;
        }
      }
    }
  }

  freeStack(PilElt);

  return 1;
}


int msh_delete_cavity(Mesh* Msh, int* cav, int nCav)
{
  int k, iTri;

  if (!Msh || !cav) return 0;

  for (k = 0; k < nCav; ++k) {
    iTri = cav[k];

    if (iTri < 1 || iTri > Msh->NbrTri)
      return 0;

    if (Msh->Tri[iTri][0] > 0)
      Msh->Tri[iTri][0] = -Msh->Tri[iTri][0];
  }

  return 1;
}

int msh_star_cavity(Mesh* Msh, int iP,
                    BoundaryEdge* PilEdg, int nPilEdg,
                    int* cav, int nCav)
{
  int k, l, iTri;
  int nNewTri;
  int *newTri, *outTri, *outEdg;
  int iEdg, jEdg;
  int common[2], nCommon;
  int va[3], vb[3];

  if (!Msh || !PilEdg || !cav) return 0;
  if (iP < 1 || iP > Msh->NbrVer) return 0;

  nNewTri = nPilEdg;
  if (nNewTri < 1) return 0;

  if (Msh->NbrTri + (nNewTri - nCav) > Msh->NbrTriMax) {
    fprintf(stderr, "## ERROR: not enough space in triangle arrays\n");
    return 0;
  }

  newTri = malloc(nNewTri * sizeof(int));
  outTri = malloc(nNewTri * sizeof(int));
  outEdg = malloc(nNewTri * sizeof(int));
  if (!newTri || !outTri || !outEdg) {
    fprintf(stderr, "## ERROR: not enough space in local neighbor update arrays\n");
    free(newTri);
    free(outTri);
    free(outEdg);
    return 0;
  }

  for (k = 0; k < nNewTri; ++k) {
    if (k < nCav)
      newTri[k] = cav[k];
    else
      newTri[k] = Msh->NbrTri + (k - nCav) + 1;

    outTri[k] = Msh->TriVoi[PilEdg[k].iTri][PilEdg[k].iEdg];
    if (outTri[k] > 0)
      outEdg[k] = tri_find_edge_idx(Msh, outTri[k], PilEdg[k].a, PilEdg[k].b);
    else
      outEdg[k] = -1;
  }

  if (nNewTri > nCav)
    Msh->NbrTri += (nNewTri - nCav);

  for (k = 0; k < nNewTri; ++k) {
    iTri = newTri[k];

    Msh->Tri[iTri][0] = PilEdg[k].a;
    Msh->Tri[iTri][1] = PilEdg[k].b;
    Msh->Tri[iTri][2] = iP;

    tri_make_ccw(Msh, iTri);

    if (fabs(tri_area(Msh->Crd[Msh->Tri[iTri][0]],
                      Msh->Crd[Msh->Tri[iTri][1]],
                      Msh->Crd[Msh->Tri[iTri][2]])) < 1.e-30) {
      fprintf(stderr, "## ERROR: degenerate triangle during cavity starring\n");
      return 0;
    }

    if (Msh->TriRef)
      Msh->TriRef[iTri] = 0;
    if (Msh->TriMrk)
      Msh->TriMrk[iTri] = 0;

    if (Msh->TriVoi) {
      Msh->TriVoi[iTri][0] = 0;
      Msh->TriVoi[iTri][1] = 0;
      Msh->TriVoi[iTri][2] = 0;
    }
  }

  // Connect each new triangle to the outside mesh across its cavity-boundary edge.
  for (k = 0; k < nNewTri; ++k) {
    iTri = newTri[k];
    iEdg = tri_find_edge_idx(Msh, iTri, PilEdg[k].a, PilEdg[k].b);
    if (iEdg < 0) {
      fprintf(stderr, "## ERROR: failed to recover boundary edge in new triangle\n");
      free(newTri);
      free(outTri);
      free(outEdg);
      return 0;
    }

    Msh->TriVoi[iTri][iEdg] = outTri[k];
    if (outTri[k] > 0 && outEdg[k] >= 0)
      Msh->TriVoi[outTri[k]][outEdg[k]] = iTri;
  }

  // Stitch the new fan locally by matching shared edges.
  for (k = 0; k < nNewTri; ++k) {
    va[0] = Msh->Tri[newTri[k]][0];
    va[1] = Msh->Tri[newTri[k]][1];
    va[2] = Msh->Tri[newTri[k]][2];

    for (l = k + 1; l < nNewTri; ++l) {
      vb[0] = Msh->Tri[newTri[l]][0];
      vb[1] = Msh->Tri[newTri[l]][1];
      vb[2] = Msh->Tri[newTri[l]][2];

      nCommon = 0;
      for (int ia = 0; ia < 3; ++ia) {
        for (int ib = 0; ib < 3; ++ib) {
          if (va[ia] == vb[ib]) {
            if (nCommon < 2)
              common[nCommon] = va[ia];
            nCommon++;
          }
        }
      }

      if (nCommon != 2)
        continue;

      jEdg = tri_find_edge_idx(Msh, newTri[l], common[0], common[1]);
      iEdg = tri_find_edge_idx(Msh, newTri[k], common[0], common[1]);
      if (iEdg < 0 || jEdg < 0) {
        fprintf(stderr, "## ERROR: failed to stitch local cavity neighbors\n");
        free(newTri);
        free(outTri);
        free(outEdg);
        return 0;
      }

      Msh->TriVoi[newTri[k]][iEdg] = newTri[l];
      Msh->TriVoi[newTri[l]][jEdg] = newTri[k];
    }
  }

  free(newTri);
  free(outTri);
  free(outEdg);

  return 1;
}

int msh_insert_point_delaunay(Mesh* Msh, double x, double y)
{
  int iTriSeed;
  int iP;
  int* cav;
  int nCav;
  BoundaryEdge* PilEdg;
  int nPilEdg;
  double beta[3];
  const double eps = 1.e-12;

  if (!Msh) return 0;
  if (!Msh->Tri || !Msh->Crd) return 0;

  if (!Msh->TriVoi && !msh_neighbors(Msh)) {
    fprintf(stderr, "## ERROR: msh_neighbors failed to initialize Delaunay insertion\n");
    return 0;
  }

  if (!find_point_walk(Msh, x, y, 1, &iTriSeed, beta)) {
    if (!find_point_in_mesh(Msh, x, y, &iTriSeed, beta)) {
      fprintf(stderr, "## ERROR: point (%.6f, %.6f) is outside the mesh\n", x, y);
      return 0;
    }
  }

  // Keep the first version simple: do not insert points already on an edge/vertex.
  for (int k = 0; k < 3; ++k) {
    if (fabs(beta[k]) < eps) {
      fprintf(stderr, "## ERROR: boundary/edge point insertion not handled yet\n");
      return 0;
    }
  }

  iP = msh_add_vertex(Msh, x, y);
  if (iP == 0) {
    fprintf(stderr, "## ERROR: failed to add Delaunay insertion vertex\n");
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

  free(cav);
  free(PilEdg);

  return iP;
}
