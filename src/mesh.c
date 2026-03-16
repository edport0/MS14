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

  Msh->NbrVerMax = Msh->NbrVer;
  Msh->NbrTriMax = Msh->NbrTri;

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


int msh_neighborsQ2(Mesh* Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  memset(Msh->TriVoi, 0, (Msh->NbrTri + 1) * sizeof(int3d));

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

          if(((iVer1 == jVer1) && (iVer2 == jVer2)) || ((iVer1 = jVer2) && (iVer2 == jVer1))){
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
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  memset(Msh->TriVoi, 0, sizeof(int3d) * (Msh->NbrTri + 1)); // reset!

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


int msh_histogram(double* Qal, int NbrTri, int NbrBins, double QMax)
{
  int    i, iTri;
  int*   bins  = calloc(NbrBins, sizeof(int));
  double width = (QMax - 1.0) / NbrBins;

  for (iTri = 1; iTri <= NbrTri; iTri++) {
    int b = (int)((Qal[iTri] - 1.0) / width); // shift by 1 since min quality = 1
    if (b >= NbrBins) b = NbrBins - 1;         // clamp outliers to last bin
    bins[b]++;
  }

  printf("\n--- Quality histogram (NbrBins=%d, QMax=%.2f) ---\n", NbrBins, QMax);
  for (i = 0; i < NbrBins; i++) {
    printf("  [%.4f - %.4f] : %d\n", 1.0 + i*width, 1.0 + (i+1)*width, bins[i]);
  }

  free(bins);
  return 1;
}

