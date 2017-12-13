//*****************************************************************************
// TITLE:         Naive Triangulation and Tetrahedralization 
// DESCRIPTION:   Functions of compute Delaunay Triangulation and Tetrahedralization
// AUTHORS:       Yaohong Wu, Jin Lin
// EDITS:         11/28/2017
//*****************************************************************************


boolean marked[][];  // mark whether edge (i,j) has been created
int[] idP, idQ;  // give each vertex in P or Q a unique ID

void initMarked(int nSite){
  marked = new boolean[nSite][nSite];
  for (int i = 0; i < nSite; ++i)
    for (int j = 0; j < nSite; ++j) marked[i][j] = false;
}

void initPQID(){
  idP = new int[P.nv];
  idQ = new int[Q.nv];
  for (int i = 0; i < P.nv; ++i) idP[i] = i;
  for (int i = 0; i < Q.nv; ++i) idQ[i] = i + P.nv;
}

void triMeshNaive(pt[] v, int[] id, int n, ArrayList<Triangle> tris, ArrayList<Edge> edges){
  tris.clear();  // tris will keep growing without this line
  edges.clear();  // edges will keep growing without this line
  
  for(int i = 0; i < n-2; ++i){
    for(int j = i + 1; j < n-1; ++j){
      for(int k = j + 1; k < n; ++k){
        pt E = circumcenter3Pts(v[i], v[j], v[k]);
        float r = d(E, v[i]);
        boolean valid = true;
        for(int t = 0; t < n; ++t){
          if(t == i || t == j || t == k) continue;
          if(d(E, v[t]) < r){
            valid = false;
            break;
          }
        }
        if(valid){
          int a = id[i], b = id[j], c = id[k];
          tris.add(new Triangle(a,b,c));  
          if(marked[a][b] == false) {edges.add(new Edge(a, b)); marked[a][b] = marked[b][a] = true;}
          if(marked[a][c] == false) {edges.add(new Edge(a, c)); marked[a][c] = marked[c][a] = true;}
          if(marked[b][c] == false) {edges.add(new Edge(b, c)); marked[b][c] = marked[c][b] = true;}
        }
      }
    }
  }
}


void tetMeshNaive31(ArrayList<Triangle> tris, pt[] v, int[] id, int n, ArrayList<Edge> outEdges){
  outEdges.clear();
  for(int i = 0; i < tris.size(); ++i){
    int a = tris.get(i).a, b = tris.get(i).b, c = tris.get(i).c;
    pt A, B, C;
    A = a < P.nv ? P.G[a] : Q.G[a-P.nv];
    B = b < P.nv ? P.G[b] : Q.G[b-P.nv];
    C = c < P.nv ? P.G[c] : Q.G[c-P.nv];
    
    float minBulge = MAX_FLOAT;
    int minIdx = -1;
    for(int j = 0; j < n; ++j){
      pt D = v[j];
      float bulge = findBulge(A, B, C, D);
      if(bulge < minBulge){
        minBulge = bulge;
        minIdx = j;
      }
    }
    
    int d = id[minIdx];
    if(marked[a][d] == false) {outEdges.add(new Edge(a, d)); marked[a][d] = marked[d][a] = true;}
    if(marked[b][d] == false) {outEdges.add(new Edge(b, d)); marked[b][d] = marked[d][b] = true;}
    if(marked[c][d] == false) {outEdges.add(new Edge(c, d)); marked[c][d] = marked[d][c] = true;}
  }
}


void tetMeshNaive22(ArrayList<Edge> inEdges1, ArrayList<Edge> inEdges2, ArrayList<Edge> outEdges){
  outEdges.clear();
  pt[] v1 = P.G, v2 = Q.G;
  int n1 = P.nv, n2 = Q.nv;
  for(int i = 0; i < inEdges1.size(); ++i){  // edges on the floor
    int a = inEdges1.get(i).a, b = inEdges1.get(i).b;
    pt A, B;
    A = a < P.nv ? P.G[a] : Q.G[a-P.nv];
    B = b < P.nv ? P.G[b] : Q.G[b-P.nv];
    for(int j = 0; j < inEdges2.size(); ++j){  // edges on the ceiling
      int c = inEdges2.get(j).a, d = inEdges2.get(j).b;
      pt C, D;
      C = c < P.nv ? P.G[c] : Q.G[c-P.nv];
      D = d < P.nv ? P.G[d] : Q.G[d-P.nv];
      
      pt F = circumcenter4Pts(A, B, C, D);
      float r = d(F, A);
      
      boolean valid = true;
      for(int k = 0; k < n1; ++k){
        
        if(k == a || k == b) continue;
        if(d(v1[k], F) < r) { 
          valid = false; 
          break; 
        }
      }
      if(valid == false) continue;    
      for(int k = 0; k < n2; ++k){
        
        if(k + P.nv == c || k + P.nv == d) continue;
        if(d(v2[k], F) < r) { 
          valid = false; 
          break; 
        }
      }
      
      if(valid == false) continue;
      tetCnt++;
      if(marked[a][c] == false) {outEdges.add(new Edge(a, c)); marked[a][c] = marked[c][a] = true;}
      if(marked[a][d] == false) {outEdges.add(new Edge(a, d)); marked[a][d] = marked[d][a] = true;}
      if(marked[b][c] == false) {outEdges.add(new Edge(b, c)); marked[b][c] = marked[c][b] = true;}
      if(marked[b][d] == false) {outEdges.add(new Edge(b, d)); marked[b][d] = marked[d][b] = true;}
    }
  }
}

pt circumcenter3Pts(pt A, pt B, pt C){
  pt M = P(A, B);
  vec AB = V(A, B), AC = V(A, C), BC = V(B, C);
  vec N = cross( cross(AB, AC), AB);
  float t = dot(AC, BC) / (2 * dot(AC, N));
  pt E = P(M, t, N);
  return E;
}


pt circumcenter4Pts(pt A, pt B, pt C, pt D){
  pt E = circumcenter3Pts(A, B, C);
  vec AB = V(A, B), AC = V(A, C), AD = V(A, D), EA = V(E, A), ED = V(E, D);
  vec N = cross(AB, AC);
  float t = dot(A(EA, ED), AD) / (2 * dot(N, AD));
  pt F = P(E, t, N);
  return F;
}


float findBulge(pt A, pt B, pt C, pt D){
  pt P = circumcenter4Pts(A, B, C, D);
  float r = d(P, A);
  vec AB = V(A, B), AC = V(A, C), AD = V(A, D), AP = V(A, P);
  vec N = cross(AB, AC);
  float dAPN = dot(AP, N), dADN = dot(AD, N);
  boolean sameSide = (dAPN * dADN > 0);
  float h = dAPN / norm(N);
  if(h < 0) h = -h;
  float b = 0;
  if(sameSide) b = r + h;
  else b = r - h;
  return b;
}


void drawEdges(ArrayList<Edge> edges, float r){
  for(int i = 0; i < edges.size(); ++i){
    int a = edges.get(i).a;
    int b = edges.get(i).b;
    pt A, B;
    A = a < P.nv ? P.G[a] : Q.G[a-P.nv];
    B = b < P.nv ? P.G[b] : Q.G[b-P.nv];
    beam(A, B, r);
  }
}


void constructTetMesh(){
  triMeshNaive(P.G, idP, P.nv, trisP, edgesP);
  triMeshNaive(Q.G, idQ, Q.nv, trisQ, edgesQ);

  tetMeshNaive31(trisP, Q.G, idQ, Q.nv, edges31);
  tetMeshNaive31(trisQ, P.G, idP, P.nv, edges13);
  tetMeshNaive22(edgesP, edgesQ, edges22);  //increment tetCnt inside this function

  tetCnt += trisP.size() + trisQ.size();
  mixBeamCnt += edges31.size() + edges13.size() + edges22.size();
}