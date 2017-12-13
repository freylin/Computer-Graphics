//*****************************************************************************
// TITLE:         Playground 
// DESCRIPTION:   Functions of debugging BPA, Octree, radius computing
// AUTHORS:       Yaohong Wu, Jin Lin
// EDITS:         11/28/2017
//*****************************************************************************


// compute the incenter of triangle ABC
pt incenter(pt A, pt B, pt C){
  float a = d(B, C), b = d(A, C), c = d(A, B);
  float p = a + b + c;
  return P(a/p, A, b/p, B, c/p, C);
}

// following 7 functions are for debugging BPA

void BPAOnFloor(){
  vertices.clear();
  for(int i = 0; i < P.nv; ++i){
    vertices.add(new Vertex(P.G[i], V(0,0,1), 0, i));
  }
  if(useOctree) root = buildOctree(vertices, r);
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges);
}


void BPAOnTerrain(){  // assum that P is the point cloud on the floor and Q is the point cloud on the ceiling
  vertices.clear();
  for(int i = 0; i < P.nv; ++i) {
    vertices.add(new Vertex(P.G[i], V(0,0,1), 0, i));
  }
  for(int i = 0; i < Q.nv; ++i){
    vertices.add(new Vertex(Q.G[i], V(0,0,1), 0, i+P.nv));
  }
  
  if(useOctree) root = buildOctree(vertices, r);
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges); 
}

void BPAOnCylinder(){
  vertices.clear();
  pt A = P(1000.0, 900.0, 0.0);
  pt B = P(1000.0, 900.0, 600.0);
  float rAB = 100;
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();
  sampleBeam(A, B, rAB, G, N);
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  if(showVertices){
    fill(pink); showSamples(vertices);
  }
  
  //int i = nBottom * (nSide + 1) / 2;
  //int a = i, b = i + 1, c = i + 1 + nBottom;
  //pt Pa = G.get(a), Pb = G.get(b), Pc = G.get(c);
  //fill(green); beam(Pa, Pb, 3);
  //pt M = P(Pa, Pb);
  //fill(blue); beam(M, Pc, 3);
  
  if(useOctree) root = buildOctree(vertices, r);
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges); 
}

void BPAOnSphere(){
  vertices.clear();
  pt C = P(1200.0, 800.0, 500.0);
  float rC = 300;
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();  
  sampleBall(C, rC, G, N);
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  if(showVertices){
    fill(pink); showSamples(vertices);
  }
  
  //int i = nLong * (nLat + 1) / 2;
  //int a = i, b = i + 1, c = i + nLong;
  //pt Pa = G.get(a), Pb = G.get(b), Pc = G.get(c);
  //fill(green); beam(Pa, Pb, 3);
  //pt M = P(Pa, Pb);
  //fill(blue); beam(M, Pc, 3);
  
  if(useOctree) root = buildOctree(vertices, r);
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges);
  
  //removeTriangleHole(borderEdges, innerEdges, mesh, vertices);
}


void BPAOnVShape(){
  vertices.clear();
  
  // two points on the floor, one point on the ceiling
  pt A = P(1000.0, 1000.0, 0.0);
  pt B = P(300.0, 700.0, 0.0);
  pt C = P(1000.0, 700.0, 600.0);
  
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();  
  
  seedIdx = nLong * (nLat + 1) + nBottom * (nSide + 1) / 2;
  
  sampleBall(A, rb, G, N);  // sample A
  sampleBeam(A, C, rt, G, N);  // sample AC
  sampleBall(C, rb, G, N);  // sample C
  sampleBeam(C, B, rt, G, N); // sample CB
  sampleBall(B, rb, G, N);  //sample B
  
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  if(showVertices){
    fill(pink); showSamples(vertices);
  }
  
  if(useRadius){
    r = getRadiusOfPivotingBall();
    println("the computed r is " + r);
  }
  
  int st = millis();
  if(useOctree) {
    root = buildOctree(vertices, r);
    //fill(red); drawCell(root, 5.0, true);
  }
  
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges);
  int ed = millis();
  durationBPA = (ed - st) / 1000.0;
}


void BPAOnTet(){
  vertices.clear();
  assert P.nv >= 3 && Q.nv >= 1;
  pt A = P.G[0], B = P.G[1], C = P.G[2], D = Q.G[0];  // 3 vertices on the floor and 1 vertex on the ceiling

  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();  
  
  seedIdx = nBottom * (nSide + 1) / 2;
  
  sampleBeam(A, B, rt, G, N);
  sampleBeam(B, C, rt, G, N);
  sampleBeam(C, A, rt, G, N);
  sampleBeam(A, D, rt, G, N);
  sampleBeam(D, B, rt, G, N);
  sampleBeam(D, C, rt, G, N);
  
  sampleBall(A, rb, G, N);
  sampleBall(B, rb, G, N);
  sampleBall(C, rb, G, N);
  sampleBall(D, rb, G, N);
  
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }

  if(showVertices){
    fill(pink); showSamples(vertices);
  }
  
  if(useRadius){
    r = getRadiusOfPivotingBall();
    println("the computed r is " + r);
  }
  
  int st = millis();
  if(useOctree) {
    root = buildOctree(vertices, r);
    //fill(red); drawCell(root, 5.0, true);
  }
  
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges);
  removeTriangleHole(borderEdges, innerEdges, mesh, vertices);
  
  int ed = millis();
  durationBPA = (ed - st) / 1000.0;
}


void BPAOnTetMesh(){
  vertices.clear();
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();
  
  ArrayList<Edge> allEdges = new ArrayList<Edge>();
  
  seedIdx = nBottom * (nSide + 1) / 2;
  
  allEdges.addAll(edgesP);
  allEdges.addAll(edgesQ);
  allEdges.addAll(edges31);
  allEdges.addAll(edges13);
  allEdges.addAll(edges22); 
  sampleMultiBeams(allEdges, rt, G, N);
  
  ArrayList<pt> allCenters = new ArrayList<pt>();
  for(int i = 0; i < P.nv; ++i) allCenters.add(P.G[i]);
  for(int i = 0; i < Q.nv; ++i) allCenters.add(Q.G[i]);
  sampleMultiBalls(allCenters, rb, G, N);
  
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  if(showVertices){
    fill(pink); showSamples(vertices);
  }
  
  if(useRadius){
    r = getRadiusOfPivotingBall();
    println("the computed r is " + r);
  }  
  
  int st = millis();
  if(useOctree) root = buildOctree(vertices, r);
  
  reconstructSurfaceBPA(vertices, r, mesh, borderEdges, innerEdges);
  removeTriangleHole(borderEdges, innerEdges, mesh, vertices);
  
  int ed = millis();
  durationBPA = (ed - st) / 1000.0;
  
  tm = new TriangleMesh();
  tm.setTo(vertices, mesh);
  tm.saveTriMesh(triMeshPathS);
  println("finish saving mesh");
  println("duration of BPA: " + durationBPA);
}


// Following 5 functions are for debugging Octree

void debugOctreeFloor(){
  vertices.clear();
  for(int i = 0; i < P.nv; ++i){
    vertices.add(new Vertex(P.G[i], V(0,0,1), 0, i));
  }
  
  root = buildOctree(vertices, r);
  fill(red); drawCell(root, 5.0, true);
}

void debugOctreeTerrain(){
  vertices.clear();
  for(int i = 0; i < P.nv; ++i) {
    vertices.add(new Vertex(P.G[i], V(0,0,1), 0, i));
  }
  
  for(int i = 0; i < Q.nv; ++i){
    vertices.add(new Vertex(Q.G[i], V(0,0,1), 0, i+P.nv));
  }
  
  root = buildOctree(vertices, r);
  fill(red); drawCell(root, 5.0, true);
  
  pt p = R.Picked();
  float rn = 100;
  getNeighbors(p, root, rn);
  fill(blue, 100); show(p, rn);
  
}

void debugOctreeCylinder(){
  vertices.clear();
  pt A = P(1000.0, 900.0, 0.0);
  pt B = P(1000.0, 900.0, 300.0);
  float rAB = 100;
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();
  sampleBeam(A, B, rAB, G, N);
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  showSamples(vertices);
  root = buildOctree(vertices, r);
  fill(red); drawCell(root, 5.0, true);
}

void debugOctreeSphere(){
  vertices.clear();
  pt C = P(1200.0, 800.0, 500.0);
  float rC = 300;
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();  
  sampleBall(C, rC, G, N);
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  showSamples(vertices);
  root = buildOctree(vertices, r);
  fill(red); drawCell(root, 5.0, true);
}


void debugOctreeVShape(){
  vertices.clear();
  // two points on the floor, one point on the ceiling
  pt A = P(1000.0, 800.0, 0.0);
  pt B = P(700.0, 1200.0, 0.0);
  pt C = P(900.0, 900.0, 500.0);
  
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();  
  
  sampleBall(A, rb, G, N);  // sample A
  sampleBeam(A, C, rt, G, N);  // sample AC
  sampleBall(C, rb, G, N);  // sample C
  sampleBeam(C, B, rt, G, N); // sample CB
  sampleBall(B, rb, G, N);  //sample B
  
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  fill(pink); showSamples(vertices);
  root = buildOctree(vertices, r);
  fill(red); drawCell(root, 5.0, true);
  
  // find neighbors
  int n = vertices.size();
  int pickId = n/2;
  pt p = vertices.get(pickId).position;
  ArrayList<Vertex> neighbors = getNeighbors(p, root, r);
  
  //println("number of neighbors = " + neighbors.size());
  fill(green); showSamples(neighbors);
  
  fill(blue, 100); show(p, r);
}



// following 2 functions are for debugging radius computing

void debugRadiusVShape(){
  vertices.clear();
  // two points on the floor, one point on the ceiling
  pt A = P(1000.0, 800.0, 0.0);
  pt B = P(700.0, 1200.0, 0.0);
  pt C = P(900.0, 900.0, 500.0);
  
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();  
  
  sampleBall(A, rb, G, N);  // sample A
  sampleBeam(A, C, rt, G, N);  // sample AC
  sampleBall(C, rb, G, N);  // sample C
  sampleBeam(C, B, rt, G, N); // sample CB
  sampleBall(B, rb, G, N);  //sample B
  
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  fill(pink); showSamples(vertices);
  
  println("max beam length = " + maxBeamLength);
  float rPivot = getRadiusOfPivotingBall();
  println("r for pivoting ball = " + rPivot);
  
  // show some balls with radius rPivot centered at random vertices
  fill(blue, 100);
  randomSeed(0);
  float nv = (float)vertices.size();
  for(int i = 0; i < 10; ++i){
    int id = (int)random(nv);
    pt p = vertices.get(id).position;
    show(p, rPivot);
  }
}

void debugRadiusTetMesh(){
  vertices.clear();
  ArrayList<pt> G = new ArrayList<pt>();
  ArrayList<vec> N = new ArrayList<vec>();
  
  ArrayList<Edge> allEdges = new ArrayList<Edge>();
  allEdges.addAll(edgesP);
  allEdges.addAll(edgesQ);
  allEdges.addAll(edges31);
  allEdges.addAll(edges13);
  allEdges.addAll(edges22); 
  sampleMultiBeams(allEdges, rt, G, N);
  
  ArrayList<pt> allCenters = new ArrayList<pt>();
  for(int i = 0; i < P.nv; ++i) allCenters.add(P.G[i]);
  for(int i = 0; i < Q.nv; ++i) allCenters.add(Q.G[i]);
  sampleMultiBalls(allCenters, rb, G, N);
  
  for(int i = 0; i < G.size(); ++i){
    vertices.add(new Vertex(G.get(i), N.get(i), 0, i));
  }
  
  fill(pink); showSamples(vertices);
  
  println("max beam length = " + maxBeamLength);
  float rPivot = getRadiusOfPivotingBall();
  println("r for pivoting ball = " + rPivot);
  
  // show some balls with radius rPivot centered at random vertices
  fill(blue, 100);
  randomSeed(0);
  float nv = (float)vertices.size();
  for(int i = 0; i < 20; ++i){
    int id = (int)random(nv);
    pt p = vertices.get(id).position;
    show(p, rPivot);
  }
}