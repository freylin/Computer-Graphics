//*****************************************************************************
// TITLE:         Ball Pivoting Algorithm
// DESCRIPTION:   Functions of finding a seed triangle, pivoting a ball around an edge, expanding triangle mesh, etc.
// AUTHORS:       Yaohong Wu, Jin Lin
// EDITS:         11/28/2017
//*****************************************************************************

// Get the center of the ball touching A, B, C with radius r.
// Note that it should be on the direction of N which is the outward normal of triangle ABC.
pt findBallCenter(pt A, pt B, pt C, vec N, float r){
  if(n2(N) > 1 + 10 * EPSILON) N.normalize();
  pt E = circumcenter3Pts(A, B, C);
  float h2 = r * r - n2(V(E, A));
  if(h2 < 0) return null;
  float h = sqrt(h2);
  pt F = P(E, h, N);
  return F;
}

 //<>//
// check whether a ball is empty
boolean isEmptyBall(pt F, float r, ArrayList<Vertex> P, int i, int j, int k){
  ArrayList<Vertex> neighbors = P;  
  if(root != null){
    neighbors = getNeighbors(F, root, r);
  }
  for(Vertex p : neighbors){
    if(p.id == i || p.id == j || p.id == k) continue;
    if(d(p.position, F) < r) return false;
  }
  return true;
}


// Find a seed triangle to start with.
// P is the point cloud
// r is the radius of the ball we use to pivot
// front is a queue of front edges. After we find a seed triangle, we push three edges of it into front.
Triangle findSeedTriangle(ArrayList<Vertex> P, float r, Queue<taggedEdge> front){
  int nv = P.size();
  int st = (seedIdx >= 0 ? seedIdx : 0);
  for(int i = st; i < nv - 2; ++i){
    if(P.get(i).type != 0) continue;
    for(int j = i + 1; j < nv - 1; ++j){
      if(P.get(j).type != 0) continue;
      for(int k = j + 1; k < nv; ++k){
        if(P.get(k).type != 0) continue;
        
        pt A = P.get(i).position;
        pt B = P.get(j).position;
        pt C = P.get(k).position;
        
        vec NA = P.get(i).normal;
        vec NB = P.get(j).normal;
        vec NC = P.get(k).normal;
        
        vec NABC = U(N(V(A,B), V(A,C)));  // |NABC| = 1
        
        boolean cw = true;  // assume that when A, B, C are in clockwise order, ABxAC facing outwards
        float da = dot(NA,NABC), db = dot(NB, NABC), dc = dot(NC, NABC);
        boolean allPos = (da > 0) && (db > 0) && (dc > 0);
        boolean allNeg = (da < 0) && (db < 0) && (dc < 0);
        if( (!allPos) && (!allNeg) ) continue;  // A, B, C are not compatible
        
        if( allNeg ) {NABC.rev(); cw = false;}  // reverse the order of A, B, C and direction of normal
        
        pt F = findBallCenter(A, B, C, NABC, r);  // if the ball is too small w.r.t. triangle ABC, then F will be null
        if(F == null) continue;
         
        boolean isEmpty = isEmptyBall(F, r, P, i, j, k);  // check whether the ball is empty
        if(!isEmpty) continue;
        
        Vertex va = P.get(i), vb = P.get(j), vc = P.get(k);
        va.type = vb.type = vc.type = 1;  // vertices on front
        
        // show an arrow from the incenter of triangle ABC to F
        if(showBallPath){
          pt H = incenter(va.position, vb.position, vc.position);
          fill(magenta); arrow(H, F, rt);
        }
        
        taggedEdge e1, e2, e3;
        if(cw) {
          e1 = new taggedEdge(i, j, k, F, 0);  // 0 means front
          e2 = new taggedEdge(j, k, i, F, 0);
          e3 = new taggedEdge(k, i, j, F, 0);
          front.add(e1);
          front.add(e2);
          front.add(e3);
          va.outEdges.put(vb, e1);  // update va's adjacent edges
          vb.outEdges.put(vc, e2);
          vc.outEdges.put(va, e3);
          return new Triangle(i, j, k);
        }
        else{
          e1 = new taggedEdge(j, i, k, F, 0);
          e2 = new taggedEdge(k, j, i, F, 0);
          e3 = new taggedEdge(i, k, j, F, 0);
          front.add(e1);
          front.add(e2);
          front.add(e3);
          va.outEdges.put(vc, e3);
          vb.outEdges.put(va, e1);
          vc.outEdges.put(vb, e2);
          return new Triangle(j, i, k);
        }
      }
    }
  }
  return null;
}


// Pivot a ball around edge e
// e is the edge the ball pivots around
// P is the point cloud
// Fd is the new ball center we will compute (if any)
// return the vertex index of the cloest valid hit
int pivot(taggedEdge e, ArrayList<Vertex> P, pt Fd){
  int a = e.e.a, b = e.e.b, c = e.c, d = -1;
  
  pt A = P.get(a).position, B = P.get(b).position;
  pt M = P(A, B);
  pt F = e.F;
  
  vec MF = V(M, F);
  vec NA = P.get(a).normal, NB = P.get(b).normal;
  vec AB = V(A, B);
  float w = 7;  // store the smallest pivoting angle, initially big enough, i.e. > 2*PI
  
  float range = r + d(M, F);
  ArrayList<Vertex> neighbors = P;
  
  if(root != null){ neighbors = getNeighbors(M, root, range); }
  
  for(Vertex v: neighbors){  // only check neighbors of M
    int i = v.id;
    if(i == a || i == b || i == c) continue;  // a new point must be different from a, b, c
    
    if(v.type == 2) continue;  // avoid hitting inner vertex, very important
                               // without this rejection test, the program may run a long long time
                               // remember to update vertex type when doing mesh expansion
    
    pt D = v.position;
    vec ND = v.normal;
    vec AD = V(A, D);
    vec NABD = N(AD, AB);  // ADxAB, the order matters!
    
    float da = dot(NA,NABD), db = dot(NB, NABD), dd = dot(ND, NABD);
    if(!((da > 0) && (db > 0) && (dd > 0))) continue; // A, B, D are not compatible
    
    //if(!(((da > 0) && (db > 0) && (dd > 0)) || ((da < 0) && (db < 0) && (dd < 0)))) continue;  
    //if((da < 0) && (db < 0) && (dd < 0)) NABD.rev();
    
    NABD.normalize();
    
    pt Ftmp = findBallCenter(A, B, D, NABD, r);  // if the ball is too small w.r.t. triangle ABD, then Ftmp will be null
    if(Ftmp == null) continue;
    
    vec MFtmp = V(M, Ftmp);
    float wd = angle(MF, MFtmp);
    if(dot(N(MF, MFtmp), AB) < 0) wd = TWO_PI - wd;  // remove sign ambiguity, the ball having the smallest angle may not be empty
    
    if(wd > w) continue;
    
    // Should we test if the new ball is empty, or just pick the one with smallest pivot angle? Yes.
    // Only need to check the ball corresponding to the point with smallest pivot angle? No.

    boolean isEmpty = isEmptyBall(Ftmp, r, P, a, b, i);
    if(!isEmpty) continue;
    
    w = wd;
    d = i;
    Fd.setTo(Ftmp);  // modify Fd using setTo(), do not use assignment, i.e. Fd = Ftmp
  }
  
  return d;
}


// reconstruct surface using BPA
// P is the point cloud
// r is the radius of the pivoting ball
// mesh is the collection of triangles we will compute
// borderEdges is the collection of border edges of our triangle mesh
// innerEdges is the collection of inner edges of our triangle mesh
void reconstructSurfaceBPA(ArrayList<Vertex> P, float r, ArrayList<Triangle> mesh, ArrayList<Edge> borderEdges, ArrayList<Edge> innerEdges){
  mesh.clear(); borderEdges.clear(); innerEdges.clear();
  
  Queue<taggedEdge> frontList = new LinkedList<taggedEdge>();
  Triangle T = findSeedTriangle(P, r, frontList);
  
  if(T == null){
    println("No seed triangle found! Please try larger radius for ball pivoting!");
    return;
  }
  mesh.add(T);
  
  // show the first triangle
  //fill(blue);
  //beginShape();
  //  vertex(P.get(T.a).position); vertex(P.get(T.b).position); vertex(P.get(T.c).position);
  //endShape();
  
  long steps = 0;
  while(frontList.size() > 0){   
    taggedEdge e = frontList.poll();  // get and remove the first element
    if(e.tag != 0) continue;  // e is not a front edge
    
    pt Fd = new pt();
    int d = pivot(e, P, Fd);
    if(d < 0) {  // no hit found
      e.tag = 1;  // mark e as a boundary edge
      borderEdges.add(e.e);
      continue;
    }
    
    int a = e.e.a, b = e.e.b;
    Vertex va = P.get(a);
    Vertex vb = P.get(b);
    Vertex vd = P.get(d);
    
    mesh.add(new Triangle(a, d, b));  // create a new triangle facet
    e.tag = 2;  // the front edge becomes an inner edge
    innerEdges.add(e.e);
    va.outEdges.remove(vb);  // remove edge ab 
    vb.inEdges.remove(va);  // remove edge ab
    
    
    // show ball path
    if(showBallPath){
      pt F = e.F;
      pt Hd = incenter(P.get(a).position, P.get(b).position, P.get(d).position);
      fill(yellow); arrow(F, Fd, rt);
      fill(cyan); arrow(Hd, Fd, rt);
      
      if(steps == stepToShow){
        ballCenterToShow.setTo(Fd);
      }  
    }
    

    // Is there a front edge connecting va and vd?
    taggedEdge e1 = null;
    if(vd.outEdges.containsKey(va)) {
      e1 = vd.outEdges.get(va); 
      vd.outEdges.remove(va);
      va.inEdges.remove(vd);
    }
    else if(va.outEdges.containsKey(vd)) {  // edge ad already exists?
      //println("edge ad already exists");
      e1 = va.outEdges.get(vd);
      va.outEdges.remove(vd);  // remove edge ad and then add it back?
      vd.inEdges.remove(va);
    }
    
    if(e1 != null){  // mark e1 as inner edge
      e1.tag = 2;
      innerEdges.add(e1.e);
      
      // va may become inner
      if( va.outEdges.isEmpty() && va.inEdges.isEmpty() ){
        va.type = 2;
      }
      
    }
    else{  // mark e1 as front edge and push it to the front list
      e1 = new taggedEdge(a, d, b, Fd, 0);  // create an edge from a to d
      va.outEdges.put(vd, e1);  // update the adjacent edges starting from a
      vd.inEdges.put(va, e1);
      frontList.add(e1);  
    }
    
    // is there a front edge connecting vd and vb?
    taggedEdge e2 = null;
    if(vb.outEdges.containsKey(vd)) {
      e2 = vb.outEdges.get(vd); 
      vb.outEdges.remove(vd);
      vd.inEdges.remove(vb);
    }
    else if(vd.outEdges.containsKey(vb)) {
      //println("edge db already exists");
      e2 = vd.outEdges.get(vb); 
      vd.outEdges.remove(vb);
      vb.inEdges.remove(vd);
    }
    
    if(e2 != null){  // mark e2 as inner edge
      e2.tag = 2;
      innerEdges.add(e2.e);
      
      // vb may become inner
      if(vb.outEdges.isEmpty() && vb.inEdges.isEmpty()) vb.type = 2;
    }
    else{  // mark e2 as front edge and push it to the front list
      e2 = new taggedEdge(d, b, a, Fd, 0);  // create an edge from d to b
      vd.outEdges.put(vb, e2);  // update the adjacent edges starting from d    
      vb.inEdges.put(vd, e2);
      frontList.add(e2);
    }
    
    // vd may become border or inner
    if(vd.outEdges.isEmpty() && vd.inEdges.isEmpty()) vd.type = 2;
    else vd.type = 1;
    
    steps++;
    if( steps % 100 == 0 ) println("pivot " + steps + " times");  
  }
  maxStep = (int)steps;
  
  // show inner edges as arrows
  if(showInnerEdges){
    fill(grey);
    for(Edge e: innerEdges){
      pt A = P.get(e.a).position;
      pt B = P.get(e.b).position;
      arrow(A, B, rt/2);
    }
  }
  
  // show boundary edges as arrows
  if(showBorderEdges){
    fill(green);
    for(Edge e: borderEdges){
      pt A = P.get(e.a).position;
      pt B = P.get(e.b).position;
      arrow(A, B, rt/2);
    }
  }
  
  // show triangle mesh
  if(showMesh){
    fill(red);
    for(Triangle t: mesh){
      pt A = P.get(t.a).position;
      pt B = P.get(t.b).position;
      pt C = P.get(t.c).position;
      beginShape(TRIANGLES);
        vertex(A); vertex(B); vertex(C);
      endShape();
    }
  }
  
  if(showVertexTypes){
    for(Vertex v : P){
      if(v.type == 1){ // border
        fill(dgreen); show(v.position, 10);
      }
      else if(v.type == 2){  // inner
        fill(blue); show(v.position, 10);
      }
      else{  // orphan
        fill(black); show(v.position, 10);
      }
    }
  }
  
  if(showPivotBall){
    fill(blue, 100); show(ballCenterToShow, r);
  }
}


void removeTriangleHole(ArrayList<Edge> borderEdges, ArrayList<Edge> innerEdges, ArrayList<Triangle> mesh, ArrayList<Vertex> P){
  // store true border edges, use hashing to fast access or remove a certain edge
  HashMap<Integer, Integer> bEdges = new HashMap<Integer, Integer>();
  for(int i = 0; i < borderEdges.size(); ++i){
    Edge e = borderEdges.get(i);
    bEdges.put(e.a, e.b);
  }
  
  for(int i = 0; i < borderEdges.size(); ++i){
    Edge e = borderEdges.get(i);
    if(!bEdges.containsKey(e.a) || bEdges.get(e.a) != e.b) continue; // edge e has been removed
    Vertex va = P.get(e.a), vb = P.get(e.b);
    for(Vertex vc : vb.outEdges.keySet()){  // bc is an edge, check whether ca is an edge
      if(va.inEdges.containsKey(vc)){  // abc is a triangle
        int a = va.id, b = vb.id, c = vc.id;
        mesh.add(new Triangle(a, b, c));
        
        // show the triangle in the triangular hole we find
        fill(blue);
        beginShape();
          vertex(va.position); vertex(vb.position); vertex(vc.position);
        endShape();
        
        bEdges.remove(a);
        bEdges.remove(b);
        bEdges.remove(c);
        innerEdges.add(new Edge(a, b));
        innerEdges.add(new Edge(b, c));
        innerEdges.add(new Edge(c, a));
        va.outEdges.remove(vb); va.inEdges.remove(vc);
        vb.outEdges.remove(vc); vb.inEdges.remove(va);
        vc.outEdges.remove(va); vc.inEdges.remove(vb);
        break;
      }
    }
  }
  
  ArrayList<Edge> newBorderEdges = new ArrayList<Edge>();
  for(Integer a : bEdges.keySet()){
    newBorderEdges.add(new Edge(a, bEdges.get(a)));
  }
  
  //println("number of border edges before post-processing: "+borderEdges.size());
  borderEdges.clear();
  borderEdges.addAll(newBorderEdges);
  //println("number of border edges after post-processing: "+borderEdges.size());
  
  return;
}