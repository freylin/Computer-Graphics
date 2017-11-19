class pts // class for manipulaitng and displaying pointclouds or polyloops in 3D 
  { 
    int maxnv = 16000;                 //  max number of vertices
    pt[] G = new pt [maxnv];           // geometry table (vertices)
    char[] L = new char [maxnv];             // labels of points
    vec [] LL = new vec[ maxnv];  // displacement vectors
    Boolean loop=true;          // used to indicate closed loop 3D control polygons
    int pv =0,     // picked vertex index,
        iv=0,      //  insertion vertex index
        dv = 0,   // dancer support foot index
        nv = 0,    // number of vertices currently used in P
        pp=1; // index of picked vertex

  pts() {}
  pts declare() 
    {
    for (int i=0; i<maxnv; i++) G[i]=P(); 
    for (int i=0; i<maxnv; i++) LL[i]=V(); 
    return this;
    }     // init all point objects
  pts empty() {nv=0; pv=0; return this;}                                 // resets P so that we can start adding points
  pts addPt(pt P, char c) { G[nv].setTo(P); pv=nv; L[nv]=c; nv++;  return this;}          // appends a new point at the end
  pts addPt(pt P) { G[nv].setTo(P); pv=nv; L[nv]='f'; nv++;  return this;}          // appends a new point at the end
  pts addPt(float x,float y) { G[nv].x=x; G[nv].y=y; pv=nv; nv++; return this;} // same byt from coordinates
  pts copyFrom(pts Q) {empty(); nv=Q.nv; for (int v=0; v<nv; v++) G[v]=P(Q.G[v]); return this;} // set THIS as a clone of Q

  pts resetOnCircle(int k, float r)  // sets THIS to a polyloop with k points on a circle of radius r around origin
    {
    empty(); // resert P
    pt C = P(); // center of circle
    for (int i=0; i<k; i++) addPt(R(P(C,V(0,-r,0)),2.*PI*i/k,C)); // points on z=0 plane
    pv=0; // picked vertex ID is set to 0
    return this;
    } 
  // ********* PICK AND PROJECTIONS *******  
  int SETppToIDofVertexWithClosestScreenProjectionTo(pt M)  // sets pp to the index of the vertex that projects closest to the mouse 
    {
    pp=0; 
    for (int i=1; i<nv; i++) if (d(M,ToScreen(G[i]))<=d(M,ToScreen(G[pp]))) pp=i; 
    return pp;
    }
  pts showPicked() {show(G[pv],23); return this;}
  pt closestProjectionOf(pt M)    // Returns 3D point that is the closest to the projection but also CHANGES iv !!!!
    {
    pt C = P(G[0]); float d=d(M,C);       
    for (int i=1; i<nv; i++) if (d(M,G[i])<=d) {iv=i; C=P(G[i]); d=d(M,C); }  
    for (int i=nv-1, j=0; j<nv; i=j++) { 
       pt A = G[i], B = G[j];
       if(projectsBetween(M,A,B) && disToLine(M,A,B)<d) {d=disToLine(M,A,B); iv=i; C=projectionOnLine(M,A,B);}
       } 
    return C;    
    }

  // ********* MOVE, INSERT, DELETE *******  
  pts insertPt(pt P) { // inserts new vertex after vertex with ID iv
    for(int v=nv-1; v>iv; v--) {G[v+1].setTo(G[v]);  L[v+1]=L[v];}
     iv++; 
     G[iv].setTo(P);
     L[iv]='f';
     nv++; // increments vertex count
     return this;
     }
  pts insertClosestProjection(pt M) {  
    pt P = closestProjectionOf(M); // also sets iv
    insertPt(P);
    return this;
    }
  pts deletePicked() 
    {
    for(int i=pv; i<nv; i++) 
      {
      G[i].setTo(G[i+1]); 
      L[i]=L[i+1]; 
      }
    pv=max(0,pv-1); 
    nv--;  
    return this;
    }
  pts setPt(pt P, int i) { G[i].setTo(P); return this;}
  
  pts drawBalls(float r) {for (int v=0; v<nv; v++) show(G[v],r); return this;}
  pts showPicked(float r) {show(G[pv],r); return this;}
  pts drawClosedCurve(float r, color stubColor) 
    {
    fill(dgreen);
    for (int v=0; v<nv; v++) show(G[v],r*3);    
    fill(stubColor);
    for (int v=0; v<nv-1; v++) stub(G[v],V(G[v],G[v+1]),r,r);  
    stub(G[nv-1],V(G[nv-1],G[0]),r,r);
    pushMatrix(); //translate(0,0,1); 
    scale(1,1,0.03);  
    fill(grey);
    for (int v=0; v<nv; v++) show(G[v],r*3);    
    for (int v=0; v<nv-1; v++) stub(G[v],V(G[v],G[v+1]),r,r);  
    stub(G[nv-1],V(G[nv-1],G[0]),r,r);
    popMatrix();
    return this;
    }
  pts set_pv_to_pp() {pv=pp; return this;}
  pts movePicked(vec V) { G[pv].add(V); return this;}      // moves selected point (index p) by amount mouse moved recently
  pts setPickedTo(pt Q) { G[pv].setTo(Q); return this;}      // moves selected point (index p) by amount mouse moved recently
  pts moveAll(vec V) {for (int i=0; i<nv; i++) G[i].add(V); return this;};   
  pt Picked() {return G[pv];} 
  pt Pt(int i) {if(0<=i && i<nv) return G[i]; else return G[0];} 



  // ********* I/O FILE *******  
 void savePts(String fn) 
    {
    String [] inppts = new String [nv+1];
    int s=0;
    inppts[s++]=str(nv);
    for (int i=0; i<nv; i++) {inppts[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z)+","+L[i];}
    saveStrings(fn,inppts);
    };
  
  void loadPts(String fn) 
    {
    println("loading: "+fn); 
    String [] ss = loadStrings(fn);
    String subpts;
    int s=0;   int comma, comma1, comma2;   float x, y;   int a, b, c;
    nv = int(ss[s++]); print("nv="+nv);
    for(int k=0; k<nv; k++) 
      {
      int i=k+s; 
      //float [] xy = float(split(ss[i],",")); 
      String [] SS = split(ss[i],","); 
      G[k].setTo(float(SS[0]),float(SS[1]),float(SS[2]));
      L[k]=SS[3].charAt(0);
      }
    pv=0;
    };
 
  // Dancer
  void setPicekdLabel(char c) {L[pp]=c;}
  
  void setFifo() 
    {
    _LookAtPt.reset(G[dv],60);
    }              

  void next() {dv=n(dv);}
  int n(int v) {return (v+1)%nv;}
  int p(int v) {if(v==0) return nv-1; else return v-1;}
  
  // demo subdivision
  pts subdivideDemoInto(pts Q){
    Q.empty();
    for(int i=0; i<nv; i++){
      Q.addPt(P(G[i])); 
      Q.addPt(P(G[i],G[n(i)]));
    }
    return this;
  }
    
  // quintic B-Spline subdivision
  pts subdivideQuinticInto(pts Q){
    float t = -0.5;
    JSplineSubdivideFromPToQ(this, Q, t);
    return this;
  }
 
  // cubic B-Spline subdivision
  pts subdivideCubicInto(pts Q){
    float t = 0;
    JSplineSubdivideFromPToQ(this, Q, t);
    return this;
  }
 
  // Jarek subdivision
  pts subdivideJarekInto(pts Q){
    float t = 0.5;
    JSplineSubdivideFromPToQ(this, Q, t);
    return this;
  }
  
  // four-point subdivision
  pts subdivideFourPointInto(pts Q){
    float t = 1;
    JSplineSubdivideFromPToQ(this, Q, t);
    return this;
  }
  
  // quadratic B-Spline subdivision
  pts subdivideQuadraticInto(pts Q){
    //not J-Spline subdivision
    //repeat (refine, dual) (or (refine, refine, kill)) for 1 time
    Q.empty();
    for(int i = 0; i < nv; ++i){
      pt A = P(2/3., G[i], 1/3., G[n(i)]);  // 2/3. is a floating number (=0.666...), 2/3 is an integer (=0)
      pt B = P(1/3., G[i], 2/3., G[n(i)]);
      Q.addPt(A);
      Q.addPt(B);
    }
    return this;
  }
  
  
  // display Skater or arrow running on a trajectory
  void displaySkater() {
    if(showCurve) {fill(yellow); for (int j=0; j<nv; j++) caplet(G[j],6,G[n(j)],6); }
    pt[] B = new pt [nv];           // geometry table (vertices)
    for (int j=0; j<nv; j++) B[j]=P(G[j],V(0,0,h));  // points of the upper curve 
    if(showPath) {fill(lime); for (int j=0; j<nv; j++) caplet(B[j],6,B[n(j)],6);} 
    if(showKeys) {fill(cyan); for (int j=0; j<nv; j+=4) arrow(B[j],G[j],3);}
    
    if(animating) f++;
    if (f % df == 0) ff = (ff+1)%nv;
    if (f == 99999) f = 0;
    
    pts trajectory;
    if(useOffsetCurve) trajectory = T;
    else trajectory = this;
    
    assert trajectory.nv == this.nv;
    assert B.length == this.nv;
    
    if (showSkater){
      if(!constStepLength){
        int l = ff, r = n(ff);
        vec forwardDirection = V(trajectory.G[l], trajectory.G[r]);
        forwardDirection.normalize();
        if (ff % 2 == 0){
          l = n(ff);
          r = ff;
        }
        showDancer(trajectory.G[l], s, trajectory.G[r], forwardDirection, B[ff]); 
      }
      else{
        Step nextStep = oneStep(curStep, trajectory, stepLength);
        pt leftFoot, rightFoot, head;
        leftFoot = positionOfStep(curStep, trajectory);
        rightFoot = positionOfStep(nextStep, trajectory);
        vec forwardDirection = V(leftFoot, rightFoot);
        forwardDirection.normalize();
        if (ff%2 == 0){
          swap(leftFoot, rightFoot);
        }
        head = P(1-curStep.t, B[curStep.i], curStep.t, B[(curStep.i+1)%nv]);
        showDancer(leftFoot, s, rightFoot, forwardDirection, head);
        curStep = nextStep;
      }
    }
    else{
      fill(red);
      if(!constStepLength) {
        arrow(B[f%nv], trajectory.G[f%nv], 20);
      }
      else{
        int i = curStep.i;
        float t = curStep.t;
        pt foot = P(1-t, trajectory.G[i], t, trajectory.G[(i+1)%nv]);
        pt head = P(1-t, B[i], t, B[(i+1)%nv]);
        arrow(head, foot, 20);
        curStep = oneStep(curStep, trajectory, stepLength);
      }
    }
  }
} // end of pts class


// J-Spline subdivision skeme
// t is parameter
// t = -0.5: quintic B-Spline
// t = 0: cubic B-Spline
// t = 0.5: jarak
// t = 1: four-point
void JSplineSubdivideFromPToQ(pts P, pts Q, float t){
  Q.empty();
  pt[] G = P.G;
  int nv = P.nv;
  for(int i = 0; i < nv; ++i){
    int a = i;
    int b = n(a,nv);
    int c = n(b,nv);
    int d = n(c,nv);
    int e = n(d,nv);
    int f = n(e,nv);
    pt A = G[a], B = G[b], C = G[c], D = G[d], E = G[e], F = G[f];
    rs(A,B,C,D,E,F, Q, t, 1);  //only subdivide by one level
  }
}


// recursively subdivide (into a span which is controled by six control points)
// r is level
// t is parameter for J-Spline
// Q stores points on subdivision curve
void rs(pt _A, pt _B, pt _C, pt _D, pt _E, pt _F, pts Q, float t, int r) {
  pt A = P(_A), B = P(_B), C = P(_C), D = P(_D), E = P(_E), F = P(_F);
  if (r==0) {
    Q.addPt(C);
  } 
  else {
    pt J = P(B,C), K = P(C,D), M=P(D,E);
    vec BL = V(1./4, V(B,P(A,C)));
    vec CL = V(1./4, V(C,P(B,D)));
    vec DL = V(1./4, V(D,P(C,E)));
    vec EL = V(1./4, V(E,P(D,F)));
    
    B.add(1-t, BL);
    C.add(1-t, CL);
    D.add(1-t, DL);
    E.add(1-t, EL);
    
    J.add(-t, V(BL,CL));
    K.add(-t, V(CL,DL));
    M.add(-t, V(DL,EL));
    
    rs(B,J,C,K,D,M, Q, t, r-1);
    rs(J,C,K,D,M,E, Q, t, r-1);
  };
};


// next index of i w.r.t. nv
int n(int i, int nv){
  if(i == nv-1) return 0;
  else return i+1;
}


// previous index of i w.r.t. nv
int p(int i, int nv){
  if(i == 0) return nv-1;
  else return i-1;
}


// compute centripetal force (acceleration)
// note that what we compute is an approximation
float centripetalForce(pt A, pt B, pt C, boolean useCurvature){
  if(useCurvature){  // (v^2)/r, where v is the magnitude of velocity and r is the radius of curvation
    vec AB = V(A,B);
    vec NB = U(R(V(A,C)));
    return 2*abs(dot(AB,NB));
  }
  else{  // a, where a is acceleration
    vec BA = V(B,A), BC = V(B,C);
    vec gB = A(BA,BC);
    return n(gB)/2;
  }
}


// displace curve according to acceleration
void displaceCurve(pts P){
  T.empty();
  int nv = P.nv;
  pt[] G = P.G;
  for(int i = 0; i < nv; ++i){
    pt P1 = G[p(i,nv)], P2 = G[i], P3 = G[n(i,nv)];
    float cpf = centripetalForce(P1,P2,P3,useCurvature);
    float offset = h*cpf / grav;  //h and grav are global constants
    vec N = U(R(V(P1,P3)));
    if(dot(V(P1,P2), N) > 0) N = M(N); // N = -N
    pt newPoint = P(P2, -offset, N);
    T.addPt(newPoint);
  }
}


// retrofit
// P is the original control polygon
// R is the iterated control polygon
// s is the parameter related to J-Spline skeme
void retrofit(pts P, pts R, float s){
  int nv = P.nv;
  pts nR = new pts();
  nR.declare();
  for(int i = 0; i < nv; ++i){
    pt A = R.G[(i-2+nv)%nv];
    pt B = R.G[(i-1+nv)%nv];
    pt C = R.G[i];
    pt D = R.G[(i+1)%nv];
    pt E = R.G[(i+2)%nv];
    float a = s*(s-1), b = 2*s*(8-s), c = 72+2*(s-9)*s;
    pt AB = P(a, A, b, B);
    pt CDE = P(c, C, b, D, a, E);
    pt L = A(AB,CDE);
    L.div(12*(6+s));
    vec discrepancy = V(L, P.G[i]);
    nR.addPt(P(C, discrepancy));    
  }
  R.copyFrom(nR);
}


// Class Step
// a step is a point on a polyloop
class Step{
  int i = 0;  //index of vertex of the polyloop, such that the step is in the line segment between p[i] and p[i+1]
  float t = 0;  //parameter used to find the position of the step (position of step = (1-t)*p[i] + t*p[i+1])
  Step(){}
  Step(int _i, float _t){
    i = _i;
    t = _t;
  }
}


// compute next step given current step, trajectory and step length
Step oneStep(Step curStep, pts points, float distance){
  int nv = points.nv;
  pt[] G = points.G;
  int i = curStep.i;
  float t = curStep.t;
  int j = n(i,nv);  //next index of i
  pt P = P(1-t, G[i], t, G[j]);
  float remain = d(P, G[j]);
  Step nextStep = new Step();
  if(remain > distance){
    nextStep.i = i;
    nextStep.t = t + distance/d(G[i],G[j]);
  }
  else{
    do{
      distance -= remain;
      remain = d(G[j], G[n(j,nv)]);
      j = n(j,nv);
    }while(remain < distance);
    int k = p(j,nv);
    nextStep.i = k;
    nextStep.t = distance / remain;
  }
  return nextStep;
}


// compute the position of a step given the step and the trajectory
pt positionOfStep(Step step, pts trajectory){
  return P(1-step.t, trajectory.G[step.i], step.t, trajectory.G[n(step.i,trajectory.nv)]);
}


// swap point A and point B
void swap(pt A, pt B){
  pt C = P(B);  //copy B into C
  B.set(A);
  A.set(C);
}