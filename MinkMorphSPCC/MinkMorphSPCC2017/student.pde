//************************************************************************
// Student's names: Jin Lin
// Date last modified:  10/30/2017

//**** PCC CLASS
//************************************************************************
class PCC{
    pt[] ctrlPts = new pt[12];
    Arc[] arcs = new Arc[360];
    pt[] C = new pt[360];  //center of each arc
    pt[] M = new pt[360];  //morph vertices
    vec[] norm = new vec[360];  //normal of each vertex
    int len = 0;
    
    PCC(){}
    
    PCC(pt[] _ctrlPts){
      ctrlPts = _ctrlPts;
    }
    
    void computeBoundary(){
      arcs = computeExactBoundary(ctrlPts);
      for (int i = 0; i < 6; i++){
        C[i] = arcs[i].C;
        if (arcs[i].convex){
          if (arcs[i].cw) M[i] = arcs[i].A;
          else M[i] = arcs[i].B;
          norm[i] = U(C[i], M[i]);
        }
        else{
          if (arcs[i].cw) M[i] = arcs[i].B;
          else M[i] = arcs[i].A;
          norm[i] = U(M[i], C[i]);
        }
      }
    }
    
    PCC copy(){
      PCC NewP = new PCC();
      for (int i = 0; i < 6; i++){
        NewP.M[i] = P(M[i]);
        NewP.C[i] = P(C[i]);
        NewP.norm[i] = V(norm[i]);
      }
      return NewP;
    }
    
    void sortByNorm(int num){
      for (int i = 0; i < num; i++)
        for (int j = i+1; j < num; j++)
        if (norm[j].angle() < norm[i].angle())
          {
            vec tempv = norm[i]; norm[i] = norm[j]; norm[j] = tempv;
            pt tempp = M[i]; M[i] = M[j]; M[j] = tempp;
            tempp = C[i]; C[i] = C[j]; C[j] = tempp;
          }
    }
     
    void drawNorm(int starti, int endi){
      if (!showNorm) return;
      stroke(orange);
      for (int i = starti; i < endi; i++) edge(M[i], P(M[i], 10, norm[i]));
    } 
    
    // split this PCC by P
    void split(PCC P, int n1, int n2){
      for (int i = 0; i < n2; i++)
      {
        int k = n1 + i;
        for (int j = 0; j < n1; j++){
          if (between(norm[j], P.norm[i], norm[(j+1)%n1], true)) {
            M[k] = P(C[j], d(C[j], M[j]), P.norm[i]);
            C[k] = C[j];
            norm[k] = P.norm[i];
            break;
          }
        }
      }
    }
    
    // split this concave PCC by P
    void split(PCC P, int n1, int n2, int orient){
      for (int i = 0; i < n2; i++)
      {
        int k = n1 + i;
        for (int j = 0; j < n1; j++){
          if (between(norm[j], P.norm[i], norm[(j+1)%n1], true)) {
            M[k] = P(C[j], d(C[j], M[j]), W(orient, P.norm[i]));
            C[k] = C[j];
            norm[k] = P.norm[i];
            break;
          }
        }
        if (norm[k]==null){
            M[k] = P(M[0]);
            C[k] = P(C[0]);
            norm[k] = V(norm[0]);
        }
      }
    }
    
    // split PCC into 360 parts (iter==360)
    PCC split(int num, int iter){
      PCC P = new PCC();
      vec N = V(1, 0);
      float a = PI * 2 / iter;
      for (int i = 0; i < num; i++)
        if (between(norm[i], N, norm[(i+1)%num], true)) {
          int index = i;
          int count = 0;
          while (count < iter){
            P.M[count] = P(C[index], d(C[index], M[index]), N);
            P.norm[count] = V(N);
            N.rotateBy(a);
            if (!between(norm[index], N, norm[(index+1)%num], true))
              index = (index+1) % num;
            count++;
          }
        break;
        }
      return P;
    }
    
    // split this PCC by arc A
    void split(Arc A, int n1){
      vec[] N = new vec[2];
      if (A.cw) {N[0] = U(A.B, A.C); N[1] = U(A.A, A.C);}
      else {N[0] = U(A.A, A.C); N[1] = U(A.B, A.C);}
      
      for (int i = 0; i < 2; i++)
      {
        int k = n1 + i;
        for (int j = 0; j < n1; j++){
          if (between(norm[j], N[i], norm[(j+1)%n1], false)) {
            M[k] = P(C[j], d(C[j], M[j]), N[i]);
            C[k] = C[j];
            norm[k] = N[i];
            
            CUT[i] = j;
            break;
          }
        }
      }
    }
    
    // trim/cut out
    PCC cut(int n1, int starti, int endi, int step, boolean BETWEEN, boolean BREAK, int first, int last){
      PCC P2 = new PCC();
      P2.M[0] = P(M[first]);
      P2.C[0] = P(C[first]);
      P2.norm[0] = V(norm[first]);
      int n2 = 1;
      int i = starti;
      while (i!=endi)
      {
        int ii = (i+n1) % n1;
        boolean find = false;
        if (BETWEEN) find =  between(norm[n1+1], norm[ii], norm[n1], false);
        else find = !(dot(R(norm[n1+1]), norm[ii])>=0 && dot(R(norm[n1]), norm[ii])<=0);
        if (find){
          P2.M[n2] = P(M[ii]);
          P2.C[n2] = P(C[ii]);
          P2.norm[n2] = V(norm[ii]);
          n2++;
        }
        else if (BREAK) break;
        i = i + step;
      }
      P2.M[n2] = P(M[last]);
      P2.C[n2] = P(C[last]);
      P2.norm[n2] = V(norm[last]);
      n2++;
      P2.len = n2;
      return P2;
    }
    
    void visualize(color col1, color col2, int num){
      if (num < 360)
        for (int i = 0; i < num; i++){
          if (i % 2 == 1) stroke(col1);
          else stroke(col2);
          arcs[i].show();
        }
      else{
        beginShape();
        for(int i = 0; i < num; i++) {
          if (i % 2 == 1) stroke(col1);
          else stroke(col2);
          v(M[i]);
        }
        if (num % 2 == 1) stroke(col1);
          else stroke(col2);
        v(M[0]);
        endShape();
      }
    }
}

boolean between(vec a, vec c, vec b, boolean CLOSED){
  if (CLOSED) return (dot(R(a), c) >= -(1.0E-5) && dot(R(b), c) <= (1.0E-5));
  return (dot(R(a), c)>(1.0E-5) && dot(R(b), c)<(-1.0E-5));
}


Arc[] computeExactBoundary(pt ctrlPts[]){
    assert ctrlPts.length >= 6;
    pt P1 = ctrlPts[0], P2 = ctrlPts[1], P3 = ctrlPts[2];
    pt P4 = ctrlPts[3], P5 = ctrlPts[4], P6 = ctrlPts[5];
    vec T1 = U(P1, P2), T2 = U(P3, P4), T3 = U(P5, P6);
    Arc[] arcs = new Arc[6];
    Arc[] tmpBiArc;
    tmpBiArc = findBiArc(P1, T1, P3, T2);
    arcs[0] = tmpBiArc[0]; arcs[1] = tmpBiArc[1];
    tmpBiArc = findBiArc(P3, T2, P5, T3);
    arcs[2] = tmpBiArc[0]; arcs[3] = tmpBiArc[1];
    tmpBiArc = findBiArc(P5, T3, P1, T1);
    arcs[4] = tmpBiArc[0]; arcs[5] = tmpBiArc[1];
    
    for(int i = 0; i < 6; ++i)
      arcs[i].convex = isConvex(i, arcs);
    return arcs;
}

//************************************************************************
//**** ARC CLASS
//************************************************************************
class Arc{
  pt C = P(0, 0);
  pt A = P(1, 0);
  pt B = P(0, 1);
  boolean cw = true;
  int n = 30;
  boolean convex = true;
  float r = 1;
  float w = PI/2;
  
  Arc(){};
  Arc(pt _C, pt _A, pt _B, boolean _cw, boolean _convex, int _n) {
    C = _C; A = _A; B = _B; cw = _cw; convex = _convex; n = _n;
    r = d(C, A);
    float tmpW = angle(V(C,A), V(C,B));
    if(cw){
      if(tmpW < 0) tmpW += 2*PI;
    }
    else{
      if(tmpW > 0) tmpW -= 2*PI;
    }
    w = tmpW;
  };
  
  pt findMidPoint(){
    float a = w/2;
    pt mid = P(C, R(V(C,A), a));
    return mid;
  }
  
  void show(){
    if (w>PI*2/10*9) return;
    vec CA = V(C, A);
    float dw = w / n;
    float a = 0;
    beginShape();
    for(int i = 0; i <= n; ++i){
      v(P(C, R(CA, a)));
      a += dw;
    }
    endShape();
  }
}

//Find bi-arc given point P1, point P2, direction T1, direction T2
Arc[] findBiArc(pt P1, vec T1, pt P3, vec T2){
    vec S = V(P1, P3);
    vec T = W(T1, T2);
    float a;
    boolean bo = true;
    float Denom = 4 - n2(T);
    float DotST = dot(S, T);
    float th=0.02; // threshold
    if (abs(Denom) > th)
    {
      float d = sq(DotST) + n2(S) * Denom;
      a = (sqrt(d) - DotST) / Denom;
    } else
    {
      DotST = dot(S, T1);
      if (DotST > th)
        a = n2(S) / (4 * DotST);
      else if (DotST < -th)
      {
        bo = false;
        a = n2(S) / (4 * DotST);
      } else 
      {
        //rarely happend
        a = 100;
        // half circles
      }
    }
    pt B1 = P(P1, a, T1), B2 = P(P3, -a, T2);
    pt C = P(B1, 0.5, V(B1, B2));
    pt CC1 = findCircleCenter(P1, B1, C), CC2 = findCircleCenter(P3, B2, C);
    
    boolean cw1 = clockwise(CC1, P1, C, bo);
    boolean convex1 = true;
    Arc[] arcs = new Arc[2];
    arcs[0] = new Arc(CC1, P1, C, cw1, convex1, 30);
    
    boolean cw2 = clockwise(CC2, C, P3, bo);
    vec norm = (V(B1, B2).rotateBy(PI / 2)).normalize();
    boolean convex2 = dot(V(C, CC1), norm) * dot(V(C, CC2), norm) > 0 ? true : false;
    
    arcs[1] = new Arc(CC2, C, P3, cw2, convex2, 30);
    return arcs;
}

//Given a bunch of arcs and a index of some arc
//determine whether arcs[idx] is convex 
boolean isConvex(int idx, Arc[] arcs){
  int cnt = 0;
  pt P = arcs[idx].findMidPoint();
  vec N = U(V(arcs[idx].C, P));
  vec D = R(N);
  for(int i = 0; i < arcs.length; ++i){
    if(i == idx) continue;
    int tmp = numOfIntersectionRayArc(P, D, arcs[i]);
    cnt += tmp;
  }
  return cnt % 2 == 0;
}

//Find the number of intersection between a ray and an arc
int numOfIntersectionRayArc(pt P, vec D, Arc arc){
  // point of intersection = P+tD
  // t is unknown
  pt C = arc.C;
  float r = arc.r;
  vec CP = V(C,P);
  float a = n2(D);
  float b = 2*dot(CP, D);
  float c = n2(CP)-r*r;
  
  float[] sol = solveQuadratic(a,b,c);
  if(sol[0] == MAX_FLOAT) return 0;
  if(sol[0] <= 0) return 0;  //we want t > 0
  
  int res = 0;
  pt testP;
  
  testP = P(P, sol[0], D);
  if(pointOnArc(testP, arc, true)) res++;
  
  if(sol[1] > 0){
    testP = P(P, sol[1], D);
    if(pointOnArc(testP, arc, true)) res++;
  }
  return res;
}

boolean pointOnArc(pt P, Arc arc, boolean onCircle){
  pt A = arc.A;
  pt C = arc.C;
  
  if(!onCircle){
    if( abs(d(P,C) - arc.r) > 8*EPSILON) return false;  
  }

  boolean onArc = false;
  float a = angle(V(C, A), V(C, P));
  if(arc.cw){
    if(a < 0) a += 2*PI;
    if(a <= arc.w + 8*EPSILON) //why 8*EPSILON here? relax the restriction because of floating-point calculation
      onArc = true;
  }
  else{
    if(a > 0) a -= 2*PI;
    if(a >= arc.w - 8*EPSILON) //why 8*EPSILON here? relax the restriction because of floating-point calculation
      onArc = true;
  }
  return onArc;
}

//Solve quadratic equation in the form of ax^2+bx+c = 0
float[] solveQuadratic(float a, float b, float c){
  float[] res = new float[2];
  float delta = b*b - 4*a*c;
  if(delta < 0) {  //no solution
    res[0] = MAX_FLOAT;  //just a label
    res[1] = MIN_FLOAT;
    return res;
  }
  if(delta > EPSILON){  //two solutions
    res[0] = (-b + sqrt(delta)) / (2*a);  //the first one is bigger
    res[1] = (-b - sqrt(delta)) / (2*a);  //the second one is smaller
    return res;
  }
  else{  //one solution
    res[0] = res[1] = (-b)/(2*a);  
    return res;
  }
}

//Given a hat defined by A, H, B, find the center of 
//the circle defined by the hat
pt findCircleCenter(pt A, pt H, pt B){
  vec HA = V(H, A);
  vec HB = V(H, B);
  vec HD = W(HA, HB);
  float d = dot(HA, HA) / dot(HA, HD);
  pt C = P(H, d, HD);
  return C;
}

//Determine the direction of an arc
//Note that this function is only called in findBiArc()
boolean clockwise(pt C, pt A, pt B, boolean bo){
  if (bo){
    return angle(A, C, B) > 0 ?  true : false;
  }
  return angle(A, C, B) > 0 ?  false : true;
}