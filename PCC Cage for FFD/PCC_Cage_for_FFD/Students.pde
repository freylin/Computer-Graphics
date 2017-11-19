// Place student's code here
// Student's names:  Yaohong Wu, Jin Lin
// Date last modified:  09/19/2017
/* Functionality provided (say what works):
High-Level:
  Compute bi-arc
  Determine convexity of each arc using shooting ray method
  Find arc with supporting circle inside the stroke
  Walk along the boundary (between arcs whose supporting circles are inside the stroke) and find points on the medial axis
  Extend medial axis by extend some arcs
  Compute transversal and extended transversal 
  Discretize each transversal and get a quad matrix for a stroke
  Map between two strokes
 
Low-Level:
  intersection of line and circle (also line and arc)
  intersection of circle and circle (also arc and arc)
  determine if a point is on an arc
  solve quadratic equations
  sample a sequence of points with constant spacing on arcs
  etc.
  
For more details, see below.
*/



//************************************************************************
//**** ARC CLASS
//************************************************************************
class Arc{
  pt C = P(0, 0);  //center of the supporting circle
  pt A = P(1, 0);  //starting point
  pt B = P(0, 1);  //ending point
  boolean cw = true;  //direction
  boolean convex = true;  //convexity w.r.t. the region
  float r = 1;  //radius
  float w = PI/2;  //angle, [-2*PI,2*PI]
  int n = 30;  //number of points sampled from this arc, only used in show()
  
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
  
  //Fine middle point on the arc
  pt findMidPoint(){
    float a = w/2;
    pt mid = P(C, R(V(C,A), a));
    return mid;
  }
  
  void show(){
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


//Given point P1 on arc1, find the corresponding max disk w.r.t. arc1 and arc2
//Return the center of the max disk and a tangency point on arc2
//Note that empty list will be returned if we cannot find such max disk
ArrayList<pt> findPointOnMA(Arc arc1, Arc arc2, pt P1){
  pt M1, P2;
  pt C1 = arc1.C;
  pt C2 = arc2.C;
  float r2= arc2.r;
  vec N = U(C1, P1);
  vec C2P1 = V(C2, P1);
  float d;  // d can be positive or negative
  if(arc1.convex == arc2.convex){  // |(P1+dN)-C2| = r2 + d
    d = (dot(C2P1, C2P1) - r2*r2) / (2 * (-dot(C2P1, N) + r2));  //what if denominator == 0
  }
  else{  // |(P1+dN)-C2| = r2 - d
    d = (dot(C2P1, C2P1) - r2*r2) / (2 * (-dot(C2P1, N) - r2));  //what if denominator == 0
  }
  M1 = P(P1, d, N);
  
  if(arc2.convex){
    P2 = P(M1, abs(d), U(C2, M1));
  }
  else{
    P2 = P(M1, abs(d), U(M1, C2));
  }
  
  if(abs(d(P2, C2)-r2) > 0.01){  //if we choose P2 wrongly
    P2 = P(P2, 2, V(P2, M1));  //change our choice
  }
  
  //check if P2 is on arc2
  boolean onArc2 = pointOnArc(P2, arc2, true);
  ArrayList<pt> res = new ArrayList<pt>();
  if(onArc2){
    res.add(M1);
    res.add(P2);
  }
  return res;
}


//Given two arcs, arc1 and arc2, find n points on their medial axis
//Note that n "hat"s will be returned for finding transversals
ArrayList<pt[]> findMATwoArcs(Arc arc1, Arc arc2, int n){
  pt C1 = arc1.C;
  pt A1 = arc1.A;
  vec V = V(C1, A1);
  float dw = arc1.w / n;
  float a = 0;
  ArrayList<pt[]> res = new ArrayList<pt[]>();
  for(int i = 0; i <= n; ++i){
    pt p = P(C1, R(V, a));
    a += dw;
    ArrayList<pt> tmp = findPointOnMA(arc1, arc2, p);
    if(tmp.size() > 0){
      pt[] hat = new pt[3];
      hat[0] = p;
      hat[1] = tmp.get(0);
      hat[2] = tmp.get(1);
      res.add(hat);
    }
  }
  println("number of points on MA = ", res.size());
  return res;
}


//Given arc1, arc2 and arc3, find the medial axis between arc1 and arc2
//The medial axis will be trimmed such that it doesn't go beyond arc3
//startLen denotes the starting position in arc1
//unitLen denotes the arc-length step we use to sample on arc1
//Return a list of hats
//The tip of each hat is on MA of arc1 and arc2
//The other two vertices of each hat are on arc3
ArrayList<pt[]> findMATwoExtArcs(Arc arc1, Arc arc2, Arc arc3, float startLen, float unitLen){
  pt C1 = arc1.C;
  pt A1 = arc1.A;
  vec V = V(C1, A1);
  
  float r1 = arc1.r;
  float ang = startLen / r1;
  float dAng = unitLen / r1;
  if(arc1.cw == false){
    ang = -ang;
    dAng = -dAng;
  }
  
  pt C3 = arc3.C;
  float r3 = arc3.r;
  ArrayList<pt[]> res = new ArrayList<pt[]>();
  while(true){
    pt p = P(C1, R(V, ang));
    ang += dAng;
    ArrayList<pt> tmp = findPointOnMA(arc1, arc2, p);
    if(tmp.size() > 0){
      pt[] hat = new pt[3];
      pt P1 = p;
      pt M = tmp.get(0);
      pt P2 = tmp.get(1);
      
      //check whether hat[1] is inside supporting circle of acr3
      //if not, we stop extending 
      if(d(M, C3) > r3 + 8*EPSILON){
        break;
      }
      
      vec MP1 = V(M, P1);
      vec C3M = V(C3, M);
      float a1 = dot(MP1,MP1);
      float b1 = 2*dot(C3M, MP1);
      float c1 = dot(C3M, C3M) - r3*r3;
      float[] sol1 = solveQuadratic(a1,b1,c1);
      float t1 = sol1[0];  //we prefer the bigger solution (positive in this case)
      pt newP1 = P(M, t1, MP1);
      
      vec MP2 = V(M, P2);
      float a2 = dot(MP2, MP2);
      float b2 = 2*dot(C3M, MP2);
      float c2 = c1;
      float[] sol2 = solveQuadratic(a2,b2,c2);
      float t2 = sol2[0];  //we prefer the bigger solution (positive in this case)
      pt newP2 = P(M, t2, MP2);
      
      hat[0] = newP1;
      hat[1] = M;
      hat[2] = newP2;
      
      res.add(hat);
    }
  }
  
  return res;
}


//Given an array of arcs, a point P and its index of arc,
//find the max disk w.r.t. all the arcs
pt[] findHatOfMaxDiskAtP(pt P, int idx, Arc[] arcs){
  Arc arc = arcs[idx];
  float d = MAX_FLOAT;
  pt M = new pt();
  pt otherP = new pt();
  int n = arcs.length;
  for(int i = 0; i < n; ++i){
    if(i == idx) continue;  //ignore the same arc
    if(i == (idx-1+n)%n || i == (idx+1)%n) continue;  //ignore nearby arcs
    ArrayList<pt> tmp = findPointOnMA(arc, arcs[i], P);
    if(tmp.size() > 0){
      pt C = tmp.get(0);
      pt tmpP = tmp.get(1);
      float r = d(C, P);
      if(r < d) {
        d = r;
        M = C;
        otherP = tmpP;
      }
    }
  }
  pt[] hat = new pt[3];
  hat[0] = P;
  hat[1] = M;
  hat[2] = otherP;
  return hat;
}


//Extend an arc by angle w
//w is always positive
//If atA == true, then extend at end point A, else at end point B
Arc extendArc(Arc arc, boolean atA, float w){
  assert w >= 0;
  boolean cw = atA ? !arc.cw : arc.cw;
  if(cw==false) w = -w;
  pt P1 = atA ? arc.A : arc.B;
  pt C = arc.C;
  pt P2 = P(C, R(V(C, P1), w));
  boolean convex = arc.convex;
  Arc res = new Arc(C, P1, P2, cw, convex, 30);
  return res;  
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


//Determine the direction of an arc
//Note that this function is only called in findBiArc()
boolean clockwise(pt C, pt A, pt B, boolean bo){
  if (bo){
    return angle(A, C, B) > 0 ?  true : false;
  }
  return angle(A, C, B) > 0 ?  false : true;
}


//Find bi-arc given point P1, point P2, direction T1, direction T2
Arc[] findBiArc(pt P1, vec T1, pt P3, vec T2){
    vec S = V(P1, P3);
    vec T = W(T1, T2);
    float a;
    boolean bo = true;
    float Denom = 4 - n2(T);
    float DotST = dot(S, T);
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


//Determin if a point P is on a given arc
//If onCircle is true, we can omit checking whether P 
//is on the supporting circle of the arc
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


//Compute the intersection points of two circles (if any)
//Reference: https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
ArrayList<pt> intersectionOfTwoCircles(pt C1, float r1, pt C2, float r2){
  ArrayList<pt> res = new ArrayList<pt>();
  float d = d(C1,C2);
  if(d > r1+r2 || d < abs(r1-r2)) return res;
  assert (d==0 && r1==r2) == false;
  
  float a = (r1*r1 - r2*r2 + d*d) / (2*d);
  float h = sqrt(r1*r1 - a*a);
  
  vec D = U(V(C1,C2));
  vec N = R(D);
  
  pt M = P(C1, a, D);
  pt P1 = P(M, h, N);
  pt P2 = P(M, -h, N);
  res.add(P1);
  res.add(P2);
  return res;
}


//Compute the intersection points of two arcs
ArrayList<pt> intersectionOfTwoArcs(Arc arc1, Arc arc2){
  pt C1 = arc1.C;
  float r1 = arc1.r;
  pt C2 = arc2.C;
  float r2 = arc2.r;
  
  ArrayList<pt> intersections = intersectionOfTwoCircles(C1,r1,C2,r2);
  if(intersections.size() == 0) return intersections;
  
  ArrayList<pt> res = new ArrayList<pt>();
  for(int i = 0; i < intersections.size(); ++i){
    pt P = intersections.get(i);
    if(pointOnArc(P,arc1,true) && pointOnArc(P,arc2,true)) res.add(P);
  }
  return res;
}


//Find which arc has a supporting circles inside 
//the region bounded by the given arcs
boolean[] findInsideCircles(Arc[] arcs){
  int n = arcs.length;
  boolean[] mask = new boolean[n];
  for(int i = 0; i < n; ++i) mask[i] = false;
  for(int i = 0; i < n; ++i){
    if(arcs[i].convex == false) continue;
    Arc arc = arcs[i];
    Arc anotherArc = new Arc(arc.C, arc.A, arc.B, !arc.cw, true, 30);
    
    int j = 0;
    for(; j < n; ++j){
      if(j == i) continue;
      if(j == (i-1+n)%n || j == (i+1)%n) continue;  //ignore nearby arcs
      ArrayList<pt> tmp = intersectionOfTwoArcs(anotherArc, arcs[j]);
      if(tmp.size() > 0) {
        break;
      }
    }
    if(j == n){
      //check if nearby arc in this circle
      int k = (i+1)%n;
      Arc rightArc = arcs[k];
      pt P = rightArc.findMidPoint();
      pt C = arcs[i].C;
      float r = arcs[i].r;
      if(d(C,P) < r) mask[i] = false;  
      else mask[i] = true;
    }
  }
  return mask;
}


//Find the length of the boundary made of the given arcs
float findLengthOfBoundary(Arc[] arcs){
  float res = 0;
  for(int i = 0; i < arcs.length; ++i){
      res += findLengthOfArc(arcs[i]);
  }
  return res;
}


//Find the length of a given arc
float findLengthOfArc(Arc arc){
  return arc.r * abs(arc.w);
}


//Find the angle between CA and CB
//If cw is true, the angle should be non-negative, else non-postive
float findRotAngle(pt C, pt A, pt B, boolean cw){
  vec CA = V(C,A);
  vec CB = V(C,B);
  float w = angle(CA,CB);
  if(cw){
    if(w<0) w += 2*PI;
  }
  else{
    if(w>0) w -= 2*PI;
  }
  return w;
}


//Sample a sequence points on boundary made of arcs
//stIdxArc specifies the index of the starting arc
//insideMask tells us which arc has a supporting circle contained in the boundary such that we know when to stop
//If forMA is true, we only sample points on certain arcs, else on all the arcs
//points stores the points we sample
//idxsArc stores the corresponding index of arc for each sampled point
//unitLen specifies the arc-length step when we sample
//Return (the index of the last point used to construct medial axis) + 1
//If forMA is true, the returned index is equal to the size of points,
//else the returned index is less than the size of points (because more points will be sampled from more arcs)
int findPointsOnBoundary(Arc arcs[], int stIdxArc, boolean[] insideMask, boolean forMA, ArrayList<pt> points, ArrayList<Integer> idxsArc, float unitLen){
  float startLen = 0;
  int n = arcs.length;
  int edIdx = 0;
  boolean stop = false;
  for(int i = 0; i < n; ++i){
    int j = (stIdxArc + i) % n;
    
    if(stop == false && insideMask[j]){
      edIdx = points.size();
      stop = true;
      if(forMA) return edIdx;
    }
   
    Arc curArc = arcs[j];
    float r = curArc.r;
    float absw = abs(curArc.w);
    float ang = startLen / r;
    float dAng = unitLen / r;
    
    if(ang > absw){
      startLen -= r*absw;
      continue;
    } 
    int maxk = floor((absw - ang) / dAng);
    if(curArc.cw == false) {
      ang = -ang;
      dAng = -dAng;
    }
    pt C = curArc.C;
    pt A = curArc.A;
    vec CA = V(C,A);
    for(int k = 0;k <= maxk; ++k){
      pt P = P(C, R(CA, ang));
      points.add(P);
      idxsArc.add(j);
      ang += dAng;
    }
    
    float remainAng = absw - abs(ang-dAng);
    startLen = unitLen - remainAng * r;
  }
  
  return edIdx;
}


//Find medial axis for the region bounded by arcs
//points are the points we sampled from some arcs
//idxsArc are the corresponding indices of arc for points
//arcs specifies the region
//startIdx specifies the index of the starting point
//endIdx specifies the index of the ending point
ArrayList<pt[]> findHatsMedialAxis(ArrayList<pt> points, ArrayList<Integer> idxsArc, Arc[] arcs, int startIdx, int endIdx){
  ArrayList<pt[]> hatsMA = new ArrayList<pt[]>();
  
  for(int i = startIdx; i < endIdx; ++i){
    pt P = points.get(i);
    int idx = idxsArc.get(i);
    pt[] hat = findHatOfMaxDiskAtP(P, idx, arcs);
    pt M = hat[1];
    if(M.x == 0 && M.y == 0) continue;
    hatsMA.add(hat);
  }

  return hatsMA;
}


//Find the true positions for a boolean mask
ArrayList<Integer> findTruePos(boolean[] mask){
  ArrayList<Integer> res = new ArrayList<Integer>();
  for(int i = 0; i < mask.length; ++i){
    if(mask[i]) res.add(i);
  }
  return res;
}


//Discretize an arc defined by "hat" ABC
//n specifies the number of points we sample from this arc
pt[] discretizeArcDefinedByHat(pt A, pt B, pt C, int n){
  vec BA = V(B,A), BC = V(B,C);
  float d = dot(BC,BC) / dot(BC,W(BA,BC));
  pt X = P(B,d,W(BA,BC));  // X is the center of the circle decided by A, B, C
  
  vec XA = V(X,A), XC = V(X,C);
  float a = angle(XA,XC), da=a/(n-1);
  pt[] res = new pt[n];
  float w = 0;
  for (int i = 0; i < n; ++i) {
    res[i] = P(X,R(XA,w));
    w += da;
  }
  return res;
}


//Compute exact boundary, i.e. six arcs, from six control points
Arc[] computeExactBoundary(pt ctrlPts[]){
    //println("ctrlPts.length = ", ctrlPts.length);
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
      
    boolean[] isDiff = new boolean[3];
    isDiff[0] = arcs[0].convex != arcs[1].convex;
    isDiff[1] = arcs[2].convex != arcs[3].convex;
    isDiff[2] = arcs[4].convex != arcs[5].convex;
      
    for(int i = 0; i < 3; ++i){
      int j = 2*i;
      if(isConvex(j, arcs)) arcs[j].convex = true;
      else arcs[j].convex = false;
        
      if(isDiff[i]) arcs[j+1].convex = !arcs[j].convex;
      else arcs[j+1].convex = arcs[j].convex;
    }
    return arcs;
}


//Find extended medial axis for the two arcs neighboring arcs[idx]
//walkAtNextArc specifies the arc we sample
//startLen specifies the starting position on the arc from which we sample
//unitLen specifies the arc-length step when we sample
ArrayList<pt[]> extendHatsMA(Arc[] arcs, int idx, boolean walkAtNextArc, float startLen,  float unitLen){
  int numArcs = arcs.length;
  int idxN = (idx + 1) % numArcs;
  int idxP = (idx - 1 + numArcs) % numArcs;
  Arc extNext = extendArc(arcs[idxN], true, PI-(abs(arcs[idxN].w)));
  Arc extPrev = extendArc(arcs[idxP], false, PI-(abs(arcs[idxP].w)));
  
  ArrayList<pt[]> extHats;
  if(walkAtNextArc){
    extHats = findMATwoExtArcs(extNext, extPrev, arcs[idx], startLen, unitLen);
  }
  else{
    extHats = findMATwoExtArcs(extPrev, extNext, arcs[idx], startLen, unitLen);
  }
  return extHats;
}


//Find remaining arc length for a point P and an arc
//i.e. the arc-length of PB if B is the ending point of the arc
float findRemainArcLen(pt P, Arc arc){
  float w = angle(V(arc.C, P), V(arc.C, arc.B));
  if(arc.cw){
    if(w<0) w += 2*PI;
  }
  else{
    if(w>0) w -= 2*PI;
  }
  return abs(w)*arc.r;
}


//Find two extended medial axes (we can get two because we can extend the natural medial axis from two end points)
//idxs specifies the index of starting arc and the index of ending arc
//endPoints specifies the starting point and ending point we use to construct the natural medial axis
//unitLen specifies the arc-length step when we sample
ArrayList<ArrayList<pt[]>> extendHatsMAFull(Arc[] arcs, int[] idxs, pt[] endPoints, float unitLen){
  pt firstPt = endPoints[0];
  pt lastPt = endPoints[1];
  
  int idxArcFirstPt = (idxs[0]+1) % (arcs.length);
  int idxArcLastPt = (idxs[1]-1) % (arcs.length);
  
  float startLenL = unitLen - findRemainArcLen(firstPt, arcs[idxArcFirstPt]);
  float startLenR = unitLen - findRemainArcLen(lastPt, arcs[idxArcLastPt]);
  startLenL = max(startLenL, 0);
  startLenR = max(startLenR, 0);
 
  ArrayList<ArrayList<pt[]>> res = new ArrayList<ArrayList<pt[]>>();
  ArrayList<pt[]> extL = extendHatsMA(arcs, idxs[0], true, startLenL, unitLen);
  ArrayList<pt[]> extR = extendHatsMA(arcs, idxs[1], false, startLenR, unitLen);
  res.add(extL);
  res.add(extR);
  return res;
}


//Find quads given natural medial axis and two extended medial axes
//nvTrans specifies the number of points we sample for each transversal
ArrayList<pt[]> findQuads(ArrayList<pt[]> hatsMA, ArrayList<ArrayList<pt[]>> hatsExt, int nvTrans){
  ArrayList<pt[]> extHatsL = hatsExt.get(0);
  ArrayList<pt[]> extHatsR = hatsExt.get(1);
  ArrayList<pt[]> verticesQuad = new ArrayList<pt[]>();
  for(int i = extHatsL.size()-1; i >= 0; --i){
    pt[] vertices = discretizeArcDefinedByHat(extHatsL.get(i)[0], extHatsL.get(i)[1], extHatsL.get(i)[2], nvTrans);
    verticesQuad.add(vertices);
  }
  for(int i = 0; i < hatsMA.size(); ++i){
    pt[] vertices = discretizeArcDefinedByHat(hatsMA.get(i)[0], hatsMA.get(i)[1], hatsMA.get(i)[2], nvTrans);
    verticesQuad.add(vertices);
  }
  for(int i = 0; i < extHatsR.size(); ++i){
    pt[] vertices = discretizeArcDefinedByHat(extHatsR.get(i)[0], extHatsR.get(i)[1], extHatsR.get(i)[2], nvTrans);
    verticesQuad.add(vertices);
  }
  return verticesQuad;
}


//************************************************************************
//**** STROKE CLASS
//************************************************************************
class Stroke{
    pt[] ctrlPts = new pt[6];  //six control points
    Arc[] arcs = new Arc[6];  //six arcs
    boolean[] insideMask = new boolean[6];  //a boolean mask specifying which arc has an supporting circle inside the stroke
    ArrayList<Integer> idxsArcLR = new ArrayList<Integer>();  //two indices, one for starting arc and one for ending arc
    ArrayList<pt> ptsBoundary = new ArrayList<pt>();  //points on boundary
    ArrayList<Integer> idxsArcBoundary = new ArrayList<Integer>(); //indices of arcs for boundary points
    boolean forMA; // if it is true, then ptsBoundary are only from some part of the boundary, not the entire boundary
    int stIdx = 0, edIdx;  //index of starting point and index of ending point (we only use the points in the range [start, end-1] to construct the medial axis)
    float unitLen;  //the arc-length step we use in sampling
    int nvTrans;  //the number of samples from a transversal
    ArrayList<pt[]> hatsMA = new ArrayList<pt[]>();  //"hat"s on medial axis
    ArrayList<ArrayList<pt[]>> hatsExt = new ArrayList<ArrayList<pt[]>>();  //"hat"s on extended medial axis
    ArrayList<pt[]> quads = new ArrayList<pt[]>();  //quads we get from discretization of the stroke
    
    Stroke(pt[] _ctrlPts, boolean _forMA, float _unitLen, int _nvTrans){
      ctrlPts = _ctrlPts;
      forMA = _forMA;
      unitLen = _unitLen;
      nvTrans = _nvTrans;
    }
    
    void computeBoundary(){
      arcs = computeExactBoundary(ctrlPts);
    }
    
    void computeInsideMask(){
      insideMask = findInsideCircles(arcs);
    }
    
    void computeIdxsArcLR(){
      idxsArcLR = findTruePos(insideMask);
    }
    
    void computeHatsMA(){
      int idxStArc = (idxsArcLR.get(0)+1)%6;
      edIdx = findPointsOnBoundary(arcs, idxStArc, insideMask, forMA, ptsBoundary, idxsArcBoundary, unitLen);
      hatsMA = findHatsMedialAxis(ptsBoundary, idxsArcBoundary, arcs, stIdx, edIdx);
    }
    
    void extendHatsMA(){
      pt[] twoEndPts = {ptsBoundary.get(stIdx), ptsBoundary.get(edIdx-1)};
      int[] idxsLR = {idxsArcLR.get(0), idxsArcLR.get(1)};
      hatsExt = extendHatsMAFull(arcs, idxsLR, twoEndPts, unitLen);
    }
    
    void computeQuads(){
      quads = findQuads(hatsMA, hatsExt, nvTrans);
    }
    
    boolean isValid(){
      int cnt = 0;
      for(int i = 0; i < 6; ++i){
        if(insideMask[i]) cnt++;
      }
      return cnt == 2;
    }
    
    void check(){
      if(isValid()) return;
      else{
        showMessageDialog(null, "Medial Axis is not branch-free!", "Alert", ERROR_MESSAGE);
        exit();
      }
    }
    
    void visulize(){
      if(showBoundary){
        for (int i = 0; i < 6; ++i){
          if (i % 2 == 0) stroke(red);
          else stroke(blue);
          arcs[i].show();
        }
      }
      
      if(showMA){
        stroke(green);
        beginShape();
        for(int i = 0; i < hatsMA.size(); ++i){
          v(hatsMA.get(i)[1]);
        }
        endShape();
      }
      
      if(showExtMA){
        stroke(grey);
        beginShape();
        for(int i = 0; i < hatsExt.get(0).size(); ++i) v(hatsExt.get(0).get(i)[1]);
        endShape();
        beginShape();
        for(int i = 0; i < hatsExt.get(1).size(); ++i) v(hatsExt.get(1).get(i)[1]);
        endShape();
        if(hatsExt.get(0).size() > 0) edge(hatsMA.get(0)[1], hatsExt.get(0).get(0)[1]);
        if(hatsExt.get(1).size() > 0) edge(hatsMA.get(hatsMA.size()-1)[1], hatsExt.get(1).get(0)[1]);
      }
      
      if(showTrans){
        stroke(sand);
        for(int i = 0; i < hatsMA.size(); ++i) 
          drawCircleArcInHat(hatsMA.get(i)[0], hatsMA.get(i)[1], hatsMA.get(i)[2]);
      }
      
      if(showExtTrans){
        stroke(pink);
        for(int i = 0; i < hatsExt.get(0).size(); ++i) 
          drawCircleArcInHat(hatsExt.get(0).get(i)[0], hatsExt.get(0).get(i)[1], hatsExt.get(0).get(i)[2]);
        for(int i = 0; i < hatsExt.get(1).size(); ++i) 
          drawCircleArcInHat(hatsExt.get(1).get(i)[0], hatsExt.get(1).get(i)[1], hatsExt.get(1).get(i)[2]);
      }
      
      if(showQuadVerts){
        stroke(220);
        for(int i = 0; i < quads.size(); ++i){
          for(int j = 0; j < nvTrans; ++j){
             quads.get(i)[j].show(1);
          }
        }
      }
      
      if(showQuads){
        stroke(220);
        beginShape(QUADS);
        for(int i = 0; i < quads.size()-1; ++i){
          beginShape(QUAD_STRIP);
          for(int j = 0; j < nvTrans; ++j){
            v(quads.get(i)[j]); v(quads.get(i+1)[j]);
          }
          endShape(QUAD_STRIP);
        }
        endShape();
      }
      
      if(showTwoEndPoints){
        stroke(magenta); ptsBoundary.get(stIdx).show(8);
        stroke(cyan); ptsBoundary.get(edIdx-1).show(8);
      }
      
      if(showInsideCircles){
        stroke(brown);
        for(int i = 0; i < 6; ++i){
          if(insideMask[i]){
            pt C = arcs[i].C;
            float r = arcs[i].r;
            C.show(r);
          }
        }
      }
    }
}