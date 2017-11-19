// Place student's code here
// Student's names:  Jin Lin
// Date last edited:  09/04/2017
/* Functionality provided (say what works):
Find the middle vertex of a potential "hat".
Find two tangent points given a point and a circle.
Choose the good one between the two tangent points.
Find circle center given a "hat".
Get a sequence of points of an arc.
Draw a circular arc given a center and two points on the arc. Note that this arc is a minor arc.
A PMC class for problem 2.
For more details, see comments before each function.
*/


/***
Given point L on a circle with center S, find a point H such that 
LH is tangent to circle S and the tangent line segment to circle E 
starting at H has the same length as |LH|.
***/
pt findHatMiddlePoint(pt S, pt E, float e, pt L){
  vec SL = V(S, L);
  vec T = U(R(SL));
  vec EL = V(E, L);
  
  float tmp = dot(EL, T);
  float numerator = e*e - dot(EL, EL);
  float x = 0;
  
  /*
  if(abs(numerator) >= EPSILON ){
    if(abs(tmp) < EPSILON) {  // special case : EL is orthogonal to T
      x = (numerator > 0 ? MAX_FLOAT : MIN_FLOAT);
    }
    else{
      x = numerator / (2*tmp) ;  
    }
  }
  */
  
  x = numerator / (2*tmp);
  pt H = P(L, x, T);
  return H;
}


/***
Given point A outside of a circle with center C and radius r,
find two points P1 and P2 such that AP1 and AP2 are both tangent to circle C.  
***/
pt[] findTangentPoints(pt A, pt C, float r){
  pt[] ps = new pt[2];
  float x = d(A, C);  // obviously, x != 0
  float w = acos(r/x);
  vec dir = U(V(C,A));
  ps[0] = P(C, r, R(dir, w));
  ps[1] = P(C, r, R(dir, -w));
  return ps;
}


/***
Choose one point from two candidate tangent points.
Here "good" means the point we find is the one on the 
differnet side from S w.r.t. line LE.
***/
pt findGoodOne(pt S, pt E, pt L, pt[] ps){
  vec LP1 = V(L, ps[0]);
  //vec LP2 = V(L, ps[1]);
  vec LE = V(L, E);
  vec LS = V(L, S);
  vec LK = R(LE);
  float dot0 = dot(LK, LS);  //  dot0 == 0 <==> S, E, L on the same line
  float dot1 = dot(LK, LP1);  //  dot1 == 0 <==> S, E, L on the same line
  //float dot2 = dot(LK, LP2);  //  dot2 == 0 <==> S, E, L on the same line
  if(dot1 * dot0 < 0) return ps[0];
  else return ps[1];
}


/***
Given a hat defined by A, H, B,
find the center of the circle defined by the hat.
***/
pt findCircleCenter(pt A, pt H, pt B){
  vec HA = V(H, A);
  vec HB = V(H, B);
  vec HD = W(HA, HB);
  float d = dot(HA, HA) / dot(HA, HD);
  pt C = P(H, d, HD);
  return C;
}


/***
Given a circle with center C and two points A, B on it,
find a sequence of points on arc AB.
Note that arc AB is a minor arc.
When A, C, B on the same line, there are two possible ways to draw arc AB.
In this case, if inverse == true, draw arc AB in cw direction, otherwise in ccw direction.
***/
void getPointsOnMinorArc(pts points, pt C, pt A, pt B, boolean inverse){
  vec CA = V(C, A);
  vec CB = V(C, B);
  float w = angle(CA, CB);
  if(inverse) w = -w;
  //print("angle is ", w, "\n");
  float dw = w / 30;
  float a = 0;
  for(int i = 0; i < 30; ++i){
    points.addPt(P(C, R(CA, a)));
    a += dw;
  }
}

/***
Given a circle with center S and two points A, B on it 
and another circle with center E, find a sequence of points
on arc AB such at it doesn't intersect with segment SE.
Note that arc AB could be minor or major.
***/
void getPointsOnSEArc(pts points, pt S, pt A, pt B, pt E){
  vec SE = V(S, E);
  vec SA = V(S, A);
  vec SB = V(S, B);
  float w = angle(SA, SB);
  //print("angle w is ", w, "\n");
  if(dot(SE, W(SA, SB)) > 0){
    w = w > 0 ? w - TWO_PI : w + TWO_PI;
  }
  float dw = w / 30;
  float a = 0;
  for(int i = 0; i < 30; ++i){
    points.addPt(P(S, R(SA, a)));
    a += dw;
  }
}


/***
Given a circle with center C and two points A, B on it,
draw arc AB.
Note that the arc we draw is a minor arc.
***/
void drawCircleMinorArc(pt C, pt A, pt B){
  vec CA = V(C, A);
  vec CB = V(C, B);
  float w = angle(CA, CB);
  //print("angle is ", w, "\n");
  float dw = w / 60;
  float a = 0;
  beginShape();
  for(int i = 0; i < 60; ++i){
    v(P(C, R(CA, a)));
    a += dw;
  }
  endShape();
}


//************************************************************************
//**** POINT MEMBERWISE CLASSIFIER CLASS
//************************************************************************
class PMC{
  pt[] Cs = new pt[4];
  pt[] Ps = new pt[4];
  float[] Rs = new float[4];
  boolean[] bs = new boolean[4];
  float dp1, dp2;
  
  PMC(pt[] _Cs, pt[] _Ps, float[] _Rs){
    Cs = _Cs;
    Ps = _Ps;
    Rs = _Rs;
    bs[0] = bs[1] = true;
    bs[2] = dot(V(Cs[0], Ps[0]), V(Ps[0], Ps[2])) < 0;
    bs[3] = dot(V(Cs[1], Ps[1]), V(Ps[1], Ps[3])) < 0; 
    
    dp1 = dot(V(Ps[0], Cs[1]), R(V(Ps[0], Ps[3])));  //  what if dp1 == 0
    dp2 = dot(V(Ps[1], Cs[0]), R(V(Ps[1], Ps[2])));  //  what if dp2 == 0
  }
  
  boolean insideRegion(pt P){
    if(d(P, Cs[0]) <= Rs[0] || d(P, Cs[1]) <= Rs[1]) return true;
    
    float dpt1 = dot(V(Ps[0], P), R(V(Ps[0], Ps[3])));
    if(dpt1 * dp1 < 0) return false;
    
    float dpt2 = dot(V(Ps[1], P), R(V(Ps[1], Ps[2]))); 
    if(dpt2 * dp2 < 0) return false;
    
    boolean tmp1;
    tmp1 = bs[2] ? d(P, Cs[2]) <= Rs[2] : d(P, Cs[2]) > Rs[2];
    if(tmp1 == false) return false;
    
    boolean tmp2;
    tmp2 = bs[3] ? d(P, Cs[3]) <= Rs[3] : d(P, Cs[3]) > Rs[3];
    if(tmp2 == false) return false;    
    
    return true;
  }
}