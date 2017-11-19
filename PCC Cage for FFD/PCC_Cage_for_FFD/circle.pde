//************************************************************************
//**** CIRCLES
//************************************************************************
CIRCLE C(pt P, float r) {return new CIRCLE(P,r);} // shortcut for creating circle

class CIRCLE 
  { 
  pt C=P(200,200); // center
  float r=200;     // radius
  CIRCLE (pt P, float s) {C.setTo(P); r=s;}; // CREATE  circle
  CIRCLE () {} 
  void drawCirc() {show(C,r);} // draw
  pt ProjectionOf(pt P) {return T(C,r,P);} // closest projection of point P onto Circle, T(C,r,P) returns point at distance r from C towards P
  pt PtOnCircle(float a) {return P(r*cos(a)+C.x , r*sin(a)+C.y);}   // Computes point on this circle (for animation or drawing)
  } //*********** END CIRCLE CLASS


void drawCircleArcInHat(pt PA, pt B, pt PC) // draws circular arc from PA to PC that starts tangent to B-PA and ends tangent to PC-B
  {
  float e = (d(B,PC)+d(B,PA))/2;
  pt A = P(B,e,U(B,PA));
  pt C = P(B,e,U(B,PC));
  vec BA = V(B,A), BC = V(B,C);
  float d = dot(BC,BC) / dot(BC,W(BA,BC));
  pt X = P(B,d,W(BA,BC));  // X is the center of the circle decided by A, B, C
  //float r=abs(det(V(B,X),U(BA)));
  vec XA = V(X,A), XC = V(X,C); 
  float a = angle(XA,XC), da=a/60;
  beginShape();
   if(a>0) for (float w=0; w<=a; w+=da) v(P(X,R(XA,w))); 
   else for (float w=0; w>=a; w+=da) v(P(X,R(XA,w)));
  endShape();
  }
  
  

  