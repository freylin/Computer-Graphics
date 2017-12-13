//*****************************************************************************
// TITLE:         Vertices Sampling
// DESCRIPTION:   Sample vertices from surfaces of cylinders and spheres
// AUTHORS:       Yaohong Wu, Jin Lin
// EDITS:         11/28/2017
//*****************************************************************************


// parameters for sampling a cylinder, or a beam
// generate nBottom * (nSide + 1) vertices for a beam
int nBottom = 10;
int nSide = 40;

// parameters for sampling a sphere, or a ball
// generate nLong * (nLat + 1) vertices for a ball
int nLong = 20;
int nLat = 20;

// maximum beam length, used for computing radius of the pivoting ball
float maxBeamLength = 0;

void sampleBeam(pt P, pt Q, float r, ArrayList<pt> G, ArrayList<vec> N){
  vec V = V(P,Q);
  vec I = U(Normal(V));
  vec J = U(N(I,V));
  
  float da = TWO_PI / nBottom, dh = 1.0 / nSide;
  float a = 0, h = 0;
  for(int i = 0; i <= nSide; ++i){
    pt C = P(P, h, V);
    a = 0;
    if (i % 2 > 0) a = PI / nBottom;  // this will improve the performance of BPA
    for(int j = 0; j < nBottom; ++j){
      pt tmp = P(C, r * cos(a), I, r * sin(a), J);
      a += da;
      G.add(tmp);
      N.add(V(C, tmp).normalize());
    }
    h += dh;
  }
  
  maxBeamLength = max(maxBeamLength, d(P, Q));
}

void sampleMultiBeams(ArrayList<Edge> edges, float r, ArrayList<pt> G, ArrayList<vec> N){
  for(Edge e : edges){
    int a = e.a, b = e.b;
    pt A = a < P.nv ? P.G[a] : Q.G[a-P.nv];  // P and Q are global
    pt B = b < P.nv ? P.G[b] : Q.G[b-P.nv];
    sampleBeam(A, B, r, G, N);
  }
}

void sampleBall(pt C, float r, ArrayList<pt> G, ArrayList<vec> N){
  float a = -PI/2, b = 0;
  float da = PI / nLat, db = TWO_PI / nLong;
  a += da;
  vec X = V(1, 0, 0), Y = V(0, 1, 0), Z = V(0, 0, 1);
  G.add(P(C, -r, Z));
  N.add(V(0, 0, -1));
  for(int i = 1; i < nLat; ++i){
    float h = r * sin(a);
    float c = r * cos(a);
    pt D = P(C, h, Z);
    b = 0;
    if(i % 2 == 0) b = db / 2;  // this will improve the performance of BPA
    
    for(int j = 0; j < nLong; ++j){
      pt tmpP = P(D, c*cos(b), X, c*sin(b), Y);
      vec tmpN = V(cos(b), sin(b), 0);
      G.add(tmpP);
      N.add(tmpN);
      b += db;
    }
    a += da;
  }
  G.add(P(C, r, Z));
  N.add(Z);
}


void sampleMultiBalls(ArrayList<pt> centers, float r, ArrayList<pt> G, ArrayList<vec> N){
  for(pt c:centers){
    sampleBall(c, r, G, N);
  }
}


void showSamples(ArrayList<Vertex> P){
  for(Vertex v: P){show(v.position, 4);}
}


// compute a reasonable radius of the pivoting ball
// according to sizes of beams and balls
// and "sampling density" on beams and balls
float getRadiusOfPivotingBall(){
  float a = (TWO_PI / nBottom) * rt;
  float b = maxBeamLength / (nSide-1);  // assume that we have known the maximum beam length
  float c = sqrt(a*a / 4 + b*b);
  
  float d = (TWO_PI / nLong) * rb;
  float e = (PI / nLat) * rb;
  float f = sqrt(d*d / 4 + e*e);
  
  return (max(a, c) + max(d, f)) / 2;
}