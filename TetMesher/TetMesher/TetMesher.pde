//  ******************* LITM: Layer-Interpolating Tet Mesh, 2017 *********************** //<>//

import java.util.*;  // use data structures, like Queue, HashMap

boolean
  animating=true,
  PickedFocus=false, 
  center=true, 
  
  showBalls = true,
  showTubes = true, 
  flipped = false;

float h_floor = 0, h_ceiling = 600, h = h_floor;
float rb = 30, rt = rb/2; // radius of the balls and tubes
float r = 300;  // radius for ball pivoting

pts P = new pts(); // polyloop in 3D
pts Q = new pts(); // second polyloop in 3D
pts R, S;

boolean useTetMesh = false;  // whether or not generate the edges of the tetrahedron mesh for part 1
boolean useOctree = true;  // whether or not use Octree
boolean useRadius = true;  // whether or not compute an appropriate radius of pivoting ball

boolean debugBPA = true;  // set debugBPA when you want to run BPA on toy examples, there are 7 examples in total, see line 145 in this file
boolean debugOctree = false;
boolean debugRadius = false;

boolean showMesh = true;  // show the reconstructed mesh
boolean showBallPath = false;  
boolean showInnerEdges = false;
boolean showBorderEdges = false;
boolean showVertices = false;  // show the sampled vertices 
boolean showVertexTypes = false;
boolean showPivotBall = false;

int seedIdx = -1;  // used to get a good start triangle

int stepToShow = 0;
int maxStep = 0;
pt ballCenterToShow = new pt();

int tetCnt = 0;
int mixBeamCnt = 0;
int siteCnt = 0;
int triCnt = 0;
int avgNumNeighbors = 0;
float duration = 0;
float durationBPA = 0;

ArrayList<Triangle> trisP = new ArrayList<Triangle>();
ArrayList<Edge> edgesP = new ArrayList<Edge>();
ArrayList<Triangle> trisQ = new ArrayList<Triangle>();
ArrayList<Edge> edgesQ = new ArrayList<Edge>();
ArrayList<Edge> edges31 = new ArrayList<Edge>();
ArrayList<Edge> edges13 = new ArrayList<Edge>();
ArrayList<Edge> edges22 = new ArrayList<Edge>();

ArrayList<Vertex> vertices = new ArrayList<Vertex>();
ArrayList<Triangle> mesh = new ArrayList<Triangle>();
ArrayList<Edge> borderEdges = new ArrayList<Edge>();
ArrayList<Edge> innerEdges = new ArrayList<Edge>();
TriangleMesh tm = new TriangleMesh();
String triMeshPathS = "data/tri_mesh_unnamed";  // the path to save the triangle mesh we reconstruct
String triMeshPathL = "data/tri_mesh_unnamed";  // the path of the triangle mesh we want to load

Cell root = null;  // the root of the octree
int rootDepth = 2;  // default value, the octree should has depth at least rootDepth
float rootSize = 0;  // default value


void setup() {
  YaohongPix = loadImage("data/Yaohong.jpg");
  JinPix = loadImage("data/Jin.jpg");
  textureMode(NORMAL);
  size(900, 900, P3D); // P3D means that we will do 3D graphics
  
  P.declare(); 
  Q.declare();
  //P.resetOnCircle(6,100); Q.copyFrom(P); // use this to get started if no model exists on file: move points, save to file, comment this line
  
  P.loadPts("data/pts_P"); // loads saved models from file (comment out if they do not exist yet)
  Q.loadPts("data/pts_Q");
  siteCnt = P.nv + Q.nv;
  
  noSmooth();
  frameRate(30);
  R=P; 
  S=Q;
}

void draw() {
  int startTime = millis();

  background(255);
  hint(ENABLE_DEPTH_TEST); 
  pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
  setView();  // see pick tab
  showFloor(h); // draws dance floor as yellow mat
  doPick(); // sets Of and axes for 3D GUI (see pick Tab)
  R.SETppToIDofVertexWithClosestScreenProjectionTo(Mouse()); // for picking (does not set P.pv)

  initMarked(siteCnt);
  initPQID();
  tetCnt = 0;
  mixBeamCnt = 0;

  if(useTetMesh) constructTetMesh();

  if(showBalls) {
    fill(orange); P.drawBalls(rb);
    fill(green); Q.drawBalls(rb);
    fill(red, 100); R.showPicked(rb+5);
  }

  if(showTubes) {
    fill(orange); drawEdges(edgesP, rt);
    fill(green); drawEdges(edgesQ, rt);

    fill(grey); drawEdges(edges31, rt);
    fill(grey); drawEdges(edges13, rt);
    fill(grey); drawEdges(edges22, rt);
  }

  if(debugRadius){
    int sel = 0;
    switch(sel){
      case 0: debugRadiusVShape(); break;
      case 1: debugRadiusTetMesh(); break;
    }
  }

  if(debugOctree){
    int sel = 0;
    switch(sel){
      case 0: debugOctreeFloor(); break;
      case 1: debugOctreeTerrain(); break;
      case 2: debugOctreeCylinder(); break;
      case 3: debugOctreeSphere(); break;
      case 4: debugOctreeVShape(); break;
    }
  }
  
  if (debugBPA) {
    int sel = 5;
    switch(sel) {
    case 0:
      BPAOnFloor();
      break; // when samples are on the floor
    case 1: 
      BPAOnTerrain(); 
      break; // when samples are on a terrain, i.e. two planes
    case 2: 
      BPAOnCylinder(); 
      break; // when samples are on a beam
    case 3: 
      BPAOnSphere(); 
      break; // when samples are on a ball
    case 4:
      BPAOnVShape(); 
      break; // when samples are on the "V" shape, it can take about 0.2 second if radius is small (around 15)
    case 5:
      BPAOnTet();  // when samples are on a tetrahedron (3 vertices on the floor and 1 vertex on the ceiling)
      break;
    case 6:
      BPAOnTetMesh();  // when samples are on surface of beams and balls of a tet mesh, it takes less than 20 seconds, about 60k triangles!
      break;
    }
  }
  
  if(debugBPA == false && showMesh == true){
    tm.loadTriMesh(triMeshPathL);
    fill(red); tm.show();
  }

  if(root != null && r > 0){
    avgNumNeighbors = avgNumberNeighbors(root, vertices, r);
  }

  int endTime = millis();
  duration = (endTime - startTime) / 1000.0;

  popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
  hint(DISABLE_DEPTH_TEST); // no z-buffer test to ensure that help text is visible

  //*** TEAM: please fix these so that they provice the correct counts
  scribeHeader("Site count: " + P.nv +" floor + " + Q.nv + " ceiling", 1);
  scribeHeader("Beam count: " + edgesP.size() + " floor + " + edgesQ.size() +" ceiling +" + mixBeamCnt + " mixed", 2);
  scribeHeader("Tet count: " + tetCnt, 3);

  int nVertices = vertices.size() > 0 ? vertices.size() : tm.nv;
  int nTriangles = mesh.size() > 0 ? mesh.size() : tm.nt;
  scribeHeader("Sample count: " + nVertices, 4);
  scribeHeader("Triangle count: " + nTriangles, 5);
  scribeHeader("Second per frame: " + duration, 6);
  scribeHeader("Radius of rolling ball: " + r, 7);
  
  if(useOctree) scribeHeader("use Octree", 8); else scribeHeader("not use Octree", 8);
  if(root != null) scribeHeader("root depth: " + root.depth + ", root size: " + root.size, 9);
  scribeHeader("average number of neighboring vertices: " + avgNumNeighbors, 10);
  
  scribeHeader("duration for BPA: " + durationBPA, 11);
  scribeHeader("number of border edges: " + borderEdges.size(), 12);
  scribeHeader("number of inner edges: " + innerEdges.size(), 13);
  
  // used for demos to show red circle when mouse/key is pressed and what key (disk may be hidden by the 3D model)
  if (mousePressed) {
    stroke(cyan); 
    strokeWeight(3); 
    noFill(); 
    ellipse(mouseX, mouseY, 20, 20); 
    strokeWeight(1);
  }
  if (keyPressed) {
    stroke(red); 
    fill(white); 
    ellipse(mouseX+14, mouseY+20, 26, 26); 
    fill(red); 
    text(key, mouseX-5+14, mouseY+4+20); 
    strokeWeight(1);
  }
  if (scribeText) {
    fill(black); 
    displayHeader();
  } // dispalys header on canvas, including my face
  if (scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if (filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++, 4)+".tif");  // save next frame to make a movie
  change=true;
}