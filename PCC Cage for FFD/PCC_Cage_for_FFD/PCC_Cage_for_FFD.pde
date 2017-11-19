// LecturesInGraphics: vector interpolation
// Template Author: Jarek ROSSIGNAC
PImage YaohongPix; // picture of author's face, should be: data/pic.jpg in sketch folder
PImage JinPix; // picture of author's face, should be: data/pic.jpg in sketch folder
PImage tex1;  // picture of texture
import static javax.swing.JOptionPane.*;  // to pop up a message window when things go wrong

//**************************** global variables ****************************
pts P = new pts();
pts Q = new pts();
float t=0, f=0;
Boolean animate=true, linear=true, circular=true, beautiful=true;
boolean b1=true, b2=true, b3=true, b4=false;
float len=30; // length of arrows
float th=0.02; // threshold, used in bi-arc computation

float unitLenA = 4;  // used in sampling from boundary
int nvTransA = 10;  // used in parameterization, which means the number of points we sample from each transversal
float unitLenB = 2;
int nvTransB = 10;

boolean showBoundary = true;
boolean showMA = false;
boolean showExtMA = false;
boolean showTrans = false;
boolean showExtTrans = false;
boolean showQuadVerts = false;
boolean showQuads = false;
boolean showTwoEndPoints = false;
boolean showInsideCircles = false;

boolean showLabelDisksAtCtrlPts = false;
boolean showArrowsAtCtrlPts = false;

boolean showTime = false;
boolean chkValid = false;

float velocity = 0.005;  //used in animation

//**************************** initialization ****************************
void setup() {               // executed once at the begining 
  size(1200, 800, P2D);      // window size
  frameRate(30);             // render 30 frames per second
  smooth();                  // turn on antialiasing
  YaohongPix = loadImage("data/Yaohong.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  JinPix = loadImage("data/Jin.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  tex1 = loadImage("data/image003.png");  // load texture image
  P.declare().resetOnCircle(6);
  P.loadPts("data/pts_A");
  Q.declare().resetOnCircle(6);
  Q.loadPts("data/pts_003");
  }

//**************************** display current frame ****************************
void draw() {      // executed at each frame
  background(white); // clear screen and paints white background
  if(snapPic) beginRecord(PDF,PicturesOutputPath+"/P"+nf(pictureCounter++,3)+".pdf"); // start recording for PDF image capture
  if(animating) {t+=0.005; if(t>=0.7){ t = 0.0; animating=false; }}
   //<>//
  pt P1 = P.G[0], P2 = P.G[1], P3 = P.G[2], P4 = P.G[3], P5 = P.G[4], P6 = P.G[5];
  
  if(animating){
    P1.add(3*velocity, V(10,1));
    if(t < 0.35) P2.rotate(-2*velocity, P1); else P2.rotate(2*velocity, P1);
    
    if(t < 0.3) P3.add(10*velocity, V(2,1)); else P3.add(10*velocity, V(-2,1));
    if(t < 0.3) P4.rotate(1*velocity, P3); else P4.rotate(-2*velocity, P3);
    
    if(t < 0.4) P5.add(10*velocity, V(1,1)); else P5.add(5*velocity, V(10,1));
    if(t < 0.35) P6.rotate(2*velocity, P5); else P6.rotate(0, P5);  
  }

  ArrayList<pt[]> verticesQuadA = new ArrayList<pt[]>();
  ArrayList<pt[]> verticesQuadB = new ArrayList<pt[]>();
  
  Stroke SA = new Stroke(P.G, true, unitLenA, nvTransA);
  Stroke SB = new Stroke(Q.G, true, unitLenB, nvTransB);
  
  strokeWeight(2);
  
  if(b1){
    image(tex1, 0, 0);
  }
  
  if(b2){
      /*
      // Compute six arcs
      
      // Find arcs whose supporting circles are inside the stroke
      
      // Sample points from boundary
      // Arc length between two points should be constant
      
      // Constructing the medial axis while storing tangency points
      // Tangency points, together with points on the medial axis, are used to define a bunch of hats which can generate circular arcs
      
      // Extend the arcs near starting arc and ending arc
      // Walk on these extended arcs and extend the medial axis we just computed
      
      // Compute vertices of the quad matrix
      */
      
    SA.computeBoundary();
    SA.computeInsideMask();
    if(chkValid) SA.check();
    SA.computeIdxsArcLR();
    SA.computeHatsMA();
    SA.extendHatsMA();
    SA.computeQuads();
    SA.visulize();  //play with '5', '6', ..., '9', '0' and 'b'
    verticesQuadA = SA.quads;
   }  
  
  if(b3){
    SB.computeBoundary();
    SB.computeInsideMask();
    if(chkValid) SB.check();
    SB.computeIdxsArcLR();
    SB.computeHatsMA();
    SB.extendHatsMA();
    SB.computeQuads();
    SB.visulize();
    verticesQuadB = SB.quads;
  }
  
  if(b4){  
    //Texture mapping
    stroke(210);
    int rowA = verticesQuadA.size();
    int colA = nvTransA;
    int rowB = verticesQuadB.size();
    int colB = nvTransB;
    textureMode(IMAGE);
    beginShape(QUADS);
    for(int i = 0; i < rowA - 1; ++i){
      beginShape(QUAD_STRIP); texture(tex1);
      int ii = floor((i*rowB)/rowA);
      for(int j = 0; j < colA; ++j){
        int jj = floor((j*colB)/colA);
        v(verticesQuadA.get(i)[j], verticesQuadB.get(ii)[jj].x, verticesQuadB.get(ii)[jj].y);
        v(verticesQuadA.get(i+1)[j], verticesQuadB.get(ii+1)[jj].x, verticesQuadB.get(ii+1)[jj].y);
      }
      endShape();
    }
    endShape();
  }

  if(showTime) scribeFooter("time is " + t, 26);
  
  if(showLabelDisksAtCtrlPts){
    noFill();
    stroke(black); P.draw(white); // paint empty disks around each control point
    fill(black); 
    label(P1,V(-1,-2),"P1"); label(P2,V(-1,-2),"P2"); label(P3,V(-1,-2),"P3"); 
    label(P4,V(-1,-2),"P4"); label(P5,V(-1,-2),"P5"); label(P6,V(-1,-2),"P6"); 
  }
  
  if(showArrowsAtCtrlPts){
    strokeWeight(3);
    stroke(dgreen);
    vec arrow1 = S(len, U(P1,P2));
    vec arrow2 = S(len, U(P3,P4));
    vec arrow3 = S(len, U(P5,P6));
    arrow1.showArrowAt(P1);
    arrow2.showArrowAt(P3);
    arrow3.showArrowAt(P5);
  }
  
  if(snapPic) {endRecord(); snapPic=false;} // end saving a .pdf of the screen
  
  fill(black); displayHeader();
  if(scribeText && !filming) displayFooter(); // shows title, menu, and my face & name 
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif"); // saves a movie frame 
  change=false; // to avoid capturing movie frames when nothing happens
  }  // end of draw()

//**************************** user actions ****************************
void keyPressed() { // executed each time a key is pressed: sets the "keyPressed" and "key" state variables, 
                    // till it is released or another key is pressed or released
  if(key=='?') scribeText=!scribeText; // toggle display of help text and authors picture
  if(key=='!') snapPicture(); // make a picture of the canvas and saves as .jpg image
  if(key=='`') snapPic=true; // to snap an image of the canvas and save as zoomable a PDF
  if(key=='~') {filming=!filming; } // filming on/off capture frames into folder FRAMES 
  if(key=='a') {animating=true; f=0; t=0;}  
  if(key=='s') P.savePts("data/pts");
  if(key=='l') P.loadPts("data/pts");
  if(key=='1') b1=!b1;
  if(key=='2') b2=!b2;
  if(key=='3') b3=!b3;
  if(key=='4') b4=!b4;
  
  if(key=='5') showMA = !showMA;
  if(key=='6') showExtMA = !showExtMA;
  if(key=='7') showTrans = !showTrans;
  if(key=='8') showExtTrans = !showExtTrans;
  if(key=='9') showQuadVerts = !showQuadVerts;
  if(key=='0') showQuads = !showQuads;
  if(key=='b') showBoundary = !showBoundary;
  if(key=='e') showTwoEndPoints = !showTwoEndPoints;
  if(key=='i') showInsideCircles = !showInsideCircles;
  if(key=='d') showLabelDisksAtCtrlPts = !showLabelDisksAtCtrlPts;
  if(key=='v') showArrowsAtCtrlPts = !showArrowsAtCtrlPts;
  
  if(key=='Q') exit();  // quit application
  change=true; // to make sure that we save a movie frame each time something changes
  }

void mousePressed() {  // executed when the mouse is pressed
  P.pickClosest(Mouse()); // used to pick the closest vertex of C to the mouse
  change=true;
  }

void mouseDragged() {
  if (!keyPressed || (key=='a')) P.dragPicked();   // drag selected point with mouse
  if (keyPressed) {
      if (key=='.') t+=2.*float(mouseX-pmouseX)/width;  // adjust current frame   
      if (key=='t') P.dragAll(); // move all vertices
      if (key=='r') P.rotateAllAroundCentroid(Mouse(),Pmouse()); // turn all vertices around their center of mass
      if (key=='z') P.scaleAllAroundCentroid(Mouse(),Pmouse()); // scale all vertices with respect to their center of mass
      }
  change=true;
  }

//**************************** text for name, title and help  ****************************
String title ="6491 2017 P2: PCC Cage for FFD", 
       name ="Student: Jin Lin,   Yaohong Wu",
       menu="?:(show/hide) help, s/l:save/load control points, a: animate, `:snap picture, ~:(start/stop) recording movie frames, Q:quit",
       guide="click and drag to edit"; // help info

float timeWarp(float f) {return sq(sin(f*PI/2));}