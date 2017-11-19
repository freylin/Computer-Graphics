//  ******************* Tango dancer 3D 2017 ***********************
//  Yaohong Wu, Jin Lin
//  Last updated: 10/02/2017
Boolean 
  animating=true, 
  PickedFocus=false, 
  center=true, 
  track=false, 
  showViewer=false, 
  showBalls=false, 
  showControl=true, 
  showCurve=true, 
  showPath=true, 
  showKeys=true, 
  showSkater=false,  //if true, show dancer, otherwise show arrow
  scene1=false,
  solidBalls=false,
  showCorrectedKeys=true,
  showQuads=true,
  showVecs=true,
  showTube=false,
  useOffsetCurve=true,  //if true, use offset curve, otherwise subdivision curve
  constStepLength=false,  //if true, step length is a constant
  doRetrofitting=false;  //if true, do retrofitting before subdivising
  
float 
  t=0, 
  s=0;
int
  f=0, maxf=2*30, level=3, method=4, retrofitLevel=2;
int ff = 0;
int df = 5;
String SDA = "angle";
float defectAngle=0;

float h = 150;  //height of the upper curve
float grav = 10;  //gravitational constant
boolean useCurvature = true;  //whether or not use curvature to compute centripetal force
float[] ts = {0, 1, 0.5, 0, -0.5};
float stepLength = 20;

pts P = new pts(); // polyloop in 3D
pts Q = new pts(); // subdivision polyloop in 3D
pts R = new pts(); // subdivision polyloop (from Q) in 3D
pts T = new pts(); // trajectory on the floor (move R out a little bit)

Step curStep = new Step();  // (i=0, t=0) when initialized


void setup() {
  myFace = loadImage("data/YaohongAndJin.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  textureMode(NORMAL);          
  size(900, 900, P3D); // P3D means that we will do 3D graphics
  //size(600, 600, P3D); // P3D means that we will do 3D graphics
  P.declare(); Q.declare(); R.declare(); T.declare(); // P is a polyloop in 3D: declared in pts
  //P.resetOnCircle(6,100); Q.copyFrom(P); // use this to get started if no model exists on file: move points, save to file, comment this line
  P.loadPts("data/pts");  Q.loadPts("data/pts2"); // loads saved models from file (comment out if they do not exist yet)
  noSmooth();
  frameRate(30);
}

void draw() {
  background(255);
  hint(ENABLE_DEPTH_TEST); 
  pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
  setView();  // see pick tab
  showFloor(); // draws dance floor as yellow mat
  doPick(); // sets Of and axes for 3D GUI (see pick Tab)
  P.SETppToIDofVertexWithClosestScreenProjectionTo(Mouse()); // for picking (does not set P.pv)
  
  R.copyFrom(P);
  
  if(doRetrofitting && method > 0 && method < 5){
    float t = ts[method];
    float s = 1 - t;
    for(int i = 0; i < retrofitLevel; ++i){
      retrofit(P,R,s);
    }
  }
  
  for(int i=0; i<level; i++) {
    Q.copyFrom(R); 
    if(method==5) {Q.subdivideDemoInto(R);}
    if(method==4) {Q.subdivideQuinticInto(R);}
    if(method==3) {Q.subdivideCubicInto(R); }
    if(method==2) {Q.subdivideJarekInto(R); }
    if(method==1) {Q.subdivideFourPointInto(R);}
    if(method==0) {Q.subdivideQuadraticInto(R); }
  }
  
  if(useOffsetCurve){
    displaceCurve(R);
    T.drawClosedCurve(3, blue);
  }
  
  R.displaySkater();
  
  if(showCurve) Q.drawClosedCurve(3, magenta);  //draw Q, the previous state of R (Q->R using one subdivision step)
  if(showControl) P.drawClosedCurve(3, red); // draw control polygon P 
  fill(yellow,100); P.showPicked(); // draw the vertex of P we pick
  
  
  //if(animating)
  //  {
  //  f++; // advance frame counter
  //  if (f>maxf) // if end of step
  //    {
  //    P.next();     // advance dv in P to next vertex
 ////     animating=true;  
  //    f=0;
  //    }
  //  }
  //t=(1.-cos(PI*f/maxf))/2; //t=(float)f/maxf;

  //if(track) F=_LookAtPt.move(X(t)); // lookAt point tracks point X(t) filtering for smooth camera motion (press'.' to activate)
 
  popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
  hint(DISABLE_DEPTH_TEST); // no z-buffer test to ensure that help text is visible
  
  if(method==4) scribeHeader("Quintic UBS",2);
  if(method==3) scribeHeader("Cubic UBS",2);
  if(method==2) scribeHeader("Jarek J-spline",2);
  if(method==1) scribeHeader("Four Points",2);
  if(method==0) scribeHeader("Quadratic UBS",2);
  if(doRetrofitting) scribeHeader("retrofitting " + retrofitLevel + " times (cannot be applied to quadratic)" , 3);
  if(useOffsetCurve) scribeHeader("run on the offset curve whose computation is based on discrete " + (useCurvature?"curvature":"acceleration"), 4);
  if(constStepLength) scribeHeader("run at a constant speed: " + stepLength, 5);
  
  // used for demos to show red circle when mouse/key is pressed and what key (disk may be hidden by the 3D model)
  if(mousePressed) {stroke(cyan); strokeWeight(3); noFill(); ellipse(mouseX,mouseY,20,20); strokeWeight(1);}
  if(keyPressed) {stroke(red); fill(white); ellipse(mouseX+14,mouseY+20,26,26); fill(red); text(key,mouseX-5+14,mouseY+4+20); strokeWeight(1); }
  if(scribeText) {fill(black); displayHeader();} // dispalys header on canvas, including my face
  if(scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif");  // save next frame to make a movie
  change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  change=true;
}