// LecturesInGraphics: vector interpolation
// Template Author: Jarek ROSSIGNAC
PImage YaohongPix; // picture of author's face, should be: data/pic.jpg in sketch folder
PImage JinPix; // picture of author's face, should be: data/pic.jpg in sketch folder

//**************************** global variables ****************************
pts P = new pts();
float t=0.5, f=0;
Boolean animate=true, linear=true, circular=true, beautiful=true;
boolean b1=true, b2=true, b3=true, b4=true;
float len=200; // length of arrows

//**************************** initialization ****************************
void setup() {               // executed once at the begining 
  size(800, 800, P2D);            // window size
  frameRate(30);             // render 30 frames per second
  smooth();                  // turn on antialiasing
  YaohongPix = loadImage("data/Yaohong.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  JinPix = loadImage("data/Jin.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  P.declare().resetOnCircle(4);
  P.loadPts("data/pts_0");
  }

//**************************** display current frame ****************************
void draw() {      // executed at each frame
  background(white); // clear screen and paints white background
  if(snapPic) beginRecord(PDF,PicturesOutputPath+"/P"+nf(pictureCounter++,3)+".pdf"); // start recording for PDF image capture

  if(animating) {t+=0.01; if(t>=1) {t=1; animating=false;}} 

  strokeWeight(2);
  pt S=P.G[0], E=P.G[1], L=P.G[2], R=P.G[3]; // named points defined for convenience
  stroke(black); edge(S,L); edge(E,R);
  float s=d(S,L), e=d(E,R); // radii of control circles computged from distances
  CIRCLE Cs = C(S,s), Ce = C(E,e); // declares circles
  stroke(dgreen); Cs.drawCirc(); stroke(red); Ce.drawCirc(); // draws both circles in green and red
 

  strokeWeight(5);
  pt cL = new pt(), cR = new pt();
  pt pL = new pt(), pR = new pt();
  
  
  if(b1)
    {
    // your code for part 1
    pts points = new pts();
    points.declare();  
    
    stroke(green);
    if(!parallel(V(S, E), V(S, L))){ 
      pt Hs = findHatMiddlePoint(S, E, e, L);
      pt[] ps = findTangentPoints(Hs, E, e);
      pL = findGoodOne(S, E, L, ps);
      cL = findCircleCenter(L, Hs, pL);
      drawCircleMinorArc(cL, L, pL);
      getPointsOnMinorArc(points, cL, L, pL, false);
      
      //debug
      /*
      stroke(grey); edge(Hs, L);
      stroke(grey); edge(Hs, pL);
      stroke(grey); edge(E, pL);
      stroke(grey); edge(E, Hs);
      stroke(cyan); drawCircleArcInHat(L, Hs, ps[0]);
      stroke(magenta); drawCircleArcInHat(L, Hs, ps[1]);
      //stroke(blue); edge(S, E);
      stroke(green); edge(L, E);
      
      stroke(black); strokeWeight(1); 
      show(Hs, 13); show(pL, 13);
      label(Hs, V(13, -13),"H"); label(pL, V(13, -13), "L'");
      strokeWeight(5);
      */
      /*
      stroke(grey); 
      edge(cL, L); edge(cL, pL); edge(cL, Hs); edge(Hs, L); edge(Hs, pL);
      stroke(black); strokeWeight(1); 
      show(Hs, 13); show(pL, 13); show(cL, 13);
      label(Hs, V(13, -13),"H"); label(pL, V(13, -13), "L'"); label(cL, V(-13, -13), "U");
      strokeWeight(5);
      */
    }
    else{
      scribeHeader("S,E,L on the same line!", 3);
      pL = P(E, e, U(V(S,E)));
      cL = P(L, pL);  //in this case, the center is (L+pL)/2
      if(dot(R(V(cL, pL)), V(cL, R)) >= 0) {
        drawCircleMinorArc(cL, L, pL);
        getPointsOnMinorArc(points, cL, L, pL, false);
      }
      else {
        drawCircleMinorArc(cL, pL, L);
        getPointsOnMinorArc(points, cL, L, pL, true);
      }
      
      //debug
      /*
      pt ppL = P(E, -e, U(V(S,E)));
      pt ccL = P(L, ppL);
      stroke(cyan); drawCircleMinorArc(cL, L, pL);
      stroke(magenta); drawCircleMinorArc(cL, pL, L);
      stroke(green); drawCircleMinorArc(ccL, L, ppL);
      stroke(yellow); drawCircleMinorArc(ccL, ppL, L);
      
      
      stroke(black); strokeWeight(1); 
      show(pL, 13); show(ppL, 13);
      label(pL, V(13, -13), "L'"); label(ppL, V(13, -13), "L''");
      strokeWeight(5);
      */
    }
    
    getPointsOnSEArc(points, E, pL, R, S);
    
    
    stroke(red);
    if(!parallel(V(E, S), V(E, R))){
      pt He = findHatMiddlePoint(E, S, s, R);
      pt[] pe = findTangentPoints(He, S, s);
      pR = findGoodOne(E, S, R, pe);
      cR = findCircleCenter(R, He, pR);
      drawCircleMinorArc(cR, R, pR);
      getPointsOnMinorArc(points, cR, R, pR, false);
    }
    else{
      scribeHeader("E,S,R on the same line!", 3);
      pR = P(S, s, U(V(E,S)));
      cR = P(R, pR);  //in this case, the center is (R+pR)/2
      if(dot(R(V(cR, pR)), V(cR, L)) >= 0) {
        drawCircleMinorArc(cR, R, pR);
        getPointsOnMinorArc(points, cR, R, pR, false);
      }
      else {
        drawCircleMinorArc(cR, pR, R);
        getPointsOnMinorArc(points, cR, R, pR, true);
      }
    }
    
    getPointsOnSEArc(points, S, pR, L, E);
    
    // debug
    /*
    stroke(cyan); drawCircleMinorArc(S, pR, L);
    stroke(black); strokeWeight(1);
    show(pR, 13); label(pR, V(13,13), "R'");
    strokeWeight(5);
    */
    
    /*
    pts points1 = new pts();
    pts points2 = new pts();
    pts points3 = new pts();
    pts points4 = new pts();
    points1.declare();
    points2.declare();
    points3.declare();
    points4.declare();
    
    getPointsOnSEArc(points1, S, pR, L, E); //<>//
    getPointsOnMinorArc(points2, cL, L, pL, false); //<>//
    getPointsOnSEArc(points3, E, pL, R, S); //<>//
    getPointsOnMinorArc(points4, cR, R, pR, false); //<>//
    points1.drawCurve(cyan);
    points2.drawCurve(red);
    points3.drawCurve(magenta);
    points4.drawCurve(green);
    */
    
    stroke(black);
    /*
    getPointsOnSEArc(points, S, pR, L, E);
    getPointsOnMinorArc(points, cL, L, pL, false);
    getPointsOnSEArc(points, E, pL, R, S);
    getPointsOnMinorArc(points, cR, R, pR, false);
    */
    points.drawCurve(pink);
    
  }
    
  if(b2)
    {
      // your code for part 2
      pt[] Centers = {S, E, cL, cR};
      pt[] Points = {L, R, pL, pR};
      float[] Radii = {s, e, d(pL, cL), d(pR, cR)};
      
      PMC classifier = new PMC(Centers, Points, Radii);
      pt p = Mouse();
      if(classifier.insideRegion(p)){
        stroke(red);
      }
      else{
        stroke(blue);
      }
      p.show();
      
      // debug
      /*
      stroke(grey);
      edge(L, pR); edge(R, pL);
      edge(cL, L); edge(cL, pL);
      edge(cR, R); edge(cR, pR);
      stroke(black); strokeWeight(1);
      show(pL, 13); label(pL, V(13,-13), "L'");
      show(pR, 13); label(pR, V(13,13), "R'");
      show(cL, 13); label(cL, V(13,-13), "U");
      show(cR, 13); label(cR, V(13,13), "D");
      strokeWeight(5);
      */
      
    }

   stroke(black);   strokeWeight(1);

   if(b3)
     {
     /*
     fill(black); scribeHeader("t="+nf(t,1,2),2); noFill();
     // your code for part 4
      strokeWeight(3); stroke(blue); 
     //    drawCircleInHat(Mr,M,Ml);  
     */
     }
   strokeWeight(1);
  
  noFill(); stroke(black); P.draw(white); // paint empty disks around each control point
  fill(black); label(S,V(-1,-2),"S"); label(E,V(-1,-2),"E"); label(L,V(-1,-2),"L"); label(R,V(-1,-2),"R"); noFill(); // fill them with labels
  
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
String title ="6491 2017 P1: Caplets", 
       name ="Student: Jin Lin,   Yaohong Wu",
       menu="?:(show/hide) help, s/l:save/load control points, a: animate, `:snap picture, ~:(start/stop) recording movie frames, Q:quit",
       guide="click and drag to edit, press '1' or '2' to toggle LINEAR/CIRCULAR,"; // help info


  
float timeWarp(float f) {return sq(sin(f*PI/2));}