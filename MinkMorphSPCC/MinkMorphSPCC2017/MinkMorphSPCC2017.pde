//*****************************************************************************
// Minkowski Morph, Jarek Rossignac, Oct 2006, revised Oct 2017 for P4 in CS6491
// Minkowski Morph of SPCC, Jin Lin, Oct 2017
//*****************************************************************************
import processing.pdf.*;    // to save screen shots as PDFs
boolean snapPic=false;
int pictureCounter=0;
String PicturesOutputPath="data/PDFimages";
void snapPicture() {saveFrame("data/JPGimages/P"+nf(pictureCounter++,3)+".jpg"); }
Boolean scribeText=true; // toggle for displaying of help text
PImage myFace; // picture of author's face, should be: data/pic.jpg in sketch folder
Boolean debug=true;
Boolean approximateSolution = false;
Boolean showNorm = true;
Boolean showConnection = true; //using when morph a concave SPCC to a convex SPCC
Boolean Movie = false;

int cap=6;                    // max number of vertices
pt[] P = new pt [cap];          // vertices of the editable curve
pt[] Q = new pt [cap];          // vertices of the other curve
pt[] PQ = new pt [cap];          // vertices of the other curve
//int nP=4, nQ=12, nPQ;           // numbers of vertices in the two curves
int nP=6, nQ=6, nPQ;           // numbers of vertices in the two curves
int bi=-1;                      // index of selected mouse-vertex, -1 if none selected
pt Mouse = new pt(0,0);         // current mouse position
pt Last = new pt(0,0);          // last point drawn
float t=0.50, dt=0.01;           // current time amd time increment
color red = color(200, 10, 10), blue = color(10, 10, 200),  green = color(0, 150, 0), orange = color(250, 150, 5), black=#000000;// COLORS
boolean showDots=true, animating=false, showFrames=false;
int pic=0;
PCC A, B;           // piecewise circular curve
PCC mm;            //MinkMorphPCC
int caNum, cbNum;
int[] ca = new int[cap];
int[] cb = new int[cap];
int[] CUT = new int[2];
pt firstP, lastP;
vec firstN, lastN;

void setup() 
  {   
  size(1200, 600, P2D);  // open window
  for (int i=0; i<cap; i++) {P[i]=new pt(0,0); Q[i]=new pt(0,0); PQ[i]=new pt(0,0); };  // create all vertices, just in case
  for (int i=0; i<nP; i++) {P[i].setTo(-sin((i+0.29)*TWO_PI/nP)*(height*0.28)+height/3.0,cos((i+0.39)*TWO_PI/nP)*(height*0.28)+height/2.0);};
  A = new PCC(P);
  A.computeBoundary();
  for (int i=0; i<nQ; i++) {Q[i].setTo(-sin((i)*TWO_PI/nQ)*(height*0.1)+(width-height/3.0),cos((i)*TWO_PI/nQ)*(height*0.2)+height/2.0);};
  Q[5].x = Q[5].x + 100.0;
  B = new PCC(Q);
  B.computeBoundary();
  PFont font = loadFont("ArialMT-24.vlw"); textFont(font, 16);      // load font for writing on the canvas
  myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  smooth();  strokeCap(ROUND);
  }   
 
void draw() 
  { 
  background(255); 
  if(snapPic) beginRecord(PDF,PicturesOutputPath+"/P"+nf(pictureCounter++,3)+".pdf"); // start recording for PDF image capture
  fill(0); if(scribeText) displayHeader();
  if (bi!= -1) {    // snap selected vertex to mouse position during dragging
    if (bi<nP) { P[bi].setToMouse(); A.computeBoundary(); }
    else { Q[bi - nP].setToMouse(); B.computeBoundary(); }
  };                    
  strokeWeight(2);
  //stroke(green); drawCurve(P,nP);                    // draw P curve and its vertices in green
  //stroke(red); drawCurve(Q,nQ);                      // draw Q curve in red
  stroke(green); A.visualize(red, blue, nP);               // draw PCC A and its vertices in green
  stroke(red); B.visualize(red, blue, nQ);                 // draw PCC B curve in red
  if (showDots) drawDots();                          // draw control points of P
  if (animating) {  t+=dt; if ((t>=1)||(t<=0)) dt=-dt; } // change time during animation to go back & forth
  
  caNum = 0; cbNum = 0;
  // judge whether existing concave arc before morphing
  boolean existConcave = false;
  for (int i = 0; i < nP; i++ ){
    if (A.arcs[i].convex == false){
      existConcave = true;
      ca[caNum++] = i;
    }
    if (B.arcs[i].convex == false){
      existConcave = true;
      cb[cbNum++] = i;
    }
  }
  if (existConcave){
      if (caNum>1 || cbNum>0) {
        scribeFooter("Too many concave arcs! 1 in A and 0 in B most!",1); 
        return;
      }
      mm = MM2(A, t, B);
  }
  else mm = MM(A, t, B); // compute morphing curve
  if (showFrames) // show a series of 9 frames of the animation
    {                  
    float dtt=1.0/8;
    for (float tt=dtt; tt<1; tt+=dtt) {
      if (!existConcave) mm = MM(A, tt, B);  // draw each frame using both colors
      else mm = MM2(A, tt, B);
    }
    }                  
  if(snapPic) {endRecord(); snapPic=false;} // end saving a .pdf of the screen
  if(Movie) snapPicture();
  if(scribeText) displayFooter(); // shows title, menu, and my face & name   
  }
  
void drawCurve(pt[] P, int nP) {beginShape();  for (int i=0; i<nP; i++) {P[i].v(); };  endShape(CLOSE); };    
 
void drawDots() {
  stroke(green);
  for (int i=0; i<nP; i++) {P[i].show(5);};
  for (int i=0; i<nQ; i++) {Q[i].show(5);};
};    

void drawTri(pt A, pt B, pt C) {beginShape(); A.v(); B.v(); C.v(); endShape(CLOSE);};    

//**** Minkowski Morph with convex SPCCS
//************************************************************************
PCC MM(PCC A, float t, PCC B){
  PCC mm = new PCC();
  PCC AA = A.copy();
  PCC BB = B.copy();
  AA.sortByNorm(nP);    //sort vertices by norm
  BB.sortByNorm(nQ);
  int nM;      // numbers of vertices in the morph PCC
  
  if (approximateSolution) {      //iterate over 360 orientations
    int iter = 360;
    nM = iter;
    PCC newA = AA.split(nP, iter);  // split A into 360 parts
    PCC newB = BB.split(nQ, iter);  // split B into 360 parts
    newA.drawNorm(0, iter);
    newB.drawNorm(0, iter);
    // morph
    for(int i = 0; i < iter; i++){
      mm.M[i] = L(newA.M[i], newB.M[i], t);
      mm.norm[i] = newA.norm[i];
    }
  }
  else {                  //split the arcs on A and on B
    nM = nP+nQ;
    AA.split(B, nP, nQ);      // split AA by B
    BB.split(A, nQ, nP);      // split BB by A
    AA.sortByNorm(nM);
    BB.sortByNorm(nM);
    AA.drawNorm(0, nM);
    BB.drawNorm(0 ,nM);
    // morph
    for (int i = 0; i < nM; i++) {
      mm.M[i] = L(AA.M[i], BB.M[i], t);
      mm.C[i] = L(AA.C[i], BB.C[i], t);
      mm.norm[i] = AA.norm[i];
      if (i > 0)
        mm.arcs[i-1] = new Arc(mm.C[i-1], mm.M[i-1], mm.M[i], true, true, 30);
    }
    mm.arcs[nM-1] = new Arc(mm.C[nM-1], mm.M[nM-1], mm.M[0], true, true, 30);
  }
  mm.visualize(orange, green, nM);
  mm.drawNorm(0, nM);
  return mm;
}

//**** Minkowski Morph with non-convex SPCC
//************************************************************************
PCC MM2(PCC A, float t, PCC B){    
  PCC mm = new PCC();
  PCC AA = A.copy();
  PCC BB = B.copy();
  Arc[] connectBiArc, connectBiArc2;
  pt p1, p2, p3, p4;
  vec T1, T2, T3, T4;
  boolean mConcave;
  
  // split AA by concave arcs
  for (int i = 0; i < caNum; i++)
    AA.split(A.arcs[ca[i]], nP);   
  int nA = nP + caNum * 2;
  AA.drawNorm(0, nA);
  // cut AA into 4 parts
  PCC A1 = new PCC(), A2 = new PCC(), A3 = new PCC(), A4 = new PCC();
  A1 = AA.cut(nP, ca[0]+2, nP+ca[0], 1, false, false, nP, nP+1);
  A2 = AA.cut(nP, ca[0]-1, ca[0]-nP, -1, true, true, nP+1, ca[0]);
  A2.sortByNorm(A2.len);
  A3 = AA.cut(nP, ca[0], ca[0], 1, true, true, ca[0], (ca[0]+1)%nP);
  A4 = AA.cut(nP, ca[0]+2, ca[0]+nP, 1, true, true, (ca[0]+1)%nP, nP);
  
  // split BB by concave arcs
  BB.sortByNorm(nQ);
  for (int i = 0; i < caNum; i++)
    BB.split(A.arcs[ca[i]], nQ);
  int nB = nQ + caNum * 2;
  BB.drawNorm(0, nB);
  // cut BB into 2 parts
  PCC B1 = new PCC(), B2 = new PCC();
  B1 = BB.cut(nQ, CUT[0]+1, CUT[0]+nQ, 1, false, false, nQ, nQ+1);
  B2 = BB.cut(nQ, CUT[1]+1, CUT[0]+1, 1, true, true, nQ+1, nQ);  
  
  // part1: morph A1 and B1
  int nM = A1.len + B1.len;
  A1.split(B1, A1.len, B1.len);      // split A1 by B1
  B1.split(A1, B1.len, A1.len);      // split B1 by A1
  A1.sortByNorm(nM);
  B1.sortByNorm(nM);
  A1.drawNorm(0, nM);
  B1.drawNorm(0, nM);
  int numM = 0;
  numM = partMorph(mm, t, A1, B1, nM, AA.norm[nP+1], numM, -1, true, true);
  mm.drawNorm(0, numM);
  
  // part2: morph A2 and B2
  numM = partMorph(mm, t, A2, BB, A2.len, AA.norm[ca[0]], numM, nQ+1, true, true);
  vec CA = V(mm.arcs[numM-1].C, mm.arcs[numM-1].A);
  p1 = P(mm.arcs[numM-1].C, R(CA, mm.arcs[numM-1].w));
  T1 = R(U(mm.arcs[numM-1].C, p1));
  
  // part3: morph A3 and B2
  A3.C[1] = P(A3.C[0]);
  A3.sortByNorm(A3.len);
  A3.split(B2, A3.len, B2.len, -1);      // split A3 by B2
  B2.split(A3, B2.len, A3.len, 1);      // split B2 by A3
  A3.len = A3.len + B2.len;
  for (int i = 0; i < A3.len-1; i++) {
    A3.M[i] = A3.M[i+1]; A3.C[i] = A3.C[i+1]; A3.norm[i] = A3.norm[i+1]; 
    B2.M[i] = B2.M[i+1]; B2.C[i] = B2.C[i+1]; B2.norm[i] = B2.norm[i+1]; 
  }
  A3.len = A3.len - 2;
  B2.len = A3.len;
  A3.sortByNorm(A3.len);
  B2.sortByNorm(B2.len);
  A3.drawNorm(0, A3.len);
  B2.drawNorm(0, B2.len);
  numM = partMorph(mm, t, A3, B2, A3.len, AA.norm[nP], numM, -1, true, false);
  
  // connect part2 and part3
  p2 = P(lastP);
  T2 = R(lastN);
  connectBiArc = findBiArc(p1, T1, p2, T2);
  p3 = P(firstP);
  T3 = R(firstN);
  mConcave = d(p1, p2) < d(p1, p3);
  
  // part4: morph A4 and B2
  numM = partMorph(mm, t, A4, BB, A4.len, AA.norm[nP], numM, nQ, false, true);
  // connect part3 and part4
  CA = V(mm.arcs[numM-1].C, mm.arcs[numM-1].A);
  p4 = P(mm.arcs[numM-1].C, R(CA, mm.arcs[numM-1].w));
  T4 = R(U(mm.arcs[numM-1].C, p4));
  connectBiArc2 = findBiArc(p3, T3, p4, T4);
  
  // visualize
  if (!mConcave){
    connectBiArc = findBiArc(p1, T1, p3, T3);
    connectBiArc2 = findBiArc(p2, T2, p4, T4);
  }
  if (showConnection){
    stroke(black);
    connectBiArc[0].show(); connectBiArc[1].show();
    connectBiArc2[0].show(); connectBiArc2[1].show();
  }
  
  mm.visualize(orange, green, numM);
  return mm;
}

// morph part of PCC, return index
int partMorph(PCC mm, float t, PCC P1, PCC P2, int nm, vec bet, int numM, int i2, boolean CW, boolean CONVEX){
  int ind = numM;
  for (int i = 0; i < nm; i++) 
  if (!between(bet, P1.norm[i], bet, true))
  {
    int ii = i;
    if (i2 > 0) ii = i2;
    mm.M[numM] = L(P1.M[i], P2.M[ii], t);
    mm.C[numM] = L(P1.C[i], P2.C[ii], t);
    mm.norm[numM] = P1.norm[i];
    if (i2 < 0) ii = (i+1)%nm;
    pt next = L(P1.M[(i+1)%nm], P2.M[ii], t);
    if (!equal(mm.M[numM], next)) 
    {
      if (CW) mm.arcs[numM] = new Arc(mm.C[numM], mm.M[numM], next, CW, CONVEX, 30);
      else mm.arcs[numM] = new Arc(mm.C[numM], next, mm.M[numM], CW, CONVEX, 30);
      lastP = P(next);
      lastN = V(P1.norm[(i+1)%nm]);
      numM++;
    }
  }
  firstP = P(mm.M[ind]);
  firstN = V(mm.norm[ind]);
  return numM;
}

boolean equal(pt a, pt b){
  return (abs(a.x - b.x)<=(1.0E-1) && abs(a.y - b.y)<=(1.0E-1));
}

void morph(pt pA, pt pB, pt C, float t) // draws edge (A,B) scaled by t towards C
  {
  pt A = P(pA, t, V(pA, C));
  pt B = P(pB, t, V(pB, C));
  edge(A, B);
  }   
   
void copyPtoQ () {for (int i=0; i<nP; i++) {Q[i].setTo(P[i]);}; nQ=nP;}; // in case the user wants to copy P to Q by pressing c
void swapPandQ () {
  for (int i=0; i<nP; i++) {PQ[i].setTo(P[i]);}; nPQ=nP;
  for (int i=0; i<nQ; i++) {P[i].setTo(Q[i]);}; nP=nQ;
  for (int i=0; i<nPQ; i++) {Q[i].setTo(PQ[i]);}; nQ=nPQ;
 }; 

//**************************************
// **** GUI FOR EDITING THE CONTROL CURVE
//**************************************
void keyPressed() {  
  if(key=='?') scribeText=!scribeText; // toggle display of help text and authors picture
  if(key=='~') snapPic=true; // to snap an image of the canvas and save as zoomable a PDF
  if(key=='!') snapPicture(); // make a picture of the canvas and saves as JPG image]
  if (key=='d') showDots=!showDots;
  if (key=='a') {animating=!animating; };    // reset time
  if (key=='f') {showFrames=!showFrames; if(!showFrames) t=0.5;};    // reset time
  if (key=='m') {t=0.5;};    // freeze time
  if (key==CODED) {if(keyCode==UP) dt*=2; if(keyCode==DOWN) dt/=2;}; // accelerate decelerate animation
  if (key=='s') {approximateSolution = !approximateSolution;};
  if (key=='n') {showNorm = !showNorm;};
  if (key=='c') {showConnection = !showConnection;};
  if (key=='v') {Movie = !Movie;};
  };       

void mouseDragged()                                                    // to do when mouse is pressed  
  {
  if(keyPressed && key=='t') {
    float bd=9999999; 
    int j = 0;
    for (int i=0; i<nP; i++) { 
      if (d2(P[i]) < bd) {                                                // select closest vertex
        bd = d2(P[i]); 
        j = i;  
      };
    };   
    for (int i=0; i<nQ; i++) { 
      if (d2(Q[i]) < bd) {                                                // select closest vertex
        bd = d2(Q[i]); 
        j = i + nP;  
      };
    }; 
    if (j < nP)
    {
      for (int i=0; i<nP; i++) A.ctrlPts[i].add(V(mouseX-pmouseX,mouseY-pmouseY));
      A.computeBoundary();
    }
    else 
    {
      for (int i=0; i<nQ; i++) B.ctrlPts[i].add(V(mouseX-pmouseX,mouseY-pmouseY));
      B.computeBoundary();
    }
  }
  }
  
void mousePressed()                                                    // to do when mouse is pressed  
  {
  Mouse.setToMouse();                                                   // save current mouse location
  float bd=900;                                                     // init square of smallest distance to selected point
  for (int i=0; i<nP; i++) { 
    if (d2(P[i]) < bd) {                                                // select closest vertex
      bd = d2(P[i]); 
      bi = i;  
    };
  };   
  for (int i=0; i<nQ; i++) { 
    if (d2(Q[i]) < bd) {                                                // select closest vertex
      bd = d2(Q[i]); 
      bi = i + nP;  
    };
  }; 
  if (bd>10)                                                         // if closest vertex is too far
    {                                                         
    //bd=600;                                                       // reinitilize distance squared
    //for (int i=0; i<nP; i++) {if (dm(i)<bd) {bd=dm(i); bi=i;};};   // closest mid-edge point
    //if (bd<20) 
    //  {  
    //  for (int i=nP-1; i>bi; i--) P[i+1].setTo(P[i]);                // shift down the rest 
    //  bi++;  P[bi].setTo(Mouse); nP++;                                  // insert new vertex at mouse position
    //  }
    //else bi=-1;     // nothing selected
    bi=-1;
    };
  }

//float d(int j) {return sqrt(d2(j));};                     //  squared distance from mouse to vertex P[j]
float d2(pt v) {return (d(Mouse, v));};                     //  squared distance from mouse to vertex v
//float dm(int j) {return (d(Mouse, P(P[j],P[in(j)])));};   // squared distance from mouse to mid-edge point

void mouseReleased() // do this when mouse released
  {                                      
//  if ( (bi!=-1) &&  P[bi].isOut() )  // if outside of port
//    {                 
//    for (int i=bi; i<nP; i++) P[i].setTo(P[in(i)]);        // shift up to delete selected vertex
//    nP--; 
//    println("deleted vertex "+bi);
//    };                                           // reduce vertex vn
  bi= -1;   
  };
 
                                                          

//**************************************
//**** next and previous functions for polyloops P and Q
//**************************************
int in(int j) {  if (j==nP-1) {return (0);}  else {return(j+1);}  };  // next vertex in control loop
int ip(int j) {  if (j==0) {return (nP-1);}  else {return(j-1);}  };  // previous vertex in control loop
int jn(int j) {  if (j==nQ-1) {return (0);}  else {return(j+1);}  };  // next vertex in control loop
int jp(int j) {  if (j==0) {return (nQ-1);}  else {return(j-1);}  };  // previous vertex in control loop

//**************************************
//**** Title and help text
//**************************************
String title ="6491 2017 P4: Minkowski Morph of Smooth Piecewise-Circular Loops", 
       name ="Student: Jin Lin",
       menu="?:(show/hide) help, a: animate, f:show frames, s:use approximate solution, n:show norms, c:show connection~/!:snap PDF/JPG",
       guide="click and drag to edit green vertices"; // help info
void scribe(String S, float x, float y) {fill(0); text(S,x,y); noFill();} // writes on screen at (x,y) with current fill color
void scribeHeader(String S, int i) { text(S,10,20+i*20); noFill();} // writes black at line i
void scribeHeaderRight(String S) {fill(0); text(S,width-8*S.length()-10,20); noFill();} // writes black on screen top, right-aligned
void scribeFooter(String S, int i) {fill(0); text(S,10,height-10-i*20); noFill();} // writes black on screen at line i from bottom
void scribeAtMouse(String S) {fill(0); text(S,mouseX,mouseY); noFill();} // writes on screen near mouse
void displayHeader() { // Displays title and authors face on screen
    scribeHeader(title,0); scribeHeaderRight(name); 
    image(myFace, width-myFace.width/2,25,myFace.width/2,myFace.height/2); 
    }
void displayFooter() { // Displays help text at the bottom
    scribeFooter(guide,1); 
    scribeFooter(menu,0); 
    }