// Student's should use this to render their model


void showDancer(pt LeftFoot, float transfer, pt RightFoot, vec Forward, pt HeadForward)
  {
  float footRadius=3, kneeRadius = 6,  hipRadius=12 ; // radius of foot, knee, hip
  float hipSpread = hipRadius; // half-displacement between hips
  float bodyHeight = 100; // height of body center B
  float ankleBackward=10, ankleInward=4, ankleUp=6, ankleRadius=4; // ankle position with respect to footFront and size
  float pelvisHeight=10, pelvisForward=hipRadius/2, pelvisRadius=hipRadius*1.3; // vertical distance form BodyCenter to Pelvis 
  float LeftKneeForward = 20; // arbitrary knee offset for mid (B,H)

  float chestRadius=hipRadius*1.5, chestHeight = 50;
  float headRadius=hipRadius*1.1, headHeight = 30;
  float shoulderRadius=chestRadius*0.5, shoulderHeight = chestHeight+10;
  float elbowRadius=shoulderRadius/3*2, elbowHeight = shoulderHeight/2;
  float wristRadius=elbowRadius/3*2, wristHeight = elbowHeight/6*5;
  float handRadius=wristRadius*1.2, handHeight = handRadius;

  // BODY
  pt BodyProjection = L(LeftFoot,1./3+transfer/3,RightFoot); // floor projection of B
  
  vec Up = V(BodyProjection, HeadForward); // up vector
  Up.normalize();
  
  vec Right = N(Up,Forward); // side vector pointing towards the right
  
  pt BodyCenter = P(BodyProjection,bodyHeight,Up); // Body center
  fill(blue); showShadow(BodyCenter,5); // sphere(BodyCenter,hipRadius);
  //fill(blue); arrow(BodyCenter,V(100,Forward),5); // forward arrow
  
 
 // ANKLES
  pt RightAnkle =  P(RightFoot, -ankleBackward,Forward, -ankleInward,Right, ankleUp,Up);
  fill(red);  
  capletSection(RightFoot,footRadius,RightAnkle,ankleRadius);  
  pt LeftAnkle =  P(LeftFoot, -ankleBackward,Forward, ankleInward,Right, ankleUp,Up);
  fill(green);  
  capletSection(LeftFoot,footRadius,LeftAnkle,ankleRadius);  
  fill(blue);  
  sphere(RightAnkle,ankleRadius);
  sphere(LeftAnkle,ankleRadius);
 
  // FEET (CENTERS OF THE BALLS OF THE FEET)
  fill(blue);  
  sphere(RightFoot,footRadius);
  pt RightToe =   P(RightFoot,5,Forward);
  capletSection(RightFoot,footRadius,RightToe,1);
  sphere(LeftFoot,footRadius);
  pt LeftToe =   P(LeftFoot,5,Forward);
  capletSection(LeftFoot,footRadius,LeftToe,1);

  // HIPS
  pt RightHip =  P(BodyCenter,hipSpread,Right);
  fill(red);  sphere(RightHip,hipRadius);
  pt LeftHip =  P(BodyCenter,-hipSpread,Right);
  fill(green);  sphere(LeftHip,hipRadius);

  // KNEES AND LEGs
  float RightKneeForward = 20;
  pt RightMidleg = P(RightHip,RightAnkle);
  pt RightKnee =  P(RightMidleg, RightKneeForward,Forward);
  fill(red);  
  sphere(RightKnee,kneeRadius);
  capletSection(RightHip,hipRadius,RightKnee,kneeRadius);  
  capletSection(RightKnee,kneeRadius,RightAnkle,ankleRadius);  
   
  pt LeftMidleg = P(LeftHip,LeftAnkle);
  pt LeftKnee =  P(LeftMidleg, LeftKneeForward,Forward);
  fill(green);  
  sphere(LeftKnee,kneeRadius);
  capletSection(LeftHip,hipRadius,LeftKnee,kneeRadius);  
  capletSection(LeftKnee,kneeRadius,LeftAnkle,ankleRadius);  

  // PELVIS
  pt Pelvis = P(BodyCenter,pelvisHeight,Up, pelvisForward,Forward); 
  fill(blue); sphere(Pelvis,pelvisRadius);
  capletSection(LeftHip,hipRadius,Pelvis,pelvisRadius);  
  capletSection(RightHip,hipRadius,Pelvis,pelvisRadius);  

  // CHEST
  pt Chest = P(BodyCenter, chestHeight, Up);
  fill(yellow); sphere(Chest,chestRadius);
  capletSection(Pelvis,pelvisRadius, Chest, chestRadius);  
  
  // HEAD
  pt Head = P(BodyCenter, chestHeight + headHeight, Up);
  fill(darkWood); sphere(Head,headRadius);
  
  // SHOULDER
  pt LeftShoulder = P(BodyCenter, shoulderHeight, Up, shoulderRadius*1.2, Right);
  pt RightShoulder = P(BodyCenter, shoulderHeight, Up, -shoulderRadius*1.2, Right);
  fill(orange); sphere(LeftShoulder,shoulderRadius);
  fill(metal); sphere(RightShoulder,shoulderRadius);
  
  // ELBOW
  vec x = V(RightKnee, LeftKnee);
  pt LeftElbow = P(LeftShoulder, 1, x, -elbowHeight, Up, elbowRadius*4, Right);
  pt RightElbow = P(RightShoulder, -1, x, -elbowHeight, Up, -elbowRadius*4, Right);
  fill(orange); sphere(LeftElbow, elbowRadius);
  capletSection(LeftElbow,elbowRadius,LeftShoulder,shoulderRadius);  
  fill(metal); sphere(RightElbow, elbowRadius);
  capletSection(RightElbow,elbowRadius,RightShoulder,shoulderRadius);  
  
  // WRIST
  x = V(LeftShoulder, LeftElbow);
  pt LeftWrist, RightWrist;
  if (d(x, Forward) > 0)
    LeftWrist = P(LeftElbow, 0.6, x, wristHeight, Forward, 15, Up);
  else LeftWrist = P(LeftElbow, 0.6, x);
  LeftWrist = P(LeftWrist, -5, Right);
  x = V(RightShoulder, RightElbow);
  if (d(x, Forward) > 0)
    RightWrist = P(RightElbow, 0.6, x, wristHeight, Forward, 15, Up);
  else RightWrist = P(RightElbow, 0.6, x);
  RightWrist = P(RightWrist, 5, Right);
  fill(orange); sphere(LeftWrist, wristRadius);
  capletSection(LeftElbow,elbowRadius,LeftWrist,wristRadius);  
  fill(metal); sphere(RightWrist, wristRadius);
  capletSection(RightElbow,elbowRadius,RightWrist,wristRadius); 
  
  // HAND
  x = V(LeftElbow, LeftWrist);
  pt LeftHand = P(LeftWrist, 0.2, x);
  x = V(RightElbow, RightWrist);
  pt RightHand = P(RightWrist, 0.2, x);
  fill(yellow); sphere(LeftHand, handRadius);
  fill(metal); sphere(RightHand, handRadius);
  
  }
  
void capletSection(pt A, float a, pt B, float b) { // cone section surface that is tangent to Sphere(A,a) and to Sphere(B,b)
  coneSection(A,B,a,b);
  }  