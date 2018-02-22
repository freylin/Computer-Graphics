//*****************************************************************************
// TITLE:         Octree
// DESCRIPTION:   Functions of building an Octree, getting neighborhood of a certain position given an Octree
// AUTHORS:       Jin Lin
// EDITS:         02/22/2018
// REASON:        Octree is a kind of poweful and elegent data structure, which can be used to accelerate 3D spatial query, I use it in many projects.
//*****************************************************************************



class Cell{
  Cell[] child = new Cell[8];  //  children
  
  float size; // the edge length of the cube represented by this cell
  int depth; // the depth, i.e. level of subdivision, note that the depth of a leaf is 0
  pt origin;  // origin point, i.e. the corner with min x, y, z coordinates
  
  ArrayList<Vertex> P = new ArrayList<Vertex>();  // store vertices in leaf cells
  
  boolean hasChild = false;
  Cell(){}
  Cell(float _size, int _depth, pt _origin){
    size = _size; depth = _depth; origin = _origin;
  }
}



// find origin point o and return size s
float findOrigin(ArrayList<Vertex> G, pt o){
  int nv = G.size();
  o.setTo(P(G.get(0).position));
  pt maxP = P(G.get(0).position);
  for (int i = 1; i < nv; i++){
    o.x = min(o.x, G.get(i).position.x);
    o.y = min(o.y, G.get(i).position.y);
    o.z = min(o.z, G.get(i).position.z);
    maxP.x = max(maxP.x, G.get(i).position.x);
    maxP.y = max(maxP.y, G.get(i).position.y);
    maxP.z = max(maxP.z, G.get(i).position.z);
  }
  return max(max(maxP.x-o.x, maxP.y-o.y), maxP.z-o.z)+10;
}

int getLoc(float p, float o, float s, int binsize){
  if(p - o >= s) return binsize-1;  // corner case
  if(p < o) return 0; // corner case
  double xx = Math.floor((p - o) / s * binsize);
  return (int)xx;
}

// build a octree and return its root cell
// r is the radius of ball pivoting which is used to determine the root depth
Cell buildOctree(ArrayList<Vertex> G, float r){
  pt o = P(0,0,0);
  float s = findOrigin(G, o);
  
  int depth = (int)(Math.log(s/r) / Math.log(2)) + 1;
  
  if(rootDepth > depth) depth = rootDepth;
  else rootDepth = depth;
  rootSize = s;
  
  int binsize = 1 << depth;
  Cell root = new Cell(s, depth, o);
  int nv = G.size();
  for (int i = 0; i < nv; i++){
    pt p = G.get(i).position;
    int xloc = getLoc(p.x, o.x, s, binsize);
    int yloc = getLoc(p.y, o.y, s, binsize);
    int zloc = getLoc(p.z, o.z, s, binsize);
    Cell cell = root;
    int l = depth - 1;
    while (cell.depth > 0){
      int childBranchBit = 1 << l;
      int x = (xloc & childBranchBit) >> l;
      int y = (yloc & childBranchBit) >> l;
      int z = (zloc & childBranchBit) >> l;
      int childIndex = (x<<2) + (y<<1) + z;
      if (cell.child[childIndex] == null){
        float childSize = cell.size/2;
        pt childOrigin = P(cell.origin.x+x*childSize, cell.origin.y+y*childSize, cell.origin.z+z*childSize);
        Cell child = new Cell(childSize, cell.depth - 1, childOrigin);
        cell.child[childIndex] = child;
      }
      cell.hasChild = true;
      cell = cell.child[childIndex];
      l = l - 1;
    }
    cell.P.add(G.get(i));
  }
  return root;
}


int findDepth(float size, int depth, float d){
  while (depth > 0 && size / 2 > d){
    size = size / 2;
    depth--;
  }
  return depth;
}


// find cells with depth level that contains p
Cell findCell(Cell root, pt p, int level){
  pt o = root.origin;
  float s = root.size;
  int binsize = 1 << root.depth;
  int xloc = getLoc(p.x, o.x, s, binsize);
  int yloc = getLoc(p.y, o.y, s, binsize);
  int zloc = getLoc(p.z, o.z, s, binsize);
  Cell cell = root;
  int l = root.depth - 1;  // depth of potential child
  while (l >= level){  // should be ">=" 
    int childBranchBit = 1 << l;
    int x = (xloc & childBranchBit) >> l;
    int y = (yloc & childBranchBit) >> l;
    int z = (zloc & childBranchBit) >> l;
    int childIndex = (x<<2) + (y<<1) + z;
    if (cell.child[childIndex] == null) {break;}
    cell = cell.child[childIndex];
    l = l - 1;
  }
  return cell;
}

void findLeafPoints(Cell C, ArrayList<Vertex> P){
  if (!C.hasChild){  // C is a leaf
    for (int i = 0; i < C.P.size(); i++)
      P.add(C.P.get(i));
    return;
  }
  for (int i = 0; i < 8; i++)
    if (C.child[i] != null) findLeafPoints(C.child[i], P);
}

ArrayList<Vertex> getNeighbors(pt p, Cell root, float r){
  ArrayList<Vertex> Nr = new ArrayList<Vertex>();
  
  int l = findDepth(root.size, root.depth, 2 * r);
  
  ArrayList<Cell> cells = new ArrayList<Cell>();
  
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      for (int k = -1; k <= 1; k++){
        Cell C = findCell(root, P(p.x+i*r, p.y+j*r, p.z+k*r), l);
                
        if(C.depth > l) continue;
        
        boolean newCell = true;
        for (int m = 0; m < cells.size(); m++)
          if (C == cells.get(m)){
            newCell = false;
            break;
          }
        
        if (newCell) cells.add(C);
      }
  
  for (int i = 0; i < cells.size(); i++){
    if(debugOctree){
      fill(yellow); drawCell(cells.get(i), 10, false);  //draw for debugging, comment this line when running on huge datasets!
    }
    ArrayList<Vertex> tmpP = new ArrayList<Vertex>();
    findLeafPoints(cells.get(i), tmpP);
    for (int j = 0; j < tmpP.size(); j++)
      if (d(tmpP.get(j).position, p) < r)  Nr.add(tmpP.get(j));
  }
  
  return Nr;
}


// for debugging
void drawCell(Cell C, float r, boolean toLeaf){
  pt a = C.origin;
  pt b = P(1, C.origin, C.size, P(1,0,0));
  beam(a,b,r);
  b = P(1, C.origin, C.size, P(0,1,0));
  beam(a,b,r);
  b = P(1, C.origin, C.size, P(0,0,1));
  beam(a,b,r);
  
  a = P(1, C.origin, C.size, P(0,1,1));
  b = P(1, a, C.size, P(0,0,-1));
  beam(a,b,r);
  b = P(1, a, C.size, P(0,-1,0));
  beam(a,b,r);
  b = P(1, a, C.size, P(1,0,0));
  beam(a,b,r);
  
  a = P(1, C.origin, C.size, P(1,1,0));
  b = P(1, a, C.size, P(-1,0,0));
  beam(a,b,r);
  b = P(1, a, C.size, P(0,0,1));
  beam(a,b,r);
  b = P(1, a, C.size, P(0,-1,0));
  beam(a,b,r);
  
  a = P(1, C.origin, C.size, P(1,0,1));
  b = P(1, a, C.size, P(0,0,-1));
  beam(a,b,r);
  b = P(1, a, C.size, P(0,1,0));
  beam(a,b,r);
  b = P(1, a, C.size, P(-1,0,0));
  beam(a,b,r);
  
  if(toLeaf){
    for (int i = 0; i < 8; i++)
      if (C.child[i] != null) drawCell(C.child[i], 5, toLeaf);  
  }
}


int avgNumberNeighbors(Cell root, ArrayList<Vertex> P, float r){
  if(root == null || P.size() == 0) return P.size();
  
  long total = 0;
  for(Vertex v : P){
    pt p = v.position;
    ArrayList<Vertex> tmp = getNeighbors(p, root, r);
    total += tmp.size();
  }
  
  return (int)(total / P.size());
}
