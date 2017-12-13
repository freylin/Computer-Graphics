//*****************************************************************************
// TITLE:         Types 
// DESCRIPTION:   Classes of vertex, edge, triangle
// AUTHORS:       Yaohong Wu, Jin Lin
// EDITS:         11/28/2017
//*****************************************************************************

class Triangle{
  int a, b, c;
  Triangle(int _a, int _b, int _c) { a = _a; b = _b; c = _c;}
}

class Edge{
  int a, b;
  Edge(int _a, int _b) { a = _a; b = _b;}
}

class taggedEdge{
  Edge e;
  int c; // the id of opposite vertex w.r.t. edge e, it only makes senses when e is a front or border edge
  pt F;  // center of ball touching c and e, it only makes senses when e is a front or border edge
  int tag;  // 0:front, 1:border, 2:inner
  taggedEdge(int _a, int _b, int _c, pt _F, int _tag){
    e = new Edge(_a,_b);
    c = _c; F = _F; tag = _tag;
  }
}


class Vertex{
  pt position;
  vec normal;
  int type; // 0:orphan (unvisited, isolated), 1:border 2:inner
  int id;
  
  HashMap<Vertex, taggedEdge> outEdges;  // front or border edges going from this vertex
  HashMap<Vertex, taggedEdge> inEdges;  // front or border edges going into this vertex
  
  Vertex(){}
  Vertex(pt _position, vec _normal, int _type, int _id){
    position = _position;
    normal = _normal;
    type = _type;
    id = _id;
    outEdges = new HashMap<Vertex, taggedEdge>();
    inEdges = new HashMap<Vertex, taggedEdge>();
  }
}