//*****************************************************************************
// TITLE:         Triangle Mesh
// DESCRIPTION:   Class of Triangle Mesh
// AUTHORS:       Yaohong Wu, Jin Lin
// EDITS:         11/28/2017
//*****************************************************************************

class TriangleMesh{
  ArrayList<pt> G;  // positions of vertices
  int nv;  // number of vertices
  ArrayList<Triangle> T;  // triangles, note that only vertex indices are stored
  int nt;  // number of triangles
  
  TriangleMesh(){
    nv = nt = 0;
  }
  
  TriangleMesh(int _nv){ 
    nv = _nv;
  }
  
  TriangleMesh(ArrayList<pt> _G, ArrayList<Triangle> _T){
    G = _G;
    T = _T;
    nv = G.size();
    nt = T.size();
  }
  
  void setTo(ArrayList<Vertex> _P, ArrayList<Triangle> _T){
    G = new ArrayList<pt>();
    for(Vertex v: _P) G.add(v.position);
    T = _T;
    nv = G.size();
    nt = T.size();
  }
  
  void saveTriMesh(String fn){
    String[] ss = new String [nv+nt+2];
    int s=0;
    ss[s++] = str(nv);
    for(int i = 0; i < nv; ++i){ ss[s++] = str(G.get(i).x) + "," + str(G.get(i).y) + "," + str(G.get(i).z); }
    ss[s++] = str(nt);  
    for(int i = 0; i < T.size(); ++i){ ss[s++] = str(T.get(i).a) + "," + str(T.get(i).b) + "," + str(T.get(i).c);}
    saveStrings(fn,ss);
  }
  
  void loadTriMesh(String fn){
    String[] ss = loadStrings(fn);
    int s = 0;
    nv = int(ss[s++]);
    
    G = new ArrayList<pt>();
    
    for(int i = 0; i < nv; ++i) { 
      String cur = ss[s++];
      String[] xyz = split(cur, ",");
      G.add(P(float(xyz[0]), float(xyz[1]), float(xyz[2])));
    }

    nt = int(ss[s++]);
    
    T = new ArrayList<Triangle>();
    
    for(int i = 0; i < nt; ++i){
      String cur = ss[s++];
      String[] abc = split(cur, ",");
      T.add(new Triangle(int(abc[0]), int(abc[1]), int(abc[2])));
    }
  }
  
  void show(){
    for(int i = 0; i < nt; ++i){
      beginShape();
        vertex(G.get(T.get(i).a));
        vertex(G.get(T.get(i).b));
        vertex(G.get(T.get(i).c));
      endShape();
    }
  }
}