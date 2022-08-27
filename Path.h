/****************************************************************************
 * Path.h                                                                   *
 *                                                                          *
 * Contains header for the Path class as well as other structs used in      *
 * conjunction with it.  The Path class is the main interface between the   *
 * Renderer and the path calculation engine.                                *
 ****************************************************************************/

#ifndef __PATH_H
#define __PATH_H

#include <string>
#include <list>

//#include "Leda.h"
#include "myCGAL.h"
#include "NNLocation.h"

#include <list>
#include <vector>
#include <iostream>

using namespace std;

class  Face;
class  Poly;
class  Edge;
class  SPath;
// class  sliceMngt;
class  SPathNode;
struct  history_t;

// used by GraphNode to store adjacency and distance information about the
// graph
struct AdjDist
{
    int m_adj;
    double m_dist;

    AdjDist(int a, double d) : m_adj(a), m_dist(d)
    {
        // do nothing
    };
};


// provides an order on the adjancency and distance information.  done
// this way to automatically remove duplicates from the set when
// adding new adjacent nodes
struct AdjDistLt
{
    bool operator()( AdjDist * a, AdjDist * b ) const {
        return a->m_adj < b->m_adj;
    }
};

// useful typedef
//typedef set<AdjDist *, AdjDistLt>   NodeSet;

// node in the path graph
struct GraphNode
{
    Point  on_sphere;
    Point  m_p; // location on polytope
    NNLocation  m_nnloc; // information about where the NN lies
    Array<AdjDist>   neighbors;
     // NodeSet  m_adjDist; // set of adjacent nodes
    int  m_deg; // degree of the node
    GraphNode  * m_prev; // for keeping track of the shortest path
    GraphNode  * m_next; // for keeping track of the shortest path

    int  degree()  const {
        return  (int)neighbors.size();
    }
    GraphNode() : neighbors(),  m_deg(0) {
        m_prev = 0;
        m_next = 0;
        neighbors.init( 10, 8943 );
    }

    GraphNode(const Point & P, NNLocation & nnl,
              const Point & P_sphere )
        : neighbors(), m_deg(0) {
        // reserve Deg nodes in the set for adjacent nodes
        //for (int i=0; i<m_deg; i++)
        //   neighbors.push_back( AdjDist(-i-1, -1.0) );
        on_sphere = P_sphere;
        m_p = P;
        m_nnloc = nnl;

        m_prev = 0;
        m_next = 0;
        neighbors.init( 10, 8943 );
    }

    ~GraphNode() {
    }

    // increases the degree of the node by one and adds an adjacent
    void IncDeg(int a) {
        neighbors.pushx( AdjDist( a, -1.0 ) );
    }

    void  addAdjacency( int  self_index, int  index, bool  fSkipCheck ) {
        // trang 2013/12: avoid unused variables
        self_index = self_index + 0; fSkipCheck = !fSkipCheck;

        if  ( index == GG_INVALID )
            return;
        //printf( "AA %d %d\n", self_index, index );
        assert( index >= 0 );
        //printf( "AA index: %d\n", index );
        for  ( int  ind = 0; ind < degree(); ind++ )
            if  ( neighbors[ ind ].m_adj == index )
                return;
        neighbors.pushx( AdjDist( index, -1.0 ) );
    }

};


// used in communicating information about the source and destination points
// and their faces between the Renderer and the engine
struct SourceDest
{
  SourceDest() {};
  SourceDest(double * src, double * dest, int srcFaceIdx, int destFaceIdx,
         Face * const * faces, Poly  * poly, bool isprojected = true );

  Point m_srcPoint;
  Point m_destPoint;
  Face * m_srcFace;
  Face * m_destFace;
  // trang: for processing of descending path
  deque< Seg > m_srcSlice;
  deque< int > m_srcSliceFaces;

// copy constructor Le 6-8-2022 for Invalid access to memory location error 
  SourceDest(SourceDest &srcdest) {
      m_srcPoint = srcdest.m_srcPoint;
      m_destPoint = srcdest.m_destPoint;
      m_srcFace = srcdest.m_srcFace;
      m_destFace = srcdest.m_destFace;
      m_srcSlice = srcdest.m_srcSlice;
      m_srcSliceFaces = srcdest.m_srcSliceFaces;
  };
//

};

// used to pass information to the Renderer that is used to animate parts
// of the algorithm
struct AnimateVecInfo
{
    Point m_p1;
    Point m_p2;
    enum VecType {LONGITUDE, REGION, NN, OTHER} m_vecType;

    AnimateVecInfo(Point & p1, Point & p2, VecType vt) : m_vecType(vt) {
        m_p1 = p1;
        m_p2 = p2;
    }
};

// header for the main interface to the outside world
class Path
{
private:
    Poly  * m_poly; // the polytope
    double  m_epsilon; // parameter epsilon
    //int  m_numNodes; // number of nodes in the graph
    //int  m_sizeGraph; // size of the array holding the graph
    int  m_pickNum; // no longer used
    //int  m_iiBound; // number of latitude lines
    //int  m_jjBound; // number of longitude lines
    //GraphNode  ** m_graph; // path graph
    GridGraph  grid;
    // trang (2015/08/29): change to public one for calling from rerender
    // because of computing the global SDP
    //SourceDest  * m_srcdest; // source and destination points
    GraphNode  * m_extraGraphNodes;  // list of extra graph nodes
    // used in 2-plane
    // wedge and path projection
    GraphNode  * m_currExtraNode; // pointer to current node in above list
    double  m_pathLength; // length of the path
    Plane  m_Hs; // source plane for 3-wedge
    Plane  m_Ht; // dest plane for 3-wedge

public: // trang (2015/08/29):changed to be public for calling from rerender
    SourceDest  * m_srcdest; // source and destination points

public:

    Path(double epsilon);
    ~Path();

    // reads a polytope in OFF format from a file and returns the number of faces
    int ReadPoly(const std::string &, const Face * const * &);

    // adds a node to the path graph in the specified location
    void AddNodeToGraph(GraphNode * gn, int idx);

    // calculates the approximate shortest path between the src and dest points - not used
    void CalcPath( double * src, double * dest,
           int srcFaceIdx, int destFaceIdx,
           bool bUseHershSuri,
           SPath  * & path );
    /* void algHershSuri( double  * src, double  * dest, */
    /*             int  srcFaceIdx, int  destFaceIdx, */
    /*             bool  bUseHershSuri, SPath  * & s_path ); */
    void algHershSuri( bool  bUseHershSuri, SPath  * & s_path );
    // calculates the approximate shortest path by multiple shooting method
    /* void spCalcPath( double *src, double *dest,  */
    /*           int srcFaceIdx, int destFaceIdx, */
    /*           bool bUseHershSuri, */
    /*           SPath  *& path ); */
    void spCalcPath( SourceDest *srcdest, bool bUseHershSuri, SPath  *& path );    // for clasical distance
    void spCalcSGPPath( SourceDest *srcdest, bool bUseHershSuri, SPath  *& path );   // for SGP distance
    // Le: to calculate intersect SGP of src & dest and slice polygon
    Point spCalcSGPPath( SourceDest *srcdest, deque <Seg> slices_s, Point atVertex_s);

    void multispCalcPath( double *src, double *dest,
                          int srcFaceIdx, int destFaceIdx,
                          bool bUseHershSuri,
                          SPath *& apath, std::vector<Point> *&multisp_path,
                          std::vector<Point> *&multisp_path_first,
                          std::deque< std::deque< Seg > > *&slices );

    // for running LiuWong algorithm 2022-8-28
    void liuWongCalcPath ( double *src, double *dest,
                             int srcFaceIdx, int destFaceIdx);

    // pre-processes the polytope to generate the path graph
    void ConstructPathGraph( std::list<AnimateVecInfo *> &);

    // returns a specific node in the graph
    GraphNode * ReturnNode(int i) const;
    const GridGraph  & getGridGraph() const {
        return  grid;
    }

    // returns the length of the path
    double PathLength() const;

    // returns the path graph
    const GraphNode * const * ReturnGraph() const;

    // returns the number of nodes in the graph
    int NumNodes() const;


    double   vertices_distance( int  a, int  b ) const {
        return  grid.vertices_distance( a, b );
    }

    void   graph_set_distance( int  to, int  j, double  dist ) {
        grid.set_distance( to, j, dist );
    }

    SPath  * convertToSPath( GraphNode  * path );

    /* void  algHershSuriExt( double  * src, double  * dest, */
    /*                        int  srcFaceIdx, int  destFaceIdx, */
    /*                        bool  bUseHershSuri, SPath  * & p_path ); */
    void  algHershSuriExt( bool  bUseHershSuri, SPath  * & p_path );
    std::vector< Point >* slices();

private:
    // iterates through all of the nodes in the graph and calculates the distance
    // from that node to each of its adjacent nodes
    void CalcDists();                                                             // clasical distance
    void CalcSGPDists ();                                                  // clasical distance

    // prints the graph
    void PrintGraph();

    // executes Djikstra's algorithm
    GraphNode * FindShortestPathInGraph(const int s, const int t);             // clasical distance
    GraphNode * FindShortestGentlePathInGraph(const int s, const int t);   // SGP distance

    // cleans up the path so that there are at most two points on any face
    GraphNode * ShortcutPath(GraphNode * path);
    SPath  * ShortcutPathSP( GraphNode * path );
    GraphNode   * ProjectPathSP( GraphNode  * path,
                                       Point  & ridge_point,
                                       SPath  & sp );
    Edge  * findExitEdge( NNLocation  & nnloc,
                          Plane  & plane, Plane  & hPositive,
                          history_t  & hist, Point  & exit_point,
                          GraphNode  * dest );

    // projects the path onto the polytope
    GraphNode * ProjectPath(GraphNode * path, bool b3Wedge=false);

    // adds another GraphNode to the array, growing it if necessary
    void AddExtraGraphNode(GraphNode * curr, const Point & p);
    void AddExtraGraphNodeExt( GraphNode  * curr, const Point  & p,
                               Edge  * e );

    // returns the middle plane for the 3-plane wedge
    Plane * MiddlePlane(Edge * e, const Vec & i);

    // computes the three wedge using the given plane
    RT ThreeWedge(Point * v, Plane * P);


public: // trang 2013/12/23: for computing shortest descending path

    // options for solving in step by step ...
    // solve to finish ...
    bool solve_finish;

    // solve in step by step ...
    bool show_cut_slices;
    bool show_initial_path;
    bool finish;

    // trang: avoiding unused variable
    bool used;

    // Le: Le modifies SDP->SGP to find shortest gentle paths
    void calcSGP(SourceDest *srcdest, bool bUseHershSuri, SPath *&s_path);      // SGP distance
    void calcSP(SourceDest *srcdest, bool bUseHershSuri, SPath *&s_path);  // clasical distance

    bool isDescending(SPath s_path );
    bool isDescending(SPath *s_path, SPathNode *first_inc, SPathNode *last_inc );


};

#endif