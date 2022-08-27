#ifndef __POLY_H
#define __POLY_H

/*****************************************************************************
 * Poly.h                                                                    *
 *                                                                           *
 * Contains the header for the Poly class.  Represents the polytope and does *
 * most of the heavy duty work associated with find the nearest neighbors    *
 * and constructing the graph.                                               *
 *****************************************************************************/

#include <list>
#include <vector>
#include <iostream>

#include "myCGAL.h"
#include "NNLocation.h"
#include "bbox.h"

class  Edge;
class  Face;
class  Vert;
class  Path;
class  SourceDest;
class  GridGraph;

void  vec_print( const Vec  & v );
void  point_print( const Point  & pnt );
void  point_print_input_frmt( const Point  & pnt );

struct AnimateVecInfo;
struct  walk_nn_t;

using namespace std;

class Poly
{
public:
    Poly();
    virtual ~Poly();

    // reads the polytope
    void Read(istream &);
    void ReadInner(istream &);

    // checks the direction of the normals of all of the faces
    // void Check(ostream &);
    void Check();

    // returns the bounding box for the polytope
    const RT * Bounds() const;

    // remains from Neill's code ... not implemented in this version because
    // only used in 3-plane wedge with Hershberger-Suri algorithm
    std::list<Edge *> HorizonEdges(const Vec &) const;

    // returns faces of the polytope
    Face * const * Faces() const;

    // returns vertices of the polytope
    Vert * const * Verts() const;

    // returns number of vertices
    int NumVerts() const;

    // returns number of faces
    int NumFaces() const;

    // returns number of faces
    int NumEdges() const;

    // returns the center of the polytope
    Point Center() const;

    // also from Neill's code ... not implemented
    Point RandomVertex() const;
    bool IsHorizonEdge(const Edge *, const Vec &) const;

    // finds the nearest neighbors and constructs the path graph
    virtual void ConstructPathGraph(Path *, const double,
                                    GridGraph  & grid,
                                    std::list<AnimateVecInfo *> &);

    // adds the source and destination points to the graph
    void  AddSourceDest( Path * pthP, SourceDest * srcdest,
			 const double epsilon,
                         GridGraph  & grid, const int pickNum );
    /*****Le add ****/
    void AddSourceDest( SourceDest * srcdest, const double epsilon,
                         GridGraph  & grid, const int pickNum );
    /*****Le add - end****/
    int   addPoint( GridGraph  & grid, const Point  & pnt,
                    const Point  & pnt_on_sphere,
                    Face  * pFace, double  epsilon );
    void  doFaceWalk( Face  * f, walk_nn_t  & info );
    void  doPolyDump( const Point  * p_pnt );
    Point  computeRealCenter()  const;

    // checks all faces, edges, and vertices for the NN of curpt
    void NNCheckAll(const Point & curpt,	NNLocation & nnloc, Point &);
    int  getFaceIndex( const Face  * face ) const;

private:
    // updates the bouding box while reading in the polytope
    void UpdateBounds(const Point &);


    // true if the point is closest to this face
    //bool CheckNNFace(const Face *, const Point &, array<Vec> &);
    bool CheckNNFace(const Face *, const Point &, vector<Vec> & );

    // true if the point is closest to this edge
    bool CheckNNEdge(const Edge *, const Point &, Vec &, Vec &, Vec &, Vec &);
    // true if the point is closest to this vertex
    //bool CheckNNVertex(const Vert *, const Point &, array<Vec> &);
    bool CheckNNVertex(const Vert *, const Point &, vector<Vec> &);

    void  WalkNNFace( walk_nn_t   & info,
                      list<AnimateVecInfo *> & debugPoints );
    void  WalkNNEdge( walk_nn_t   & info,
                      list<AnimateVecInfo *> & debugPoints );
    void  WalkNNVertex( walk_nn_t   & info,
                      list<AnimateVecInfo *> & debugPoints );
    void  convexify( Vert  ** arr, int  num );

    // true if an edge is a horizon edge
    bool  hEdge(const Edge * e, const Vec & u) const;
    void  handlePoint( walk_nn_t  & info,
                            double  phi, double  theta,
                            std::list<AnimateVecInfo *> & debugPoints);
    void  doPolyTest( const Point  * p_pnt );
    void  sortAdjInfo();

protected:
    int m_numVerts; // number of vertices
    int m_numFaces; // number of faces
    int m_numEdges; // number of edges
    Vert ** m_verts; // vertices of the polytope
    Face ** m_faces; // faces of the polytope
    std::list<Edge *> m_edges; // edges of the polytope

    BBox  m_bBox;
    //rational m_boundingBox[6]; // bounding box

    // radius and center of surrounding sphere for ApproxPolytope and
    // AddSourceDest
    RT m_r;
    Point m_C;

    // trang: create the array to hold the angles of Verts of each Face
    // NOTE that a face here is a triangle, i.e., it includes three vertices ...???
public:
    std::vector < double * > m_total_vert_angles;
};

#endif
