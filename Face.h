#ifndef __FACE_H
#define __FACE_H

/*****************************************************************************
 * Face.h                                                                    *
 *                                                                           *
 * Contains the class header for the Face object.  Represents a face of the  *
 * polytope and contains information about adjacent edges and vertices.      *
 *****************************************************************************/

//#include "Leda.h"
#include "myCGAL.h"
#include "FEVBase.h"

//#include <set>

class Edge;
class Vert;
struct GraphNode;
class Poly;
class NNLocation;
typedef Array<int>  ArrayInt;

class Face : public FEVBase
{
public:
    Face(Poly *, int degree);
    virtual ~Face();
    
    int Degree() const;     // returns the degree of the face
    int EdgeDegree() const; // returns the num edges adj to face
    int VertDegree() const; // returns the num verts adj to face
    Edge ** Edges() const;  // returns the edges adj to face
    Vert ** Verts() const;  // returns the verts adj to face

    double * Angles() const; // trang: return the angles of vertices
    
    virtual int FaceDegree() const; // returns 1 because object is a face
    virtual Face ** Faces() const;  // returns pointer to this

    void   calcNormal();
    void   dump() const;

    Point  project( const Point  & pnt ) const;

    //Vec exactNormal() const ;           // returns normal of this face

    virtual Vec UnitNormal() const ;           // returns normal of this face
    virtual Vec NormalDir() const;              // returns normal of face1

    virtual Face * Face1();         // returns this
    virtual bool ContainsFace(const Face *) const; // true if arg is this
    virtual bool ContainsEdge(const Edge *) const; // true if contains edge
    virtual bool Contains(const FEVBase *) const;  // true if contains arg
    virtual bool Contains(const Face *) const;     // true if contains face
    virtual bool Contains(const Edge *) const;     // true if contains edge
    virtual bool Contains(const Vert *) const;     // true if contains vert

    void AddEdge(Edge *); // adds adj edge
    void AddVert(Vert *); // adds adj face

    void AddVertAngle( double ); // trang: adds angles of vertices of the face

    virtual void SetLastInPath(GraphNode *);    // set last node in path on this
    // object
    virtual void SetLastInPathAdj(GraphNode *); // set last node in path from
    // adj object

    virtual void AddNN(int);            // keeps track of NNs on face
    const ArrayInt & ReturnNNs() const; // returns list of NNs on face

    friend class Edge;
    friend class Vert;

    int   getSideSign( const Point  & pnt ) const;
    bool  isFaceOrNeigborFaceVisible( const Point  & pnt ) const;

    bool   isOn( const Point  & curr_pnt ) const;
    bool   isFaceNN( const Point  & curr_pnt ) const;
    bool   getNearestSqrDistance( const Point  & pnt,
                                  RT  & res,
                                  NNLocation  & test_nnloc, int  walk_id );
    Vert  * getNonAdjVertex( const Edge  * e ) const;
    RT  getPointValue( const Point  & pnt ) const;

    // Check if in the neighberhood of this face the polytope is
    // convex.
    void   checkIfConvex() const;  
    bool   isVertexAdj( const Vert  * v );

    // trang: ckeck if the face is adjacent with a vertex and get the vertex
    bool   isVertexAdj( const Vert  * v, int & iVertex );

    bool  isVisiblePoint( const Point  & pnt ) const;
    Edge  * getVertAdjEdge( const Vert  * v );
    Edge  * getOtherAdjEdge( const Vert  * v, const Edge  * e );
    virtual Vec exactNormal() const;
    const Vec  & getExactNorm() const {
        if (!m_normCalcd)
            ((Face *)this)->calcNormal();    
        return  m_exact_norm;
    }
    void  findFacePlanePoint( const Point  & pnt_a,
                              const Point  & pnt_b,
                              Point  & out_pnt,
                              NNLocation  & out_loc );
    void  writeNNLocInfo( Point  & pnt,
                          NNLocation  & out_loc );


private:
//public:                // trang's comment for angles of vertices
    int m_deg;         // degree of face
    int m_eDeg;        // num of edges on face
    int m_vDeg;        // num of verts on face
    Edge ** m_edges;   // edges on face
    Vert ** m_verts;   // verts on face
    Vec m_aprx_norm;   // normal of face
    Vec m_exact_norm;  // normal of face
    bool m_normCalcd;  // true if the normal has been calculated
    Face * m_faces[1]; // used for Faces()
    Array<int> m_NNs;    // list of NNs on face

    // trang: for saving the angle at a vertex of a face
    double * m_vert_angles;
    int m_angleDeg;

};

#endif


