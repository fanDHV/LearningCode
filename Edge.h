#ifndef __EDGE_H
#define __EDGE_H

/*****************************************************************************
 * Edge.h                                                                   *
 *                                                                           *
 * Contains the class header for the Edge object.  Represents an edge of the *
 * polytope and contains information about adjacent faces and vertices       *
 *****************************************************************************/

//#include "Leda.h"
#include "myCGAL.h"
#include "FEVBase.h"

class   Vert;
class   Face;
struct  GraphNode;
class   Poly;
class   NNLocation;

class Edge : public FEVBase
{
private:
    Vert * m_v1; // first vertex
    Vert * m_v2; // second vertex
    Face * m_f1; // first face
    Face * m_f2; // second face

    Face  * faces_arr[ 2 ];

public:
    Edge(Poly *, Vert *, Vert *);
    Edge(Poly *, const Edge *);
    virtual ~Edge();

    Vert * V1() const; // returns first vertex of edge
    Vert * V2() const; // returns second vertex of edge
    Face * F1() const; // returns first adjacent face of edge
    Face * F2() const; // returns second adjacent face of edge

    Face  * getOtherFace( Face  * f ) {
        if  ( F1() == f )
            return  F2();
        if  ( F2() == f )
            return  F1();

        assert( false );

        return  NULL;
    }

    RT  SqrDistance( const Point  & pnt ) const;
    RT  LineSqrDistance( const Point  & pnt ) const;

    // methods defined in FEVBase
    virtual int FaceDegree() const; // returns 2 b/c there are 2 adjacent faces
    virtual Face ** Faces() const;  // returns the adjacent faces
    virtual Vec exactNormal() const;                    // returns normal of face1
    virtual Vec UnitNormal() const;           // returns the normal of first face
    virtual Vec NormalDir() const;           // returns the normal of first face
    virtual Face * Face1();         // returns first face
    virtual bool ContainsFace(const Face *) const; // true if edge is adj to face
    virtual bool ContainsEdge(const Edge *) const; // true if arg is this
    virtual bool Contains(const FEVBase *) const;  // true if edge is adj to arg
    virtual bool Contains(const Face *) const;     // true if edge is adj to face
    virtual bool Contains(const Edge *) const;     // true if arg is this
    virtual bool Contains(const Vert *) const;     // true if edge contains vert

    // Check wheather a point is on an edge
    bool  isOnEdge( const Point  & pnt ) const;
    Point  project_dbl( const Point  & pnt ) const;

    // project a point to the line that defines the edge
    Point   project( const Point  & pnt ) const;

    Face  * getCommonFace( const Point  & pnt )  const;

    // adds face to adjacency list
    void addFace(Face *);
    
    // true if edge contains vert
    bool hasVert(const Vert *);

    // is the point in the 3d slab induced by the edge?
    bool isInSlab( const Point & pnt ) const;

    // Create the nnlocation record for a point
    void  writeNNLocRecord( const Point  & in_pnt,
                            NNLocation  & out );
    const Point  & getOtherPoint( const Point  & pnt ); 

    // Check if in the neighberhood of this edge the polytope is
    // convex.
    void   checkIfConvex() const;  

    bool   getSqrDistance( const Point  & in_pnt, RT  & dist,
                          NNLocation  & test_nnloc, int  walk_id );

    virtual void SetLastInPath(GraphNode *);    // used in shortcutting to keep
    // track of last node on edge
    virtual void SetLastInPathAdj(GraphNode *); // called by Face or Vert to
    // keep track of last node

    virtual void AddNN(int); // keeps track of NNs on this edge

    friend class Vert;
    friend class Face;

};

#endif






