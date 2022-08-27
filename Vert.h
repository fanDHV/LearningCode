#ifndef __VERT_H
#define __VERT_H

/*****************************************************************************
 * Vert.h                                                                    *
 *                                                                           *
 * Contains the class header for the Vert object.  Represents a vertex of    *
 * the polytope and contains information about adjacent edges and faces.     *
 *****************************************************************************/

//#include "Leda.h"
#include "myCGAL.h"
#include "FEVBase.h"

class Face;
class Edge;
struct GraphNode;
class Poly;

class Vert : public FEVBase
{
public:
    Vert(Poly*, const Point & p);
    virtual ~Vert();
    
    // location of the vertex
    const Point & Location() const;

    // number of edges that contain this vertex
    int EdgeDegree() const; 
    
    // edges that contain this vertex
    Edge ** Edges() const; 
    
    virtual int FaceDegree() const; // number of adj faces
    virtual Face ** Faces() const;  // adj faces
    virtual Vec exactNormal() const;
    virtual Vec UnitNormal() const;           // normal of first face
    virtual Vec NormalDir() const;           // normal of first face

    virtual Face * Face1();         // first face
    virtual bool ContainsFace(const Face *) const; // true if adj to face
    virtual bool ContainsEdge(const Edge *) const; // true if adj to edge
    virtual bool Contains(const FEVBase *) const;  // true if adj to arg
    virtual bool Contains(const Face *) const;     // true if contains face
    virtual bool Contains(const Edge *) const;     // true if contains edge
    virtual bool Contains(const Vert *) const;     // true if contains vert

    // sort the faces, edges, vertices info around the vertex. It is
    // sorted either CW or CCW - I dont care as long as it is sorted.
    void  sortAdjInfo();

    const Point  & getEdgeEndPoint( int  ind );

    void AddFace(Face *);           // adds an adj face
    void AddEdge(Edge *);           // adds an adj edge
    Edge * Connected(const Vert *); // if this is connected to arg, returns the
    // edge
    
    virtual void SetLastInPath(GraphNode *);    // sets last node of path on this
    // vert
    virtual void SetLastInPathAdj(GraphNode *); // sets last node of path adj
    // from this vertex
    
    static int MaxEdgeDegree();                 // returns the maximum edge
    // degree of any vertex
    
    virtual void AddNN(int);                    // keeps track of NNs on the vert

    friend class Face;
    friend class Edge;
    
private:
    void Grow();     // increases the sizes of the edge and face arrays
    
    Point m_loc;     // location of the vertex
    int m_fDeg;      // number of adj faces
    int m_eDeg;      // number of adj edges
    int m_maxDeg;    // size of arrays
    Edge ** m_edges; // adj edges
    Face ** m_faces; // adj faces
    
    static int sm_maxEDeg; // maximum edge degree

    // Those functions sorts the edges and faces arrays around the
    // vertex by walking on the adjacent faces.
    void  sortFacesInfo();
    void  sortEdgesInfo();
};

#endif
