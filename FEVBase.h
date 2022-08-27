#ifndef __FEVBASE_H
#define __FEVBASE_H

/*****************************************************************************
 * FEVBase.h                                                                 *
 *                                                                           *
 * Base class for Face, Edge, and Vert.  Almost purely abstract.  Defines    *
 * interface used mostly in shortcutting and projection.                     *
 *****************************************************************************/

struct GraphNode;
class Face;
class Edge;
class Vert;
class Poly;

class FEVBase
{
public:
    FEVBase(Poly * p) : m_poly(p)                    // set m_poly
    {
        face_walk_id = 0;
    }
    virtual ~FEVBase()                               // don't need to do anything
    {}

    virtual void SetLastInPath(GraphNode *)=0;       // set the last node in the
    // path on this object to
    // be arg
    virtual Vec UnitNormal() const=0;               // returns normal of face1
    virtual Vec NormalDir() const=0;              // returns normal of face1
    virtual Face * Face1()=0;                        // returns first face
    virtual Face ** Faces() const=0;                 // returns adj faces
    virtual int FaceDegree() const=0;                // returns num adj faces
    virtual bool ContainsFace(const Face *) const=0; // true if contains face
    virtual bool ContainsEdge(const Edge *) const=0; // true if contains edge
    virtual bool Contains(const FEVBase *) const=0;  // true if contains arg
    virtual bool Contains(const Face *) const=0;     // true if contains face
    virtual bool Contains(const Edge *) const=0;     // true if contains edge
    virtual bool Contains(const Vert *) const=0;     // true if contains vert
    
    bool  isVisibleEdge( const Edge  * e ) const;
    bool  isVisiblePoint( const Point  & pnt ) const;

    GraphNode * ReturnLastInPath() const             // returns last node in the
    { return m_lastInPath; }                       // path on this object

    virtual void AddNN(int)=0;                       // keeps track of NNs on
                                                   // this object
    int  getLastFaceWalkID() const {
        return  face_walk_id;
    }
    void  setLastFaceWalkID( int  id ) {
        face_walk_id = id;
    }


protected:
    virtual void SetLastInPathAdj(GraphNode *)=0;    // set the last node in the
                                                     // from adj object
    GraphNode * m_lastInPath;                        // last node in the path
                                                     // on this object

    int  face_walk_id;  
    // poly that contains this feature
    Poly * m_poly;                                   
};

#endif



