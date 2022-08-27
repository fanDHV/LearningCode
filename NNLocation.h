#ifndef __NNLOCATION_H
#define __NNLOCATION_H

/*****************************************************************************
 * NNLocation.h                                                              *
 *                                                                           *
 * Contains the structure that holds information about each NN.              *
 *****************************************************************************/

class FEVBase;
class Vert;
class Edge;

struct NNLocation
{
    NNLocation() : m_loc(0) {}
  
    FEVBase  * m_loc;                   // Face, Edge, or Vert the NN lies on
    enum {FACE, EDGE, VERT} m_locType; // On what the NN lies

    FEVBase  * loc() {
        return  m_loc;
    }

    bool  isEmpty() const {
        return  m_loc == NULL;
    }

    void  dump() {
        printf( "[type: %d, %p ] ", m_locType, m_loc );
    }

    bool  isFace() const {
        return  m_locType == FACE;
    }

    bool  isVertex() const {
        return  m_locType == VERT;
    }

    bool  isEdge() const {
        return  m_locType == EDGE;
    }

    Vert  & vertex() {
        assert( isVertex() );
        return  *(Vert *)m_loc;
    }
    Edge  * edgePtr() {
        assert( isEdge() );
        return  (Edge *)m_loc;
    }

};

#endif

