/**************************************************************************
 * Vert.cc
 * 
 * Implementation of the class that represents a vertex on the polytope
\**************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>

#include  "sarray.h"
#include  "Vert.h"
#include  "Edge.h"
#include  "Face.h"

int
  Vert::sm_maxEDeg = 0;

Vert::Vert (Poly * P, const Point & p):
FEVBase (P)
{
  m_loc = p;
  m_maxDeg = 3;
  m_fDeg = m_eDeg = 0;
  m_edges = new Edge *[m_maxDeg];
  m_faces = new Face *[m_maxDeg];
  m_lastInPath = 0;
}

Vert::~Vert ()
{
  delete m_edges;
  delete m_faces;
}

const Point &
     Vert::Location () const
     {
       return m_loc;
     }

     int Vert::FaceDegree () const
     {
       return m_fDeg;
     }

     int Vert::EdgeDegree () const
     {
       return m_eDeg;
     }

     Edge **Vert::Edges () const
     {
       return m_edges;
     }

     Face **Vert::Faces () const
     {
       return m_faces;
     }

Vec
     Vert::UnitNormal () const
     {
       return m_faces[0]->UnitNormal ();
     }

Vec
     Vert::NormalDir () const
     {
       return m_faces[0]->NormalDir ();
     }

Vec
     Vert::exactNormal () const
     {
       return m_faces[0]->exactNormal ();
     }

     Face *Vert::Face1 ()
{
  return m_faces[0];
}

     bool Vert::ContainsFace (const Face * f) const
     {
       for (int i = 0; i < m_fDeg; i++)
         {
           if (m_faces[i] == f)
             return true;
         }
       return
         false;
     }

     bool Vert::ContainsEdge (const Edge * e) const
     {
       for (int i = 0; i < m_eDeg; i++)
         if (m_edges[i] == e)
           return true;

       return
         false;
     }

bool
     Vert::Contains (const FEVBase * fev) const
     {
       return fev->Contains (this);
     }

bool
     Vert::Contains (const Face * f) const
     {
       return ContainsFace (f);
     }

bool
     Vert::Contains (const Edge * e) const
     {
       return ContainsEdge (e);
     }

bool
     Vert::Contains (const Vert * v) const
     {
       return (this == v);
     }

     void Vert::AddFace (Face * f)
{
  if (m_fDeg == m_maxDeg)
    Grow ();

  m_faces[m_fDeg++] = f;
}

void
Vert::AddEdge (Edge * e)
{
  if (m_eDeg == m_maxDeg)
    Grow ();

  m_edges[m_eDeg++] = e;

  if (m_eDeg > sm_maxEDeg)
    sm_maxEDeg = m_eDeg;
}

Edge *
Vert::Connected (const Vert * v)
{
  for (int i = 0; i < m_eDeg; i++)
    {
      if (m_edges[i]->hasVert (v))
        return m_edges[i];
    }

  return NULL;
}

// Since we don't know ahead of time how many faces and edges are adjacent to
// this vertex, grow the vectors if necessary
void
Vert::Grow ()
{
  int i;

  m_maxDeg *= 2;

  Face **newFaces = new Face *[m_maxDeg];
  for (i = 0; i < m_fDeg; i++)
    newFaces[i] = m_faces[i];

  delete [] m_faces;
  m_faces = newFaces;

  Edge **newEdges = new Edge *[m_maxDeg];
  for (i = 0; i < m_eDeg; i++)
    newEdges[i] = m_edges[i];

  delete [] m_edges;
  m_edges = newEdges;
}

// If the information about the last node in the path on this vertex is not
// already set, set it.  Or, if gn == 0, clear all of the settings
void
Vert::SetLastInPathAdj (GraphNode * gn)
{
  // set information for use in shortcutting
  if (!gn || !m_lastInPath)
    m_lastInPath = gn;

  // set the information for all adjacent edges and verts too
  if (m_fDeg != m_eDeg)
    {
      std::
        cerr << "ERROR: Edge degree does not match face degree!" << std::endl;
      exit (1);
    }
  for (int i = 0; i < m_fDeg; i++)
    {
      if (!gn || !m_edges[i]->m_lastInPath)
        m_edges[i]->m_lastInPath = gn;
      if (!gn || !m_faces[i]->m_lastInPath)
        m_faces[i]->m_lastInPath = gn;
    }
}

// Tell all adjacent faces that this node is that last one in the path
// on them (a point on the vertex is on all adjacent faces)
void
Vert::SetLastInPath (GraphNode * gn)
{
  for (int i = 0; i < m_fDeg; i++)
    m_faces[i]->SetLastInPathAdj (gn);
}

int
Vert::MaxEdgeDegree ()
{
  return sm_maxEDeg;
}

// Tell all adjacent faces that there is a nearest neighbor on this vertex
void
Vert::AddNN (int nn)
{
  for (int i = 0; i < m_fDeg; i++)
    m_faces[i]->AddNN (nn);
}

void
Vert::sortFacesInfo ()
{
  Face **arr_faces;
  Face *first_face, *curr;
  Edge *e;
  int count, size;

  count = 0;
  arr_faces = new Face *[m_maxDeg];
  assert (arr_faces != NULL);

  size = sizeof (Face *) * m_maxDeg;
  memset (arr_faces, 0, size);

  curr = first_face = m_faces[0];
  e = first_face->getVertAdjEdge (this);
  do
    {
      arr_faces[count] = curr;
      count++;

      assert (count <= m_fDeg);

      e = curr->getOtherAdjEdge (this, e);
      curr = e->getOtherFace (curr);
    }
  while (curr != first_face);
  assert (count == m_fDeg);
  assert (m_eDeg == m_fDeg);

  memcpy (m_faces, arr_faces, size);
  delete[]arr_faces;
}


// used to regenerate the vertices adjacency info so that it is sorted
// around the vertex. This assumes that the faces adj list is already
// sorted
void
Vert::sortEdgesInfo ()
{
  int count = 0;
  Face *curr, *first_face;
  Edge *e;

  // rewrite the vertices info
  curr = first_face = m_faces[0];
  e = first_face->getVertAdjEdge (this);
  do
    {
      m_edges[count] = e;
      count++;

      assert (count <= m_eDeg);

      e = curr->getOtherAdjEdge (this, e);
      curr = e->getOtherFace (curr);
    }
  while (curr != first_face);
}


void
Vert::sortAdjInfo ()
{
  sortFacesInfo ();
  sortEdgesInfo ();
}

const Point &
Vert::getEdgeEndPoint (int ind)
{
  Edge *e = m_edges[ind];
  if (e->V1 () == this)
    return e->V2 ()->Location ();
  if (e->V2 () == this)
    return e->V1 ()->Location ();
  assert (false);

  return e->V1 ()->Location ();
}


/* Vert.C - End of File ------------------------------------------*/
