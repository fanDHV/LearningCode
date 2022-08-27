/*****************************************************************************
 * Face.cc
 *
 * Implementation of the class that represents a face on the polytope
 *****************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>

#include  "sarray.h"
#include "Face.h"
#include "Vert.h"
#include "Edge.h"
#include "Functions.h"
#include "Poly.h"

Face::Face (Poly * p, int degree):FEVBase (p)
{
  m_deg = degree;
  m_eDeg = m_vDeg = m_angleDeg = 0;

  m_edges = new Edge *[m_deg];
  m_verts = new Vert *[m_deg];

  // trang: holding the angles of vertices of the face
  m_vert_angles = new double [m_deg];

  m_normCalcd = false;
  m_lastInPath = 0;
  m_faces[0] = this;
  m_NNs.init (10, 54898);

}

Face::~Face ()
{
  delete m_edges;
  delete m_verts;
}

int
Face::Degree () const
{
  return m_deg;
}

int Face::EdgeDegree () const
{
  return m_eDeg;
}

int Face::VertDegree () const
{
  return m_vDeg;
}

int Face::FaceDegree () const
{
  return 1;
}

Face **Face::Faces () const
{
  return const_cast < Face ** >(m_faces);
}

Edge **Face::Edges () const
{
  return m_edges;
}

Vert **Face::Verts () const
{
  return m_verts;
}

// trang: return the angles of vertices of the face
double * Face::Angles() const
{
    return m_vert_angles;
}

void Face::calcNormal ()
{
  assert (m_vDeg == 3);

  /* Hoai - 6/8/2012
     m_aprx_norm = Unit(CrossProductDir(m_verts[1]->Location() -
     m_verts[0]->Location(),
     m_verts[2]->Location() -
     m_verts[0]->Location(), 7786 ));
     m_exact_norm = CrossProductDir(m_verts[1]->Location() -
     m_verts[0]->Location(),
     m_verts[2]->Location() -
     m_verts[0]->Location(), 89871);
   */
  m_aprx_norm = Unit (CrossProduct (m_verts[1]->Location () -
                                    m_verts[0]->Location (),
                                    m_verts[2]->Location () -
                                    m_verts[0]->Location ()));
  m_exact_norm = CrossProduct (m_verts[1]->Location () -
                               m_verts[0]->Location (),
                               m_verts[2]->Location () -
                               m_verts[0]->Location ());

  /*
     cout << "XX: " << DotProduct( m_exact_norm, m_verts[1]->Location() -
     m_verts[0]->Location() ) << "\n";
     cout << "XX: "
     << DotProduct( m_exact_norm, m_verts[1]->Location().to_vector() )
     - DotProduct( m_exact_norm, m_verts[0]->Location().to_vector() )
     << "\n";
     fflush( stdout );
   */
  // make sure that the normal is pointing outwards
  // this relies on the origin being inside the poly
  // also assumes that faces have degree 3
  if (sign (DotProduct (m_verts[0]->Location () - m_poly->Center (),
                        m_exact_norm)) < 0)
    {
      // if it's pointing the wrong way, flip it and change the order
      // of v1,v2 and e1,e2
      m_exact_norm = -m_exact_norm;
      m_aprx_norm = -m_aprx_norm;

      Vert *v = m_verts[1];
      m_verts[1] = m_verts[2];
      m_verts[2] = v;

      Edge *e = m_edges[0];
      m_edges[0] = m_edges[2];
      m_edges[2] = e;
      assert( false );
    }

  assert (DotProduct (m_exact_norm, m_verts[1]->Location () -
                      m_verts[0]->Location ()) == RT (0));

  assert (DotProduct (m_exact_norm, m_verts[2]->Location () -
                      m_verts[1]->Location ()) == RT (0));

  assert (DotProduct (m_exact_norm, m_verts[2]->Location () -
                      m_verts[0]->Location ()) == RT (0));

  m_normCalcd = true;
  dump();
  cout << "\tNorm vector (" << m_exact_norm << ")" << endl;
}



Face *
Face::Face1 ()
{
  return this;
}

bool Face::ContainsFace (const Face * f) const
{
  return (this == f);
}

bool
Face::ContainsEdge (const Edge * e) const
{
  for (int i = 0; i < m_eDeg; i++)
    if (m_edges[i] == e)
      return true;
  return false;
}

bool
Face::Contains (const FEVBase * fev) const
{
  return fev->Contains (this);
}

bool
Face::Contains (const Face * f) const
{
  return ContainsFace (f);
}

bool
Face::Contains (const Edge * e) const
{
  return ContainsEdge (e);
}

bool
Face::Contains (const Vert * v) const
{
  for (int i = 0; i < m_vDeg; i++)
    if (m_verts[i] == v)
      return true;
  return false;
}



void Face::AddEdge (Edge * e)
{
  if (m_eDeg == m_deg)
    {
      cerr << "ERROR: Face at " << this <<
        " cannot add another edge as a neighbor." << endl;
    }
  else
    m_edges[m_eDeg++] = e;
}


void
Face::AddVert (Vert * v)
{
  if (m_vDeg == m_deg)
    {
      cerr << "ERROR: Face at " << this <<
        " cannot add another vert as a neighbor." << endl;
    }
  else
    m_verts[m_vDeg++] = v;
}

// trang: holding the angles of vertices ----------------------------------
void
Face::AddVertAngle( double dValue )
{
    if (m_angleDeg == m_deg)
      {
        cerr << "ERROR: Face at " << this <<
          " cannot add another vert as a neighbor." << endl;
      }
    else
      m_vert_angles[m_angleDeg++] = dValue;
}
// ------------------------------------------------------------------------

// If the information about the last node in the path on this face is not
// already set, set it.  Or, if gn == 0, clear all of the settings
void
Face::SetLastInPathAdj (GraphNode * gn)
{
  // set information for use in shortcutting
  if (!gn || !m_lastInPath)
    m_lastInPath = gn;

  // set the information for all adjacent edges and verts too
  for (int i = 0; i < 3; i++)
    {
      if (!gn || !m_edges[i]->m_lastInPath)
        m_edges[i]->m_lastInPath = gn;
      if (!gn || !m_verts[i]->m_lastInPath)
        m_verts[i]->m_lastInPath = gn;
    }
}

// Don't need to notify any faces about this node since it's only on this face
void
Face::SetLastInPath (GraphNode * gn)
{
  SetLastInPathAdj (gn);
}

// Add the nearest neighbor to the list of nearest neighbors on the face
void
Face::AddNN (int nn)
{
  m_NNs.pushx (nn);
}

const ArrayInt &
Face::ReturnNNs () const
{
  return m_NNs;
}

int Face::getSideSign (const Point & pnt) const
{
  if (!m_normCalcd)
    ((Face *) this)->calcNormal ();
  RT q (DotProduct (pnt - Verts ()[0]->Location (), m_exact_norm));

  return sign (q);
}

RT
Face::getPointValue (const Point & pnt) const
{
  if (!m_normCalcd)
    ((Face *) this)->calcNormal ();
  RT q (DotProduct (pnt - Verts ()[0]->Location (), m_exact_norm));

  return q;
}


bool
Face::isFaceOrNeigborFaceVisible (const Point & pnt) const
{
  if (getSideSign (pnt) > 0)
    return true;

  Edge *e;
  int i;

  for ( i = 0; i < EdgeDegree (); ++i )
    {
      e = Edges ()[i];
      if ((e->F1 () != this) && (e->F1 ()->getSideSign (pnt) > 0))
        return true;
      if ((e->F2 () != this) && (e->F2 ()->getSideSign (pnt) > 0))
        return true;
    }

  return false;
}

bool Face::isFaceNN (const Point & curr_pnt) const
{
  //leda_array<Vec>   vInwardPerps( 3 );
  vector < Vec > vInwardPerps (3);
  int signSum = 0;
  int i, sgn;

  assert (EdgeDegree () == 3);
  if (!m_normCalcd)
    ((Face *) this)->calcNormal ();

  //printf( "side: %d\n", getSideSign( curr_pnt ) );

  //printf( "isFaceNN: " );
  for (i = 0; i < 3; i++)
    {
 
      vInwardPerps[i] = CrossProduct (m_exact_norm,
                                      Verts ()[(i +
                                                1) % 3]->Location () -
                                      Verts ()[i]->Location ());
      sgn =
        sign (DotProduct(vInwardPerps[i], curr_pnt - Verts ()[i]->Location ()));

      //printf( "%3d ",sgn );
      signSum +=  sgn;
    }
  //printf( "\n" );

  // the first clause makes sure that the current point and the face
  // are on the same side of the polytope and the second clause
  // makes sure that the current point is in the region of the face
  if ((sign (DotProduct (getExactNorm (), Vec (CGAL::ORIGIN, curr_pnt))) >= 0)
      && (abs (signSum) == EdgeDegree ()))
    return true;
  else
    return false;

  return false;
}

bool Face::isOn (const Point & pnt) const
{
  if (getSideSign (pnt) != 0)
    {
      // printf ("Face::isOn: not on face at all!\n");
      fflush (stdout);
      return false;
    }
  //printf( "on FACE!!!\n" );
  //fflush( stdout );

  for (int ind = 0; ind < EdgeDegree (); ind++)
    if (m_edges[ind]->isOnEdge (pnt))
      return true;

  //printf( "trying if isFaceNN!\n" );
  return isFaceNN (pnt);
}


bool
Face::getNearestSqrDistance (const Point & pnt,
                             RT & res, NNLocation & test_nnloc, int walk_id)
{
  //bool  f_ret;
  bool fInitDist, fRes;

  fInitDist = false;
  //printf( "face -bogi\n" );
  assert (walk_id > 0);

  //f_ret = false;
  if (getSideSign (pnt) < 0)
    return false;

  if (isFaceNN (pnt))
    {
      test_nnloc.m_loc = (FEVBase *) this;
      test_nnloc.m_locType = NNLocation::FACE;

      RT q = DotProduct (pnt - Verts ()[0]->Location (),
                         getExactNorm ());
      res = q * q / getExactNorm ().squared_length ();
      //printf( "RES: %g\n", res.to_double() );
      return true;
    }
  //Edge  * e;
  int i;

  res = 0;
  //for ( i = 0, e = Edges()[ i ]; i < EdgeDegree();
  //      e = Edges()[++i] ) {
  for (i = 0; i < EdgeDegree (); ++i)
    {
      RT dist;

      NNLocation tmp_loc;
      fRes = Edges ()[i]->getSqrDistance (pnt, dist, tmp_loc, walk_id);
      //assert( fRes );
      if (!fRes)
        continue;

      if ((!fInitDist) || (dist < res))
        {
          fInitDist = true;
          res = dist;
          test_nnloc = tmp_loc;
        }
    }

  return fInitDist;
}

void
Face::writeNNLocInfo (Point & pnt, NNLocation & out_loc)
{
  //bool  f_ret;
  //bool  fInitDist, fRes;
  Edge *e;
  int i;

  //printf( "Generating info for point: " );
  /*
     point_print( pnt );
     printf( "\n" );
     fflush( stdout );
   */
  for (i = 0, e = Edges ()[i]; i < EdgeDegree (); e = Edges ()[++i])
    {
      if (e->isOnEdge (pnt))
        {
          pnt = e->project_dbl (pnt);
          e->writeNNLocRecord (pnt, out_loc);
          //printf( "On edge...\n" );fflush( stdout );

          return;
        }
    }

  out_loc.m_loc = (FEVBase *) this;
  out_loc.m_locType = NNLocation::FACE;
}


Vert *
Face::getNonAdjVertex (const Edge * e) const
{
  bool f_edge_found = false;

  /* lets first verify that we got a real edge that belogns to this
     face! */
  for (int ind = 0; ind < EdgeDegree (); ind++)
    {
      if (m_edges[ind] == e)
        {
          f_edge_found = true;
          break;
        }
    }
  assert (f_edge_found);

  {
    for (int ind = 0; ind < VertDegree (); ind++)
      {
        if (m_verts[ind] == e->V1 ())
          continue;
        if (m_verts[ind] == e->V2 ())
          continue;

        return m_verts[ind];
      }
  }

  assert (false);

  return NULL;
}

void
Face::checkIfConvex () const
{
  for (int ind = 0; ind < EdgeDegree (); ind++)
    m_edges[ind]->checkIfConvex ();
}

void Face::dump () const
{
  printf ("/-- Face ----------------------------------------\n");

  for (int ind = 0; ind < VertDegree (); ind++)
    {
      printf ("%d: ", ind);
      point_print (m_verts[ind]->Location ());
      printf ("\n");
    }

  printf ("\\-----------------------------------------------\n");
}


bool Face::isVertexAdj (const Vert *v)
{
  for (int ind = 0; ind < m_vDeg; ind++)
    if (m_verts[ind] == v)
      return true;

  return false;
}

// trang: ckeck if the face is adjacent with a vertex and get the vertex
bool Face::isVertexAdj (const Vert *v, int & iVertex)
{
  for (int ind = 0; ind < m_vDeg; ind++)
    if (m_verts[ind] == v)
    {
        iVertex = ind;
        return true;
    }
  //
  return false;
}

//
Edge *
Face::getOtherAdjEdge (const Vert * v, const Edge * e)
{
  for (int ind = 0; ind < m_vDeg; ind++)
    {
      if (m_edges[ind] == e)
        continue;
      if (m_edges[ind]->Contains (v))
        return m_edges[ind];
    }

  assert (false);

  return NULL;
}


Edge *
Face::getVertAdjEdge (const Vert * v)
{
  for (int ind = 0; ind < m_vDeg; ind++)
    if (m_edges[ind]->Contains (v))
      return m_edges[ind];

  assert (false);

  return NULL;
}

Point Face::project (const Point & pnt) const
{
  cout << "m_normCalcd: " << endl;
  if (!m_normCalcd)
    ((Face *) this)->calcNormal ();

  RT coff (DotProduct (pnt - Verts ()[0]->Location (), m_exact_norm)
           / m_exact_norm.squared_length ());
  Point new_pnt (pnt - coff * m_exact_norm);

  assert (getSideSign (new_pnt) == 0);
  if (isOn (new_pnt))
    return new_pnt;

  //dump();

  RT dist;
  int pos = 0;

  dist = m_edges[0]->SqrDistance (new_pnt);
  pos = 0;
  for (int ind = 1; ind < EdgeDegree (); ind++)
    {
      RT
        candid_val = m_edges[ind]->SqrDistance (new_pnt);

      if (candid_val < dist)
        {
          pos = ind;
          dist = candid_val;
        }
    }

  return m_edges[pos]->project_dbl (new_pnt);
}

Vec Face::exactNormal () const
{
  return getExactNorm ();
}

bool FEVBase::isVisibleEdge (const Edge * e) const
{
  for (int ind = 0; ind < FaceDegree (); ind++)
    {
      Face *f (Faces ()[ind]);
      if (f->Contains (e))
        return true;
    }

  return false;
}

bool FEVBase::isVisiblePoint (const Point & pnt) const
{
  for (int ind = 0; ind < FaceDegree (); ind++)
    {
      Face *f (Faces ()[ind]);
      if (f->isOn (pnt))
        return
          true;
    }

  return false;
}

Vec Face::UnitNormal () const
{
  // only calculate the normal once
  if (!m_normCalcd)
    ((Face *) this)->calcNormal ();

  return m_aprx_norm;
}


Vec
Face::NormalDir () const
{
  // only calculate the normal once
  if (!m_normCalcd)
    ((Face *) this)->calcNormal ();

  return m_exact_norm;
}

void
Face::findFacePlanePoint (const Point & pnt_a,
                          const Point & pnt_b,
                          Point & out_pnt, NNLocation & out_loc)
{

  Plane pl (m_verts[0]->Location (),
            m_verts[1]->Location (), m_verts[2]->Location ());

  //int  val pl.intersection( pnt_a, pnt_b, out_pnt );

  CGAL::Object result = CGAL::intersection (Line (pnt_a, pnt_b), pl);
  if (const Point * ip = CGAL::object_cast < Point > (&result))
    {
      out_pnt = *ip;

      if (!isOn (out_pnt))
        {
          //memset (&out_loc, 0, sizeof (NNLocation));
          out_pnt = approx_pnt (out_pnt);

          return;
        }
      writeNNLocInfo (out_pnt, out_loc);
    }
  else
    assert (false);
}


/* Face.C - End of File ------------------------------------------*/
