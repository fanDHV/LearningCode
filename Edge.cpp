/**************************************************************************
 * Edge.cc
 * 
 * Implementation of the class that represents an edge on the polytope 
\**************************************************************************/
#include  <assert.h>
#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>

#include  "sarray.h"

#include  "Edge.h"
#include  "Face.h"
#include  "Vert.h"
#include  "Functions.h"
#include  "Poly.h"
#include "myCGAL.h"

Edge::Edge (Poly * p, Vert * v1, Vert * v2):FEVBase (p)
{
  m_v1 = v1;
  m_v2 = v2;
  m_f1 = m_f2 = 0;
  m_lastInPath = 0;
}

Edge::Edge (Poly * p, const Edge * e):FEVBase (p)
{
  m_v1 = e->m_v1;
  m_v2 = e->m_v2;
  m_f1 = e->m_f1;
  m_f2 = e->m_f2;
}

Edge::~Edge ()
{
}

int
Edge::FaceDegree () const
{
  return 2;
}

Face **Edge::Faces () const
{
  ((Face **) faces_arr)[0] = m_f1;
  ((Face **) faces_arr)[1] = m_f2;
  
  return const_cast < Face ** >(faces_arr);
}

Vert *Edge::V1 () const
{
  return m_v1;
}

Vert *Edge::V2 () const
{
  return m_v2;
}

Face *Edge::F1 () const
{
  return m_f1;
}

Face *Edge::F2 () const
{
  return m_f2;
}

Vec
Edge::exactNormal () const
{
  return m_f1->exactNormal ();
}

Vec
Edge::UnitNormal () const
{
  return m_f1->UnitNormal ();
}

Vec
Edge::NormalDir () const
{
  return m_f1->NormalDir ();
}

Face *Edge::Face1 ()
{
  return m_f1;
}

bool Edge::ContainsFace (const Face * f) const
{
  return ((m_f1 == f) || (m_f2 == f));
}

bool
Edge::ContainsEdge (const Edge * e) const
{
  return (this == e);
}

bool
Edge::Contains (const FEVBase * fev) const
{
  return fev->Contains (this);
}

bool
Edge::Contains (const Face * f) const
{
  return ContainsFace (f);
}

bool
Edge::Contains (const Edge * e) const
{
  return ContainsEdge (e);
}

bool
Edge::Contains (const Vert * v) const
{
  return ((m_v1 == v) || (m_v2 == v));
}

void Edge::addFace (Face * f)
{
  if (!m_f1)
    m_f1 = f;
  else if (!m_f2)
    m_f2 = f;
  else
    {
      cerr << "ERROR: Edge at " << this <<
        " cannot add third face as neighor." << endl;
      exit (-1);
    }
  if ((!f->isVertexAdj (V1 ())) || (!f->isVertexAdj (V2 ())))
    {
      printf ("something wrong with this face: \n");
      f->dump ();
      
      //assert( f->isVertexAdj( V1() ) );
      //      assert( f->isVertexAdj( V2() ) );
    }
}

bool Edge::hasVert (const Vert * v)
{
  return ((v == m_v1) || (v == m_v2));
}

// If the information about the last node in the path on this edge is not
// already set, set it.  Or, if gn == 0, clear all of the settings
void
Edge::SetLastInPathAdj (GraphNode * gn)
{
  // set information for use in shortcutting
  if (!gn || !m_lastInPath)
    m_lastInPath = gn;
  
  // set the information for all adjacent faces and verts too
  if (!gn || !m_f1->m_lastInPath)
    m_f1->m_lastInPath = gn;
  if (!gn || !m_f2->m_lastInPath)
    m_f2->m_lastInPath = gn;
  if (!gn || !m_v1->m_lastInPath)
    m_v1->m_lastInPath = gn;
  if (!gn || !m_v2->m_lastInPath)
    m_v2->m_lastInPath = gn;
}

// Tell both adjacent faces that this node is that last one in the path
// on them (a point on the edge is on both faces)
void
Edge::SetLastInPath (GraphNode * gn)
{
  m_f1->SetLastInPathAdj (gn);
  m_f2->SetLastInPathAdj (gn);
}

// Tell both adjacent faces that there is a nearest neighbor on the edge 
// (once again, the nearest neighbor is also on the faces)
void
Edge::AddNN (int nn)
{
  m_f1->AddNN (nn);
  m_f2->AddNN (nn);
}


bool Edge::isInSlab (const Point & pnt) const
{
  Vec
    vP (pnt - V1 ()->Location ());
  /*    
        cout << "isInSlab - vP: " << "(" 
        << vP.xcoord().to_double() << ", "
        << vP.ycoord().to_double() << ", "
        << vP.zcoord().to_double() << ")"
        << "\n";
  */
  // check to see if the current point is within the "skinny" faces
  // of the wedge created by the edge
  Vec
    vV1 = (V2 ()->Location () - V1 ()->Location ());
  Vec
    vV2 = -vV1;
  
  if (sign (DotProduct (vV1, vP)) < 0)
    return
      false;
  Vec
    vQ (pnt - V2 ()->Location ());
  
  if (sign (DotProduct (vV2, vQ)) < 0)
    return
      false;

  return
    true;
}

/*
  inline leda_rational point_sqr_length( const Point  & pnt ) 
  {
    return  pnt.xcoord() * pnt.xcoord()
    + pnt.ycoord() * pnt.ycoord()
        + pnt.zcoord() * pnt.zcoord();
}
*/
inline
RT
point_sqr_length (const Point & pnt)
{
  Vec
    v (CGAL::ORIGIN, pnt);
  return v.squared_length ();
}


void
Edge::checkIfConvex () const
{
  const Point & pnt_f1 (F1 ()->getNonAdjVertex (this)->Location ());
  const Point & pnt_f2 (F2 ()->getNonAdjVertex (this)->Location ());
  
  if (F1 ()->getSideSign (pnt_f2) > 0)
    {
      printf ("1 - Convexity failure at edge: %p\n", this);
      
      F1 ()->dump ();
      F2 ()->dump ();

      printf ("pnt_f2: \n");
      point_print (pnt_f2);
      printf ("\n");
      
      cout << "F1: " << F1 ()->getPointValue (pnt_f2) << "\n";
      cout << "F1: " << F1 ()->getPointValue (pnt_f2).
        to_double () << "\n";
      
      Point m_C = m_poly->computeRealCenter ();

      cout << "m_cF1: " << F1 ()->getPointValue (m_C) << "\n";
      cout << "m_cF1: " << F1 ()->getPointValue (m_C).
        to_double () << "\n";
      
      fflush (stdout);
      fflush (stderr);
      
      assert (false);
    }
  if (F2 ()->getSideSign (pnt_f1) > 0)
    {
      fflush (stdout);
      fflush (stderr);
           
      F1()->dump();
      F2()->dump();
      cout << pnt_f1 << endl;
      cout << pnt_f2 << endl;
      //cout << F1()->getExactNorm() << endl;
      cout << "2 - Convexity failure at edge: " 
           << this->V1()->Location() << " -- "
           << this->V2()->Location() << endl;
      assert (false);
    }
}


Point Edge::project (const Point & pnt) const
{
  //assert( isInSlab( pnt ) );
  if (!isInSlab (pnt))
    {
      RT
        dist1 = CGAL::squared_distance (V1 ()->Location (), pnt);
      RT
        dist2 = CGAL::squared_distance (V2 ()->Location (), pnt);

      if (dist1 < dist2)
        return
          V1 ()->
          Location ();
      else
        return
          V2 ()->
          Location ();
    }
  
  Vec
    vV1 = (V2 ()->Location () - V1 ()->Location ());
  RT
    q =
    DotProduct (vV1, pnt - V1 ()->Location ()) / vV1.squared_length ();
  
  Point
    tmp (V1 ()->Location () + q * vV1);
  
  assert ((F1 ()->getSideSign (tmp) == 0)
          && (F2 ()->getSideSign (tmp) == 0));
  //    printf( "VV %d %d\n", F1()->getSideSign( tmp ),
  //        F2()->getSideSign( tmp ) );
  //printf( "VV1 %d %d\n", F1()->getSideSign( V1()->Location() ),
  //       F2()->getSideSign( V1()->Location() ) );
  // printf( "VV2 %d %d\n", F1()->getSideSign( V2()->Location() ),
  //        F2()->getSideSign( V2()->Location() ) );
  
  
  return Point (tmp);
}


Point Edge::project_dbl (const Point & pnt) const
{
  //assert( isInSlab( pnt ) );
  if (!isInSlab (pnt))
    {
      RT
        dist1 = CGAL::squared_distance (V1 ()->Location (), pnt);
      RT
        dist2 = CGAL::squared_distance (V2 ()->Location (), pnt);
      
      if (dist1 < dist2)
        return
          V1 ()->
          Location ();
      else
        return
             V2 ()->
          Location ();
         }
  
  Vec
    vV1 = (V2 ()->Location () - V1 ()->Location ());
  RT
    qx =
         DotProduct (vV1, pnt - V1 ()->Location ()) / vV1.squared_length ();
  RT
    q (qx.to_double ());
  //qx.normalize();
  //q.normalize();
  
  //cout << "q: " << q;
  //cout << "\nqx : " << qx << "\n";
  //fflush( stdout );
  
  Point
    tmp (V1 ()->Location () + q * vV1);
  
  assert ((F1 ()->getSideSign (tmp) == 0)
          && (F2 ()->getSideSign (tmp) == 0));
  /*
    printf( "VV %d %d\n", F1()->getSideSign( tmp ),
    F2()->getSideSign( tmp ) );
    printf( "VV1 %d %d\n", F1()->getSideSign( V1()->Location() ),
    F2()->getSideSign( V1()->Location() ) );
    printf( "VV2 %d %d\n", F1()->getSideSign( V2()->Location() ),
    F2()->getSideSign( V2()->Location() ) );
  */
  
  return Point (tmp);
}


RT Edge::SqrDistance (const Point & pnt) const
{
  if (!isInSlab (pnt))
    {
      RT
        dist1 = CGAL::squared_distance (V1 ()->Location (), pnt);
      RT
        dist2 = CGAL::squared_distance (V2 ()->Location (), pnt);
      
      if (dist1 < dist2)
        return
          dist1;
      else
        return
          dist2;
    }
  
  Vec
    vV1 = (V2 ()->Location () - V1 ()->Location ());
  RT
    q =
    DotProduct (vV1, pnt - V1 ()->Location ()) / vV1.squared_length ();
  //Point  proj_diff( pnt - q*vV1 - V1()->Location());
  Vec
    proj_diff (pnt - q * vV1 - V1 ()->Location ());
  
  //return  point_sqr_length( proj_diff );
  return proj_diff.squared_length ();
}


RT Edge::LineSqrDistance (const Point & pnt) const
{
  Vec
    vV1 = (V2 ()->Location () - V1 ()->Location ());
  RT
    q =
    DotProduct (vV1, pnt - V1 ()->Location ()) / vV1.squared_length ();
  //Point  proj_diff( pnt - q*vV1 - V1()->Location());
  Vec
    proj_diff (pnt - q * vV1 - V1 ()->Location ());
  
  //return  point_sqr_length( proj_diff );
  return proj_diff.squared_length ();
}

void
Edge::writeNNLocRecord (const Point & in_pnt, NNLocation & out)
{
  if (PointEq (in_pnt, V1 ()->Location ()))
    {
      out.m_locType = NNLocation::VERT;
      out.m_loc = V1 ();
      return;
    }
  if (PointEq (in_pnt, V2 ()->Location ()))
    {
      out.m_locType = NNLocation::VERT;
      out.m_loc = V2 ();
      return;
    }
  
  if (!isInSlab (in_pnt))
    {
      RT
        dist1 = CGAL::squared_distance (V1 ()->Location (), in_pnt);
      RT
        dist2 = CGAL::squared_distance (V2 ()->Location (), in_pnt);
      
      out.m_locType = NNLocation::VERT;
      
      if (dist1 < dist2)
        out.m_loc = V1 ();
      else
        out.m_loc = V2 ();
      return;
    }
  
  out.m_locType = NNLocation::EDGE;
  out.m_loc = (FEVBase *) this;
}


bool
Edge::getSqrDistance (const Point & in_pnt,
                      RT & dist, NNLocation & test_nnloc, int walk_id)
{
  //printf( "edge -bogi\n" );

  if ((getLastFaceWalkID () > 0) && (getLastFaceWalkID () == walk_id))
    {
      //printf( "edge breaking!\n" );
      return false;
    }
  setLastFaceWalkID (walk_id);

  //printf( "edge: %p get_sqr_dist: %p\n", this, &in_pnt );
  //printf( "VV dist1: %g\n", dist1.to_double() );
  //printf( "VV dist2: %g\n", dist2.to_double() );

  if (!isInSlab (in_pnt))
    {
      RT dist1 = CGAL::squared_distance (V1 ()->Location (), in_pnt);
      RT dist2 = CGAL::squared_distance (V2 ()->Location (), in_pnt);

      test_nnloc.m_locType = NNLocation::VERT;

      if (dist1 < dist2)
        {
          dist = dist1;
          test_nnloc.m_loc = V1 ();
        }
      else
        {
          dist = dist2;
          test_nnloc.m_loc = V2 ();
        }

      //printf( "EdgeDist1: %g\n", dist.to_double() );

      return true;
    }

  test_nnloc.m_locType = NNLocation::EDGE;
  test_nnloc.m_loc = (FEVBase *) this;

  Vec vV1 = (V2 ()->Location () - V1 ()->Location ());
  RT q =
    DotProduct (vV1, in_pnt - V1 ()->Location ()) / vV1.squared_length ();

  //Point  proj_diff( in_pnt - q*vV1 - V1()->Location());
  Vec proj_diff (in_pnt - q * vV1 - V1 ()->Location ());

  //dist = point_sqr_length( proj_diff );
  dist = proj_diff.squared_length ();

  //printf( "EdgeDist2: %g\n", dist.to_double() );
  return true;
}

bool
isSamePoint (const Point & a, const Point & b)
{
  //if  ( ( a.xcoord() == b.xcoord() )
  //      &&  ( a.ycoord() == b.ycoord() )
  //      &&  ( a.zcoord() == b.zcoord() ) )
  //    return  true;
  if (a.x () == b.x () && a.y () == b.y () && a.z () == b.z ())
    return true;

  return false;
}

     bool Edge::isOnEdge (const Point & pnt) const
     {
       /*
          printf( "Edge::isOnEdge( " );
          point_print( pnt );
          point_print( V1()->Location() );
          point_print( V2()->Location() );
          printf( "\n" );
        */

       if (PointEq (pnt, V1 ()->Location ()))
         return
           true;
       if (PointEq (pnt, V2 ()->Location ()))
         return
           true;

       Point
       prj (
 (pnt));

       //printf( "CMP:\n" );

       //point_print( prj );
       //point_print( pnt );
       //printf( "\n" );
       bool
         f_res = (prj == pnt);

       if (f_res)
         return
           true;

       //printf( "%d\n", (int)f_res );
       //printf( "_%d\n", (int)isSamePoint( prj, pnt ) );
       //  printf( "%g\n", (double)prj.sqr_dist( pnt ).to_double() );
       //if  ( (double)prj.sqr_dist( pnt ).to_double() < 1e-20 )
       if ((prj - pnt).squared_length ().to_double () < 1e-20)
         return
           true;

       return
         f_res;
     }

     const
       Point &
     Edge::getOtherPoint (const Point & pnt)
{
  if (V2 ()->Location () == pnt)
    return V1 ()->Location ();
  if (V1 ()->Location () == pnt)
    return V2 ()->Location ();

  assert (false);

  return V2 ()->Location ();
}

Face *
     Edge::getCommonFace (const Point & pnt) const
     {
       if (F1 ()->isOn (pnt))
         return F1 ();
       if (F2 ()->isOn (pnt))
         return F2 ();

       assert (false);

       return NULL;
     }


/* Edge.C - End of File ------------------------------------------*/
