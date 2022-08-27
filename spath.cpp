/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * spath.C -
 *
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#include  <assert.h>
#include  <stdlib.h>
#include  <stdio.h>
#include  <limits.h>
#include  <memory.h>

#include  "sarray.h"
#include  "Functions.h"
#include  "Poly.h"
#include  "Face.h"
#include  "Vert.h"
#include  "Edge.h"
#include  "Path.h"
#include  "rgeneric.h"
#include  "Timer.h"
#include  "spath.h"
#include  "history.h"
#include  "global-settings.h"

/*--- Constants ---*/


/*--- Start of Code ---*/

SPath::SPath (SPath * ptr)
  : list_first(NULL), list_last(NULL), p_poly(NULL)
{
  //memset( this, 0, sizeof (SPath) );
  init( ptr->list_first, ptr->p_poly );
}


SPath::SPath ()
  : list_first(NULL), list_last(NULL), p_poly(NULL)
{
  //memset (this, 0, sizeof (SPath));
}


SPath::SPath (GraphNode * first, Poly * _p_poly)
{
  init (first, _p_poly);
}

SPath::SPath (SPathNode * _first, Poly * _p_poly)
{
  init (_first, _p_poly);
}


void
SPath::init (GraphNode * gn, Poly * _p_poly)
{
  SPathNode *node;

  //memset (this, 0, sizeof (SPath));
  if (gn == NULL)
    return;

  p_poly = _p_poly;
  list_first = new SPathNode (gn);
  list_last = list_first;

  gn = gn->m_next;
  while (gn != NULL)
    {
      node = new SPathNode (gn);
      list_last->append (node);
      list_last = node;
      gn = gn->m_next;
    }
}


void
SPath::init (SPathNode * gn, Poly * _p_poly)
{
  SPathNode *node;

  //memset (this, 0, sizeof (SPath));
  if (gn == NULL)
    return;

  p_poly = _p_poly;
  list_first = new SPathNode (gn);
  list_last = list_first;

  gn = gn->next ();
  while (gn != NULL)
    {
      node = new SPathNode (gn);
      list_last->append (node);
      list_last = node;
      gn = gn->next ();
    }
}


void
SPath::init_first (SPathNode * gn, Poly * _p_poly)
{
  //SPathNode  * node;

  //memset (this, 0, sizeof (SPath));
  if (gn == NULL)
    return;

  p_poly = _p_poly;
  list_first = new SPathNode (gn);
  list_last = list_first;
}


void
SPath::init_first (GraphNode * gn, Poly * _p_poly)
{
  //SPathNode  * node;

  //memset (this, 0, sizeof (SPath));
  if (gn == NULL)
    return;

  p_poly = _p_poly;
  list_first = new SPathNode (gn);
  list_last = list_first;
}


SPathNode::SPathNode (Face * f, const Point & _pnt)
  : gn(NULL), p_next(NULL), p_prev(NULL), f_in_sky(false)
{
  //memset (this, 0, sizeof (SPathNode));

  pnt = _pnt;
  nnloc.m_locType = NNLocation::FACE;
  nnloc.m_loc = f;
}


SPathNode::SPathNode (const NNLocation & _nnloc, const Point & _pnt)
  : gn(NULL), p_next(NULL), p_prev(NULL), f_in_sky(false)
{
  //memset (this, 0, sizeof (SPathNode));

  pnt = _pnt;
  nnloc = _nnloc;
}


SPathNode::SPathNode (GraphNode * _gn)
  : gn(NULL), p_next(NULL), p_prev(NULL), f_in_sky(false)
{
  //cout << "Before memset..." << endl;
  //memset (this, 0, sizeof (SPathNode));
  //cout << "After memset" << endl;

  gn = _gn;
  pnt = gn->m_p;
  nnloc = gn->m_nnloc;
}


SPathNode::SPathNode (SPathNode * _gn)
  : nnloc(_gn->nnloc),
    pnt(_gn->point()),
    gn(_gn->gn),
    p_next(NULL), p_prev(NULL), f_in_sky(false)
{
  //memset (this, 0, sizeof (SPathNode));

  //gn = _gn->gn;
  //pnt = _gn->point ();
  //nnloc = _gn->nnloc;
}


void
SPathNode::append (SPathNode * node)
{
  //printf( "appending!\n" );

  p_next = node;
  node->p_prev = this;
}


int
SPathNode::getFaceAdjDegree ()
{
  return nnloc.m_loc->FaceDegree ();
}

Face *
SPathNode::getAdjFace (int ind)
{
  Face *ptr;
  ptr = nnloc.m_loc->Faces ()[ind];

  assert (ptr != NULL);

  //printf( "p: %p\n", ptr );

  return ptr;
}


SPathNode *
SPathNode::getLastPointAdjToFace (Face * f)
{
  SPathNode *curr, *last_adj;
  int adj_deg;

  last_adj = NULL;
  curr = next ();
  while (curr != NULL)
    {
      if (curr->nnloc.isEmpty ())
        {
          curr = curr->next ();
          continue;
        }

      adj_deg = curr->getFaceAdjDegree ();

      for (int ind = 0; ind < adj_deg; ind++)
        {
          if (curr->getAdjFace (ind) == f)
            {
              if (f->getSideSign (curr->pnt) != 0)
                continue;

              last_adj = curr;
              break;
            }
        }

      curr = curr->next ();
    }

  return last_adj;
}



void
SPathNode::shortcut_on_faces ()
{
  int adj_deg;
  SPathNode *last_adj;

  if (nnloc.isEmpty ())
    return;

  adj_deg = getFaceAdjDegree ();
  for (int ind = 0; ind < adj_deg; ind++)
    {
      Face *f = getAdjFace (ind);
      if (f->getSideSign (pnt) != 0)
        {
          printf ("Strange!\n");
          continue;
        }

      last_adj = getLastPointAdjToFace (f);
      if ((last_adj == NULL) || (last_adj == this) || (last_adj == next ()))
        continue;

      printf ("this: %p  next: %p, last_adj: %p\n", this, next (), last_adj);
      p_next = last_adj;
      last_adj->p_prev = this;
    }
}

void
SPath::dump ()
{
  SPathNode *curr;

  printf ("-- P a t h -------------------------------------\n");
  curr = first ();
  while (curr != NULL)
    {
      curr->dump ();
      curr = curr->next ();
    }
  printf ("- - - - - - - - - - - - - - --------------------\n");
  fflush (stdout);
}

void
SPath::shortcut_on_faces ()
{
  SPathNode *curr;

  //dump();
  curr = first ();
  while (curr != NULL)
    {
      curr->shortcut_on_faces ();
      curr = curr->next ();
    }
  //printf( "--------\n" );
  //dump();
}


double
SPath::getAprxLength () const
{
  double sum = .0;
  SPathNode *curr;

  curr = list_first;
  while( curr != NULL && curr->next () != NULL )
    {
     // sum += SGP_distance (curr->point (), curr->next ()->  point () ) ;
    // Le test SP SGP
         if (SGP)
         sum += SGP_distance (curr->point (), curr->next ()->  point () ) ;
         else
         sum += sqrt (CGAL::squared_distance ( curr->point (), curr->next ()->  point () ).to_double ());
      curr = curr->next ();
    }

  return sum;
}
/****Le 13/10/2019 purpose of calculate the length of SPath connecting prevp to nextp point
double
SPath::getSubpathLength (const Point &prevp, const Point &nextp) const
{
    int tmp_idx = 0;
    double sum = 0;
    SPathNode *tmp, *curr;
    // assert if prevp and nextp are on SPath?
    tmp = list_first;
    while( tmp != NULL && tmp->next () != NULL )
    {
    if (tmp->point() == prevp)
    tmp_idx +=1;
    if (tmp->point() == nextp)
    tmp_idx +=1;
    tmp = tmp ->next();
    // assert if prevp and nextp are on SPath tmp =2?
    assert (tmp_idx == 2);
    }
    curr = list_first;
    while( curr != NULL && curr->next () != NULL )
    { if ( curr->point() == prevp ) break;
    else
    curr = curr->next();
    }
    while ( curr->next () != NULL && curr->point() != nextp )
    {
         sum += SGP_distance (curr->point (), curr->next ()->  point () ) ;
         curr = curr->next ();
    }

    return sum;
}
*/

// trang (2015/08/30): print the list of the path ... purpose of printing global path
void
SPath::printNodeList () const
{
  SPathNode *curr;

  curr = list_first;
  while( curr != NULL )
    {
      cout << "\t "; point_print(curr->point()); cout << endl;
      //
      curr = curr->next ();
    }
}


void
SPathNode::dump ()
{
    // trang's test: ok
    //cout << "%% Dump the shortest path by Agarwal" << endl; exit(1);

    printf ("[ ");
    if (!nnloc.isEmpty ())
      nnloc.dump ();
    printf (", ");
    point_print (pnt);

    if (!nnloc.isEmpty ())
      {
        int adj_deg = getFaceAdjDegree ();

        for (int ind = 0; ind < adj_deg; ind++)
          {
            Face *f = getAdjFace (ind);
            printf (" %2d ", f->getSideSign (pnt));
          }
      }

    printf ("]\n");
}

void
SPath::pushTop (SPathNode * gn)
{
  pushTop (gn->loc (), gn->point ());
}


void
SPath::pushTop (GraphNode * gn)
{
  pushTop (gn->m_nnloc, gn->m_p);
}


bool
is_location_better (const NNLocation & candid, const NNLocation & old)
{
  if (candid.m_locType == old.m_locType)
    return false;
  if (old.m_locType == NNLocation::VERT)
    return false;
  if (old.m_locType == NNLocation::EDGE)
    {
      if (candid.m_locType == NNLocation::VERT)
        return true;
      else
        return false;           // candid.m_locType is face
    }

  if (old.m_locType == NNLocation::FACE)
    return true;

  assert (false);

  return false;
}

void
SPath::append (SPathNode * node)
{
  list_last->append (node);
  list_last = node;
}


double
distance3 (const Point & src, const Point & mid, const Point & trg)
{
  double sum;

  // sum = SGP_distance (src, mid)  + SGP_distance (mid, trg);
 // Le test SGP va SP
 if (SGP){
 sum = SGP_distance (src, mid)  + SGP_distance (mid, trg);
 }
 else
    sum = sqrt (CGAL::squared_distance ( src, mid ).to_double ())
            +sqrt (CGAL::squared_distance ( mid, trg ).to_double ());
  return sum;
}


void
SPathNode::shortcut_on_edge_inner (Edge * e, Poly * pPoly)
{
  if ((next () == NULL) || (next ()->next () == NULL))
    return;

  //printf( "shortcut on edge inner called!\n" ); fflush( stdout );

  assert (e->isOnEdge (next ()->point ()));

  const Point & target (next ()->next ()->point ());
  /*
     Plane  h1( e->V1()->Location(),
     e->V2()->Location(),
     point() );
     Plane  h2( e->V1()->Location(),
     e->V2()->Location(),
     target );
   */

  Face *pntFace (e->getCommonFace (point ()));
  Face *trgFace (e->getCommonFace (target));

  if (pntFace == trgFace)
    {
      assert (pntFace != NULL);
      p_next = next ()->next ();
      p_next->p_prev = this;
      return;
      //assert( pntFace != trgFace );
    }

  assert (pntFace->isOn (point ()));

  Vec vv (-pntFace->NormalDir ());
  Plane h1 (point (), vv);
  Vec vv2 (-trgFace->NormalDir ());
  Plane h2 (target, vv2);
  /*
     printf( "h1.oriented_side( m_C ): %d %d\n",
     (int)h1.oriented_side( pPoly->Center() ),
     (int)h2.oriented_side( pPoly->Center() ) ) ;
   */

  /*
     m_Hs = Plane(m_srcdest->m_srcPoint,
     -(m_srcdest->m_srcFace->Normal()));
     m_Ht = Plane(m_srcdest->m_destPoint,
     -(m_srcdest->m_destFace->Normal()));
   */
  /*
     cout << "Normals: ";
     cout << -pntFace->NormalDir() << "     -        ";
     cout << pntFace->NormalDir() << "\n";
     cout << DotProduct( pntFace->NormalDir(),
     pntFace->Verts()[ 0 ]->Location() -
     pntFace->Verts()[ 1 ]->Location() );
     cout << "VV"
     << h2.oriented_side( trgFace->Verts()[ 0 ]->Location() )
     << h2.oriented_side( trgFace->Verts()[ 1 ]->Location() )
     << h2.oriented_side( trgFace->Verts()[ 2 ]->Location() )
     << "\n";
     cout << "\n";
     cout << DotProduct( pntFace->NormalDir(),
     pntFace->Verts()[ 1 ]->Location() -
     pntFace->Verts()[ 2 ]->Location() );
     cout << "\n";

     cout << "\n";

     cout << "h1 norm: " << -pntFace->NormalDir() << "\n";
     cout << "h2 norm: " << -trgFace->NormalDir() << "\n";
   */
  Point tmp_ridge_point = ExtRidgePoint (point (), h1,
                                         -pntFace->NormalDir (),
                                         target,
                                         h2, -trgFace->NormalDir (),
                                         pPoly->Center ());
  //cout << "lsd: " << e->LineSqrDistance( tmp_ridge_point );

  Point ridge_point = e->project_dbl (tmp_ridge_point);
  /*printf( "\ntmp, ridge_point\n" );
     point_print( tmp_ridge_point );
     printf( "\n" );
     point_print( ridge_point );
     printf( "\n" );
   */
  if (ridge_point == next ()->point ())
    return;

  double old_dist, tmp_dist, new_dist;

  old_dist = distance3 (point (), next ()->point (),
                        next ()->next ()->point ());
  tmp_dist = distance3 (point (), tmp_ridge_point,
                        next ()->next ()->point ());
  new_dist = distance3 (point (), ridge_point, next ()->next ()->point ());
  (void) tmp_dist;              // kill warnning
  //assert( tmp_dist <= old_dist );
  if (new_dist >= old_dist)
    {
      printf ("Faces are probably coplanar - shortcutting aborted!\n");
      return;
    }
  /*
     printf( "vvold distance: %g\n", distance3( point(),
     next()->point(),
     next()->next()->point() ) );

     printf( "new distance: %g\n", distance3( point(),
     ridge_point,
     next()->next()->point() ) );
     printf( "tmp distance: %g\n", distance3( point(),
     tmp_ridge_point,
     next()->next()->point() ) );
     printf( "e->v1         : %g\n", distance3( point(),
     e->V1()->Location(),
     next()->next()->point() ) );
     printf( "e->v2         : %g\n", distance3( point(),
     e->V2()->Location(),
     next()->next()->point() ) );
     printf( "e->V1: \n" );
     point_print( e->V1()->Location() );
     printf( "\n" );
     printf( "e->V2: \n" );
     point_print( e->V2()->Location() );
     printf( "\n" );
   */


  // We found a better point. Rewrite the record
  e->writeNNLocRecord (ridge_point, next ()->nnloc);

  SPathNode *ptr = new SPathNode (next ()->nnloc, ridge_point);

  ptr->p_next = next ()->p_next;
  ptr->p_prev = this;

  if (ptr->p_next != NULL)
    ptr->p_next->p_prev = ptr;

  p_next = ptr;
  //cout << ptr->point() << "\n";
  //Point p( ptr->point().xcoord().normalize(),
  //         ptr->point().ycoord().normalize(),
  //         ptr->point().zcoord().normalize() );
  Point p (ptr->point ().x (), ptr->point ().y (), ptr->point ().z ());
  //cout << "p: " << p << "\n";

  //fflush( stdout );
  //printf( "shortcut on edge inner done!\n" ); fflush( stdout );
}


// is mid lies on the middle of a common edge to src and trg?
// If so -> return this edge.
Edge *
getCommonEdge (NNLocation & src, NNLocation & mid, NNLocation & trg)
{
  if (mid.isFace ())
    return NULL;

  if (mid.isEdge ())
    {
      if (src.loc ()->isVisibleEdge (mid.edgePtr ())
          && trg.loc ()->isVisibleEdge (mid.edgePtr ()))
        return mid.edgePtr ();
      else
        return NULL;
    }

  // middle is a vertex. Maybe on of its adjacent edges answers our
  // requirement?
  Vert & v (mid.vertex ());
  for (int ind = 0; ind < v.EdgeDegree (); ind++)
    {
      Edge *e (v.Edges ()[ind]);

      if (src.loc ()->isVisibleEdge (e) && trg.loc ()->isVisibleEdge (e))
        return e;
    }

  return NULL;
}


bool
isOnCommonFace (NNLocation & src, NNLocation & mid, NNLocation & trg)
{
  for (int ind = 0; ind < src.loc ()->FaceDegree (); ind++)
    {
      Face *f (src.loc ()->Faces ()[ind]);
      if (mid.loc ()->Contains (f) && trg.loc ()->Contains (f))
        return true;
    }

  return false;
}


void
SPathNode::shortcut_on_edge (Poly * pPoly)
{
  if ((next () == NULL) || (next ()->next () == NULL))
    return;

  // shoudl we shortcut on a common face?
  if (isOnCommonFace (nnloc, next ()->nnloc, next ()->next ()->nnloc))
    {
      p_next = next ()->next ();
      return;
    }
  Edge *e;
  e = getCommonEdge (nnloc, next ()->nnloc, next ()->next ()->nnloc);
  //printf( "Common edge: %p\n", e );
  if (e == NULL)
    return;

  /*
     point_print( next()->point() );
     printf( "\n" );
     point_print( e->V1()->Location() );
     printf( "\n" );
     point_print( e->V2()->Location() );
     printf( "\n" );
   */

  if (e->isOnEdge (next ()->point ()))
    shortcut_on_edge_inner (e, pPoly);
}


void
SPath::shortcut_on_edge ()
{
  SPathNode *node = list_first;

  //printf( "shortcut on edge ___\n" ); fflush( stdout );

  while (node != NULL)
    {
      node->shortcut_on_edge (p_poly);
      node = node->next ();
    }
}


void
SPath::pushTop (const NNLocation & nnloc, const Point & _pnt)
{
  Point pnt_tmp = _pnt;
  //Point  pnt( pnt_tmp.xcoord().normalize(),
  //            pnt_tmp.ycoord().normalize(),
  //            pnt_tmp.zcoord().normalize() );
  Point pnt (pnt_tmp.x (), pnt_tmp.y (), pnt_tmp.z ());
  //printf( "list_last: %p\n", list_last );
  //fflush( stdout );

  //if  ( list_last != NULL ) {
  //printf( "list_last point: " );
  //point_print( list_last->point() );
  //}

  //printf( "\n" );
  //fflush( stdout );

  //cout << "DISTT: " << list_last->point().sqr_dist( pnt ) << "\n";
  if (list_last->point () == pnt)
    {
      if (is_location_better (nnloc, list_last->loc ()))
        list_last->setLocation (nnloc);
      return;
    }

  //printf( "appending point: "  );
  //point_print( pnt );
  //printf( "\n" );
  //fflush( stdout );

  SPathNode *node = new SPathNode (nnloc, pnt);

  append (node);

  //printf( "new_list_last: %p\n", list_last );
  //    fflush( stdout );
  //printf( "xlist_last point: " );
  //fflush( stdout );
  //point_print( list_last->point() );
  //printf( "\n" );
  //fflush( stdout );
}


// What is the const c, so that dir * c + src == pnt?
double
getProjectionValDbl (const Point & src, const Vec & dir, const Point & pnt)
{
  //Vec   vV1 = (V2()->Location() - V1()->Location());
  RT q = DotProduct (dir, pnt - src) / dir.squared_length ();

  return q.to_double ();
}


class Sweeper
{
private:
  Edge ** edges;
  int edges_num, max_edges;
  Plane plane;
  /* Point  direction; // this is the ray defined by the vertex */
  Vec direction;                // this is the ray defined by the vertex
  double low, high;             // [0,high] is the range of our parameter
  Point base1, base2;
  const Point & peak;
  Line line;
  bool f_init_range;

  double min_param, min_length;

public:
    Sweeper (const Point & _base1, const Point & _peak,
             const Point & _base2,
             int _max_degree):plane (_base1, _peak, _base2), peak (_peak)
  {


    assert (_max_degree > 0);

    base1 = _base1;
    base2 = _base2;
    //peak = _peak;
    //direction = CrossProduct( base1 - peak, peak - base2 );
    //line = Line( peak, peak.to_vector() + direction.to_vector() );

    low = high = 0;
    f_init_range = false;
    max_edges = _max_degree;
    edges_num = 0;

    edges = (Edge **) calloc (sizeof (Edge *), _max_degree);
    assert (edges != NULL);
  }

  void computeRangeForEdge (Edge * e)
  {
    //printf( "computeRangeForEdge( %p )\n", e );
    //fflush( stdout );

    const Point & pnt (e->getOtherPoint (peak));
    /*
       printf( "3 points : " );
       point_print( base1 );
       printf( "\n" );
       point_print( pnt );
       printf( "\n" );
       point_print( base2 );
       printf( "\n" );
       fflush( stdout );
     */

    Plane pl_tmp (base1, pnt, base2);
    Point out_pnt;
    double val;

    /*
       printf( "PEAK 3 POints : " );
       point_print( base1 );
       printf( "\n" );
       point_print( peak );
       printf( "\n" );
       point_print( base2 );
       printf( "\n" );
     */
    //int  sign_of_point = plane.oriented_side( pnt );
    //if  ( sign_of_point == 0 )
    if (plane.has_on (pnt))
      return;

    //if  ( plane.oriented_side( pnt ) < 0 ) {
    if (plane.has_on_negative_side (pnt))
      {
        printf ("sign = %d\n", plane.oriented_side (pnt));
        /* Hoai - ridiculous */
        //assert( plane.oriented_side( pnt ) > 0 );
        assert (plane.oriented_side (pnt) > 0);
      }

    /* Hoai - 09/08/2012 - only checking for single-point intersection */
    //assert( pl_tmp.intersection( peak, peak.to_vector()
    //                             + direction.to_vector(), out_pnt )
    //        == 1 );

    // the followign is subtle. Our representation is not closed,
    // and it might include negativ values. However, we are not
    // interested in such values anyway, as they will exclude one
    // of the neigbors of our vertex. So we can jsut skip such
    // edges.
    //if   ( plane.oriented_side( out_pnt ) <= 0 )
    if (plane.oriented_side (out_pnt) <= 0)
      return;
    //assert( plane.oriented_side( out_pnt ) > 0 );
    assert (plane.oriented_side (out_pnt) > 0);

    val = getProjectionValDbl (peak, direction, out_pnt);
    if (!f_init_range)
      {
        f_init_range = true;
        low = high = val;
        if (low < 0.0)
          high = 0;
        else
          low = 0;
      }
    if (val < low)
      low = val;
    if (val > high)
      high = val;
    printf ("RANGE: [%g, %g]\n", low, high);
  }

  void computeRange ()
  {
    //printf( "\n\n\n\n\nSTARTING RANGE COMPUTATION!\n" );
    //fflush( stdout );

    int ind = 0;
    //if  ( plane.oriented_side( edges[ 0 ]->getOtherPoint( peak ) ) == 0 )
    if (plane.oriented_side (edges[0]->getOtherPoint (peak)) == 0)
      ind++;
    if (ind >= edges_num)
      return;
    //if  ( plane.oriented_side( edges[ ind ]->getOtherPoint( peak ) ) == 0 ) {
    if (plane.oriented_side (edges[ind]->getOtherPoint (peak)) == 0)
      {
        printf ("Empty range - shoud not be possible\n");
        fflush (stdout);
        assert (false);
        return;
      }

    direction = edges[ind]->getOtherPoint (peak) - peak;
    line = Line (peak, peak + direction);

    {
      for (int ind = 0; ind < edges_num; ind++)
        computeRangeForEdge (edges[ind]);
    }
    if (!(((low == 0.0) && (high > 0.0)) || ((low < 0.0) && (high == 0.0))))
      {
        printf ("Range failure: [%g, %g]!\n", low, high);
        fflush (stdout);
        assert (false);
      }

  }

  double computeLengthForValue (double val)
  {
    Point new_pnt (peak + RT (val) * direction);
    Plane pl_tmp (base1, new_pnt, base2);
    double sum;
    Point prev_pnt;

    /*
       point_print( base1 );
       printf( "\n" );
       point_print( new_pnt );
       printf( "\n" );
       point_print( base2 );
       printf( "\n" );
       printf( "Direction!" );
       point_print( direction );
       printf( "\n" );
     */

    sum = 0.0;
    prev_pnt = base1;
    for (int ind = 0; ind < edges_num; ind++)
      {
        Point out_pnt;
        Edge & e (*(edges[ind]));
        if (plane.oriented_side (e.getOtherPoint (peak)) == 0)
          continue;

        /* Hoai - 09/08/2012 - only checking single-point intersection */
        //assert( pl_tmp.intersection( e.V1()->Location(),
        //                             e.V2()->Location(),
        //                             out_pnt ) == 1 );
        // Le test SGP va SP
           // sum += SGP_distance (prev_pnt, out_pnt);
           if (SGP){
           sum += SGP_distance (prev_pnt, out_pnt);
           }
           else
           sum += sqrt (CGAL::squared_distance (prev_pnt, out_pnt).to_double ());
        prev_pnt = out_pnt;
      }
    // Hoai :
    // sum += sqrt (CGAL::squared_distance (prev_pnt, base2).to_double ());
    // Le sua
    // sum += SGP_distance (prev_pnt, base2);
    // Le test SGP va SP
    if (SGP){
           sum += SGP_distance (prev_pnt, base2);
           }
    else
           sum += sqrt (CGAL::squared_distance (prev_pnt, base2).to_double ());

    return sum;
  }

  int find_first_non_neg (Vert & v, int ind)
  {
    int count = 0;

    while (plane.oriented_side (v.getEdgeEndPoint (ind)) < 0)
      {
        ind++;

        if (ind >= v.EdgeDegree ())
          ind = 0;

        assert (count <= v.EdgeDegree ());
        count++;
      }
    return ind;
  }


  int non_neg_start (Vert & v)
  {
    int sign, zero_pos;

    zero_pos = -1;
    for (int ind = 0; ind < v.EdgeDegree (); ind++)
      {
        sign = plane.oriented_side (v.getEdgeEndPoint (ind));
        if (sign < 0)
          return find_first_non_neg (v, ind);
        if (sign == 0)
          zero_pos = ind;
      }

    if (zero_pos >= 0)
      return zero_pos;

    assert (false);

    return -1;
  }


  // we store
  void register_positive_range (Vert & v)
  {
    int ind = 0;
    bool f_not_pos = false;

    ind = non_neg_start (v);

    // Now, we just spill all the positive edges into our array
    //int  edges_num;
    edges_num = 0;
    while (!f_not_pos)
      {
        int endpoint_sign;

        endpoint_sign = plane.oriented_side (v.getEdgeEndPoint (ind));
        if (endpoint_sign < 0)
          {
            f_not_pos = true;
            break;
          }

        // we are cycling over ourself.
        if ((edges_num > 0) && (edges[0] == v.Edges ()[ind]))
          break;

        printf ("inserting: edge[ %d ]\n", ind);
        edges[edges_num++] = v.Edges ()[ind];
        if ((endpoint_sign == 0) && (edges_num > 1))
          {
            f_not_pos = true;
            break;
          }

        ind++;
        if (ind >= v.EdgeDegree ())
          ind = 0;
      }
    assert (edges_num > 0);

    printf ("edges_num; %d\n", edges_num);
    {
      for (int ind = 0; ind < v.EdgeDegree (); ind++)
        {
          int endpoint_sign = plane.oriented_side (v.getEdgeEndPoint (ind));
          printf ("%d: [%d]\n", ind, endpoint_sign);
        }
    }

    // The only problem that we might encouter now is that the
    // edges are in inverse ordering to what we want them to be.
    // I would like edges[ 0 ] to adjacent to the point base1. If
    // it is not, then we are in trouble. So lets be sure it does.
    assert (edges[0]->isVisiblePoint (base1)
            || edges[edges_num - 1]->isVisiblePoint (base1));
    assert (edges[0]->isVisiblePoint (base2)
            || edges[edges_num - 1]->isVisiblePoint (base2));

    if (!edges[0]->isVisiblePoint (base1))
      {
        for (int ind = 0; ind < edges_num / 2; ind++)
          {
            Edge *tmp;
            tmp = edges[ind];
            edges[ind] = edges[edges_num - 1 - ind];
            edges[edges_num - 1 - ind] = tmp;
          }
        //swap( edges[ ind ], edges[ edges_num - 1 - ind ] );
      }
    assert (edges[0]->isVisiblePoint (base1)
            && edges[edges_num - 1]->isVisiblePoint (base2));

    computeRange ();
  }

  void testValue (double weight)
  {
    double pos, res;

    //printf( "testValue( %g )\n", weight ); fflush( stdout );
    pos = low * weight + high * (1 - weight);
    assert ((low <= pos) && (pos <= high));

    res = computeLengthForValue (pos);
    //printf( "testValue: res: %g  min_length: %g  pos: %g\n",
    //        res, min_length, pos );
    if (res < min_length)
      {
        min_length = res;
        min_param = pos;
      }
  }

  double find_best (double &out_length)
  {
    min_param = 0.0;
    min_length = computeLengthForValue (0.0);

    printf ("min_length = %g\n", min_length);

    testValue (0.0001);
    testValue (0.9999);
    testValue (1.0);

    for (double val = 0.1; val < 1.0; val += 0.1)
      testValue (val);

    printf ("improved min_length = %g\n", min_length);

    out_length = min_length;

    return min_param;
  }

  SPathNode *appendIntoPath (SPathNode * node, double val, bool f_reverse)
  {
    int true_ind;

    Plane pl_tmp (base1, peak + RT (val) * direction, base2);

    for (int ind = 0; ind < edges_num; ind++)
      {
        Point out_pnt;
        if (f_reverse)
          true_ind = edges_num - ind - 1;
        else
          true_ind = ind;

        Edge & e (*(edges[true_ind]));

        /* Hoai - 09/08/2012 - only checking for single-point intersection */
        //assert( pl_tmp.intersection( e.V1()->Location(),
        //                             e.V2()->Location(),
        //                             out_pnt ) == 1 );
        Point e_pnt (e.project_dbl (out_pnt));
        NNLocation nnloc;

        e.writeNNLocRecord (e_pnt, nnloc);
        SPathNode *new_node;

        new_node = new SPathNode (nnloc, e_pnt);
        node->append (new_node);
        node = new_node;
      }

    return node;
  }
};


bool SPathNode::isOnVertex () const
{
  return
    nnloc.
    isVertex ();
}


void
SPathNode::rewrite_path (Sweeper & sw, double val, bool f_reverse)
{
  //printf( "val: %g\n", val );
  //fflush( stdout );

  assert (val >= 0.0);

  SPathNode *
  tail, *
    node;

  tail = next ()->next ();

  p_next = NULL;

  node = sw.appendIntoPath (this, val, f_reverse);

  // note that we had manipulated an interior interval of the
  // connected list of the path, so we do not need to fix the last
  // pointer in SPath
  node->p_next = tail;
  tail->p_prev = node;
}


void
SPathNode::shortcut_on_vertex ()
{
  if ((next () == NULL) || (next ()->next () == NULL))
    return;

  // should we shortcut on a common face?
  if (isOnCommonFace (nnloc, next ()->nnloc, next ()->next ()->nnloc))
    {
      p_next = next ()->next ();
      return;
    }

  if (!next ()->isOnVertex ())
    return;

  Vert & v ((next ()->nnloc.vertex ()));
  Sweeper sw_pos (point (), next ()->point (), next ()->next ()->point (),
                  v.FaceDegree ());
  Sweeper sw_neg (next ()->next ()->point (), next ()->point (), point (),
                  v.FaceDegree ());

  sw_pos.register_positive_range (v);
  sw_neg.register_positive_range (v);

  double val_pos, val_neg;
  double length_pos, length_neg;

  val_pos = sw_pos.find_best (length_pos);
  val_neg = sw_neg.find_best (length_neg);

  if (length_pos < length_neg)
    rewrite_path (sw_pos, val_pos, false);
  else
    rewrite_path (sw_neg, val_neg, true);
}


void
SPath::shortcut_on_vertex ()
{
  SPathNode *node = list_first;

  //printf( "shortcut on vertex` ___\n" ); fflush( stdout );

  while (node != NULL)
    {
      node->shortcut_on_vertex ();
      node = node->next ();
    }
}



// This shortcuts the path, finds the two-plane wedges, and projects the path
// onto the polytope.  Also removes duplicates
SPath *
SPath::shortcut ()
{
  SPath *sp;
  SPathNode *curr;
  Vec n1, n2;

  optimize ();

  if (list_first == NULL)
    return NULL;

  sp = new SPath ();

  assert (p_poly != NULL);

  curr = list_first;

  sp->init_first (curr, p_poly);
  while (curr->next () != NULL)
    {
      n1 = curr->loc ().m_loc->NormalDir ();
      n2 = curr->next ()->loc ().m_loc->NormalDir ();

      sp->pushTop (curr);
      debug ("n1 = (" << n1.x() << "," << n1.y() << "," << n1.z() << ")" );
      debug ("n2 = (" << n2.x() << "," << n2.y() << "," << n2.z() << ")" );
      if (n1 != n2)
        lift_and_shortcut (curr, curr->next (), *sp);

      curr = curr->next ();
    }
  sp->pushTop (curr);

  return sp;
}


// Shortcut betweeh start and end using a lifting, and then projecting
// down.
void
SPath::lift_and_shortcut (SPathNode * start, SPathNode * end, SPath & sp)
{
  Vec n1, n2;

  //printf( "SPath::lift_and_shrotcut!\n" ); fflush( stdout );

  assert (start->loc ().m_loc != NULL);
  assert (end->loc ().m_loc != NULL);

  n1 = start->loc ().m_loc->NormalDir ();
  n2 = end->loc ().m_loc->NormalDir ();

  //printf( "normals!\n" );  fflush( stdout );

  Plane Hs (start->point (), -n1);
  Plane Ht (end->point (), -n2);

  Point rp = ExtRidgePoint (start->point (), Hs, -n1,
                            end->point (), Ht, -n2, p_poly->Center ());
  /*
     printf( "center:\n" );
     point_print( p_poly->Center() );
     printf( "\n" );

     printf( "the new path before projection!\n" );
     point_print( start->point() );
     printf( "\n" );
     point_print( rp );
     printf( "\n" );
     point_print( end->point() );
     printf( "%d %d\n",
     ( start->point() == rp ),
     ( rp == end->point() ) );

     printf( "\n" ); fflush( stdout );
   */

  if (PointEq (start->point (), rp) || PointEq (rp, end->point ()))
    {
      sp.pushTop (start);
      sp.pushTop (end);
    }
  else
    sp.projectWedgePath (start, rp, end);
  //printf( "lift_and_shrotcut done!\n" ); fflush( stdout );
}


double
SPath::get_length (SPathNode * start, SPathNode * end)
{
  double sum = 0;
  //SPathNode  * curr;

  while ((start != NULL) && (start->next () != NULL) && (start != end))
    {
     // Hoai
     //    sum += sqrt (CGAL::squared_distance (start->point (),
      //                  start->next ()-> point ()).to_double ());
     // Le
     // sum += SGP_distance (start->point (),  start->next ()-> point ());
     // Le test SGP va SP
      if (SGP){
           sum += SGP_distance (start->point (),  start->next ()-> point ());
           }
           else
           sum += sqrt (CGAL::squared_distance (start->point (),
                        start->next ()-> point ()).to_double ());

      start = start->next ();
    }

  return sum;
}


double
SPath::lift_length (SPathNode * start, SPathNode * end)
{
  Vec n1, n2;
  double sum;

  n1 = start->loc ().m_loc->NormalDir ();
  n2 = end->loc ().m_loc->NormalDir ();

  Plane Hs (start->point (), -n1);
  Plane Ht (end->point (), -n2);

  Point rp = ExtRidgePoint (start->point (), Hs, -n1,
                            end->point (), Ht, -n2,
                            p_poly->Center ());

 // sum = SGP_distance (start->point (), rp);
  //sum += SGP_distance (rp, end->point ());
  // Le test SGP va SP
    if (SGP){
        sum = SGP_distance (start->point (), rp);
        sum += SGP_distance (rp, end->point ());
    }
    else
    {

        sum = sqrt (CGAL::squared_distance (start->point (), rp).to_double ());
        sum += sqrt (CGAL::squared_distance (rp, end->point ()).to_double ());
    }
  return sum;
}


void
SPath::shortcut_by_lifting (SPathNode * start, SPathNode * end)
{
  SPathNode *old_list_last, *tail, *ptr;
  SPath *new_sp;

  //printf( "shortcut_by_lifting\n" ); fflush( stdout );
  //dump();
  //printf( "\n\n\n\n" );

  new_sp = new SPath (new SPathNode (start), p_poly);

  new_sp->lift_and_shortcut (start, end, *new_sp);
  //new_sp->pushTop( start );
  //    new_sp->pushTop( end );

  if (get_length (start, end) <= new_sp->getAprxLength ())
    {
      delete new_sp;
      return;
    }

  //printf( "found a shortcut!\n" );
  //fflush( stdout );


  //dump();

  old_list_last = list_last;

  tail = end->next ();
  list_last = start;
  start->p_next = NULL;

  ptr = new_sp->list_first;
  while (ptr != NULL)
    {
      pushTop (ptr);
      ptr = ptr->next ();
    }

  //printf( "YYY\n" );
  //fflush( stdout );

  //dump();
  if (tail == NULL)
    return;

  list_last->p_next = tail;
  tail->p_prev = list_last;

  list_last = old_list_last;
  //printf( "XXX\n" );

  //dump();
}

void
SPath::shortcut_by_relifting_with_dist (int dist)
{
  SPathNode *start, *end;
  static int count = 0;

  optimize ();
  printf ("shortcut with dist: %d\n", dist);
  //dump();
  fflush (stdout);

  if (dist <= 1)
    return;

  start = list_first;
  end = start;
  for (int ind = 0; ind < dist; ind++)
    {
      end = end->next ();
      if (end == NULL)
        return;
    }

  while (end != NULL)
    {
      shortcut_by_lifting (start, end);
      //printf( "after shortcut by lifting!\n" );
      //dump();

      count++;
      if ((count & 0xf) == 0)
        {
          printf ("[%10d]\n", count++);
          fflush (stdout);
        }

      end = end->next ();
      start = start->next ();
    }
}


void
SPath::shortcut_by_relifting ()
{
  SPathNode *start, *end;

  optimize ();
  start = list_first;
  while (start != NULL)
    {
      end = start->next ();

      while (end != NULL)
        {
          shortcut_by_lifting (start, end);

          end = end->next ();
        }
      start = start->next ();
    }
}


void
SPath::shortcut_by_relifting_jump ()
{
  int len, bound;

  shortcut_by_relifting_with_dist (2);
  shortcut_by_relifting_with_dist (3);
  shortcut_by_relifting_with_dist (list_first->length () - 1);

  len = list_first->length ();
  bound = 4;
  while (bound < len)
    {
      shortcut_by_relifting_with_dist (bound);
      bound *= 2;
    }
}


// We compute the plane that passes through v0 and v1, and is
// penpendicular to the plane spanned by this triangle. Morover, v1 is
// positive in this plane.
Plane
getPositivePlane (Point & v0, Point & v1, Point & v2)
{
  //Point i = CrossProductDir(approx_vec( v1 - v0 ),
  //                          approx_vec( v2 - v1 ), 975678123 );
  //Point i = CrossProduct( approx_vec( v1 - v0 ),
  //approx_vec( v2 - v1 ) );
  Vec i = CrossProduct (approx_vec (v1 - v0),
                        approx_vec (v2 - v1));
  //Point  new_pnt( i + v0 );
  Point new_pnt (v0 + i);

  // v0, v2, new_pnt lives in the plane that I am interested in
  Plane h (v0, v2, new_pnt);
  if (h.oriented_side (v1) > 0)
    return h;
  Plane h1 (v2, v0, new_pnt);
  if (h1.oriented_side (v1) > 0)
    return h1;

  assert (false);

  return h1;
}

void
SPath::projectWedgePath (SPathNode * start, Point & mid, SPathNode * end)
{
  Edge *e1;
  //Face  * f;
  NNLocation loc;

  pushTop (start);

  Point v0 (start->point ());
  Point & v1 (mid);
  Point v2 (end->point ());
  //Point  vec_dir( v2 - v0 );
  Point exit_point (v0);

  // if any of these points are equal then we have no work to do
  if (PointEq (v0, v1) || PointEq (v1, v2) || PointEq (v2, v0))
    {
      //dump();
      //printf( "Nothing was DONE!\n" );
      //assert( false );
      return;
    }

  history_t hist;

  hist.init (v0);

  // all the interesting things happen on the side where
  // hPositive is, well, positive
  Plane hPositive (getPositivePlane (v0, v1, v2));

  //Plane   pathPlane(v0, CrossProductDir( approx_vec( v1-v0 ),
  //                                       approx_vec( v2-v0 ), 3215));
  Plane pathPlane (v0, CrossProduct (approx_vec (v1 - v0),
                                     approx_vec (v2 - v0)));

  e1 = NULL;
  //f = NULL;
  loc = start->loc ();
  //printf( "Current point: " );
  //point_print( start->point() );
  do
    {
      //printf( "-!  " );
      //loc().dump();
      //printf( "-!\n" );

      e1 = findExitEdge (loc, pathPlane, hPositive,
                         hist, exit_point, end->loc ());
      if (e1 == NULL)
        break;
      //sp.pushTop( e1->intersect( pathPlane ), e1 );

      loc.m_loc = e1;
      loc.m_locType = NNLocation::EDGE;
      if ((pathPlane.oriented_side (e1->V1 ()->Location ()) == 0)
          && (pathPlane.oriented_side (e1->V2 ()->Location ()) == 0))
        {
          assert (false);
        }
      if (pathPlane.oriented_side (e1->V1 ()->Location ()) == 0)
        {
          loc.m_loc = e1->V1 ();
          loc.m_locType = NNLocation::VERT;
        }
      if (pathPlane.oriented_side (e1->V2 ()->Location ()) == 0)
        {
          loc.m_loc = e1->V2 ();
          loc.m_locType = NNLocation::VERT;
        }
      //printf( "exit_point: \n" );
      //point_print( exit_point );
      //fflush( stdout );

      pushTop (loc, exit_point);
    }
  while (e1 != NULL);

  pushTop (end);
}


Edge *
SPath::findExitEdge (NNLocation & loc,
                     Plane & pathPlane, Plane & hPositive,
                     history_t & hist, Point & exit_point,
                     NNLocation & dest_loc)
{
  //bool  f_defined;
  RT val;
  Edge *out_edge = NULL;
  Face *f;
  Edge *e;
  Point q;
  bool contains;

  //f_defined = false;
  for (int j = 0; j < loc.m_loc->FaceDegree (); j++)
    {
      f = loc.m_loc->Faces ()[j];

      // we dont need to do any smart thing...
      if (dest_loc.m_loc->ContainsFace (f))
        return NULL;
      for (int i = 0; i < f->EdgeDegree (); i++)
        {
          e = f->Edges ()[i];

          if ((pathPlane.oriented_side (e->V1 ()->Location ()) != 0)
              && (pathPlane.oriented_side (e->V1 ()->Location ()) ==
                  pathPlane.oriented_side (e->V2 ()->Location ())))
            continue;
          q = IntersectLinePlane (e->V1 ()->Location (),
                                  e->V2 ()->Location (), pathPlane, contains);
          if (hist.isIn (q))
            continue;
          if (hPositive.oriented_side (q) < 0)
            continue;

          if (e->isInSlab (q))
            {
              exit_point = q;
              hist.register_pt (q);

              out_edge = e;
            }
        }
    }

  return out_edge;
}

int
     SPathNode::length () const
     {
       const SPathNode *curr;
       int len = 0;

         curr = this;
       while (curr != NULL)
         {
           curr = curr->p_next;
           len++;
         }

       return len;
     }


void
SPath::shortcutify_by_face (Face * f, SPathNode * curr, SPathNode * last)
{
  Point first_pnt, last_pnt;
  NNLocation first_loc, last_loc;

  printf ("shortcutify_by_face( %p, %p )\n", curr, last);
  //dump();

  assert (curr->prev () != NULL);
  //printf(
  /*
     printf( "[%d %d]",
     f->getSideSign( curr->prev()->point() ),
     f->getSideSign( curr->point() ) );
     printf( "[%d %d]",
     f->getSideSign( last->point() ),
     f->getSideSign( last->next()->point() ) );
     printf( "\n" );
     fflush( stdout );
   */

  f->findFacePlanePoint (curr->prev ()->point (), curr->point (),
                         first_pnt, first_loc);
  f->findFacePlanePoint (last->point (),
                         last->next ()->point (), last_pnt, last_loc);

  /*
     printf( "------------------@@\n" );
     point_print( first_pnt );
     first_loc.dump();
     printf( "\n" );
     point_print( last_pnt );
     printf( "\n" );
     fflush( stdout );
   */

  SPathNode *new_first, *new_last;

  new_first = new SPathNode (first_loc, first_pnt);
  new_last = new SPathNode (last_loc, last_pnt);

  new_last->p_prev = new_first;
  new_last->p_next = last->next ();

  new_first->p_next = new_last;
  new_first->p_prev = curr->prev ();

  curr->prev ()->p_next = new_first;
  last->next ()->p_prev = new_last;

  //dump();
}


void
SPath::shortcut_by_face (Face * f)
{
  SPathNode *curr, *last, *temp;

  //printf( "SPath::shortcut_by_face!\n" );
  //dump();
  //fflush( stdout );

  assert (f->getSideSign (list_last->point ()) <= 0);
  assert (f->getSideSign (list_first->point ()) <= 0);

  curr = list_first;
  while (curr != list_last)
    {
      if (f->getSideSign (curr->point ()) > 0)
        {
          last = curr;
          temp = curr->next ();

          while (temp != NULL)
            {
              if (f->getSideSign (temp->point ()) > 0)
                last = temp;

              //printf( "iteraiton\n" ); fflush( stdout );
              temp = temp->next ();
            }

          //SPathNode * ump;

          //ump = curr->p_prev;
          shortcutify_by_face (f, curr, last);
          //curr = ump->next();
          return;
        }
      curr = curr->next ();
    }
  //printf( "D O N E !!\n" ); fflush( stdout );
}


void
SPath::optimize ()
{
  SPathNode *curr;

  curr = list_first;
  while ( curr != NULL
          && curr->next () != NULL
          && curr->next ()->next () != NULL )
    {

      if (curr->next ()->point () == curr->next ()->next ()->point ())
        {
          curr->p_next = curr->next ()->next ();
          curr->next ()->p_prev = curr;
          continue;
        }
      if (curr->point () == curr->next ()->point ())
        {
          curr->nnloc = curr->next ()->nnloc;
          curr->p_next = curr->next ()->next ();
          curr->next ()->p_prev = curr;
          continue;
        }

      curr = curr->next ();
    }
}



void
SPath::shortcut_brute_force ()
{
  //dump();
  //printf( "gogn\n" );
  //fflush( stdout );

  optimize ();

  for (int ind = 0; ind < p_poly->NumFaces (); ind++)
    shortcut_by_face (p_poly->Faces ()[ind]);

  optimize ();
}




/* spath.C - End of File ------------------------------------------*/
