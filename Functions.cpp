/*****************************************************************************
 * Functions.cc
 *
 * Contains various useful functions.  Mostly related to math and vectors
 *****************************************************************************/

#include  <stdio.h>
#include  <memory.h>
#include  <assert.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "myCGAL.h"

#include  "sarray.h"

#include  "Functions.h"
#include  "Path.h"
#include  "FEVBase.h"
#include  "Face.h"
#include  "global-settings.h"

// used to determine if a number "is" zero
#define DELTA 0.000001

using namespace std;

Point
approx_pnt (const Point & p)
{
  return Point (RT (p.x ().to_double ()),
                RT (p.y ().to_double ()), RT (p.z ().to_double ()));
}

Vec
approx_vec (const Vec & v)
{
  return Vec (RT (v.x ().to_double ()),
              RT (v.y ().to_double ()), RT (v.z ().to_double ()));
}


Vec
ApproxUnit (const Vec & v)
{
  return Unit (v);
  
}



Vec
CrossProduct (const Vec & u, const Vec & v)
{
  return CGAL::cross_product (u, v);
 
}



Vec
Unit (const Vec & _v)
{
  RT shrink = (abs (_v[0]) + abs (_v[1]) + abs (_v[2])) / RT (6);
  Vec v (approx_vec (Vec (_v[0] / shrink, _v[1] / shrink, _v[2] / shrink)));
  double length =
    sqrt ((v.hx () * v.hx () + v.hy () * v.hy () +
           v.hz () * v.hz ()).to_double ());
  RT l (length);
  IT nl = l.numerator ();
  IT dl = l.denominator ();
  if (nl == 0)
    return Vec (0, 0, 0, 1);
  //Vec uv = Vec (dl * v.hx (), dl * v.hy (), dl * v.hz (), nl);
  return Vec (dl * v.hx (), dl * v.hy (), dl * v.hz (), nl);

  
}



// simplify the homogenous coordinates of the vector
void
Simplify (Vec & v)
{
  RT a = v.homogeneous (0);
  RT b = v.homogeneous (1);
  RT c = v.homogeneous (2);
  RT d = v.homogeneous (3);
  IT an = a.numerator (), ad = a.denominator ();
  IT bn = b.numerator (), bd = b.denominator ();
  IT cn = c.numerator (), cd = a.denominator ();
  IT dn = d.numerator (), dd = a.denominator ();
  IT gn = gcd (an, gcd (bn, gcd (cn, dn)));
  IT gd = gcd (ad, gcd (bd, gcd (cd, dd)));
  if (gn != IT (1))
    {
      an /= gn;
      bn /= gn;
      cn /= gn;
      dn /= gn;
    }
  if (gd != IT (1))
    {
      ad /= gd;
      bd /= gd;
      cd /= gd;
      dd /= gd;
    }

  v = Vec (RT (an, ad), RT (bn, bd), RT (cn, cd), RT (dn, dd));
  
}


// Simplify the vector so it points in the same direciton
void
SimplifyCrossProduct (Vec & v)
{
  RT a = v.homogeneous (0);
  RT b = v.homogeneous (1);
  RT c = v.homogeneous (2);
  RT d = v.homogeneous (3);
  IT an = a.numerator (), ad = a.denominator ();
  IT bn = b.numerator (), bd = b.denominator ();
  IT cn = c.numerator (), cd = a.denominator ();
  IT dn = d.numerator (), dd = a.denominator ();
  IT gn = gcd (an, gcd (bn, cn));
  IT gd = gcd (ad, gcd (bd, cd));
  if (gn != IT (1))
    {
      an /= gn;
      bn /= gn;
      cn /= gn;
    }
  if (gd != IT (1))
    {
      ad /= gd;
      bd /= gd;
      cd /= gd;
    }

  if (d > RT (0))
    v = Vec (RT (an, ad), RT (bn, bd), RT (cn, cd), RT (1));
  else
    v = Vec (RT (an, ad), RT (bn, bd), RT (cn, cd), RT (-1));

  
}

RT
DotProduct (const Vec & u, const Vec & v)
{
  return u * v;
  /*
     return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
   */
}



// calculates the theta parameter of the spherical coordinate representation
// of x0
double
CalcTheta (const Point & x0, const RT & r)
{
  double sinphi = sqrt (1.0 - ((x0.z () * x0.z ()) / (r * r)).to_double ());
  double x0theta =
    std::acos ((x0.x () / (r.to_double () * sinphi)).to_double ());
  // this is a hack, why is acos(-1) == 0??? (got that result once)
  if ((x0theta == 0) && (x0.x () < 0))
    x0theta = CGAL_PI;
  if (x0.y () < 0)
    x0theta = 2 * CGAL_PI - x0theta;
  return x0theta;

  
}

// calculates the phi parameter of the spherical coordinate representation
// of x0
double
CalcPhi (const Point & x0, const RT & r)
{
  return std::acos ((x0.z () / r).to_double ());
  /*
     return acos((x0.zcoord()/r).to_double());
   */
}

// projects V onto sphere centered at origin with radius r in direction of
// N
//Point ProjectPointOntoSphere(const Point & V, const Vec & N, const RT r)
Point
ProjectPointOntoSphere (const Point & V, const Vec & N, const RT r)
{
  Vec v (CGAL::ORIGIN, V);
  Vec n (Unit (N));

  RT A = DotProduct (v, n);
  RT d = -A +
    RT (sqrt
        ((r * r +
          v.squared_length () * (A * A / v.squared_length () -
                                 1)).to_double ()));
  return V + d * n;

  
}


class AdjGraph
{
private:
  int iiBound, jjBound;

public:
    AdjGraph (int _iiBound, int _jjBound)
  {
    iiBound = _iiBound;
    jjBound = _jjBound;
  }
};


// creates a new node to add to the path graph and adds the initial adjacency
// information
GraphNode *
CreateGraphNode (GridGraph & grid,
                 Point & p, NNLocation & nnloc,
                 int i, int j, Point & pnt_on_sphere)
{
  printf ("CreateGraphNode( %d, %d)\n", i, j);
  //int deg=8;
  //bool bNoTop=false;
  //bool bNoBottom=false;
  //int x=0;
  GraphNode *gn;

  // i == -1      - north pole
  // i == 0       - just below the north pole
  // i == iiBound - south Pole

  gn = new GraphNode (p, nnloc, pnt_on_sphere);
  if (i == -1)
    {
      for (int si = 0; si < grid.getThetaBound (); si++)
        gn->addAdjacency (-1, grid.calcIndex (0, si), true);
      return gn;
    }

  // is It the south pole ?
  if (i == grid.getPhiBound ())
    {
      for (int si = 0; si < grid.getThetaBound (); si++)
        gn->addAdjacency (999999, grid.calcIndex (i - 1, si), true);
      return gn;
    }
  int tmp = grid.calcIndex (i, j);

  gn->addAdjacency (tmp, grid.calcIndex (i - 1, j - 1), false);
  gn->addAdjacency (tmp, grid.calcIndex (i - 1, j), false);
  gn->addAdjacency (tmp, grid.calcIndex (i - 1, j + 1), false);
  gn->addAdjacency (tmp, grid.calcIndex (i, j - 1), false);
  gn->addAdjacency (tmp, grid.calcIndex (i, j + 1), false);
  gn->addAdjacency (tmp, grid.calcIndex (i + 1, j - 1), false);
  gn->addAdjacency (tmp, grid.calcIndex (i + 1, j), false);
  gn->addAdjacency (tmp, grid.calcIndex (i + 1, j + 1), false);

  for (int ind = 0; ind < nnloc.m_loc->FaceDegree (); ind++)
    if (nnloc.m_loc->Faces ()[ind]->getSideSign (p) != 0)
      {
        nnloc.dump ();
        printf ("\npointers: %p %p\n",
                nnloc.m_loc, nnloc.m_loc->Faces ()[ind]);
        printf ("what?\n");
        assert (false);
      }

  return gn;

#ifdef  OLD_JUNK


  // if we're just below the north pole, we won't have as many connections
  if (i == 0)
    {
      deg -= 2;
      bNoTop = true;
    }

  // if we're just above the south pole, we won't have as many connections
  if (i + 1 == iiBound)
    {
      deg -= 2;
      bNoBottom = true;
    }

  // handle the north pole
  if (i == -1)
    {
      deg = jjBound;

      //      for (x=0; x<deg; x++)
      for (int si = 0; si < gn->degree (); si++)
        {
          //printf( "si; %d\n", si );
          //fflush( stdout );
          gn->neighbors[si].m_adj = x * iiBound;
        }

      //printf( "---return gn\n" ); fflush( stdout );
      return gn;
    }

  // handle the south pole
  if (i == iiBound)
    {
      deg = jjBound;
      gn = new GraphNode (p, nnloc, deg);

      //      for (x=0; x<deg; x++)
      for (int si = 0; si < gn->degree (); si++)
        gn->neighbors[si].m_adj = i - 1 + x * iiBound;

      return gn;
    }

  gn = new GraphNode (p, nnloc, deg);
  int si = 0;
  //gn->m_adjDist.begin();

  gn->neighbors[si].m_adj = i + ((j == 0) ? jjBound - 1 : j - 1) * iiBound;
  // left
  si++;
  gn->neighbors[si].m_adj = i + ((j + 1 == jjBound) ? 0 : j + 1) * iiBound;
  // right
  si++;

  if (!bNoTop)
    {
      gn->neighbors[si].m_adj = i - 1 + j * iiBound;  // up
      si++;
      gn->neighbors[si].m_adj = i - 1 + ((j == 0) ? jjBound - 1 : j - 1) * iiBound; // up-left
      si++;
      gn->neighbors[si].m_adj = i - 1 + ((j + 1 == jjBound) ? 0 : j + 1) * iiBound; // up-right
      si++;
    }
  else
    {
      gn->neighbors[si].m_adj = jjBound * iiBound + 1;  // north_pole
      si++;
    }

  if (!bNoBottom)
    {
      gn->neighbors[si].m_adj = i + 1 + j * iiBound;  // down
      si++;
      gn->neighbors[si].m_adj = i + 1 + ((j == 0) ? jjBound - 1 : j - 1) * iiBound; // down-left
      si++;
      gn->neighbors[si].m_adj = i + 1 + ((j + 1 == jjBound) ? 0 : j + 1) * iiBound; // down-right
      si++;
    }
  else
    {
      gn->neighbors[si].m_adj = jjBound * iiBound;  // south_pole
      si++;
    }

  //printf( "return gn\n" ); fflush( stdout );

  return gn;
#endif // OLD_JUNK
}

// PointOn was written by Neill Occhiogrosso but the fact that there could
// be division by zero was never taken into account.  This function is
// used to permute the parameters before calculating an intersection point
// so that we will never divide by zero
void
PermuteForIntersect (RT & a, RT & b, RT & c, int p)
{
  RT i, j;

  switch (p)
    {
    case 1:
      {
        break;
      }
    case 2:
      {
        i = b;
        b = c;
        c = i;
        break;
      }
    case 3:
      {
        i = a;
        j = b;
        a = c;
        b = i;
        c = j;
        break;
      }
    case 4:
      {
        i = a;
        a = b;
        b = i;
        break;
      }
    case 5:
      {
        i = a;
        j = b;
        b = c;
        a = j;
        c = i;
        break;
      }
    case 6:
      {
        i = a;
        a = c;
        c = i;
        break;
      }
    }
}

// this code was written by Neill Ochiogrosso
Point
PointOn (const Plane & hp, const Plane & hq)
{
  // permutations are numbered as follows:
  //   1: A, B, C
  //   2: A, C, B
  //   3: C, A, B
  //   4: B, A, C
  //   5: B, C, A
  //   6: C, B, A
  int perm = 0;

  Point p = hp.point ();
  Point q = hq.point ();

  Vec a = normalVector (hp);
  Vec b = normalVector (hq);

  RT A0 = a[0];
  RT B0 = a[1];
  RT C0 = a[2];

  RT A1 = b[0];
  RT B1 = b[1];
  RT C1 = b[2];

  RT x0 = p[0];
  RT y0 = p[1];
  RT z0 = p[2];

  RT x1 = q[0];
  RT y1 = q[1];
  RT z1 = q[2];

  RT x, y;
  RT i, j;

  if (A1 == 0)
    {
      if (B1 == 0)
        {
          if (C1 * B0 - C0 * B1 == 0)
            {
              // C, A, B
              perm = 3;
            }
          else
            {
              // C, B, A
              perm = 6;
            }
        }
      else
        {
          if (B1 * A0 - B0 * A1 == 0)
            {
              // B, C, A
              perm = 5;
            }
          else
            {
              // B, A, C
              perm = 4;
            }
        }
    }
  else
    {
      if (A1 * B0 - A0 * B1 == 0)
        {
          // A, C, B
          perm = 2;
        }
      else
        {
          // A, B, C
          perm = 1;
        }
    }

  // see above description for PermuteForIntersect
  PermuteForIntersect (A0, B0, C0, perm);
  PermuteForIntersect (A1, B1, C1, perm);
  PermuteForIntersect (x0, y0, z0, perm);
  PermuteForIntersect (x1, y1, z1, perm);


  y =
    (-A0 * (C1 * z1 + B1 * y1 + A1 * (x1 - x0)) +
     A1 * (C0 * z0 + B0 * y0)) / (A1 * B0 - A0 * B1);
  x = x1 + C1 * z1 / A1 + (y1 - y) * B1 / A1;

  IT n1 = x.numerator ();
  IT d1 = x.denominator ();

  IT n2 = y.numerator ();
  IT d2 = y.denominator ();

  switch (perm)
    {
    case 1:
      return Point (n1 * d2, n2 * d1, 0, d1 * d2);
    case 2:
      return Point (n1 * d2, 0, n2 * d1, d1 * d2);
    case 3:
      return Point (0, n1 * d2, n2 * d1, d1 * d2);
    case 4:
      return Point (n2 * d1, n1 * d2, 0, d1 * d2);
    case 5:
      return Point (n2 * d1, 0, n1 * d2, d1 * d2);
    case 6:
      return Point (0, n2 * d1, n1 * d2, d1 * d2);
    default:
      assert (false);
      return Point (0, 0, 0, 0);
    }
}

// rotates a around r by angle radians
Vec
Rotate (const Vec & a, const Vec & r, double angle)
{
  RT c_angle (cos (angle));
  RT s_angle (sin (angle));

  Vec r1 =
    c_angle * r +
    s_angle * RT (sqrt (r.squared_length ().to_double ())) *
    (Unit (CrossProduct (a, r))) + DotProduct (a, r) * (1 - c_angle) * a;
  Simplify (r1);

  return r1;
}

// unfolds p on h1 or h2 onto the other
Point
Unfold (const Point & p, const Plane & h1, const Plane & h2)
{
  // get a point on the intersection of the two planes
  Point O = PointOn (h1, h2);

  // create a vector from the above point to the point to unfold
  Vec r = p - O;
  Vec n1 = normalVector (h1);
  Vec n2 = normalVector (h2);
  // find a vector that lies along the line of intersection of the two planes
  Vec a = Unit (CrossProduct (n1, n2));

  // calculate the angle to rotate by
  double angle = angleBetween (n1, n2);

  // perform the rotation
  Vec r1 = Rotate (a, r, angle);

  // compute the location of the new point
  Point q = O + r1;

  Vec v = Vec (CGAL::ORIGIN, q);
  Simplify (v);
  q = Point (v.x (), v.y (), v.z ());


  return q;
}

// returns the intersection of a line and a plane
Point
IntersectLinePlane (const Point & q1, const Point & q2,
                    const Plane & P, bool & contains)
{
  Line l (q1, q2);
  contains = false;
  CGAL::Object result = CGAL::intersection (l, P);
  if (const Point * ip = CGAL::object_cast < Point > (&result))
    {
      // intersect at unique point
      return *ip;
    }
  else if (CGAL::object_cast < Line > (&result))
    {
      // line is on plane
      contains = true;
      return q1;
    }
  else
    {
      // line is parallel with plane
      contains = true;
      return Point (0, 0, 0);
    }

  
}


// returns the ridgepoint of the hershberger-suri 2-plane wedge
Point
RidgePoint (const Point & p1, const Point & p2,
            const Plane & h1, const Plane & h2)
{
  bool temp;                    // for Inter, but value not used here
  Point q = Unfold (p1, h1, h2);
  Point rp = IntersectLinePlane (p2, q, h1, temp);

  //assert( h1.contains( rp ) );
  //assert( h2.contains( rp ) );

  return rp;
}

int
my_side_of (const Plane & pl, const Point & pnt)
{
  int sgn;

  sgn = sign (DotProduct (pnt - pl.point (), normalVector (pl)));
  return sgn;
}

void
flip (Plane & pl)
{
  //pl = Plane( pl.point2(), pl.point1(), pl.point3() );
  pl = pl.opposite ();
}


void vec_print (const Vec & v);

Point
ExtUnfold (const Point & center,
           Plane target_plane, Plane src_plane, const Point & pnt)
{
  assert (src_plane.has_on (pnt));

  Vec src_norm (normalVector (src_plane));
  Vec target_norm (normalVector (target_plane));

  if (target_plane.has_on (pnt))
    return pnt;

  Point O = PointOn (target_plane, src_plane);

  assert (src_plane.has_on (O));
  assert (target_plane.has_on (O));

  if (DotProduct (src_norm, center - O) >= 0)
    src_norm = -src_norm;
  if (DotProduct (target_norm, center - O) >= 0)
    target_norm = -target_norm;

  if (DotProduct (src_norm, center - O) >= 0)
    assert (false);
  if (DotProduct (target_norm, center - O) >= 0)
    assert (false);

  src_norm = Unit (src_norm);
  target_norm = Unit (target_norm);
  /*
     printf( "normalized: " );
     vec_print( target_norm );
     printf( "\n" );
   */

  Vec r = pnt - O;

  Vec pre_vec_along_line (CrossProduct (target_norm, src_norm));
  Vec vec_along_line (Unit (pre_vec_along_line));
  Vec vec_on_trg = Unit (CrossProduct (target_norm,
                                       pre_vec_along_line));
  Vec vec_on_src = Unit (CrossProduct (src_norm,
                                       pre_vec_along_line));
  RT c_along_line, c_on_src;
  RT c_on_src_norm;

  c_along_line = DotProduct (vec_along_line, r);
  c_on_src = DotProduct (vec_on_src, r);
  c_on_src_norm = DotProduct (src_norm, r);

  /* Hoai - 23/08/2012 */
  assert (fabs (c_on_src_norm.to_double ()) <= 0.01);

  //Point  new_p2 = approx_pnt( c_along_line * vec_along_line
  //    + c_on_src * vec_on_trg + Vec( CGAL::ORIGIN, O ) );
  Point new_p2 = approx_pnt (O + c_along_line * vec_along_line
                             + c_on_src * vec_on_trg);

  //Vec v = target_plane.normal_project( new_p2 );
  Vec v = (target_plane.perpendicular_line (new_p2)).to_vector ();
  //cout << "v : " << v << "\n";
  if (sqrt (v.squared_length ().to_double ()) > 0.01)
    {
      printf ("sqrt distance: %g\n", sqrt (v.squared_length ().to_double ()));
      printf ("\n");
      point_print (new_p2);
      printf ("\n");
      cout << c_along_line << "\n";
      cout << c_on_src << "\n";
      cout << c_on_src_norm << "\n";
      cout << c_along_line.to_double () << "\n";
      cout << c_on_src.to_double () << "\n";
      cout << c_on_src_norm.to_double () << "\n";
      printf ("\n");
      printf ("\n");
      vec_print (src_norm);
      printf ("\n");
      vec_print (target_norm);
      printf ("\n");
      cout << "SRC NORM:" << src_norm;
      printf ("\n");
      cout << target_norm;
      printf ("\n");
      printf ("Origin:");
      point_print (O);

      cout << v << "\n";
      fflush (stdout);
      assert (sqrt (v.squared_length ().to_double ()) <= 0.01);
    }

  Point image_on_trg (new_p2 + v);

  return image_on_trg;
}


Point
ExactExtRidgePoint (const Point & _p1, const Plane h1,
                    const Vec & _h1_norm,
                    const Point & _p2, const Plane h2,
                    const Vec & _h2_norm, const Point & center)
{
    // trang: avoiding unused varibale
    Point tmp(CGAL::ORIGIN); tmp = center;

    //int  h2_sign_p1, h1_sign_p2;
    Vec h1_norm (_h1_norm);
    Vec h2_norm (_h2_norm);


    //if ((!h1.has_on (_p1) != 0) || (!h2.has_on (_p2) != 0))
    if ( !h1.has_on (_p1) || !h2.has_on (_p2) )
      {
        assert (h1.has_on (_p1));
        assert (h2.has_on (_p2));
      }

    //h2_sign_p1 = h2.side_of( _p1 );
    //if  ( h2_sign_p1 == 0 ) {
    if (h2.has_on (_p1))
      {
        // h2 contains both p1 and p2.  and p1 is on the intersecting
        // line - it is the required ridgem point,
        //printf( "NNAN?\n" );
        return _p1;
      }

    //h1_sign_p2 = h1.side_of( _p2 );
    //if  ( h1_sign_p2 == 0 ) {
    if (h1.has_on (_p2))
      {
        //printf( "NNAN2?\n" );
        // h1 contains both p1 and p2.  and p2 is on the intersecting
        // line - it is the required ridge point,
        return _p2;
      }

    const Point & p1 (_p1);
    const Point & p2 (_p2);

    Point O = PointOn (h1, h2);
    assert (h1.has_on (O));
    assert (h2.has_on (O));

    /* Hoai - 23/08/2012 - redundant */
    assert (h2.has_on (p2));
    assert (h1.has_on (p1));

    Vec r = p2 - O;
 
    if (DotProduct (h2_norm, p1 - O) >= 0)
      h2_norm = -h2_norm;
    if (DotProduct (h1_norm, p2 - O) >= 0)
      h1_norm = -h1_norm;

    Vec pre_vec_along_line (CrossProduct (h1_norm, h2_norm));
    assert (h1.has_on (O + pre_vec_along_line));
    assert (h2.has_on (O + pre_vec_along_line));

    Vec vec_along_line (Unit (pre_vec_along_line));
    Vec vec_on_h2 = Unit (CrossProduct (h2_norm, pre_vec_along_line));
    Vec vec_on_h1 = Unit (CrossProduct (h1_norm, pre_vec_along_line));
  
    RT c_along_line, c_on_h2;
    RT c_on_h2_norm;

    c_along_line = DotProduct (vec_along_line, r);
    c_on_h2 = DotProduct (vec_on_h2, r);
    c_on_h2_norm = DotProduct (h2_norm, r);

 
    Point new_p2 = O + c_along_line * vec_along_line + c_on_h2 * vec_on_h1;
 
    Point old_p2 = O + c_along_line * vec_along_line + c_on_h2 * vec_on_h2;

    Point image_p2_on_h1 = h1.projection( new_p2 );

    assert( CGAL::squared_distance( image_p2_on_h1, new_p2 ).to_double() <= .01 );

    CGAL::Oriented_side a, b;
    a = h2.oriented_side (image_p2_on_h1);
    b = h2.oriented_side (p1);
    //cout << "a = " << a << ", b = " << b << endl;
    if (a == CGAL::ON_ORIENTED_BOUNDARY)
      a = -b;
    if (b == CGAL::ON_ORIENTED_BOUNDARY)
      b = -a;
 
    assert (a == -b);

    bool temp;

    Point rp = IntersectLinePlane (p1, image_p2_on_h1, h2, temp);
    return rp;
 
}




Point
ExtRidgePoint (const Point & _p1, const Plane h1,
               Vec h1_norm,
               const Point & _p2, const Plane h2,
               Vec h2_norm, const Point & center)
{
  //int  h2_sign_p1, h1_sign_p2;

  //printf( "ExtRigePoint called!\n" );

  if ((!h1.has_on (_p1)) || (!h2.has_on (_p2)))
    {
      assert (h1.has_on (_p1));
      assert (h2.has_on (_p2));
    }

  //h2_sign_p1 = h2.side_of( _p1 );
  //if  ( h2_sign_p1 == 0 ) {
  if (h2.has_on (_p1))
    {
      // h2 contains both p1 and p2.  and p1 is on the intersecting
      // line - it is the required ridge point,
      //printf( "NNAN?\n" );
      return _p1;
    }

  //h1_sign_p2 = h1.side_of( _p2 );
  //if  ( h1_sign_p2 == 0 ) {
  if (h1.has_on (_p2))
    {
      //printf( "NNAN2?\n" );
      // h1 contains both p1 and p2.  and p2 is on the intersecting
      // line - it is the required ridge point,
      return _p2;
    }

  Point p1 (approx_pnt (_p1));
  Point p2 (approx_pnt (_p2));

  Point O = approx_pnt (PointOn (h1, h2));
  //Vec  h2_norm = Unit( h2.normal() );
  //Vec  h1_norm = Unit( h1.normal() );

  if (DotProduct (h2_norm, p1 - O) >= 0)
    h2_norm = -h2_norm;
  if (DotProduct (h1_norm, p2 - O) >= 0)
    h1_norm = -h1_norm;

  Vec r = p2 - O;
  Vec vec_along_line (ApproxUnit (CrossProduct (h1_norm, h2_norm)));
  Vec vec_on_h2 = ApproxUnit (CrossProduct (h2_norm, vec_along_line));
  Vec vec_on_h1 = ApproxUnit (CrossProduct (h1_norm, vec_along_line));
 
  RT c_along_line, c_on_h2;

  c_along_line = DotProduct (vec_along_line, r);
  c_on_h2 = DotProduct (vec_on_h2, r);

  Point new_r = O + c_along_line * vec_along_line + c_on_h2 * vec_on_h1;

  Vec v = (h1.perpendicular_line (new_r)).to_vector ();

  Point image_p2_on_h1 (new_r + v);

  //if  ( h2.side_of( image_p2_on_h1 ) != - h2.side_of( p1 ) ) {
  if (h2.oriented_side (image_p2_on_h1) != -h2.oriented_side (p1))
    {
      //printf( "calling ExactExtRidgePoint!\n" ); fflush( stdout );
      Point out (ExactExtRidgePoint (_p1, h1, h1_norm, _p2, h2, h2_norm,
                                     center));
      //printf( "     done ExactExtRidgePoint!\n" ); fflush( stdout );
      return out;
    }


  //assert( h2.side_of( image_p2_on_h1 ) == - h2.side_of( p1 ) );

  bool temp;

  Point rp = IntersectLinePlane (p1, image_p2_on_h1, h2, temp);

  Point a_center (approx_pnt (center));
  Vec vv ((approx_pnt (rp) - a_center));

  Point new_rp (approx_pnt (a_center + RT (1.0000001) * vv));

  return new_rp;
}


// returns true if the square distance between two points at most 1/100
bool
PointEq (const Point & p1, const Point & p2)
{
  return (CGAL::squared_distance (p1, p2) <= RT (1, 100));
}


double
GridGraph::vertices_SGP_distance (int a, int b) const
{

        return SGP_distance (arr[a]->m_p, arr[b]->m_p);
}

double
GridGraph::vertices_distance (int a, int b) const
{

    return sqrt (CGAL::squared_distance (arr[a]->m_p, arr[b]->m_p).to_double ());
    
}



void GridGraph::set_distance (int to, int j, double dist)
{
  GraphNode & gn (*(arr[to]));

  for (int sk = 0; sk < (int) gn.neighbors.size (); sk++)
    {
      if (gn.neighbors[sk].m_adj == j)
        {
          gn.neighbors[sk].m_dist = dist;
          break;
        }
    }
}



void
GridGraph::addNodeToGraph (GraphNode * gn, int idx)
{
  if (idx != entriesNum)
    {
      printf ("idx: %d  limit: %d  entriesNum: %d\n", idx, limit, entriesNum);

      fflush (stdout);

      assert (idx == entriesNum);
    }

  if (idx >= limit)
    {
      int new_limit, new_size;
      GraphNode **new_arr;

      new_limit = limit * 2 + 1;
      new_size = new_limit * sizeof (GraphNode *);

      new_arr = (GraphNode **) malloc (new_size);
      assert (new_arr != NULL);
      memset (new_arr, 0, new_size);

      memcpy (new_arr, arr, limit * sizeof (GraphNode *));
      free (arr);
      arr = new_arr;
      limit = new_limit;
    }

  //printf( "addNodeToGraph - idx: %d\n", idx );
  arr[idx] = gn;
  entriesNum++;
}

void
GridGraph::connectByIndices (int src, int dst)
{
  assert ((0 <= src) && (src < entriesNum));
  assert ((0 <= dst) && (dst < entriesNum));

  //printf( "UU %d %d\n", src, dst );
  //assert( false );
  arr[src]->addAdjacency (src, dst, false);
  arr[dst]->addAdjacency (dst, src, false);
}


void
GridGraph::connect (int i, int j, int ind_new)
{
  int tmp;

  tmp = calcIndex (i, j);
  if (tmp == GG_INVALID)
    return;
  assert ((0 <= tmp) && (tmp < entriesNum));
  assert ((0 <= ind_new) && (ind_new < entriesNum));

  //printf( "VVL %d %d\n", tmp, ind_new );
  arr[tmp]->addAdjacency (tmp, ind_new, false);
  arr[ind_new]->addAdjacency (ind_new, tmp, false);
}


int
GridGraph::AppendVertex (GraphNode * gn)
{
  addNodeToGraph (gn, entriesNum);

  return entriesNum - 1;
}


void
GridGraph::PrintGraph ()
{
  //NodeSet::iterator si;
  for (int j = 0; j < entriesNum; j++)
    {
      cerr << "[" << j << "]: ";
      //      for (int i=0; i<m_graph[j]->m_deg; i++)
      for (int si = 0; si < (int) arr[j]->neighbors.size (); si++)
        {

          cerr << arr[j]->neighbors[si].m_adj << " ";
        }
      cerr << endl;
    }
}


// for clasical_distance
void
GridGraph::CalcDists ()
{
  printf( "GridGraph::CalcDists\n" );

  //printf( "( %d == %d )\n", entriesNum, getNorthPoleIndex() );
  //fflush( stdout );
  if (entriesNum == getNorthPoleIndex ())
    addNodeToGraph (pNorthPole, getNorthPoleIndex ());

  //    NodeSet::iterator si;

  for (int j = 0; j < entriesNum; j++)
    {
      debug ("j=" << j);

      for (int si = 0; si < (int) arr[j]->neighbors.size (); si++)
        {
          if (arr[j]->neighbors[si].m_dist == -1)
            {
              int to = arr[j]->neighbors[si].m_adj;

              //printf( "j: %d         to: %d\n", j, to );
              //fflush( stdout );

              if ((arr[j] == NULL) || (arr[to] == NULL))
                {
                  printf ("j = %d, to = %d\n", j, to);
                  printf ("entriesNum : %d\n", entriesNum);
                  printf ("north pole: %d  south pole: %d\n",
                          getNorthPoleIndex (), getSouthPoleIndex ());
                  fflush (stdout);

                  assert (arr[j] != NULL);
                  assert (arr[to] != NULL);
                }

              arr[j]->neighbors[si].m_dist = vertices_distance (j, to);

              set_distance (to, j, arr[j]->neighbors[si].m_dist);
            }
        }
    }
  printf( "after calc dist!\n" );
}

// for SGP_distance
void
GridGraph::CalcSGPDists ()
{
  //printf( "GridGraph::CalcDists\n" );

  //printf( "( %d == %d )\n", entriesNum, getNorthPoleIndex() );
  //fflush( stdout );
  if (entriesNum == getNorthPoleIndex ())
    addNodeToGraph (pNorthPole, getNorthPoleIndex ());

  //    NodeSet::iterator si;

  for (int j = 0; j < entriesNum; j++)
    {
      debug ("j=" << j);

      for (int si = 0; si < (int) arr[j]->neighbors.size (); si++)
        {
          if (arr[j]->neighbors[si].m_dist == -1)
            {
              int to = arr[j]->neighbors[si].m_adj;

              //printf( "j: %d         to: %d\n", j, to );
              //fflush( stdout );

              if ((arr[j] == NULL) || (arr[to] == NULL))
                {
                  printf ("j = %d, to = %d\n", j, to);
                  printf ("entriesNum : %d\n", entriesNum);
                  printf ("north pole: %d  south pole: %d\n",
                          getNorthPoleIndex (), getSouthPoleIndex ());
                  fflush (stdout);

                  assert (arr[j] != NULL);
                  assert (arr[to] != NULL);
                }

              arr[j]->neighbors[si].m_dist = vertices_SGP_distance (j, to);

              set_distance (to, j, arr[j]->neighbors[si].m_dist);
            }
        }
    }
  //printf( "after calc dist!\n" );
}



struct PQ_node
{
PQ_node (double len = 0.0, int ind = -1):m_len (len), m_ind (ind)
  {
  }
  double m_len;
  int m_ind;
};

class Compare_PQ_node
{
public:
  bool operator  () (PQ_node & n1, PQ_node & n2)
  {
    return (n1.m_len > n2.m_len);
  }
};

void
pushHeap (const PQ_node & node, vector < PQ_node > &H)
{
  H.push_back (node);
  push_heap (H.begin (), H.end (), Compare_PQ_node ());
}

void
popHeap (vector < PQ_node > &H)
{
  pop_heap (H.begin (), H.end (), Compare_PQ_node ());
  H.pop_back();
}




/*----------------------------- IMPORTANCE FUNC ---------------------------------
// --return the original version find shortest paths in graph for clasical distance
-------------------------------------------------------------------------------*/
GraphNode * GridGraph::FindShortestPathInGraph (const int src, const int trg)
{
  PQ_node curr_node;
  vector < PQ_node > PQ;
  double *dist_arr = new double[entriesNum];
  bool *visited = new bool[entriesNum];
  int j;
  GraphNode *node;
  double new_dist;
  FEVBase *dest = arr[trg]->m_nnloc.m_loc;

  // clear out array
  for (int ind = 0; ind < entriesNum; ind++)
    {
      dist_arr[ind] = -1;
      visited[ind] = false;
      arr[ind]->m_prev = NULL;
      arr[ind]->m_next = NULL;
      arr[ind]->m_nnloc.m_loc->SetLastInPath (NULL);
    }

  // init the priority queue with source node
  pushHeap (PQ_node (0.0, src), PQ);
  dist_arr[src] = 0.0;

  //cout << "Source: " << src << ", destionation: " << trg << endl;
  // main loop
  while (!PQ.empty ())
  {
    // pop the first (smallest) elements
    curr_node = PQ.front ();
    //cout << "Before queue size: " << PQ.size() << endl;
    popHeap (PQ);

    node = arr[curr_node.m_ind];
    visited[curr_node.m_ind] = true;

    // if using the current point makes you land on the same face as the
    // destination point, "add" that edge to the graph
    if ( dest->Contains (node->m_nnloc.m_loc)
         //-----------------------------------
         // Hoai - 16/3/2013
         //&& node->m_p.z() >= arr[trg]->m_p.z()
         //-----------------------------------
       )
      {
        new_dist = dist_arr[curr_node.m_ind] +
                 sqrt (CGAL::squared_distance  (arr[trg]->m_p, node->m_p).to_double ());

         if (dist_arr[trg] < 0 || new_dist < dist_arr[trg])
          {
            pushHeap ( PQ_node(new_dist, trg), PQ);
            dist_arr[trg] = new_dist;
            arr[trg]->m_prev = node;
          }
      }

    for (int si = 0; si < (int) node->neighbors.size (); si++)
    {
        AdjDist & info (node->neighbors[si]);
        new_dist = dist_arr[curr_node.m_ind] + info.m_dist;
        j = info.m_adj;
        if (visited[j])
          continue;

        //----------------------------------------------
        if (dist_arr[j] < 0 || new_dist < dist_arr[j])
        {
            pushHeap (PQ_node (new_dist, j), PQ);
            dist_arr[j] = new_dist;
            arr[j]->m_prev = node;
        }
    }
  }

  cerr << "Distance after Dijkstra: " << dist_arr[trg] << endl;

  // release memory
  delete[]dist_arr;
  delete[]visited;

  return arr[trg];
}

// --find shortest paths in graph for SGP distance

GraphNode * GridGraph::FindShortestGentlePathInGraph (const int src, const int trg)
{
  PQ_node curr_node;
  vector < PQ_node > PQ;
  double *dist_arr = new double[entriesNum];
  bool *visited = new bool[entriesNum];
  int j;
  GraphNode *node;
  double new_dist;
  FEVBase *dest = arr[trg]->m_nnloc.m_loc;

  // clear out array
  for (int ind = 0; ind < entriesNum; ind++)
    {
      dist_arr[ind] = -1;
      visited[ind] = false;
      arr[ind]->m_prev = NULL;
      arr[ind]->m_next = NULL;
      arr[ind]->m_nnloc.m_loc->SetLastInPath (NULL);
    }

  // init the priority queue with source node
  pushHeap (PQ_node (0.0, src), PQ);
  dist_arr[src] = 0.0;

  //cout << "Source: " << src << ", destionation: " << trg << endl;
  // main loop
  while (!PQ.empty ())
  {
    // pop the first (smallest) elements
    curr_node = PQ.front ();
    //cout << "Before queue size: " << PQ.size() << endl;
    popHeap (PQ);

    node = arr[curr_node.m_ind];
    visited[curr_node.m_ind] = true;

    // if using the current point makes you land on the same face as the
    // destination point, "add" that edge to the graph
    if ( dest->Contains (node->m_nnloc.m_loc)
        //-----------------------------------
        // Hoai - 16/3/2013
        //&& node->m_p.z() >= arr[trg]->m_p.z()
        //-----------------------------------
       )
      {
        new_dist = dist_arr[curr_node.m_ind] +
          SGP_distance (arr[trg]->m_p, node->m_p); 

        // new_dist = dist_arr[curr_node.m_ind] +
        //             sqrt (CGAL::squared_distance  (arr[trg]->m_p, node->m_p).to_double ());
    

        if (dist_arr[trg] < 0 || new_dist < dist_arr[trg])
          {
            pushHeap ( PQ_node(new_dist, trg), PQ);
            dist_arr[trg] = new_dist;
            arr[trg]->m_prev = node;
          }
      }

    for (int si = 0; si < (int) node->neighbors.size (); si++)
    {
        AdjDist & info (node->neighbors[si]);
        new_dist = dist_arr[curr_node.m_ind] + info.m_dist;
        j = info.m_adj;
        if (visited[j])
          continue;

        
        if ( node->m_p.z() < arr[ j ]->m_p.z() )
        {
            //cout << " %% z of the path is decrseasing, continued ...." << endl;
            continue;
        }
        else // trang's debug
        {
            //cout << " %% z of the path is increasing" << endl;
            //exit(0);
        }
        

        //----------------------------------------------
        if (dist_arr[j] < 0 || new_dist < dist_arr[j])
        {
            pushHeap (PQ_node (new_dist, j), PQ);
            dist_arr[j] = new_dist;
            arr[j]->m_prev = node;
        }
    }
  }

  cerr << "SGP Distance after Dijkstra: " << dist_arr[trg] << endl;

  // release memory
  delete[]dist_arr;
  delete[]visited;

  return arr[trg];
}


double
     GridGraph::getVertex_DoubleXCoord (int i) const
     {
       assert ((0 <= i) && (i < entriesNum));

       //return  (double)arr[ i ]->m_p.xcoordD();
       return arr[i]->m_p.x ().to_double ();
     }

     double GridGraph::getVertex_DoubleYCoord (int i) const
     {
       assert ((0 <= i) && (i < entriesNum));

       //return  (double)arr[ i ]->m_p.ycoordD();
       return arr[i]->m_p.y ().to_double ();
     }

     double GridGraph::getVertex_DoubleZCoord (int i) const
     {
       assert ((0 <= i) && (i < entriesNum));

       //return  (double)arr[ i ]->m_p.zcoordD();
       return arr[i]->m_p.z ().to_double ();
     }

     int GridGraph::calcIndex (int phi_ind, int theta_ind)
{
  int res;

  if (phi_ind < -1)
    res = GG_INVALID;
  else if (phi_ind == -1)
    res = getNorthPoleIndex ();
  else if (phi_ind > getPhiBound ())
    res = GG_INVALID;
  else if (phi_ind == getPhiBound ())
    res = getSouthPoleIndex ();
  else
    {
      if (theta_ind == -1)
        theta_ind = getThetaBound () - 1;
      if (theta_ind >= getThetaBound ())
        theta_ind = theta_ind - getThetaBound ();

      res = phi_ind + theta_ind * getPhiBound ();
    }

  return res;
}

/* Hoai 07/08/2012 */
// angle between two vectors
// pre-condition checking code is not considered yet
double angleBetween (const Vec & v1, const Vec & v2)
{
  if ( v1 * v2 == 0 ){
    // cerr << v1 << endl;
    // cerr << v2 << endl;
    return CGAL_PI / 2.0;
  }

  Vec _v1 = v1 / sqrt ((v1 * v1).to_double ());
  Vec _v2 = v2 / sqrt ((v2 * v2).to_double ());
  double mycos = (_v1 * _v2).to_double ();
  if ( mycos < -1.0 ) mycos = -1.0;
  if ( mycos > 1.0 ) mycos = 1.0;
  return acos( mycos );
}

double angleBetween_2(const Vector_2 & v1, const Vector_2 & v2)
{
  if ( v1 * v2 == 0 ){
    // cerr << v1 << endl;
    // cerr << v2 << endl;
    return CGAL_PI / 2.0;
  }

  Vector_2 _v1 = v1 / sqrt ((v1 * v1).to_double ());
  Vector_2 _v2 = v2 / sqrt ((v2 * v2).to_double ());
  return acos ((_v1 * _v2).to_double ());
}

// normal vector to a plane
Vec
normalVector (const Plane & pl)
{
  Vec ov = pl.orthogonal_vector ();
  return ov / sqrt ((ov * ov).to_double ());
}



/* --------------------- trang: projecting using triangle as a face -----------------------*/
Point myLineIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2 )
{
  Plane f1( p1, l1, l2 );
  //Plane h2( p2, l2, l1 );
  Line l( l1, l2 );

  // project p2 on plane(l1,l2,p1)
  Point p2ronf1;
  if ( f1.has_on( p2 ) )
  {
      p2ronf1 = p2;
  } else
  {
      Projection_3 myproj;
      //
      p2ronf1 = myproj( f1, p2 );
  }

  if ( Seg( p1, p2ronf1 ).is_degenerate() )
  {
      err_point_print( p1 );
      err_point_print( p2ronf1 );
      err_point_print( p2 );
  }
  assert( !Seg( p1, p2ronf1 ).is_degenerate() );
  assert( !l.is_degenerate() );

  CGAL::Object obj = CGAL::intersection( Seg( p1, p2ronf1 ), l );
  //
  if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
  {
      return *pobj;
  }
  else
      return CGAL::ORIGIN;
}


/*-------------------trang: find intersection in which there some faces between-----------------*/
Point myLineIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2,
                       const std::vector< Triangle > TriBetween, unsigned int onFace)
{
    Line l( l1, l2 );

    cout << "\t\t\t- find intersection using sequence of between faces ... "; //PressEnterToContinue();

    Plane shootingface = TriBetween[ onFace ].supporting_plane();

    // 1. unfold face sequences from p1 onto the shooting face
    std::vector< Point > p1ron_shootingface;
    p1ron_shootingface.push_back(p1);

    if ( !shootingface.has_on( p1 ) )
    {
        // trang's debug
        cout << "\t\t\t- unfold face sequences from p1 onto the shooting face" << endl;
        //
        for ( unsigned int i = 1; i <= onFace; i++ )
        {
            Projection_3 myproj;

            Point bar_p1 = myproj( TriBetween[i].supporting_plane(), p1ron_shootingface.back() );
            p1ron_shootingface.push_back( bar_p1 );

            // trang's debug:
            cout << "\t\t\t\t - unfold from " << i-1 << " to " << i << " ... "; //PressEnterToContinue();
        }
    }

    // 2. unfold face sequences from p2 onto the shooting face
    std::vector< Point > p2ron_shootingface;
    p2ron_shootingface.push_back(p2);

    if ( !shootingface.has_on( p2 ) )
    {
        // trang's debug
        cout << "\t\t\t- unfold face sequences from p2 onto the shooting face" << endl;
        //
        for ( unsigned int i = TriBetween.size() - 2; i >= onFace; i-- )
        {
            Projection_3 myproj;

            Point bar_p2 = myproj( TriBetween[i].supporting_plane(), p2ron_shootingface.back() );
            p2ron_shootingface.push_back( bar_p2 );

            // trang's debug:
            cout << "\t\t\t\t - unfold from " << i+1 << " to " << i << " ... "; //PressEnterToContinue();
        }
    }

    // 3. determine next shooting point by intersecting line through [ bar_p1, bar_p2 ] and l
    Seg bar_p1p2(p1ron_shootingface.back(), p2ron_shootingface.back());
    //
    if ( bar_p1p2.is_degenerate() )
    {
        err_point_print( p1 );
        err_point_print( p1ron_shootingface.back() );
        err_point_print( p2ron_shootingface.back() );
        err_point_print( p2 );
    }
    assert( !bar_p1p2.is_degenerate() );
    assert( !l.is_degenerate() );

    // find the next shooting point as the intersectrion ...
    CGAL::Object obj = CGAL::intersection( bar_p1p2, l );
    //
    if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
    {
        return *pobj;
    }
    else
        return CGAL::ORIGIN;
}

/*---------------- trang: find intersection in which there some segments between --------------*/
Point myLineIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2,
                       const std::vector< Seg > SegBetween)
{
    Plane h1( p1, l1, l2 );
    Line l( l1, l2 );

    // trang's debug
    cout << "\t\t\t- find intersection using sequence of between edges ... "; //PressEnterToContinue();

    // project p2 on plane(l1,l2,p1)
    std::vector< Point > p2ronhi;
    //
    p2ronhi.push_back(p2);

    // unfolfing p2 into plane h1: (p1, l1, l2)
    if ( !h1.has_on( p2 ) )
    {
        for ( unsigned int i = 0; i < SegBetween.size(); i++ )
        {
            Line li(SegBetween[i]);

            Projection_3 myproj;
            //Vec n1 = h1.orthogonal_vector();
            //Vec n2 = h2.orthogonal_vector();

            Point p2onl = myproj( li, p2ronhi.back() );
            Point p1onl = myproj( li, p1 );

            double d1 = sqrt( ( p1 - p1onl ).squared_length().to_double() );
            double d2 = sqrt( ( p2 - p2onl ).squared_length().to_double() );

            if ( d1 == 0 )
              return p1;

            p2ronhi.push_back( p2onl + ( p1onl - p1 ) * d2 / d1 );
        }
    }

    //
    if ( Seg( p1, p2ronhi.back() ).is_degenerate() )
    {
        err_point_print( p1 );
        err_point_print( p2ronhi.back() );
        err_point_print( p2 );
    }
    assert( !Seg( p1, p2ronhi.back() ).is_degenerate() );
    assert( !l.is_degenerate() );

    CGAL::Object obj = CGAL::intersection( Seg( p1, p2ronhi.back() ), l );
    //
    if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
    {
        return *pobj;
    }
    else
        return CGAL::ORIGIN;
}

/*-----------------------------Sua SGP_distance???----*/
Point myCircularIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2,
                           double alphabar )
{
  Projection_3 myproj;                 // project a point on line
  Line l( l1, l2 );
  Plane h1( p1, l1, l2 );
  Plane h2( p2, l2, l1 );
  Vec n1 = h1.orthogonal_vector();
  Vec n2 = h2.orthogonal_vector();

  // unfold p1 on h2
  Point p1onl = myproj( l, p1 );
  Point p2onl = myproj( l, p2 );
  double d1 = sqrt( ( p1 - p1onl ).squared_length().to_double() );
  double d2 = sqrt( ( p2 - p2onl ).squared_length().to_double() );
  if ( d2 == 0 ){   // p1 is on intersection line
    cerr << "myCircularIntersect: Cannot happen" << endl;
    exit( -1 );
  }
  Point p1onh2 = p1onl + ( p2onl - p2 ) * d1 / d2;

  Transformation T( 1, 0, 0, - p1onh2.x(),
                    0, 1, 0, - p1onh2.y(),
                    0, 0, 1, - p1onh2.z() );
  Direction d = n2.direction();
  double demxy = sqrt( pow( d.dx().to_double(), 2.0 )
                       + pow( d.dy().to_double(), 2.0 ) );
  Transformation Txz( d.dx()/demxy, d.dy()/demxy, 0, 0,
                      -d.dy()/demxy, d.dx()/demxy, 0, 0,
                      0, 0, 1, 0 );
  double demxyz = sqrt( pow( d.dx().to_double(), 2.0 )
                        + pow( d.dy().to_double(), 2.0 )
                        + pow( d.dz().to_double(), 2.0 ) );
  Transformation Tz( d.dz()/demxyz, 0, - RT( demxy /demxyz ), 0,
                     0, 1, 0, 0,
                     RT( demxy /demxyz ), 0, d.dz()/demxyz, 0 );
  //Transformation Rz( cos( alphabar ), -sin( alphabar ), 0, 0,
  //                   sin( alphabar ), cos( alphabar ), 0, 0,
  //                   0, 0, 1, 0 );
  Transformation Tf = Tz * Txz * T;
  Transformation Tb = T.inverse() * Txz.inverse() * Tz.inverse();

  Point p22d = Tf( p2 );
  Point l12d = Tf( l1 );
  Point l22d = Tf( l2 );

  // assert( fabs( p22d.z().to_double() ) < 1e-6 );
  // assert( fabs( l12d.z().to_double() ) < 1e-6 );
  // assert( fabs( l22d.z().to_double() ) < 1e-6 );

  // project all points to plane h2 for 2D algorithms
  Point_2 A( 0, 0 );
  Point_2 B( p22d.x(), p22d.y() );
  Point_2 X( l12d.x(), l12d.y() );
  Point_2 Y( l22d.x(), l22d.y() );
  Line_2 L( X, Y );


  // compute perpendicular bisector of A and B
  Point_2 M = CGAL::midpoint( A, B );
  Line_2 centerline = Line_2( A, B ).perpendicular( M );

  // compute perpendicular line of the tangent of circle
  Transformation_2 R1( cos( -alphabar ), - sin( -alphabar ), 0,
                       sin( -alphabar ), cos( -alphabar ), 0 );
  Transformation_2 R2( cos( 2*CGAL_PI-alphabar ), - sin( 2*CGAL_PI-alphabar ), 0,
                       sin( 2*CGAL_PI-alphabar ), cos( 2*CGAL_PI-alphabar ), 0 );
  Line_2 sideline = ( alphabar < CGAL_PI )? Line_2( A, R1(B) ).perpendicular( A ):
    Line_2( A, R2(B) ).perpendicular( A );

  // compute the center of the circle
  CGAL::Object o1 = CGAL::intersection( centerline, sideline );
  const Point_2 *C = CGAL::object_cast< Point_2 >( &o1 );
  if ( C == NULL ){
    cerr << "Cannot happen" << endl;
    exit( -1 );
  }

  // compute the intersections
  double radius2 = CGAL::squared_distance( *C, A ).to_double();
  Vector_2 Delta = X - *C;
  Vector_2 D = Y - X;
  double delta = pow( ( D * Delta ).to_double(), 2.0 )
    - D.squared_length().to_double() * ( Delta.squared_length().to_double() - radius2 );
  if ( delta + 1e-6 < 0.0 ){
    cerr << "Cannot happen ! No intersection ?!?!?!" << endl;
    exit( -1 );
  }
  delta = fabs( delta );
  double t1 = ( - ( D * Delta ).to_double() + sqrt( delta ) )/D.squared_length().to_double();
  double t2 = ( - ( D * Delta ).to_double() - sqrt( delta ) )/D.squared_length().to_double();
  Point_2 I1 = X + t1 * D;
  Point_2 I2 = X + t2 * D;
  // cerr << "\t\t\talpha_2(1) = " << angleBetween_2( A - I1, B - I1 )*180/CGAL_PI << endl;
  // cerr << "\t\t\talpha_2(2) = " << angleBetween_2( A - I2, B - I2 )*180/CGAL_PI << endl;

  // lift back to 3D
  Point I1in3d = Tb( Point( I1.x(), I1.y(), 0 ) );
  Point I2in3d = Tb( Point( I2.x(), I2.y(), 0 ) );
  // project point on [l1,l2] (although it's quite near)
  Point N1 = myproj( l, I1in3d );
  Point N2 = myproj( l, I2in3d );

 
  return ( alphabar < CGAL_PI )? N1: N2;
}

/*-----------------------------------------------------------------------------------------*/
// save a path to file
void savePath( const char *fname, int iteration, vector<Point> *path )
{
  // create file name
  char fullname[256];
  sprintf( fullname, "%s%d.vect", fname, iteration );

  // open file
  FILE *fp;
  fp = fopen( fullname, "w" );

  // write the path
  fprintf( fp, "VECT\n" );
  fprintf( fp, "1 %d 1\n", (unsigned int)path->size() );
  fprintf( fp, "%d\n", (unsigned int)path->size() );
  fprintf( fp, "1\n\n" );

  for( unsigned int i = 0; i < path->size(); i ++ ){
    fprintf( fp, "%.5f %.5f %.5f  ",
             (*path)[i].x().to_double(),
             (*path)[i].y().to_double(),
             (*path)[i].z().to_double() );
  }
  fprintf( fp, "\n\n0 0.6 1 1" );

  // close file
  fclose( fp );
}

// trang: pausing the program ... ------------------------------------------
void PressEnterToContinue()
{
    std::cout << "Press ENTER to continue... " << flush;
    std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
}

/* -----------------------------  Le _ Khoang cách SGP ------------------------
 * -----------------------------------------------------------------------*/

  double SGP_distance (const Point p1, const Point p2)
  {
  double t = CGAL::abs(p1.z().to_double() - p2.z().to_double())/SIN_THETA;
  return
          CGAL::max(sqrt (CGAL::squared_distance (p1, p2).to_double ()), t);
  }

/*---------------------------------------------------------------------------*/




