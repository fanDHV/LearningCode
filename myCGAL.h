#ifndef __MYCGAL_H
#define __MYCGAL_H

#include <CGAL/Cartesian.h>
#include <CGAL/array.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Aff_transformation_2.h>

#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>

// rational number
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq                       RT;
typedef CGAL::Gmpz                       IT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::Quotient<CGAL::MP_Float>   RT;
typedef CGAL::Quotient<CGAL::MP_Integer> IT;
#endif

// typing shortcuts
//typedef CGAL::Exact_predicates_exact_constructions_kernel<RT> K;
typedef CGAL::Cartesian<RT>             K;
typedef K::Point_3                      Point;
typedef K::Plane_3                      Plane;
typedef K::Vector_3                     Vec;
typedef K::Direction_3                  Direction;
typedef K::Line_3                       Line;
typedef K::Segment_3                    Seg;
typedef K::Triangle_3                   Triangle;
typedef CGAL::Polyhedron_3<K>           KPoly;
typedef K::Aff_transformation_3         Transformation;
typedef K::Construct_projected_point_3  Projection_3;
typedef K::Construct_supporting_plane_3 Supporting_plane;

typedef K::Point_2                     Point_2;
typedef K::Aff_transformation_2        Transformation_2;
typedef K::Line_2                      Line_2;
typedef K::Vector_2                    Vector_2;
typedef K::Direction_2                 Direction_2;


void  point_print( const Point  & pnt );
void err_point_print( const Point & pnt );

// Surpresses output of debugging information if pre-processor flag is not set
#ifdef DEBUG
#define debug(s) {cerr << "DBG:\t" << s << endl;}
#else
#define debug(s)
#endif

#endif
