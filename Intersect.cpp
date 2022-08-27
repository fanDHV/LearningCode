#include <math.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
//#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include "myCGAL.h"
#include "Functions.h"

typedef CGAL::AABB_polyhedron_triangle_primitive<K,KPoly> Primitive;
//typedef CGAL::AABB_face_graph_triangle_primitive<K,KPoly> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Primitive_id Primitive_id;
typedef KPoly::Facet_iterator Facet_iterator;
typedef KPoly::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef K::Triangle_3 Triangle;

Point myLSP( Point l1, Point l2, Point p1, Point p2 );

Point ArbitraryRotate1(Point p,double theta,Point p1,Point p2)
{
  //Point u,q1,q2;
  Vec u, q1, q2;
  double d;

  /* Step 1 */
  q1 = p - p1;
  // q1.x = p.x - p1.x;
  // q1.y = p.y - p1.y;
  // q1.z = p.z - p1.z;
  
  u = p2 - p1;
  // u.x = p2.x - p1.x;
  // u.y = p2.y - p1.y;
  // u.z = p2.z - p1.z;
  //Normalise(&u);
  //d = sqrt(u.y*u.y + u.z*u.z);
  d = sqrt( u.squared_length().to_double() );

  /* Step 2 */
  if (d != 0) {
    q2 = Vec( q1.x(), 
              q1.y() * u.z() / d - q1.z() * u.y() / d, 
              q1.y() * u.y() / d + q1.z() * u.z() / d );
    // q2.x = q1.x;
    // q2.y = q1.y * u.z / d - q1.z * u.y / d;
    // q2.z = q1.y * u.y / d + q1.z * u.z / d;
  } else {
    q2 = q1;
  }

  /* Step 3 */
  q1 = Vec( q2.x() * d - q2.z() * u.x(), 
            q2.y(),
            q2.x() * u.x() + q2.z() * d );
  // q1.x = q2.x * d - q2.z * u.x;
  // q1.y = q2.y;
  // q1.z = q2.x * u.x + q2.z * d;

   /* Step 4 */
  q2 = Vec( q1.x() * cos(theta) - q1.y() * sin(theta),
            q1.x() * sin(theta) + q1.y() * cos(theta),
            q1.z() );
  // q2.x = q1.x * cos(theta) - q1.y * sin(theta);
  // q2.y = q1.x * sin(theta) + q1.y * cos(theta);
  // q2.z = q1.z;

  /* Inverse of step 3 */
  q1 = Vec( q2.x() * d + q2.z() * u.x(),
            q2.y(),
            - q2.x() * u.x() + q2.z() * d );
  // q1.x =   q2.x * d + q2.z * u.x;
  // q1.y =   q2.y;
  // q1.z = - q2.x * u.x + q2.z * d;

  /* Inverse of step 2 */
  if (d != 0) {
    q2 = Vec( q1.x(),
              q1.y() * u.z() / d + q1.z() * u.y() / d,
              - q1.y() * u.y() / d + q1.z() * u.z() / d );
    // q2.x =   q1.x;
    // q2.y =   q1.y * u.z / d + q1.z * u.y / d;
    // q2.z = - q1.y * u.y / d + q1.z * u.z / d;
  } else {
    q2 = q1;
  }

   /* Inverse of step 1 */
  Point pp = p1 + q2;
  // q1.x = q2.x + p1.x;
  // q1.y = q2.y + p1.y;
  // q1.z = q2.z + p1.z;
  return(pp);
}

Point ArbitraryRotate2(Point p,double theta,Point p1,Point p2)
{
   Point q(0.0,0.0,0.0);
   double costheta,sintheta;
   Vec r( p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z() );

   // r.x = p2.x - p1.x;
   // r.y = p2.y - p1.y;
   // r.z = p2.z - p1.z;
   p = p - Vec( p1.x(), p1.y(), p1.z() ); 
   // p.x -= p1.x;
   // p.y -= p1.y;
   // p.z -= p1.z;
   //Normalise(&r);

   costheta = cos(theta);
   sintheta = sin(theta);

   q = q + Vec( (costheta + (1 - costheta) * r.x() * r.x()) * p.x() +
                ((1 - costheta) * r.x() * r.y() - r.z() * sintheta) * p.y() +
                ((1 - costheta) * r.x() * r.z() + r.y() * sintheta) * p.z(),
                ((1 - costheta) * r.x() * r.y() + r.z() * sintheta) * p.x() +
                (costheta + (1 - costheta) * r.y() * r.y()) * p.y() + 
                ((1 - costheta) * r.y() * r.z() - r.x() * sintheta) * p.z(),
                ((1 - costheta) * r.x() * r.z() - r.y() * sintheta) * p.x() + 
                ((1 - costheta) * r.y() * r.z() + r.x() * sintheta) * p.y() +
                (costheta + (1 - costheta) * r.z() * r.z() ) * p.z() );


   q = q + Vec( p1.x(), p1.y(), p1.z() );
   // q.x += p1.x;
   // q.y += p1.y;
   // q.z += p1.z;
   return(q);
}


typedef K::Construct_projected_point_3 projection_3;

Point myLSP( Point l1, Point l2, Point p1, Point p2 )
{
  projection_3 myproj;
  Line l( l1, l2 );
  Plane h1( p1, l1, l2 );
  Plane h2( p2, l2, l1 );
  Vec n1 = h1.orthogonal_vector();
  Vec n2 = h2.orthogonal_vector();
  Point p2onl = myproj( l, p2 );
  std::cout << "p2onl = " << std::flush; point_print( p2onl ); std::cout << std::endl;
  Point p1onl = myproj( l, p1 );
  std::cout << "p1onl = " << std::flush; point_print( p1onl ); std::cout << std::endl;
  double d1 = sqrt( ( p1 - p1onl ).squared_length().to_double() );
  double d2 = sqrt( ( p2 - p2onl ).squared_length().to_double() );
  std::cout << "d1 = " << d1 << ", d2 = " << d2 << std::endl;
  if ( d1 == 0 )
    return p1;
  Point p2ronh1 = p2onl + ( p1onl - p1 ) * d2 / d1;
  std::cout << "p2ronh1 = " << std::flush; point_print( p2ronh1 ); std::cout << std::endl;

  // Point p2onh1 = myproj( h1, p2 );
  // std::cout << "p2onh1 = " << std::flush; point_print( p2onh1 ); std::cout << std::endl;
  // double d2 = sqrt( ( p2 - p2onl ).squared_length().to_double() );
  // double d2onh1 = d2 * cos( angleBetween( n2, n1 ) );
  // Vec v2 = ( p2onh1 - p2onl ) * d2 / d2onh1;
  // Point p2ronh1 = p2onl + v2;

  assert( !Seg( p1, p2ronh1 ).is_degenerate() );
  assert( !l.is_degenerate() );
  CGAL::Object obj = CGAL::intersection( Seg( p1, p2ronh1 ), l );
  if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) ){
    return *pobj;
  } else
    return Point( 0,0,0 );
} 
