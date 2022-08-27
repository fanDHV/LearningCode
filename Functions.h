#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

/*****************************************************************************
 * Functions.h                                                               *
 *                                                                           *
 * Contains prototypes for various useful functions.  Mostly math related    *
 *****************************************************************************/

//#include "Leda.h"
#include <vector>
#include "myCGAL.h"
#include "NNLocation.h"

struct GraphNode;

//#define  GG_NORTH_POLE_IND  -10
#define  GG_INVALID         -999

//  phi - goes from the north pole to the south pole -
//  theta - longtitued 0 ... 2pi
class  GridGraph {
private:
    int  iiBound,
        jjBound,
        limit;
    int  entriesNum, src_index, dst_index;
    GraphNode   ** arr;
    GraphNode  * pNorthPole;

public:
    GridGraph() {
        //memset( this, 0, sizeof( GridGraph ) );
    }

    int   AppendVertex( GraphNode  * gn );
    void  connect( int  i, int  j, int  ind_new );
    void  connectByIndices( int  src, int  dst );

    //void  resetShortestPathPointers();

    void  addNorthPole( GraphNode  * gn ) {
        pNorthPole = gn;
    }

    int  calcIndex( int  phi_ind, int  theta_ind );

    void  dump()  const {
        printf( "GRID: ----------------------------------\n" );
        printf( "GRID: entriesNum: %d  iiBound: %d,  jjBound: %d\n",
                entriesNum, iiBound, jjBound );
        printf( "GRID: ----------------------------------\n" );

    }

    int  getVerticesNum() const {
        return  entriesNum;
    }
    int  getPhiBound() const {
        return  iiBound;
    }
    int  getThetaBound() const {
        return  jjBound;
    }

    int   getNorthPoleIndex() const {
        return  iiBound * jjBound + 1;
    }

    int   getSouthPoleIndex() const {
        return  iiBound * jjBound;
    }

    GraphNode  * getNode( int  ind ) const {
        //if  ( ind == GG_NORTH_POLE_IND )
        //    return  pNorthPole;

        assert( ( 0 <= ind )  &&  ( ind < entriesNum ) );

        return  arr[ ind ];
    }

    void  init( int  _iiBound, int   _jjBound ) {
        int  size;

        //memset( this, 0, sizeof( GridGraph ) );

        iiBound = _iiBound;
        jjBound = _jjBound;
        limit = iiBound * jjBound + 10;

        size = limit * sizeof( GraphNode * );

        arr = (GraphNode **)malloc( size );
        memset( arr, 0, size );
        assert( arr != NULL );

        entriesNum = 0;
    }

    void  term() {
        //memset( this, 0, sizeof( GridGraph ) );
    }

    double  vertices_distance( int  a, int  b )  const;               // classical distance
    double  vertices_SGP_distance( int  a, int  b )  const;           // SGP distance

    double  getVertex_DoubleXCoord( int  i ) const;
    double  getVertex_DoubleYCoord( int  i ) const;
    double  getVertex_DoubleZCoord( int  i ) const;

    void  set_distance( int  to, int  j, double  dist );

    void  addNodeToGraph( GraphNode  * gn, int  idx );
    void  PrintGraph();
    void  CalcDists();                                                                  // classical distance
    void  CalcSGPDists();                                                               // SGP distance                                              

    void  setSourceVertex( int  _source_ind ) {
        src_index = _source_ind;
    }

    void  setDestVertex( int  _dest_ind ) {
        dst_index = _dest_ind;
    }


    GraphNode * FindShortestPathInGraph( const int  src, const int  trg );
    GraphNode * FindShortestGentlePathInGraph( const int  src, const int  trg );

};


// returns the cross product of the two vectors
Vec CrossProduct( const Vec &, const Vec & );
//Vec CrossProductExact(const Vec &, const Vec &);
//Vec CrossProductDir(const Vec & u, const Vec & v, int  place );

// returns the unit vector in the direction of the given vector
Vec Unit(const Vec &);

// reduces the vector
Vec Reduce(const Vec &);

// simplifies the rationals in the vector
void Simplify(Vec &);

void SimplifyCrossProduct(Vec & v);

// returns the dot product of the two vectors
RT DotProduct(const Vec &, const Vec &);

// returns true if the first angle ( in radian) is in between the next two
//bool AngleBetween(double, double, double);
double angleBetween( const Vec &v1, const Vec &v2 );
double angleBetween_2( const Vector_2 &v1, const Vector_2 &v2 );

// calculates the theta parameter of the spherical coordinates of the given
// point
double CalcTheta(const Point &x0, const RT &r);

// calculates the phi parameter of the spherical coordinates of the given
// point
double CalcPhi(const Point &x0, const RT &r);

// projects a point inside of a sphere along a vector onto the sphere
Point ProjectPointOntoSphere(const Point &V, const Vec &N, const RT r);

// creates a new node with the appropriate adjacent nodes connected
//GraphNode * CreateGraphNode(Point &, NNLocation &, int, int, int, int);
GraphNode * CreateGraphNode( GridGraph  & grid,
                             Point & p, NNLocation & nnloc,
                             int i, int j,
                             Point  & pnt_on_sphere );

// returns a point on the intersection of two planes
Point PointOn(const Plane &, const Plane &);

// rotates a vector around another vector by a given number of degrees and
// returns the result
Vec Rotate(const Vec &, const Vec &, double);

// "unfolds" a point on a plane so that it lies in the second plane
Point Unfold(const Point &, const Plane &, const Plane &);


// returns the intersection point of a line given by two points and a plane
// also returns a parameter which is true if the line lies on the plane
Point IntersectLinePlane(const Point &, const Point &, const Plane &, bool &);


// Creates a vector from rationals, simplifies it, and converts it to a point
Point VecToPoint( const RT  & x,
                  const RT  & y,
                  const RT  & z );

// returns the point of intersection of an "unfolded" line given by the points
// and the line of intersection of the two planes
Point RidgePoint(const Point &, const Point &, const Plane &, const Plane &);

Point ExactExtRidgePoint(const Point & _p1, const Plane  h1,
                         const Vec   & h1_norm,
                         const Point & _p2, const Plane  h2,
                         const Vec   & h2_norm,
                         const Point  & center );
Point ExtRidgePoint(const Point & p1, const Plane  h1,
                    Vec    h1_norm,
                    const Point & p2, const Plane  h2,
                    Vec    h2_norm,
                    const Point  & center );
Point ExtUnfold( const Point  & center,
                 Plane   target_plane,
                 Plane   src_plane,
                 const Point  & pnt );


// returns true if the square distance between the two points is
// less than or equal to 10^-6
bool PointEq(const Point & p1, const Point & p2);
Point  approx_pnt( const Point  & p );
Vec  approx_vec( const Vec  & p );

// return normal vector to a plane
Vec normalVector( const Plane & );

Point myLineIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2 );

Point myCircularIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2,
                           double alphabar );

// save a path to file
void savePath( const char *fname, int iteration, std::vector<Point> *path );

// trang: to pause program
void PressEnterToContinue();

// trang: find intersection in which there some faces between
Point myLineIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2,
                       const std::vector< Triangle > TriBetween, unsigned int onFace);

// trang: find intersection in which there some segments between
Point myLineIntersect( const Point &l1, const Point &l2, const Point &p1, const Point &p2,
                       const std::vector< Seg > SegBetween );

// norm_SGP
double SGP_distance (const Point p1, const Point p2);



#endif
