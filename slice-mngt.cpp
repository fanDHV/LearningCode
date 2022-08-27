/*
  Authors: P.T. An, T.V. Hoai, N.N. Hai, Trang
  Creation date: 6 Jan 2013

  Modified by Trang for shortest descending path problems
  Starting date: ... Dec 2013

  Modified by Le for shortest gentle path problems
  Starting date: ... Dec 2019

 use taking specify prevp and nextp

 */
#include <algorithm>
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <cassert>
#include <cstdlib>
#include <limits.h>

#include "myCGAL.h"
#include "sarray.h"
#include "Functions.h"
#include "Path.h"
#include "Face.h"
#include "Vert.h"
#include "Edge.h"
#include "Poly.h"

#include "global-settings.h"
#include "slice-mngt.h"

using namespace std;


// adding to wite to file
#include <fstream>
// for number decimals
#include <iomanip>
// for generate random number throuhout time
#include <time.h>


typedef struct {
  RT z_;
  const Vert *v_;
  int id_;
} szbasedVert;
szbasedVert * zsortVerts;

//to sort vestices by descreasing z-coords.
static bool zOrder( szbasedVert v1, szbasedVert v2 )
{
  return v1.z_ > v2.z_;
}



double alpha_max = 0, alpha_min = CGAL_PI;

/*------------------------------------------------------------------------
  original one: no including source point.
--------------------------------------------------------------------------*/
sliceMngt::sliceMngt( Poly *poly )
  : poly_( poly ), atVertex_( NULL ), onSegment_( NULL ), spaths_( NULL ),
    spathlens_( NULL ), movedirs_( NULL ) //, movedirs_SGP_( NULL )
{
  int i, j, k;
  int nverts = poly->NumVerts();
  Vert * const *verts = poly->Verts();

  // debug message
  cout << "sliceMngt::constructor(): start" << endl;

  // create a temporary array for vertices
  szbasedVert * zsortVerts = new szbasedVert[nverts];
  for( i=0; i<nverts; i++ ){
    zsortVerts[i].z_ = verts[i]->Location().z();
    zsortVerts[i].v_ = verts[i];
    zsortVerts[i].id_ = i;
  }

  // sort vertices w.r.t. z-coordinate
  stable_sort( zsortVerts, zsortVerts+nverts, zOrder );

  // create slices paralel to Oxy, passing thru given points
  Vec norm0xy( 0.0, 0.0, 1.0 );
  vector<Plane> planes;
  for( i=0; i<nverts; i++ ){
    // do not add planes with the same z-coordinate
    if ( !planes.empty() && zsortVerts[i].v_->Location().z() == planes.back().point().z() )
      continue;
    planes.push_back( Plane( zsortVerts[i].v_->Location(), norm0xy ) );
  }

  // trang's comment: release memory for sorting
  delete [] zsortVerts;

  // intersect planes with polyhedron
  Face * const *faces = poly->Faces();
  //
  for( i=0; i<(int)planes.size(); i++ )
  {
    cout << "------------- Plane " << i << " --------------" << endl;
    list< Seg > seglist;
    list< int > facelist;

    list< pair< Edge*, Edge* > > ledgelist; // Le's: list of lateral edges
    list< pair< bool, bool > > islaterallist;  //Le's: check list of lateral edges

    for( j=0; j<poly->NumFaces(); j++ ){
      deque< Point > pp;
      deque< Edge* > ledges;
      assert( faces[ j ]->EdgeDegree() == 3 );
      for( k = 0; k < faces[ j ]->EdgeDegree(); k ++ ){
        // create edge of triangle
        Seg triedge( faces[ j ]->Edges()[ k ]->V1()->Location(),
                     faces[ j ]->Edges()[ k ]->V2()->Location() );
        // intersect with plane
        CGAL::Object obj = CGAL::intersection( planes[ i ], triedge );
        // two possibilities for intersection
        // - point
        //   + vertex of polytope
        //   + interior of polytope edge
        // - segment: the vertex set must be z-coordinate monotone
        deque< Point > ipoints;
        if ( const Seg *sobj = CGAL::object_cast< Seg >( &obj ) ){
          ipoints.push_back( sobj->source() );
          ipoints.push_back( sobj->target() );
        } else if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) ){
          ipoints.push_back( *pobj );
        }
        for( int yyy = 0; yyy < (int)ipoints.size(); yyy ++ ){
          int m = 0;
          for( ; m < (int)pp.size(); m ++ )
            if ( pp[ m ] == ipoints[ yyy ] )
              break;
          if ( m == (int)pp.size() ){
            pp.push_back( ipoints[ yyy ] );
            // check if not lateral edge
            if ( ipoints[ yyy ] == faces[ j ]->Edges()[ k ]->V1()->Location() ||
                 ipoints[ yyy ] == faces[ j ]->Edges()[ k ]->V2()->Location() )
              ledges.push_back( NULL );
            else
              ledges.push_back( faces[ j ]->Edges()[ k ] );
          }
        }
      }

      // check there is intersection with the plane
      //assert( pp.empty() || pp.size() == 2 );
      if ( pp.size() == 2 ){
        seglist.push_back( Seg( pp[ 0 ], pp[ 1 ] ) );
        ledgelist.push_back( pair< Edge*, Edge* >( ledges[ 0 ], ledges[ 1 ] ) );
        facelist.push_back( j );
      }

    }

    // construct the slicing polygon (if #segments > 0)
    cout << "\tcontructing polygon" << endl;
    deque< Seg > ccwpolygon;
    deque< int > ccwfaceidx;
    deque< pair< Edge*, Edge* > > ccwledges;
    while( ! seglist.empty() ){
      list< Seg >::iterator si;
      list< int >::iterator fi;
      list< pair< Edge*, Edge* > >::iterator lei;
      cout << "\t number of segments in the list of segment |L|=" << seglist.size() << endl;
      for( si = seglist.begin(), fi = facelist.begin(), lei = ledgelist.begin();
           si != seglist.end();
           si ++, fi ++, lei ++ )
        if ( ccwpolygon.empty() ){
          ccwpolygon.push_back( *si );
          ccwfaceidx.push_back( *fi );
          ccwledges.push_back( *lei );
          seglist.erase( si );
          facelist.erase( fi );
          ledgelist.erase( lei );
          break;
        }
        else if ( si->source() == ccwpolygon.front().source() ){
          ccwpolygon.push_front( Seg( si->target(), si->source() ) );
          ccwfaceidx.push_front( *fi );
          ccwledges.push_front( pair<Edge*,Edge*>( lei->second, lei->first ) );
          seglist.erase( si );
          facelist.erase( fi );
          ledgelist.erase( lei );
          break;
        }
        else if ( si->target() == ccwpolygon.front().source() ){
          ccwpolygon.push_front( *si );
          ccwfaceidx.push_front( *fi );
          ccwledges.push_front( *lei );
          seglist.erase( si );
          facelist.erase( fi );
          ledgelist.erase( lei );
          break;
        }
        else if ( si->source() == ccwpolygon.back().target() ){
          ccwpolygon.push_back( *si );
          ccwfaceidx.push_back( *fi );
          ccwledges.push_back( *lei );
          seglist.erase( si );
          facelist.erase( fi );
          ledgelist.erase( lei );
          break;
        }
        else if ( si->target() == ccwpolygon.back().target() ){
          ccwpolygon.push_back( Seg( si->target(), si->source() ) );
          ccwfaceidx.push_back( *fi );
          ccwledges.push_back( pair<Edge*,Edge*>( lei->second, lei->first ) );
          seglist.erase( si );
          facelist.erase( fi );
          ledgelist.erase( lei );
          break;
        }
      // if cannot add more, it is due to duplication of intersection segments
      if ( !seglist.empty() && si == seglist.end() ){
        seglist.clear();
      }
    }

    if ( !ccwpolygon.empty() ){
      // reverse the order to counterclockwise looking from z-infinitive
      Vec orthvec = CrossProduct( ccwpolygon[1].source() - ccwpolygon[0].source(),
				  ccwpolygon[2].source() - ccwpolygon[1].source() );

      if ( orthvec.z() < 0 ){
        for( j=0; j < (int)ccwpolygon.size() / 2; j++ ){
          Seg tmp = ccwpolygon[ j ].opposite();
          ccwpolygon[ j ] = ccwpolygon[ ccwpolygon.size() - j - 1 ].opposite();
          ccwpolygon[ ccwpolygon.size() - j - 1 ] = tmp;
          int itmp = ccwfaceidx[ j ];
          ccwfaceidx[ j ] = ccwfaceidx[ ccwfaceidx.size() - j - 1 ];
          ccwfaceidx[ ccwfaceidx.size() - j - 1 ] = itmp;
          pair< Edge*, Edge* > ltmp = ccwledges[ j ];
          ccwledges[ j ] = ccwledges[ ccwledges.size() - j - 1 ];
          ccwledges[ ccwledges.size() - j - 1 ] = ltmp;
        }
        if ( ccwpolygon.size() % 2 == 1 ){
          ccwpolygon[ ccwpolygon.size() / 2 ] = ccwpolygon[ ccwpolygon.size() / 2 ].opposite();
          ccwledges[ ccwledges.size() / 2 ] = pair<Edge*,Edge*>( ccwledges[ ccwledges.size() / 2 ].second,
                                                                 ccwledges[ ccwledges.size() / 2 ].first );
        }
      }
    }

    // dump the polygon of slice
    // cout << "\t" << flush;
    // for( j = 0; j < (int)ccwpolygon.size(); j++ ){
    //   point_print( ccwpolygon[j].source() );
    //   cout<< "--";
    // }
    // cout << endl;

    // check the face
    Face * const * faces = poly_->Faces();
    for( int x = 0; x < (int)ccwpolygon.size(); x ++ )
      assert( faces[ ccwfaceidx[ x ] ]->isOn( ccwpolygon[ x ].source() ) );

    // save for later use
    slices_.push_back( ccwpolygon );
    slicefaces_.push_back( ccwfaceidx );
    slicelaterals_.push_back( ccwledges );

    cout << "\tdone" << endl;
  }

  // create a data structure to store path intersection with slices
  atVertex_ = new Point[ slices_.size() ];
  onSegment_ = new int[ slices_.size() ];

  // // XXX (start): should not be here (possibly)
  // spaths_ = new deque<Point>*[ slices_.size() + 1 ];
  // for( i = 0; i < (int)slices_.size() + 1; i++ )
  //   spaths_[i] = NULL;
  // // XXX (end):

  // set effective planes to all
  fromSlice_ = 0;
  toSlice_ = slices_.size();

  cout << "sliceMngt::constructor(): end" << endl;
}

/*-------------------------------------------------------------------*/
sliceMngt::~sliceMngt()
{
  if ( atVertex_ != NULL )
    delete [] atVertex_;
  if ( onSegment_ != NULL )
    delete [] onSegment_;
  if ( spaths_ != NULL ){
    for( int i = 0; i<(int)slices_.size() + 1; i++ )
      if ( spaths_[i] != NULL )
        delete spaths_[i];
    delete [] spaths_;
  }
  if ( spathlens_ != NULL )
    delete [] spathlens_;
  if ( movedirs_ != NULL )
    delete [] movedirs_;
//      if ( movedirs_SGP_ != NULL )
 //   delete [] movedirs_SGP_;
}

/*-------------------------- trang: new constructor ---------------------------
 including the slice through source point for shortest descending path problems
-----------------------------------------------------------------------------*/
sliceMngt::sliceMngt(Poly *poly, SourceDest *source)
    : poly_( poly ), atVertex_( NULL ), onSegment_( NULL ), spaths_( NULL ),
      spathlens_( NULL ), movedirs_( NULL ) //, movedirs_SGP_( NULL )
{
    int i, j, k;
    int nverts = poly->NumVerts();
    Vert * const *verts = poly->Verts();

    // debug message
    cout << "sliceMngt::constructor(): start" << endl;

    // original one: create a temporary array for vertices
    /*
    szbasedVert *zsortVerts = new szbasedVert[nverts];
    for( i=0; i<nverts; i++ ){
      zsortVerts[i].z_ = verts[i]->Location().z();
      zsortVerts[i].v_ = verts[i];
      zsortVerts[i].id_ = i;
    }
    */

    // trang (2013/12/22): including source point ...
    nverts++;
    //
    v_source = new Vert(poly, source->m_srcPoint);
    // trang's debug
    cout << "%% v_source: "
         << v_source->Location().x() << " "
         << v_source->Location().y() << " "
         << v_source->Location().z() << " " << endl;

    zsortVerts = new szbasedVert[nverts];
    zsortVerts[0].z_  = source->m_srcPoint.z();
    zsortVerts[0].v_  = v_source;
    zsortVerts[0].id_ = 0;

    // trang's debug
    cout << "%% first item is ok" << endl;

    //
    for( i = 1; i < nverts; i++ ){
      zsortVerts[i].z_ = verts[i-1]->Location().z();
      zsortVerts[i].v_ = verts[i-1];
      zsortVerts[i].id_ = i;
    }

    // trang's debug
    cout << "%% Sorting the vertices ..." << endl;

    // sort vertices w.r.t. z-coordinate
    stable_sort( zsortVerts, zsortVerts+nverts, zOrder );

    // trang's debug
    /*
    cout << "\t- After sorting the vertices ..." << endl;
    for ( i = 0; i < nverts; i++ )
    {
        cout << "\t\t vertex " << i << ": ";
        point_print( zsortVerts[i].v_->Location() );
        cout << endl;
    }
    //if ( isDebug_) PressEnterToContinue();
    */

    // trang: hold ing id_ of vertices after sorting by z-coordiate
    /*
    sorted_verts_ = new Vert * [ nverts ];
    for ( i = 0; i < nverts + 1; i++ )
    {
        sorted_verts_[i] = zsortVerts[i].v_;
    }
    */

    // trang: create slices paralel to Oxy through all vertices of the terrain
    cout << "%% start to creating the cutting slices: " << nverts << endl;
    Vec norm0xy( 0.0, 0.0, 1.0 );
    vector<Plane> planes;
    //
    for( i=0; i<nverts; i++ )
    {
      // do not add planes with the same z-coordinate
      if ( !planes.empty() && zsortVerts[i].v_->Location().z() == planes.back().point().z() )
        continue;
      planes.push_back( Plane( zsortVerts[i].v_->Location(), norm0xy ) );
    }

    // trang's comment: release memory for sorting
    // delete [] zsortVerts;

    // intersect planes with polyhedron
    Face * const * faces = poly->Faces();
    //
    for( i=0; i<(int)planes.size(); i++ )
    {
      cout << "------------- Plane " << i << " --------------" << endl;
      list< Seg > seglist;
      list< int > facelist;

      list< pair< Edge*, Edge* > > ledgelist;
      //list< pair< bool, bool > > islaterallist;

      // trang: saving faces sharing current vertex for using straightest geodesic ...
      double total_angle_vi = 0;

      //
      for( j = 0; j < poly->NumFaces(); j++ )
      {
          // trang: add adjacent face of current vertex for using straightest geodesic ...
          if ( faces[ j ]->isVertexAdj( zsortVerts[i].v_, k ) )
          {
              double angle;
              int vDeg_fi = faces[ j ]->VertDegree();
              //
              if ( k==0 )
              {
                  angle = angleBetween
                  (
                    faces[ j ]->Verts()[ k + 1 ]->Location() - faces[ j ]->Verts()[ k ]->Location(),
                    faces[ j ]->Verts()[ vDeg_fi - 1 ]->Location() - faces[ j ]->Verts()[k]->Location()
                  );
              }
              else
                  if ( k == vDeg_fi - 1 )
                  {
                      angle = angleBetween
                      (
                        faces[ j ]->Verts()[ 0 ]->Location() - faces[ j ]->Verts()[ k ]->Location(),
                        faces[ j ]->Verts()[ k - 1 ]->Location() - faces[ j ]->Verts()[k]->Location()
                      );
                  }
              else
              {
                  angle = angleBetween
                  (
                              faces[ j ]->Verts()[ k+1 ]->Location() - faces[ j ]->Verts()[ k ]->Location(),
                              faces[ j ]->Verts()[ k-1 ]->Location() - faces[ j ]->Verts()[ k ]->Location()
                  );
              }

              // trang's debug
              cout << "\t - face " << j << " is adjacent, angle at vertex "
                   << k << ": " <<  angle << " - ";
              //if ( isDebug_) PressEnterToContinue();

              // get the angle of vertex k-th of the face
              total_angle_vi += angle;

              /* --- should do this by the following way -------------------
              // trang's debug
              cout << "\t - face " << j << " is adjacent, angle at vertex "
                   << k << ": " <<  faces[ j ]->Angles()[k] << " - ";
              //if ( isDebug_) PressEnterToContinue();

              // get the angle of vertex k-th of the face
              total_angle_vi += faces[ j ]->Angles()[k];
               ----------------------------------------------------------- */
          }

          // original code ----------------------------------------------------------
          deque< Point > pp;
          deque< Edge* > ledges;

          assert( faces[ j ]->EdgeDegree() == 3 );
          for( k = 0; k < faces[ j ]->EdgeDegree(); k ++ )
          {
            // create edge of triangle
            Seg triedge( faces[ j ]->Edges()[ k ]->V1()->Location(),
                         faces[ j ]->Edges()[ k ]->V2()->Location() );
            // intersect with plane
            CGAL::Object obj = CGAL::intersection( planes[ i ], triedge );

            // two possibilities for intersection
            // - point
            //   + vertex of polytope: how to recognize this?
            //   + interior of polytope edge
            // - segment: the vertex set must be z-coordinate monotone
            deque< Point > ipoints;

            if ( const Seg *sobj = CGAL::object_cast< Seg >( &obj ) )
            {
              ipoints.push_back( sobj->source() );
              ipoints.push_back( sobj->target() );
            }
            else if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
            {
              ipoints.push_back( *pobj );
            }

            for( int yyy = 0; yyy < (int)ipoints.size(); yyy ++ )
            {
              int m = 0;
              for( ; m < (int)pp.size(); m ++ )
                if ( pp[ m ] == ipoints[ yyy ] )
                  break;

              if ( m == (int)pp.size() )
              {
                pp.push_back( ipoints[ yyy ] );
                // check if not lateral edge
                if ( ipoints[ yyy ] == faces[ j ]->Edges()[ k ]->V1()->Location() ||
                     ipoints[ yyy ] == faces[ j ]->Edges()[ k ]->V2()->Location() )
                  ledges.push_back( NULL );
                else
                  ledges.push_back( faces[ j ]->Edges()[ k ] );
              }
            }
          }

          // check there is intersection with the plane
          //assert( pp.empty() || pp.size() == 2 );
          if ( pp.size() == 2 ){
            seglist.push_back( Seg( pp[ 0 ], pp[ 1 ] ) );
            ledgelist.push_back( pair< Edge*, Edge* >( ledges[ 0 ], ledges[ 1 ] ) );
            facelist.push_back( j );
          }
      }
      // trang's debug:
      cout << " ==> total vertex angle at " << i << ": " << total_angle_vi << " --- ";
      //if ( total_angle_vi > 2 * CGAL_PI )
      if ( total_angle_vi > 2 * CGAL_PI +1E-4 )
      {
          cout << " -- greater than 2*pi, the program stopped! " << endl; exit(0);
      }
      //if ( isDebug_) PressEnterToContinue();

      // trang: holding the total vertex angle of current vertex of the terrain ...
      vertex_total_angles_.push_back( total_angle_vi );

      // construct the slicing polygon (if #segments > 0)
      cout << "\tcontructing polygon" << endl;
      deque< Seg > ccwpolygon;
      deque< int > ccwfaceidx;
      deque< pair< Edge*, Edge* > > ccwledges;
      while( ! seglist.empty() )
      {
        list< Seg >::iterator si;
        list< int >::iterator fi;
        list< pair< Edge*, Edge* > >::iterator lei;
        cout << "\t the number of segments of polygon is |L|=" << seglist.size() << endl;

        for( si = seglist.begin(), fi = facelist.begin(), lei = ledgelist.begin();
             si != seglist.end();
             si ++, fi ++, lei ++ )

            if ( ccwpolygon.empty() ){
            ccwpolygon.push_back( *si );
            ccwfaceidx.push_back( *fi );
            ccwledges.push_back( *lei );
            seglist.erase( si );
            facelist.erase( fi );
            ledgelist.erase( lei );
            break;
          }
          else if ( si->source() == ccwpolygon.front().source() ){
            ccwpolygon.push_front( Seg( si->target(), si->source() ) );
            ccwfaceidx.push_front( *fi );
            ccwledges.push_front( pair<Edge*,Edge*>( lei->second, lei->first ) );
            seglist.erase( si );
            facelist.erase( fi );
            ledgelist.erase( lei );
            break;
          }
          else if ( si->target() == ccwpolygon.front().source() ){
            ccwpolygon.push_front( *si );
            ccwfaceidx.push_front( *fi );
            ccwledges.push_front( *lei );
            seglist.erase( si );
            facelist.erase( fi );
            ledgelist.erase( lei );
            break;
          }
          else if ( si->source() == ccwpolygon.back().target() ){
            ccwpolygon.push_back( *si );
            ccwfaceidx.push_back( *fi );
            ccwledges.push_back( *lei );
            seglist.erase( si );
            facelist.erase( fi );
            ledgelist.erase( lei );
            break;
          }
          else if ( si->target() == ccwpolygon.back().target() ){
            ccwpolygon.push_back( Seg( si->target(), si->source() ) );
            ccwfaceidx.push_back( *fi );
            ccwledges.push_back( pair<Edge*,Edge*>( lei->second, lei->first ) );
            seglist.erase( si );
            facelist.erase( fi );
            ledgelist.erase( lei );
            break;
          }
        // if cannot add more, it is due to duplication of intersection segments
        if ( !seglist.empty() && si == seglist.end() ){
          seglist.clear();
        }
      }

      if ( !ccwpolygon.empty() ){
        // reverse the order to counterclockwise looking from z-infinitive
        Vec orthvec = CrossProduct( ccwpolygon[1].source() - ccwpolygon[0].source(),
                    ccwpolygon[2].source() - ccwpolygon[1].source() );

        if ( orthvec.z() < 0 ){
          for( j=0; j < (int)ccwpolygon.size() / 2; j++ ){
            Seg tmp = ccwpolygon[ j ].opposite();
            ccwpolygon[ j ] = ccwpolygon[ ccwpolygon.size() - j - 1 ].opposite();
            ccwpolygon[ ccwpolygon.size() - j - 1 ] = tmp;
            int itmp = ccwfaceidx[ j ];
            ccwfaceidx[ j ] = ccwfaceidx[ ccwfaceidx.size() - j - 1 ];
            ccwfaceidx[ ccwfaceidx.size() - j - 1 ] = itmp;
            pair< Edge*, Edge* > ltmp = ccwledges[ j ];
            ccwledges[ j ] = ccwledges[ ccwledges.size() - j - 1 ];
            ccwledges[ ccwledges.size() - j - 1 ] = ltmp;
          }
          if ( ccwpolygon.size() % 2 == 1 ){
            ccwpolygon[ ccwpolygon.size() / 2 ] = ccwpolygon[ ccwpolygon.size() / 2 ].opposite();
            ccwledges[ ccwledges.size() / 2 ] = pair<Edge*,Edge*>( ccwledges[ ccwledges.size() / 2 ].second,
                                                                   ccwledges[ ccwledges.size() / 2 ].first );
          }
        }
      }

      // dump the polygon of slice
      // cout << "\t" << flush;
      // for( j = 0; j < (int)ccwpolygon.size(); j++ ){
      //   point_print( ccwpolygon[j].source() );
      //   cout<< "--";
      // }
      // cout << endl;

      // check the face
      Face * const * faces = poly_->Faces();
      for( int x = 0; x < (int)ccwpolygon.size(); x ++ )
        assert( faces[ ccwfaceidx[ x ] ]->isOn( ccwpolygon[ x ].source() ) );

      // save for later use
      slices_.push_back( ccwpolygon );
      slicefaces_.push_back( ccwfaceidx );
      slicelaterals_.push_back( ccwledges );

      cout << "\tdone" << endl;
    }


    // trang: do not delete this list before calculating the total vertex angles
    // at vertices of the terrain
    // delete [] zsortVerts;

    // create a data structure to store path intersection with slices
    atVertex_ = new Point[ slices_.size() ];
    onSegment_ = new int[ slices_.size() ];

    // // XXX (start): should not be here (possibly)
    // spaths_ = new deque<Point>*[ slices_.size() + 1 ];
    // for( i = 0; i < (int)slices_.size() + 1; i++ )
    //   spaths_[i] = NULL;
    // // XXX (end):

    // set effective planes to all
    fromSlice_ = 0;
    toSlice_ = slices_.size();

    cout << "sliceMngt::constructor(): end" << endl;
}

/*---------------------------------------------------------------------
  determining the plane between source and destination points by means
  of z-coornate
---------------------------------------------------------------------*/
void sliceMngt::detectEffectivePlanes( SourceDest *srcdest )
{
    int i;

    // debug message
    cout << "sliceMngt::detect EffectivePlanes: start" << endl;

    // only using unique z-coordinates, and
    // between z-coordinate of source and destination
    for( i = 0; i < (int)slices_.size(); i++ )
      if ( ! slices_[i].empty() &&
           // slices_[i].front().source().z() < srcdest->m_srcPoint.z() ){ // old one
           // trang: to include the plane through source point
           slices_[i].front().source().z() == srcdest->m_srcPoint.z() )
      {
          // trang: include the plane through source point
          fromSliceSource_ = i;
          fromSlice_ = i + 1;
          break;
      }

    for( i++;
         i<(int)slices_.size() && slices_[i].front().source().z() >=
                                  srcdest->m_destPoint.z(); i++ ) ;
    toSlice_ = i;

    // debug message
    cout << "sliceMngt::detect EffectivePlanes: end" << endl;
}

/*------------------------------------------------------------------------*/
void sliceMngt::detectSlicePaths( SPath *path )
{
  SPathNode *ptr = path->first();
  int s = fromSlice_, k, i;
  deque<Point> *slicepath = new deque<Point>;

  // debug message
  cout << "sliceMngt::detect slice path: start" << endl;

  // release and allocate memory to save slice paths
  if ( spaths_ != NULL )
  {
    cout << "... spaths_ is not NULL..." << endl;

    for( i = 0; i < (int)slices_.size() + 1; i++ )
      if ( spaths_[ i ] != NULL )
        delete spaths_[ i ];
    delete [] spaths_;
  }
  cout << "... slices_.size() + 1 =..., toSlice_ =...." << slices_.size() + 1<< " ; " << toSlice_ << endl;

  spaths_ = new deque<Point>*[ slices_.size() + 1 ];
  for( i = 0; i < (int)slices_.size() + 1; i++ )
    spaths_[i] = NULL;
  spathlens_ = new double[ slices_.size() + 1 ];
  movedirs_ = new double[ slices_.size() + 1 ];
  //movedirs_SGP_ = new double[ slices_.size() + 1 ];
  //if ( isDebug_3) PressEnterToContinue();


  // compute paths between slices
  while( s < toSlice_ && ptr != NULL )
  {
    //if ( isDebug_3) PressEnterToContinue();
    if ( ptr->point().z() > slices_[s].front().source().z() )
    {
      // save point to slice path
      slicepath->push_back( ptr->point() );
      // increase pointer
      ptr = ptr->next();
      // continue the loop
      continue;
    }

    if ( ptr->point().z() <= slices_[s].front().source().z() )
    {
      // find intersection of the slice and path (only current segment)
      Seg pseg( slicepath->back(), ptr->point() );
      k = (int)slices_[ s ].size();
      for( i = 0; i < k; i++ ){
        CGAL::Object obj = CGAL::intersection( pseg, slices_[ s ][ i ] );
        //----------------------------------------------------
        // XXX: it could be segment intersection, not solved yet
        //----------------------------------------------------
        if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
        {
          if ( *pobj != slices_[s][i].target() )
          {
            atVertex_[s] = *pobj;
            onSegment_[s] = i;

            // check the face
            Face * const * faces = poly_->Faces();
            assert( faces[ slicefaces_[ s ][ i ] ]->isOn( atVertex_[s] ) );

            break;
          }
        }
      }
      //if ( isDebug_3) PressEnterToContinue();
      // save intersection as last point of previous path
      slicepath->push_back( atVertex_[s] );

      assert( slicepath->size() >= 2 );
      // save previous path and advance to next slice
      spaths_[ s - 1 ] = slicepath;
      s ++;
      slicepath = new deque<Point>;

      // save intersection as first point of next slice
      slicepath->push_back( atVertex_[ s - 1 ] );

      // if ptr->z == intersection.z, advance ptr to next point
      if ( ptr->point().z() == atVertex_[ s - 1 ].z() )
        ptr = ptr->next();
    }
  }
  assert( s == toSlice_ );

  //if ( isDebug_3) PressEnterToContinue();

  // create the last slice path by remaining points
  while( ptr != NULL ){
    slicepath->push_back( ptr->point() );
    ptr = ptr->next();
  }
  assert( slicepath->size() >= 2 );
  spaths_[ s - 1 ] = slicepath;

  // compute length of subpaths
  if ( isDebug_3) PressEnterToContinue();
  for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ )
  {
    spathlens_[ s - 1 ] = 0;
    //updatePathLens_[ s - 1 ] = 0;
    for( int i = 1; i < (int)spaths_[ s - 1 ]->size(); i ++ )
    {

        spathlens_[ s - 1 ] +=  SGP_distance ( ( *spaths_[ s - 1 ])[ i - 1 ],
                                                           (*spaths_[ s - 1 ])[ i ] ) ;
    }
  }

  // debug message
  cout << "sliceMngt::detect slice path: end" << endl;
}




/*------------------------------------------------------------------------*/
void sliceMngt::vibrateSlicePaths( SourceDest *srcdest,
                                   int srcFaceIdx, int destFaceIdx,
                                   Path *spalg,
                                   double RATIO, int nedges )
{
  int s;
  Point src, dest;
  int srcface, destface;
  SPath *slicepath;
  SPathNode *ptr;
  Face * const * faces = poly_->Faces();
  Point spoint;
  //Point newVertex;
  //int newSegment;

  // Le 6-8-2022
  cout << "\t start vibrates path to get initial path..." << endl;
  srand((int)time(0));
  // Le 6-8-2022

  for( s = fromSlice_; s <= toSlice_; s ++ ){
    // vibrate point on slice
    if ( s < toSlice_ &&
         atVertex_[ s ] != slices_[ s ][ onSegment_[ s ] ].source() ){
      if ( rand() % 2 == 0 ){ // go left
        // shift left a nedges edges
        onSegment_[ s ] = ( onSegment_[ s ] + slices_[ s ].size() - nedges ) % slices_[ s ].size();
        // compute the start point
        if ( nedges )
          spoint = slices_[ s ][ onSegment_[ s ] ].target();
        else
          spoint = atVertex_[ s ];
        // get to point w.r.t. RATIO
        // atVertex_[ s ] = atVertex_[ s ] +
        //   RATIO * ( slices_[ s ][ onSegment_[ s ] ].source() - atVertex_[ s ] );
        atVertex_[ s ] = spoint +
          RATIO * ( slices_[ s ][ onSegment_[ s ] ].source() - spoint );
      } else {                // go right
        // shift right a nedges edges
        onSegment_[ s ] = ( onSegment_[ s ] + nedges ) % slices_[ s ].size();
        // compute the start point
        if ( nedges )
          spoint = slices_[ s ][ onSegment_[ s ] ].source();
        else
          spoint = atVertex_[ s ];
        // get to point w.r.t. RATIO
        // atVertex_[ s ] = atVertex_[ s ] +
        //   RATIO * ( slices_[ s ][ onSegment_[ s ] ].target() - atVertex_[ s ] );
        atVertex_[ s ] = spoint +
          RATIO * ( slices_[ s ][ onSegment_[ s ] ].target() - spoint );
      }
    }

    // recompute LSP
    cerr << "From: " << s - 1 << " --> " << s << endl;
    SourceDest sd;
    src = ( s == fromSlice_ )? srcdest->m_srcPoint: atVertex_[ s - 1 ];
    srcface = ( s == fromSlice_ )? srcFaceIdx: slicefaces_[ s - 1 ][ onSegment_[ s - 1 ] ];
    dest = ( s == toSlice_ )? srcdest->m_destPoint: atVertex_[ s ];
    destface = ( s == toSlice_ )? destFaceIdx: slicefaces_[ s ][ onSegment_[ s ] ];
    sd.m_srcPoint = src;
    sd.m_destPoint = dest;
    sd.m_srcFace = faces[ srcface ];
    sd.m_destFace = faces[ destface ];

    // trang: for processing of descending path
    sd.m_srcSlice = slices_[s-1];
    sd.m_srcSliceFaces = slicefaces_[s-1];

    // clear old path for saving new one
    spaths_[ s - 1 ]->clear();
    // check if on the same face
    if ( srcface == destface )
    {
        spaths_[ s - 1 ]->push_back( src );
        spaths_[ s - 1 ]->push_back( dest );
    }
    else
    {
        // original...
        // spalg->spCalcPath( &sd, false, slicepath );

        // trang 2013/12/22: calculate shortest descending path instead
        // spalg->calcSDP( &sd, false, slicepath );

        // Le 12/10: calculate shortest gentle path instead
        spalg->calcSGP( &sd, false, slicepath );

        // save the new slice path
        for( ptr = slicepath->first(); ptr != NULL; ptr = ptr->next() ){
          spaths_[ s - 1 ]->push_back( ptr->point() );
      }
      delete slicepath;
    }
  }
  cout << "\t end vibrates path to get initial path..." << endl;
}

/*-----------------------------trang: initial paths -------------------------*/
void sliceMngt::initializePath( SourceDest *srcdest,
                        int srcFaceIdx, int destFaceIdx,
                        Path *spalg,
                        double RATIO, int nedges )
{
    int s;
    Point src, dest;
    int srcface, destface;
    SPath *slicepath;
    SPathNode *ptr;
    Face * const * faces = poly_->Faces();
    Point spoint;
    //Point newVertex;
    //int newSegment;

    // Le 6-8-2022
    srand((int)time(0));
    // Le 6-8-2022

    for( s = fromSlice_; s <= toSlice_; s ++ )
    {
      // vibrate point on slice
      if ( s < toSlice_ && atVertex_[ s ] != slices_[ s ][ onSegment_[ s ] ].source() )
      {
        if ( rand() % 2 == 0 ){ // go left
          // shift left a nedges edges
          onSegment_[ s ] = ( onSegment_[ s ] + slices_[ s ].size() - nedges ) % slices_[ s ].size();
          // compute the start point
          if ( nedges )
            spoint = slices_[ s ][ onSegment_[ s ] ].target();
          else
            spoint = atVertex_[ s ];
          // get to point w.r.t. RATIO
          // atVertex_[ s ] = atVertex_[ s ] +
          //   RATIO * ( slices_[ s ][ onSegment_[ s ] ].source() - atVertex_[ s ] );
          atVertex_[ s ] = spoint +
            RATIO * ( slices_[ s ][ onSegment_[ s ] ].source() - spoint );
        }
        else {                // go right
          // shift right a nedges edges
          onSegment_[ s ] = ( onSegment_[ s ] + nedges ) % slices_[ s ].size();
          // compute the start point
          if ( nedges )
            spoint = slices_[ s ][ onSegment_[ s ] ].source();
          else
            spoint = atVertex_[ s ];
          // get to point w.r.t. RATIO
          // atVertex_[ s ] = atVertex_[ s ] +
          //   RATIO * ( slices_[ s ][ onSegment_[ s ] ].target() - atVertex_[ s ] );
          atVertex_[ s ] = spoint + RATIO * ( slices_[ s ][ onSegment_[ s ] ].target() - spoint );
        }
      }

      // recompute LSP
      cerr << "From: " << s - 1 << " --> " << s << endl;
      SourceDest sd;
      src = ( s == fromSlice_ )? srcdest->m_srcPoint: atVertex_[ s - 1 ];
      srcface = ( s == fromSlice_ )? srcFaceIdx: slicefaces_[ s - 1 ][ onSegment_[ s - 1 ] ];
      dest = ( s == toSlice_ )? srcdest->m_destPoint: atVertex_[ s ];
      destface = ( s == toSlice_ )? destFaceIdx: slicefaces_[ s ][ onSegment_[ s ] ];
      sd.m_srcPoint = src;
      sd.m_destPoint = dest;
      sd.m_srcFace = faces[ srcface ];
      sd.m_destFace = faces[ destface ];

      // trang: for processing of descending path
      sd.m_srcSlice = slices_[s-1];
      sd.m_srcSliceFaces = slicefaces_[s-1];

      // clear old path for saving new one
      spaths_[ s - 1 ]->clear();
      // check if on the same face
      if ( srcface == destface )
      {
          spaths_[ s - 1 ]->push_back( src );
          spaths_[ s - 1 ]->push_back( dest );
      }
      else
      {
          // original...
          // spalg->spCalcPath( &sd, false, slicepath );

          // trang 2013/12/22: calculate shortest descending path instead
          spalg->calcSGP( &sd, false, slicepath );

          // save the new slice path
          for( ptr = slicepath->first(); ptr != NULL; ptr = ptr->next() ){
            spaths_[ s - 1 ]->push_back( ptr->point() );
        }
        delete slicepath;
      }
    }
}

/*--------------------------------------------------------------------------*/
void sliceMngt::detectPoorSlicePaths( SPath *path )
{
  SPathNode *ptr = path->first();
  int s = fromSlice_, k, i;
  deque<Point> *slicepath = new deque<Point>;

  // debug message
  cout << "sliceMngt::detect slice path: start" << endl;

  // release and allocate memory to save slice paths
  if ( spaths_ != NULL ){
    for( i = 0; i < (int)slices_.size() + 1; i++ )
      if ( spaths_[ i ] != NULL )
        delete spaths_[ i ];
    delete [] spaths_;
  }
  spaths_ = new deque<Point>*[ slices_.size() + 1 ];
  for( i = 0; i < (int)slices_.size() + 1; i++ )
    spaths_[i] = NULL;
  spathlens_ = new double[ slices_.size() + 1 ];
  movedirs_ = new double[ slices_.size() + 1 ];
//   movedirs_SGP_ = new double[ slices_.size() + 1 ];

  // computer paths between slices
  while( s < toSlice_ && ptr != NULL ){
    //--------------------------------------------------------
    // if ptr->z > slice[s].z
    //--------------------------------------------------------
    if ( ptr->point().z() > slices_[s].front().source().z() ){
      // save point to slice path
      slicepath->push_back( ptr->point() );
      // increase pointer
      ptr = ptr->next();
      // continue the loop
      continue;
    }

    //--------------------------------------------------------
    // if ptr->z <= slice[s].z
    //--------------------------------------------------------
    if ( ptr->point().z() <= slices_[s].front().source().z() ){
      // find intersection of the slice and path (only current segment)
      Seg pseg( slicepath->back(), ptr->point() );
      k = (int)slices_[s].size();
      for( i = 0; i < k; i++ ){
        CGAL::Object obj = CGAL::intersection( pseg, slices_[s][i] );
        if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) ){
          if ( *pobj != slices_[s][i].target() ){
            //atVertex_[s] = *pobj;
            atVertex_[s] = slices_[s][i].source();
            onSegment_[s] = i;

            // check the face
            Face * const * faces = poly_->Faces();
            assert( faces[ slicefaces_[ s ][ i ] ]->isOn( atVertex_[s] ) );

            break;
          }
        }
      }

      // save intersection as last point of previous path
      slicepath->push_back( atVertex_[s] );

      // save previous path and advance to next slice
      spaths_[ s - 1 ] = slicepath;
      s ++;
      slicepath = new deque<Point>;

      // save intersection as first point of next slice
      slicepath->push_back( atVertex_[ s - 1 ] );

      // if ptr->z == intersection.z, advance ptr to next point
      if ( ptr->point().z() == atVertex_[s].z() )
        ptr = ptr->next();
    }
  }
  assert( s == toSlice_ );

  // create the last slice path by remaining points
  while( ptr != NULL ){
    slicepath->push_back( ptr->point() );
    ptr = ptr->next();
  }
  spaths_[ s - 1 ] = slicepath;

  // compute length of subpaths
  for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ ){
    spathlens_[ s - 1 ] = 0;
    for( int i = 1; i < (int)spaths_[ s - 1 ]->size(); i ++ ){

        spathlens_[ s - 1 ] +=  SGP_distance ( ( *spaths_[ s - 1 ])[ i - 1 ],
                                                           (*spaths_[ s - 1 ])[ i ] ) ;
    
    }
  }

  // debug message
  cout << "sliceMngt::detect slice path: end" << endl;
}





/*-------------------------multiple shooting process (Le's modification for sgp)----------*/
vector< Point >* sliceMngt::multiShooting( SourceDest *srcdest,
                                           int srcFaceIdx, int destFaceIdx,
                                           Path *spalg )
{
  int s;
  Point prevp, nextp;
  //int prev, next;
  Point src, dest;
  int srcface, destface;
  SPath *slicepath;
  SPathNode *ptr;
  Face * const * faces = poly_->Faces();
  Point newVertex[ toSlice_ - fromSlice_ ];
  int newSegment[ toSlice_ - fromSlice_ ];
  bool allstays = false;

  // debug message
  cout << "sliceMngt::multiple shooting: start" << endl;

  // print the Agarwal's length
  double initiallen = getAprxLength();
  // sliceMngt::getAprxLength(); not SPath::getAprxLength();
  cout << "\t- My Agarwal's length: " << initiallen << endl;
  cout << "\t- Number of slices: " << toSlice_ - fromSlice_
       << " (" << fromSlice_ << ", " << toSlice_ << ")" << endl;
  double length_first, length_last;
  
  length_last = initiallen;

  // %%%%%%%%%%% main loop
  int count = 0;
  for( int j = 0; j < MAX_ITERATIONS; j ++ ) 
  {
    cout << "NOW UP" << j ;
    count += 1;
    length_first = length_last;
    // trang: test the maximum iterations
    //cout << "%% max iteration: " << MAX_ITERATIONS << endl; exit(0);

    /////////////////////////////////////////////////////////////////
    // log incumbent shortest path
    /////////////////////////////////////////////////////////////////

    if ( PATH_LOG )
    {
      vector<Point> *logpath = getShootingPath();
      savePath( "log/sp", j, logpath );
      delete logpath;
    }

    //Le:new type of stopping condition of shortest gentle path
     allstays = checkSGPCollinearConditionInterior( srcdest, spalg, srcFaceIdx, destFaceIdx);


    /////////////////////////////////////////////////////////////////
    // exit the loop because alg. is convergent
    /////////////////////////////////////////////////////////////////
    if ( allstays )
    {
      cerr << "Co-linear condition for SGP matched at iteration: "<< count << endl;
      getSmoothPath();
      //PressEnterToContinue();
      break;
    }

    /////////////////////////////////////////////////////////////////
    // otherwise, compute new points of the slices accordingly
    /////////////////////////////////////////////////////////////////
    // trang: update for horizontal paths
    if ( existHorizontalPath_ ) 
    {
        // trang's debug
        cout << "\t- existing horizontal path ... updating ..." << endl;
        //if ( isDebug_) PressEnterToContinue();

        for( s = fromSlice_; s < toSlice_; s++ ) 
        {
            // trang's debug
            cout << "\t\t- processing on slice " << s << ": " << endl;
            
            // check the between two slices if it contains a horizontal: NAIVE way?!
            // then assign next shooting point ...
            if ( (*spaths_[ s ])[0].z() ==  (*spaths_[ s ])[1].z() )
            {
                // trang's debug
                cout << "\t\t\t- next shooting point is determined" << endl;

                // determine next shooting point
                newVertex[ s - fromSlice_ ] = (*spaths_[ s ])[1];

                // detect the segment containing next shooting : NAIVE way?!
                int i = 0;
                while ( i == onSegment_[s] || !slices_[s][i].has_on( (*spaths_[ s ])[1] ) )
                {
                    i++;
                    if ( (unsigned int)i == slices_[s].size() )
                    {
                        cout << " NOT found ... stopped" << endl; exit(0);
                    }
                }

                // trang's debug
                cout << "\t\t\t- next segment is determined: " << i << endl;

                // determine the face containing the next shooting point
                newSegment[ s - fromSlice_ ] = i;
            }
            else
            {
                // trang's debug
                cout << "\t\t\t- no horizontal on this slice" << endl;

                // dont move shooting point
                {
                    newVertex[ s - fromSlice_ ] = atVertex_[ s ];
                    newSegment[ s - fromSlice_ ] = onSegment_[ s ];
                }
            }
        }

        // mark that the horizontal paths are processed.
        existHorizontalPath_ = false;

        cout << "\t\t- end of updating horizontal paths." << endl;
    }
    else // hoai's original one: refine normally by means of angle-based update
    {
      // trang's debug
      cout << "\t- NO existing horizontal path ..." << endl;

      assert( srcdest->m_srcPoint == (*spaths_[fromSlice_-1])[0]);
      assert( srcdest->m_destPoint == (*spaths_[toSlice_-1])[ (int)spaths_[toSlice_-1]->size()-1] );
      // assert( srcdest->m_destPoint == (*spaths_[toSlice_])[0] ); -> error

      for( s = fromSlice_; s < toSlice_; s++ )
      {
        //
        int prevpFaceIdx, nextpFaceIdx;
        prevpFaceIdx = nextpFaceIdx = -1;

        if ( (int)spaths_[ s - 1 ]->size() == 2)
          {
            prevp = CGAL::midpoint( (*spaths_[ s - 1 ])[0], (*spaths_[ s - 1 ])[1] );
            //prevpFaceIdx = slicefaces_[ s-1 ][ onSegment_[ s-1 ] ];
            //assert (faces[prevpFaceIdx]->isOn(prevp));
          }
          
          else 
          {
            // prevp = this->findColinearPrevpNextp( s, prevpFaceIdx, true);
            prevp = (*spaths_[s-1])[0];
            //prevpFaceIdx = slicefaces_[ s-1 ][ onSegment_[ s-1 ] ];
            //assert (faces[prevpFaceIdx]->isOn(prevp));
          }
        //
          if (spaths_[ s ]->size() == 2)
          {
            nextp = CGAL::midpoint( (*spaths_[ s ])[0], (*spaths_[ s ])[1] );
            //nextpFaceIdx = slicefaces_[ s+1 ][ onSegment_[ s+1 ] ];
            //assert (faces[nextpFaceIdx]->isOn(nextp));
          }
          
          else 
          {
            // nextp = this->findColinearPrevpNextp( s, nextpFaceIdx, false);
            nextp = (*spaths_[s])[(int) spaths_[s]->size()-1];
            //nextpFaceIdx = slicefaces_[ s+1 ][ onSegment_[ s+1 ] ];
            //assert (faces[nextpFaceIdx]->isOn(nextp));          
          }
        // } 
        // %%%%%%%%%%%%% end of getting temporary previous and next points for current slice

        // %%%%%%%%%%%%% compute INDEX
        //prev = ( onSegment_[ s ] + (int)slices_[ s ].size() - 1 ) % (int)slices_[ s ].size();
        //next = ( onSegment_[ s ] + 1 ) % (int)slices_[ s ].size();
        // if the current point is interior to some edge of the slice
        // check whether if prevp and nextp are on a face or not?
        int commonIdx = -1;
        list<int> x = findFaceIndx( prevp );
        list<int> y = findFaceIndx( nextp );
        bool check_empty = empty_intersection( x, y, commonIdx);
        if ( check_empty )
        {
               assert( commonIdx < 0 );
               prevpFaceIdx = x.front();
               nextpFaceIdx = y.front();
               assert (faces[prevpFaceIdx]->isOn(prevp));
               assert (faces[nextpFaceIdx]->isOn(nextp));
        }
        else
        {
                assert( commonIdx >= 0 );
               prevpFaceIdx = nextpFaceIdx = commonIdx;
        }
        cout << "\n list of face containing prevpFaceIdx is .." <<  endl;
        for_each(x.begin(), x.end(), [](const auto &e){ cout << e << endl;});
        cout << "\n prevpFaceIdx is .." << prevpFaceIdx << endl;
        cout << "\nlist of face containing nextpFaceIdx is .." << endl;
        for_each(y.begin(), y.end(), [](const auto &e){ cout << e << endl;});
        cout << "\nnextpFaceIdx is .." << nextpFaceIdx << endl;
        assert (find(x.begin(), x.end(), prevpFaceIdx) != x.end());
        assert (find(y.begin(), y.end(), nextpFaceIdx) != y.end());
        // %%%%%%%%%%%%% end of computting INDEX

        // %%%%%%%%%%%%% upadating
        if ( movedirs_[ s ] != 0) {
          // trang's debug
          cout << "\t\t- slice: " << s - fromSlice_ << " - updating point" << " ... ";
          //if ( isDebug_) PressEnterToContinue();

          computeUpdateSGPInterior( prevp, nextp, prevpFaceIdx, nextpFaceIdx,
                                s, j,
                                newVertex[ s - fromSlice_ ],
                                newSegment[ s - fromSlice_ ], spalg);
        }
        // if no need to move
        else if (nextpFaceIdx == prevpFaceIdx)
        { //if prevp and nextp are on a face...
          assert (prevp.z() >= atVertex_[s].z() && atVertex_[s].z() >= nextp.z());
          //cout << "\t\t\t prevp and nextp are on a face...";
          if ( isDebug_2) PressEnterToContinue();
        
          Seg pseg(prevp, nextp);
          //find intersection of the seg and slices_[s]
          unsigned int k = (int)slices_[ s ].size();
          for (unsigned int i = 0; i < k; i++)
          {
            CGAL::Object obj = CGAL::intersection( pseg, slices_[s][i] );
            if (const Point *pobj = CGAL::object_cast< Point >( &obj ) )
            {
              newVertex[ s - fromSlice_ ] = *pobj;
              newSegment[ s - fromSlice_ ] = i;
            }
          }
        } // end else if (nextpFaceIdx == prevpFaceIdx)
        else 
        { //prevp and nextp are NOT on a face...
          // trang's debug
          cout << "\t\t- slice: " << s - fromSlice_ << " - DON'T' move" << " ... ";
          //if ( isDebug_) PressEnterToContinue();

          newVertex[ s - fromSlice_ ] = atVertex_[ s ];
          newSegment[ s - fromSlice_ ] = onSegment_[ s ];
        }
        //%%%%%%%%%%%%% end of upadating

       
      }
      // end of for-loop for( s = fromSlice_; s < toSlice_; s++ )
    }
    // end of test existHorizontalPath_ 

    // %%%%%%%% recompute shortest paths between slices 
    cout << "\t- recompute shortest paths between slices ..." << endl;
    for( s = fromSlice_; s < toSlice_ +1 ; s ++ )       // return original one b/c from $s-1$ not $s$
    {
      // for( s = fromSlice_; s < toSlice_ +1 ; s ++ ) {  // Le removes it 6-8-2022
      // trang's debug
      cout << "\t\t- computing SGP from slice " << s-1 << " to slice " << s << endl;
      // compute shortest path between consecutive slices (by any shortest
      // path algorithm).
      // Currently, we use Agarwal's algorithm
      // XXX: how to know srcFaceIdx, destFaceIdx ????
      slicepath = NULL;

      // trang's debug
      cout << "\t\t\t- create source/dest ... " << endl;
      //if ( isDebug_) PressEnterToContinue();

      // create source-destination
      src = ( s == fromSlice_ )? srcdest->m_srcPoint: newVertex[ s - fromSlice_ - 1 ];
      //src = ( s == fromSlice_ )? spalg->m_srcdest->m_srcPoint: newVertex[ s - fromSlice_ - 1 ]; // OK like above
      srcface = ( s == fromSlice_ )? srcFaceIdx: slicefaces_[ s - 1 ][ newSegment[ s - fromSlice_ - 1 ] ];
      // trang's debug
      cout << "\t\t\t\t- determining source: ok" << endl;
      if ( isDebug_3) PressEnterToContinue();

      dest = ( s == toSlice_ )? srcdest->m_destPoint: newVertex[ s - fromSlice_ ];
      //dest = ( s == toSlice_ )? spalg->m_srcdest->m_destPoint: newVertex[ s - fromSlice_ ]; // OK like above
      destface = ( s == toSlice_ )? destFaceIdx: slicefaces_[ s ][ newSegment[ s - fromSlice_ ] ];
      // trang's debug
      cout << "\t\t\t\t- determining destination: ok" << endl;
      //if ( isDebug_) PressEnterToContinue();

      SourceDest sd;
      sd.m_srcPoint = src;
      sd.m_destPoint = dest;
      sd.m_srcFace = faces[ srcface ];
      sd.m_destFace = faces[ destface ];

      // trang's debug
      cout << "\t\t\t- completed" << endl;
      //if ( isDebug_) PressEnterToContinue();

      // trang: for processing of descending path
      sd.m_srcSlice = slices_[s-1];
      sd.m_srcSliceFaces = slicefaces_[s-1];

      // clear old path for saving new one
      spaths_[ s - 1 ]->clear();

      // check if on the same face
      if ( srcface == destface )
      {
          // trang's debug
          cout << "\t\t\t- source/dest are on the same face ...";
          //if ( isDebug_) PressEnterToContinue();

          spaths_[ s - 1 ]->push_back( src );
          spaths_[ s - 1 ]->push_back( dest );
          // save length of sub-path
        //          spathlens_[ s - 1 ] = SGP_distance( src, dest );

        // Le sua test SGP va SP
        if(SGP)
        {
        spathlens_[ s - 1 ] += SGP_distance( src, dest );
        }
        else
        {
        spathlens_[ s - 1 ] +=  sqrt (CGAL::squared_distance ( src, dest ).to_double ());
        }
          cout << "completed" << endl;
      }
      else
      {
          // trang's debug
          cout << "\t\t\t- source/dest are NOT on the same face ...";
          //if ( isDebug_) PressEnterToContinue();

          // ----- importance thing: shortest path between consecutive slices ----- //
          AGARWAL_ORIGINAL = false;
          //
          // original: sp on polytope
          // spalg->spCalcPath( &sd, false, slicepath );
          //
          // trang: compute shortest gentle path between consesutive slices
          spalg->calcSGP(&sd, false, slicepath);
          cout << "\t\t\t\t- calling path::calcSGP() ... completed" << endl;
          //
          AGARWAL_ORIGINAL = true;
          // ---------------------------------------------------------------------- //

          // save the new slice path
          // slicepath->dump();

          for( ptr = slicepath->first(); ptr != NULL; ptr = ptr->next() )
          {
            spaths_[ s - 1 ]->push_back( ptr->point() );
          }

          // save length of sub-path
          spathlens_[ s - 1 ] = slicepath->getAprxLength();

          // delete the slice path to save memory
          delete slicepath;

          // trang's debug
          cout << "completed" << endl;
      }


      // move the new points on the slices
      // trang's debug
      cout << "\t\t\t- save new vertex and new edge for next shooting ...";
      if ( s < toSlice_ )
      {
        atVertex_[ s ]  = newVertex [ s - fromSlice_ ];
        onSegment_[ s ] = newSegment[ s - fromSlice_ ];
      }
      cout << " completed" << endl;

    }
    // %%%%%%%% end of recompute shortest paths between slices 
    
    //Le2022: stop because the change of the lengths of path are small
    getSmoothPath();
    length_last = getAprxLength();

    cout << "length_last is "<< length_last << "length_first" << length_first << endl;
    //PressEnterToContinue();
    if (fabs(length_last - length_first) < 1E-8)
    { 
      cout << "...OUT OF LOOP..." << endl;
      //PressEnterToContinue();
      break; // break out of for-loop
    }
    // end of stopping by the change of the lengths of path are small

    // trang's debug
    cout << "\t- end of recomputing shortest paths between slices ..." << endl;
    
  }
  // end of main for-loop for( i = 0; i < MAX_ITERATIONS; i ++ )...
  
  // print the length
  cerr << "sliceMngr::multiShooting - Iteration: " << count
       << " % Length: " << getAprxLength();
    
  
  // create the final approximate path
  vector< Point > *finalpath = getShootingPath();

  double finallen = getAprxLength();
  cerr << "Final length: " << finallen << endl;

  /////////////////////////////////////////////////////////////////
  // log the final shortest path
  /////////////////////////////////////////////////////////////////
  if ( PATH_LOG )
  {
    vector<Point> *logpath = getShootingPath();
    savePath( "log/sp", count, logpath );
    delete logpath;
  }

  // debug message
  cout << "sliceMngt::multiple shooting: end" << endl << endl;

  // final information
  cout << "==================================== Output Information ==================================" << endl;
  cout << "%% Source and destination: " << endl;
  cout << "\t "; point_print( srcdest->m_srcPoint  ); cout << endl;
  cout << "\t "; point_print( srcdest->m_destPoint ); cout << endl;
  cout << "%% Iterations: " << count << endl;
  cout << "%% Cutting slices number: " << toSlice_ - fromSlice_ << endl;
  cout << "%% Lengths of paths: " << endl;
  cout << "\t - initial path:" << initiallen << endl;
  //
  cout << "\t - final path:" << endl;
  // trang (2015/08/30): print the list of path obtained by multiple shooting
  for (unsigned int i = 0; i < finalpath->size(); i++)  
  {
      cout << "\t\t "; point_print( finalpath->at(i) ); cout << endl;
  }
  //
  cout << "\t\t - approximate length: " << finallen << endl;

  //%%%%%%%%% Ghi vao file 2022-4-6
  //Mở file bằng ofstream 

  ofstream ofs("result.txt", ios::app);
  //Kiểm tra file đã mở thành công hay chưa
  if(!ofs)
  {
    cerr << "Error: file not opened." << endl;
    exit(1);
  }
    
  //Ghi vào file 
  ofs << "\n \n New test---";
  ofs << "\n  Source : (" << srcdest->m_srcPoint.x().to_double() << ","
    << srcdest->m_srcPoint.y().to_double() << ","<< srcdest->m_srcPoint.z().to_double() << ")" ;
  ofs <<"\n  Destination: (" << srcdest->m_destPoint.x().to_double() << ","
    << srcdest->m_srcPoint.y().to_double() << ","<< srcdest->m_destPoint.z().to_double() << ")";
  ofs << "\n Index of source and dest faces is " << srcFaceIdx << " and " << destFaceIdx ;
  ofs << "\n ===== Iterations: " <<  count ;
  ofs << "\n checkCollinearCondition is : " <<  allstays ;
  ofs << "\n Cutting slices number: " << toSlice_ - fromSlice_ ;
  ofs << "\n Lengths of initial path: " << setprecision(5) << fixed << initiallen ;
  ofs << "\n Lengths of SGP final path: " << setprecision(5) << fixed << finallen ;

  //Đóng file
  ofs.close();
  //%%%%%%%%%%%%%%%%

  return finalpath;

}
/*------------------------- end of multishooting for SGP-----------------------*/


vector<Point> *sliceMngt::getShootingPath()
{
  vector<Point> *sp = new vector<Point>;
  for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ ){
    for( unsigned int i = 0; i < spaths_[ s - 1 ]->size() - 1; i ++ )
      sp->push_back( (*spaths_[ s - 1 ])[ i ] );
    if ( s == toSlice_ )
      sp->push_back( spaths_[ s - 1 ]->back() );
  }
  return sp;
}


/*-------------------------Le: 19_12_2019: make path smooth-----------------------*/


void sliceMngt::getSmoothPath()
{
    Face * const * faces = poly_->Faces();
    vector<Point> *sp = new vector<Point>;
    vector<Point> *pp = new vector<Point>;
  // copy spaths_ to sp
  for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ )
  {
    for( unsigned int i = 0; i < spaths_[ s - 1 ]->size(); i ++ )
      sp->push_back( (*spaths_[ s - 1 ])[ i ] );
    if ( s == toSlice_ )
      sp->push_back( spaths_[ s - 1 ]->back() );
  }

    // create pp spath
    pp -> push_back( sp->at(0));
    //

    for (unsigned int i = 0; i < sp->size()-2; i++){
    list <int> prevFaceIdx = findFaceIndx( sp->at(i));
    cout << "prevFaceIdx.size(): " << prevFaceIdx.size() << endl;
    // sassert ( prevFaceIdx.size() >=0 );
    list <int> nextFaceIdx = findFaceIndx( sp->at(i+2));
    cout << "nextFaceIdx.size(): " << nextFaceIdx.size() << endl;
    // assert ( nextFaceIdx.size() >=0 );

    int m = -1;
      if ( empty_intersection( prevFaceIdx, nextFaceIdx, m) ){
        pp -> push_back( sp->at(i+1) );
      }
    } // end for
    pp -> push_back( sp->at(sp->size()-1));


    // start to detect new spaths_
    if ( spaths_ != NULL ){
    for( int i = 0; i < (int)slices_.size() + 1; i++ )
      if ( spaths_[ i ] != NULL )
        delete spaths_[ i ];
    delete [] spaths_;
  }
  spaths_ = new deque<Point>*[ slices_.size() + 1 ];
  for( int i = 0; i < (int)slices_.size() + 1; i++ )
    spaths_[i] = NULL;
  spathlens_ = new double[ slices_.size() + 1 ];
  movedirs_ = new double[ slices_.size() + 1 ];
 int s = fromSlice_;
 unsigned int k = 0;
    deque<Point> *slicepath = new deque<Point>;
 while( s < toSlice_ && k < pp->size() )
  {
    //--------------------------------------------------------
    // if iter_f.z > slice[s].z
    //--------------------------------------------------------
    if ( pp->at(k).z() > slices_[s].front().source().z() )
    {
      // save point to slice path
      slicepath->push_back( pp->at(k) );
      // increase pointer
      k++;
      // continue the loop
      continue;
    }

    //--------------------------------------------------------
    // if iter_f.z <= slice[s].z
    //--------------------------------------------------------
    if ( pp->at(k).z() <= slices_[s].front().source().z() ){
      // find intersection of the slice and path (only current segment)
      Seg pseg( slicepath->back(), pp->at(k) );
      for( int i = 0; i < (int)slices_[ s ].size(); i++ ){
        CGAL::Object obj = CGAL::intersection( pseg, slices_[ s ][ i ] );
        //----------------------------------------------------
        // XXX: it could be segment intersection, not solved yet
        //----------------------------------------------------
        if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
        {
          if ( *pobj != slices_[s][i].target() )
          {
            atVertex_[s] = *pobj;
            onSegment_[s] = i;

            // check the face
            // Face * const * faces = poly_->Faces();
            assert( faces[ slicefaces_[ s ][ i ] ]->isOn( atVertex_[s] ) );

            break;
          }
        }
      }

// save intersection as last point of previous path
      slicepath->push_back( atVertex_[s] );

      // save previous path and advance to next slice
      spaths_[ s - 1 ] = slicepath;
      s ++;
      slicepath = new deque<Point>;

      // save intersection as first point of next slice
      slicepath->push_back( atVertex_[ s - 1 ] );

      // if ptr->z == intersection.z, advance ptr to next point
      if ( pp->at(k).z() == atVertex_[s].z() )
        k++;
    }
  }
  assert( s == toSlice_ );

  // create the last slice path by remaining points
  while( k < pp->size() ){
    slicepath->push_back( pp->at(k) );
    k++;
  }
  spaths_[ s - 1 ] = slicepath;

  // compute length of subpaths
  for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ ){
    spathlens_[ s - 1 ] = 0;
    for( int i = 1; i < (int)spaths_[ s - 1 ]->size(); i ++ ){

        spathlens_[ s - 1 ] +=  SGP_distance ( ( *spaths_[ s - 1 ])[ i - 1 ],
                                                           (*spaths_[ s - 1 ])[ i ] ) ;

    }
  }

  // debug message
  cout << "sliceMngt::make path smooth: end" << endl;
  // if (isDebug_2) PressEnterToContinue();

}


/*----------------------------------------------------------------------*/
double sliceMngt::convF( int n, int m, double epsilon )
{
  return epsilon / (double)( n + m );
}

/*----------------------------------------------------------------------*/
double sliceMngt::getAprxLength() const
{
  double len = .0;
  bool isfirst = true;
  Point prevPoint;


     // Le edits to calculate SGP va SP
 if(SGP)
 {
    for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ )
    {
      for( int i = 0; i < (int)spaths_[ s - 1 ]->size() - 1; i ++ )
      {
        if ( isfirst ) 
                isfirst = false;
        else
          len += SGP_distance( prevPoint, (*spaths_[ s - 1 ])[ i ] );
        prevPoint = (*spaths_[ s - 1 ])[ i ];
      }
      if ( s == toSlice_ )
        len += SGP_distance( prevPoint,  spaths_[ s - 1 ]->back() );
    }
    return len;
 }
 else
 {
    for( int s = fromSlice_; s < toSlice_ + 1 ; s ++ ){
    for( int i = 0; i < (int)spaths_[ s - 1 ]->size() - 1; i ++ )
    {
      if ( isfirst )
        isfirst = false;
      else
        len += sqrt (CGAL::squared_distance ( prevPoint, (*spaths_[ s - 1 ])[ i ] ).to_double ());
      prevPoint = (*spaths_[ s - 1 ])[ i ];
    }
    if ( s == toSlice_ )
      len += sqrt (CGAL::squared_distance ( prevPoint,  spaths_[ s - 1 ]->back() ).to_double ());
        }
    return len;
 }
}



/*----------------------------------------------------------------------
  get slices: just for slice showing purpose ...
----------------------------------------------------------------------*/
void sliceMngt::getSlices( deque< deque< Seg > > *&myslices )
{
    int s;

    cout << "\t get slices: start..." << endl;

    if ( myslices != NULL )
      delete myslices;
    myslices = new deque< deque< Seg > >;
    //for( s = fromSlice_; s < toSlice_; s++ ){
    for( s = fromSliceSource_; s < toSlice_; s++ ) // trang: showing slice through source
    {
        myslices->push_back( slices_[ s ] );
    }
    cout << "\t get slices: end..." << endl;
}



/******** Le: new type of stopping condition of shortest gentle path 1/11/2019************/
bool sliceMngt::checkSGPCollinearConditionInterior(SourceDest *srcdest, Path *spalg, int srcFaceIdx, int destFaceIdx)
{

  // trang's debug
  cout << "\t- sliceMngt::checkSGPColinearConditionInterior with new way: start ... " << endl;

  Point prevp, nextp;
  // Point upoint, lpoint;
  int s; // prev, next;
  // double anglein, angleout, alphaL, alphaR;
  Face * const *faces = poly_->Faces();
  // Le 6-8-2022: avoid unused variables
  assert ( faces[srcFaceIdx]->isOn(srcdest->m_srcPoint) );
  assert ( faces[destFaceIdx]->isOn(srcdest->m_destPoint) );
  // Le 6-8-2022: avoid unused variables

  bool allstay = true;
  // double a_max = 0, a_min = CGAL_PI;
  // trang's debug: checking the sorted vertices of terrain w.r.t z-coordinate ...
  // print to check atVertex_[s]
  cerr<< " \t points atVertex_[s] are: ";
  for( s = fromSlice_; s < toSlice_; s ++ )
  {
    cerr<< " atVertex_[" << s << "] = "; cerr << flush;
    err_point_print( atVertex_[s] );
    cerr<< "-->";
  }
  //if ( isDebug_4) PressEnterToContinue();
  cerr << endl;

  for( s = fromSlice_; s < toSlice_; s ++ )
  {
    //--------------------------------------------------------------
    // get previous and next points for current slice
    int prevpFaceIdx, nextpFaceIdx;
    prevpFaceIdx = nextpFaceIdx = -1;

    if ( (int)spaths_[ s - 1 ]->size() == 2)
      {
        prevp = CGAL::midpoint( (*spaths_[ s - 1 ])[0], (*spaths_[ s - 1 ])[1] );
        //prevpFaceIdx = slicefaces_[ s-1 ][ onSegment_[ s-1 ] ];
        //assert (faces[prevpFaceIdx]->isOn(prevp));
      }
      
      else 
      {
        // prevp = this->findColinearPrevpNextp( s, prevpFaceIdx, true);
        prevp = (*spaths_[s-1])[0];
        //prevpFaceIdx = slicefaces_[ s-1 ][ onSegment_[ s-1 ] ];
        //assert (faces[prevpFaceIdx]->isOn(prevp));
      }
    //
      if (spaths_[ s ]->size() == 2)
      {
        nextp = CGAL::midpoint( (*spaths_[ s ])[0], (*spaths_[ s ])[1] );
        //nextpFaceIdx = slicefaces_[ s+1 ][ onSegment_[ s+1 ] ];
        //assert (faces[nextpFaceIdx]->isOn(nextp));
      }
      
      else 
      {
        // nextp = this->findColinearPrevpNextp( s, nextpFaceIdx, false);
        nextp = (*spaths_[s])[(int) spaths_[s]->size()-1];
        //nextpFaceIdx = slicefaces_[ s+1 ][ onSegment_[ s+1 ] ];
        //assert (faces[nextpFaceIdx]->isOn(nextp));          
      }
    // } 
    // %%%%%%%%%%%%% end of getting temporary previous and next points for current slice

    // %%%%%%%%%%%%% compute INDEX
    // if the current point is interior to some edge of the slice
    // check whether if prevp and nextp are on a face or not?
    int commonIdx = -1;
    list<int> x = findFaceIndx( prevp );
    list<int> y = findFaceIndx( nextp );
    bool check_empty = empty_intersection( x, y, commonIdx);
    if ( check_empty )
    {
           assert( commonIdx < 0 );
           prevpFaceIdx = x.front();
           nextpFaceIdx = y.front();
           assert (faces[prevpFaceIdx]->isOn(prevp));
           assert (faces[nextpFaceIdx]->isOn(nextp));
    }
    else
    {
            assert( commonIdx >= 0 );
           prevpFaceIdx = nextpFaceIdx = commonIdx;
    }
    
    //--------create SourcsDest-------
    SourceDest srcdest_pn;
    srcdest_pn.m_srcPoint = prevp;
    srcdest_pn.m_destPoint = nextp;
    //srcdest_pn.m_srcSlice = slices_[s-1];
    //srcdest_pn.m_srcSliceFaces = slicefaces_[s-1];

    srcdest_pn.m_srcFace = faces[prevpFaceIdx];
    srcdest_pn.m_destFace = faces[nextpFaceIdx];
    cout << "\t\t\t\t- determining prevp and nextp points: ok" << endl;

    // Distance between prevp to nextp through out atVertex_[s] is:
    double o_dist_prevp_nextp = getSubpathLength( prevp, nextp, s);
    cout << "\t\t old distance between prevp to nextp through out atVertex_[s] is: ";
    cout << o_dist_prevp_nextp << endl;

    //-------------check prevp and nextp are on a face---------
    //if ( isDebug_2) PressEnterToContinue();
    if (prevpFaceIdx == nextpFaceIdx)
    {
      assert (prevp.z() > atVertex_[s].z() && atVertex_[s].z() > nextp.z());
      cout << "\t\t\t prevp and nextp are on a face...bb-\n";
 
       Seg pseg(prevp, nextp);
      //find intersection of the seg and slices_[s]
      unsigned int k = (int)slices_[ s ].size();
      for (unsigned int i = 0; i < k; i++)
      {
        CGAL::Object obj = CGAL::intersection( pseg, slices_[s][i] );
        if (const Point *pobj = CGAL::object_cast< Point >( &obj ) )
        {
          if ( fabs(o_dist_prevp_nextp - SGP_distance( prevp, *pobj)
                              - SGP_distance( *pobj, nextp) )<1E-3 )
            movedirs_ [ s ] = 0;
          else
            movedirs_ [ s ] = 1;
          
        }
      }
    }
    else // prevp and nextp are not on a face---------
    {
      cout << "\t\t\t prevp and nextp are NOT on a face...\n";

      //----- find SGP from prevp to nextp-------
      SPath *upSubSpaths_;

      // ----- importance thing: shortest path between consecutive slices ----- //
      AGARWAL_ORIGINAL = false;
      // pay attention

      spalg->calcSGP(&srcdest_pn, false, upSubSpaths_);
      // pay attention
      // AGARWAL_ORIGINAL = true;
      double len_new_prevp_next = upSubSpaths_->getAprxLength();
      cout<< "\t \t calculating SGP from prevp to nextp is completed"<< endl;
      if (fabs( len_new_prevp_next - o_dist_prevp_nextp) < 1E-3)
      {
        movedirs_ [ s ] = 0;
      }
      else if (len_new_prevp_next < o_dist_prevp_nextp )
      {
        movedirs_[ s ] = 1;
      }
      else // i.e. upSubSpaths_->getAprxLength() > getSubpathLength( prevp, nextp, s) 
      {
        //cout << "\t\t...it is strange b/c the graph is not dense..." << endl;
        //cout << "\n new path ... upSubSpaths_->getAprxLength() = " << len_new_prevp_next;
        //cout << "\n old path ... getSubpathLength( prevp, nextp, s) = "<< o_dist_prevp_nextp;
        //PressEnterToContinue();
        movedirs_[ s ] = 0;
      }
    }

    // if one moves, the co-linear condition is not matched
    allstay &= ( movedirs_[ s ] == 0 );
  }   // end for loop
  // trang's debug
  cout << "\t- sliceMngt::checkSGPColinearCondition: Ending ... ";

  return allstay;
}




/**** Le: 10/10/2019*****/
Point sliceMngt::Intersection_spath_slices( SPath *spaths_, int s, int &newUpSeg, bool isSeg)
{
    SPathNode *ptr = spaths_->first();
    deque<Point> *edit_path = new deque<Point>;
    Point newUpPoint = CGAL::ORIGIN;
    Face * const * faces = poly_->Faces();
    Seg pseg;
    assert (isSeg == false);

    assert ( ptr->point().z() > slices_[ s ].front().source().z() );

    while( ptr != NULL )
    {
      if ( fabs ( ptr->point().z().to_double() - slices_[ s ].front().source().z().to_double() ) < 1E-3 )
      {
        return ptr->point();  // it is indeed the intersection
      }
      else if ( ptr->point().z() > slices_[ s ].front().source().z() )
      {
        // save point to edit_path
        edit_path->push_back( ptr->point() );
        // increase pointer
        ptr = ptr->next();
        // continue the loop
        continue;
      }
      else  // if ( ptr->point().z() < slices_[ s ].front().source().z() )
      {
        // find intersection of the slice and path (only current segment)
        pseg = Seg( edit_path->back(), ptr->point() );
        break;
      }
    }  // end while
    int k = (int)slices_[ s ].size();
        
    for( int i = 0; i < k; i++ )
    {
        CGAL::Object obj = CGAL::intersection( pseg, slices_[ s ][ i ] );
        if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
        {
            // intersetion is a point ...
            newUpPoint = *pobj;
            newUpSeg = i;
            // check the face
            assert( faces[ slicefaces_[ s ][ i ] ]->isOn( newUpPoint ) );
            break;
        }
        else
        { // intersetion is a segment
            isSeg = true;
            if ( const Seg *sobj = CGAL::object_cast< Seg >( &obj )  )
            {
              newUpPoint = sobj->target();
              newUpSeg = i;
              // check the face
              assert( faces[ slicefaces_[ s ][ i ] ]->isOn( newUpPoint ) );
              break;
            }
        }
    }
  return newUpPoint;
}


/*** Le: find Faces' index containing pnt ***13/10/2019 ***/
// return a array of indices of all faces containing a given point
list<int> sliceMngt::findFaceIndx(const Point &pnt)
{
  list<int> FaceIndx;
  Face * const * faces = poly_->Faces();
  int numFaces_ = poly_->NumFaces();
  for (int i = 0; i < numFaces_; i++)
  {
    if (faces[ i ]->isOn( pnt ))
    {   
        FaceIndx.push_back( i );
        cout << "i= " << i << " findFaceidx continues";
    // if ( isDebug_2) PressEnterToContinue();
    }
  }
  cout << " \n size of FaceIndx =... " << FaceIndx.size();
  //assert (FaceIndx.size() > 0);
  return FaceIndx;
}

bool sliceMngt::empty_intersection(const list<int>& x, const list<int>& y, int & commonIdx)
{
    list<int>::const_iterator i = x.begin();
    list<int>::const_iterator j = y.begin();
    while (i != x.end() && j != y.end())
    {
      if (*i == *j)
      {
        commonIdx =*i;
        return false;
      }
      else if (*i < *j)
        ++i;
      else
        ++j;
    }
    return true;
}


/***************************Le (12/10/2019) update for all cases***************/
void sliceMngt::computeUpdateSGPInterior( Point &prevp, Point &nextp,         // previous, next points
                                                                          // on current path
                                       int &prevpFaceIdx, int &nextpFaceIdx,                 // previous, next segments
                                                                          // on the slice
                                      int s,                              // considered slice
                                      int iteration,                      // iteration of alg.
                                      Point &newPoint, int &newSeg, Path *spalg)
{
  Face * const * faces = poly_->Faces();
  // trang: avoiding unused variable
  iteration = iteration + 0;

  if ( movedirs_[ s ] == 0 )  // not update
  {
      newPoint = atVertex_[ s ];
      newSeg = onSegment_[ s ];
  }
  else  // update
  {
    //----- find SGP from prevp to nextp-------
    //--------create SourcsDest-------
    SourceDest srcdest_pn;
    srcdest_pn.m_srcPoint = prevp;
    srcdest_pn.m_destPoint = nextp;
    // srcdest_pn.m_srcSlice = NULL;
    // srcdest_pn.m_srcSliceFaces = NULL;
    //--------------find face indices of prevp and nextp-------------
    srcdest_pn.m_srcFace = faces[ prevpFaceIdx ];
    srcdest_pn.m_destFace = faces[ nextpFaceIdx ];
    //--------end of ---create SourcsDest-------
    // find SGP from prevp to nextp
    // ----- importance thing: shortest path between consecutive slices ----- //
    //    AGARWAL_ORIGINAL = false;
    // pay attention
    // assert (prevp.z() >= atVertex_[s].z() && atVertex_[s].z() >= nextp.z());
    if (prevpFaceIdx == nextpFaceIdx)
    { // prevp and nextp are on a face
      Seg pseg(prevp, nextp);
      //find intersection of the seg and slices_[s]
      unsigned int k = (int)slices_[ s ].size();
      for (unsigned int i = 0; i < k; i++)
      {
        CGAL::Object obj = CGAL::intersection( pseg, slices_[s][i] );
        if (const Point *pobj = CGAL::object_cast< Point >( &obj) )
        {
          // check if nescesary to update
          if ( fabs(getSubpathLength( prevp, nextp, s)- SGP_distance( prevp, *pobj)
                                - SGP_distance( *pobj, nextp) )<1E-4 )
          {
            cout << "why need to update.., Colinear condition has problem" << endl;
            PressEnterToContinue();
          }
            else
          {
            newPoint = *pobj;
            newSeg = i;
          }
        }
          // else intersection is an empty set
          // there is no case in which intersection is a segment
      }
    }
    else  // prevp and nextp are not on a face
    {
      AGARWAL_ORIGINAL = false;
      SPath *upSubSpaths_(NULL);
       spalg->calcSGP(&srcdest_pn, false, upSubSpaths_);
       // pay attention
      AGARWAL_ORIGINAL = true;
      // if (fabs(upSubSpaths_.getAprxLength() - getSubpathLength( prevp, nextp, s)) < 1E-4)
      // assert (  upSubSpaths_->getAprxLength() <= getSubpathLength( prevp, nextp, s) + 1E-4 ); // Le: 24/12/2019
       //----- find intersection between SGP from prevp to nextp and slices_[s]-------
      bool isSeg_ = false;
      newPoint = Intersection_spath_slices( upSubSpaths_, s, newSeg, isSeg_);
    }

  }
}

 /***************************Le (6/10/2019) end ***************/



/****Le 13/10/2019 purpose of calculate the length of SPath connecting prevp to nextp point */
// prevp and nextp are endpoints of segments of spath_s.
double
sliceMngt::getSubpathLength(const Point &prevp, const Point &nextp, int s) const
{
  double sum_1, sum_2;
  sum_1 = sum_2 = 0;

  if (spaths_[ s - 1 ]->size() == 2)
  {
    assert ( PointEq( prevp, CGAL::midpoint ( (*spaths_[ s - 1 ])[0],
                    (*spaths_[ s - 1 ])[1] )) );
    sum_1 = ( SGP_distance( (*spaths_[ s - 1 ])[0], (*spaths_[ s - 1 ])[1] ) ) *0.5;
  }
  else
  {
    assert ( PointEq( prevp, (*spaths_[ s - 1 ])[0] ) ); // Le: 6-8-2022
    double spathlength_1 = 0;
    for( int i = 1; i < (int)spaths_[ s - 1 ]->size(); i ++ ){
    spathlength_1 +=  SGP_distance ( ( *spaths_[ s - 1 ])[ i - 1 ],
                                                       (*spaths_[ s - 1 ])[ i ] );
    }
      sum_1 = spathlength_1;
  }
  cout << "\t\t calc dist_1 done...sum_1= " << sum_1 << endl;
  if (isDebug_3) PressEnterToContinue();

  if (spaths_[ s ]->size() == 2)
  {
    assert ( PointEq( nextp, CGAL::midpoint  ((*spaths_[ s ])[0], (*spaths_[ s ])[1] )) );
    sum_2 = ( SGP_distance( (*spaths_[ s ])[0], (*spaths_[ s ])[1] ) ) * 0.5;
  }
  else
  {
    assert ( PointEq( nextp, (*spaths_[ s ])[(int)spaths_[ s ]->size()-1] ) );  // Le: 6-8-2022

    cout << "nextp and (*spaths_[ s ])[(int)spaths_[ s ]->size()-1] at " << s << " is: " << endl;
    point_print(nextp);
    point_print((*spaths_[ s ])[0]);
    point_print((*spaths_[ s ])[1]);
    point_print ((*spaths_[ s ])[(int)spaths_[ s ]->size()-1]);
    if ( isDebug_3) PressEnterToContinue();
    //
    double spathlength_2 = 0;
    for( int i = 1; i < (int)spaths_[ s ]->size(); i ++ ){
    spathlength_2 +=  SGP_distance ( ( *spaths_[ s ])[ i - 1 ],
                                                       (*spaths_[ s ])[ i ] );
    }
    //
    sum_2 = spathlength_2;
  }
  cout << "\t\t calc dist_2 done...sum_2= " << sum_2 << endl;
  if (isDebug_3) PressEnterToContinue();
  
  return sum_1 + sum_2;
}



// Le: 20/10/2019 find prevp and nextp points to update and check traightness condition
// if prevp_nextp == true -> output is prevp point
// else output is nextp point

Point sliceMngt::findPrevpNextp(int s, int &prep_nextpFaceIdx, bool prevp_nextp){
Point prevp_or_nextp;
Face * const * faces = poly_->Faces();
bool bool_1, bool_2;
 bool_1 = bool_2 = false;
for (unsigned int i = 0; i < spaths_[s]->size(); i++){
    for ( unsigned int j = 0; j< slicefaces_[s-1].size(); j++){
        if ( faces[ slicefaces_[s-1][j] ]->isOn((*spaths_[s])[i]) ){
        bool_1 = true;
        break;
        }
    }
//
    for ( unsigned int j = 0; j< slicefaces_[s].size(); j++){
        if ( faces[ slicefaces_[s][j] ]->isOn((*spaths_[s])[i]) ){
        bool_2 = true;
        break;
        }
    }
    if (bool_1 == bool_2 && bool_1 ){
        prevp_or_nextp = (*spaths_[s])[i];
        // debug message
        cout << "\t \t prevp and nextp, prevpFaceIdx and nextpFaceIdx are determined..." << endl;
        // if (isDebug_2) PressEnterToContinue();
        break;
    }
}
    if (prevp_nextp){
    for ( unsigned int j = 0; j< slicefaces_[s-1].size()-1; j++){
        if ( faces[ slicefaces_[s-1][j] ]->isOn( prevp_or_nextp ) ){
        prep_nextpFaceIdx = slicefaces_[s-1][j];
        break;
        }
    }
    }
    else{
    for ( unsigned int j = 0; j< slicefaces_[s].size()-1; j++){
        if ( faces[ slicefaces_[s][j] ]->isOn( prevp_or_nextp ) ){
        prep_nextpFaceIdx = slicefaces_[s][j];
        break;
        }
    }
    }
   // assert ( prevp_or_nextp.z() <= slices_[s-1].front().source().z()
     //           && prevp_or_nextp.z() >= slices_[s].front().source().z() );

return prevp_or_nextp;
}

// Le: 1/11/2019 find prevp and nextp points to check traightness condition
// b_prevp_nextp =true -> output is colinear prevp point
// else output is nextp point
Point sliceMngt::findColinearPrevpNextp(int s, int &co_prevp_nextpFaceIdx, bool b_prevp_nextp){
// Face * const * faces = poly_->Faces();
Point co_prevp_or_nextp;
//
if (b_prevp_nextp){
    assert ( atVertex_[s-1] == (*spaths_[s-1])[0] );
    co_prevp_or_nextp = CGAL::midpoint( atVertex_[s-1], (*spaths_[s-1])[1] );
    co_prevp_nextpFaceIdx = slicefaces_[s-1][onSegment_[s-1]];
}
else
{
    assert ( atVertex_[s+1] == (*spaths_[s])[ (int)spaths_[s]->size()-1]);
    co_prevp_or_nextp = CGAL::midpoint( (*spaths_[s])[ (int)spaths_[s]->size()-2],
                            atVertex_[s+1] );
    co_prevp_nextpFaceIdx = slicefaces_[s+1][onSegment_[s+1]];
}
// assert ( co_prevp_or_nextp.z() <= slices_[s-1].front().source().z()
         // && co_prevp_or_nextp.z() >= slices_[s+1].front().source().z() );

return co_prevp_or_nextp;
}


/*----------------------------------------------------------------------------------*/
