/*
  Authors: P.T. An, T.V. Hoai, N.N. Hai, Trang
  Creation date: 6 Jan 2013
 */
#ifndef __SLICE_MNGT_H_
#define __SLIDE_MNGT_H_

#include "Poly.h"
#include "Path.h"
#include "spath.h"

class sliceMngt {
 public:
  sliceMngt( Poly *poly );
  sliceMngt( Poly *poly, SourceDest *source );
  ~sliceMngt();

  // detect planes within range between [source.z,destination.z]
  void detectEffectivePlanes( SourceDest *srcdest );
  // detect the intersection between Agarwal's path with slices
  void detectSlicePaths( SPath* sp );
  void detectPoorSlicePaths( SPath* sp );
  void vibrateSlicePaths( SourceDest *srcdest,
                          int srcFaceIdx, int destFaceIdx,
                          Path *spalg,
                          double ratio = 0.5, int nedges = 0 );

  // improve the path by multiple shooting method
  vector<Point>* multiShooting( SourceDest *srcdest, int srcFaceIdx, int destFaceIdx,
				Path *spalg );

  // return shortest path
  vector<Point> *getShootingPath();
  // Le: 19_12_2019: make path smooth 
  void getSmoothPath();
  // get slices
  void getSlices( deque< deque< Seg > > *&myslices );

  //Point myLineIntersect( Point prevp, Point nextp, int s);      // previous, next points
                                                                       // on current path

  // trang: initialize a path through slices
  void initializePath( SourceDest *srcdest,
                          int srcFaceIdx, int destFaceIdx,
                          Path *spalg,
                          double RATIO = 0.5, int nedges = 0 );

  // trang: check if there exists horizontal paths
  // bool existHorizontalPath_; // moved to global setting ...

 private:
  // convergent control function
  double convF( int n, int m, double epsilon );
  // get approximate length
  double getAprxLength() const;
  
  

 /***************************Le (12/10/2019) update for all cases***************/
  void computeUpdateSGPInterior( Point &prevp, Point &nextp,         // previous, next points
                                                                          // on current path
                                     int &prevpFaceIdx, int &nextpFaceIdx, // indices of faces containing prevp, nextp
                                                                          // on the slice
                                      int s,                              // considered slice
                                      int iteration,                      // iteration of alg.
                                      Point &newPoint, int &newSeg, Path *spalg);


  // Le: for shortest gentle paths with new way to check SGP - used
  bool checkSGPCollinearConditionInterior(SourceDest *srcdest, Path *spalg, int srcFaceIdx, int destFaceIdx);
  // Le: for finding intersection between spath and slices[s]
  Point Intersection_spath_slices( SPath *spaths_, int s, int &newUpSeg, bool isSeg);

  /*** Le: find Faces index containing pnt ***13/10/2019 ***/
  // return a array of indices of all faces containing a given point
  list<int> findFaceIndx(const Point &pnt);
  bool empty_intersection(const list<int>& x, const list<int>& y, int & commonIdx);
  /****Le 13/10/2019 purpose of calculate the length of SPath connecting prevp to nextp point*/
  double getSubpathLength (const Point &prevp, const Point &nextp, int s) const;
  // Le: 20/10/2019 find prevp and nextp points to update and check traightness condition
  Point findPrevpNextp(int s, int &prep_nextpFaceIdx, bool prevp_nextp);
  // Le: 1/11/2019 find prevp and nextp points to check traightness condition
  Point findColinearPrevpNextp(int s, int &co_prevp_nextpFaceIdx, bool b_prevp_nextp);

  Vert * v_source;

 private:
  Poly *poly_;
  deque< deque< Seg > > slices_;
  deque< deque< int > > slicefaces_;
  // Le: add to check list of Face's Indices of vertices of polygon
  // deque< deque< list < int > > > listIdx_slicefaces_;
   //
  deque< deque< pair< Edge*, Edge* > > > slicelaterals_;
  Point * atVertex_;
  int *onSegment_;
  deque< Point > ** spaths_;
  double *spathlens_;
  double *movedirs_;

  //
  int fromSlice_, toSlice_, fromSliceSource_;
  // fromSlice_ -1 passing through the source point $p$
  // toSlice_  passing through  the dest point $q$


  // ----------------------------------------------------------
  // trang: for using straightest geodesic ...
  // ----------------------------------------------------------
  deque < double > vertex_total_angles_;
  Vert * const * sorted_verts_;
  std::vector< Triangle > triBetween_right_;
  std::vector< Triangle > triBetween_left_;
  bool right_angle_;
};

void createSlices( Poly *poly, SourceDest *srcdest );

#endif
