/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * spath.h -
 *     Holds the shortest path.
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#ifndef  __SPATH__H
#define  __SPATH__H

class Sweeper;

class  SPathNode {
private:
    NNLocation  nnloc;
    Point  pnt;
    GraphNode  * gn;

    SPathNode  * p_next, * p_prev;
    bool  f_in_sky;

    void  rewrite_path( Sweeper  & sw, double  val,
                        bool  f_reverse );

public:
    SPathNode( GraphNode  * _gn );
    SPathNode( const NNLocation  & _nnloc, const Point   & pnt );
    SPathNode( SPathNode  * _gn );

    SPathNode( Face  * f, const Point  & pnt );

    NNLocation  & loc() {
        return  nnloc;
    }
    void   setLocation( const NNLocation  & loc ) {
        nnloc = loc;
    }
    int  length() const;
    SPathNode  * next() {
        return  p_next;
    }
    SPathNode  * prev() {
        return  p_prev;
    }

    void  append( SPathNode  * node );

    const Point  & point() {
        return  pnt;
    }
    void  dump();

    void  shortcut_on_faces();
    int   getFaceAdjDegree();
    Face  * getAdjFace( int  ind );
    SPathNode  * getLastPointAdjToFace( Face  * f );
    void  shortcut_on_edge_inner( Edge  * e, Poly  * pPoly );
    void  shortcut_on_edge( Poly  * pPoly );
    void  shortcut_on_vertex();
    bool  isOnVertex()  const;

    bool  isInSky() const {
        return  f_in_sky;
    }
    void  setInSky( bool  _f_in_sky ) {
        f_in_sky = true;

        // trang: avoiding unused variable
        _f_in_sky = !_f_in_sky;
    }
    friend class SPath;
};


class  SPath {
private:
    SPathNode  * list_first, * list_last;
    Poly  * p_poly;
// Le need to exchange private -> public
    double   get_length( SPathNode  * start,
                         SPathNode  * end );
    double  lift_length( SPathNode  * start,
                         SPathNode  * end );
    void  shortcut_by_lifting( SPathNode  * start,
                               SPathNode  * end );

public:
    SPath();
    SPath( SPath  * ptr );
    SPath( GraphNode  * first, Poly  * _p_poly );
    SPath( SPathNode  * first, Poly  * _p_poly );

    SPathNode  * first() {
        return  list_first;
    }
    /* Le exchange preivate -> public
    double   get_length( SPathNode  * start,
                         SPathNode  * end ); */

    void  pushTop( GraphNode  * gn );
    void  pushTop( const NNLocation  & nnloc, const Point  & pnt );
    void  pushTop( SPathNode  * gn );

    void  optimize();

    void  init( SPathNode  * gn, Poly  * _p_poly );
    void  init( GraphNode  * gn, Poly  * _p_poly );
    void  init_first( GraphNode  * gn, Poly  * _p_poly );
    void  init_first( SPathNode  * gn, Poly  * _p_poly );
    void  dump();
    void  append( SPathNode  * node );
    double  getAprxLength() const;
    /********Le 13/10/2019 purpose of calculate the length of SPath connecting prevp point to next point
    double getSubpathLength (const Point &prevp, const Point &nextp) const; */

    // trang (2015/08/30): print the list of the path ... purpose of printing global path
    void  printNodeList () const;

    void  shortcut_on_faces();
    void  shortcut_on_edge();
    void  shortcut_on_vertex();
    void  shortcut_by_relifting();
    void  shortcut_by_relifting_with_dist( int  dist );
    void  shortcut_by_relifting_jump();
    void shortcutify_by_face( Face  * f,
                              SPathNode  * curr,
                              SPathNode  * last );
    void  shortcut_by_face( Face  * f );
    void  shortcut_brute_force();


    void  lift_and_shortcut( SPathNode  * start,
                             SPathNode  * end,
                             SPath  & sp );
    SPath  * convertToSPath( GraphNode  * path );

    void  projectWedgePath( SPathNode  * start,
                            Point  & mid,
                            SPathNode  * end );
    Edge  * findExitEdge( NNLocation  & nnloc,
                          Plane  & pathPlane, Plane  & hPositive,
                          history_t  & hist, Point  & exit_point,
                          NNLocation  & dest_loc );

    SPath  * shortcut();
};

Plane  getPositivePlane( Point  & v0, Point  & v1,
                         Point  & v2 );

#else   /* __SPATH__H */
#error  Header file spath.h included twice
#endif  /* __SPATH__H */

/* spath.h - End of File ------------------------------------------*/
