/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * bbox.h -
 *     Defines an axis parallel bounding box.
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#ifndef  __BBOX__H
#define  __BBOX__H

#include  <limits.h>

class   BBox {
private:
    RT  min_coords[ 3 ];
    RT  max_coords[ 3 ];

public:
    void  init() { 
        min_coords[ 0 ] = 
            min_coords[ 1 ] = 
            min_coords[ 2 ] = RT(IT(INT_MAX));

        max_coords[ 0 ] = 
            max_coords[ 1 ] = 
            max_coords[ 2 ] = RT(IT(INT_MIN));
    }

    const RT  & minx() {
        return  min_coords[ 0 ];
    }
    const RT  & miny() {
        return  min_coords[ 1 ];
    }
    const RT  & minz() {
        return  min_coords[ 2 ];
    }

    const RT  & maxx() {
        return  max_coords[ 0 ];
    }
    const RT  & maxy() {
        return  max_coords[ 1 ];
    }
    const RT  & maxz() {
        return  max_coords[ 2 ];
    }

    void  bound( const Point  & pnt ) {
        //cout << "bounding: " << pnt << "\n";
        for  ( int  ind = 0; ind < 3; ind++ ) {
            if  ( pnt[ ind ] < min_coords[ ind ] )
                min_coords[ ind ] = pnt[ ind ];
            if  ( pnt[ ind ] > max_coords[ ind ] )
                max_coords[ ind ] = pnt[ ind ];
        }
    }
};


#else   /* __BBOX__H */
#error  Header file bbox.h included twice
#endif  /* __BBOX__H */

/* bbox.h - End of File ------------------------------------------*/
