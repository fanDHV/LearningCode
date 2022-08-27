/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * history.h -
 *     
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#ifndef  __HISTORY__H
#define  __HISTORY__H

#define  HIST_SIZE  5

struct  history_t 
{
    Point   hist[ HIST_SIZE ];
    int  count;

    void  init( Point  & pnt ) {
        hist[ 0 ] = pnt;
        count = 1;
    }
    void register_pt( Point  & pnt ) {
        if  ( count < HIST_SIZE ) {
            hist[ count ] = pnt;
            count++;
            return;
        }
        for  ( int  ind = 0; ind < ( HIST_SIZE - 1 ); ind++ ) 
            hist[ ind ] = hist[ ind + 1 ];
        hist[ HIST_SIZE - 1 ] = pnt;
    }

    bool  isIn( const Point  & pnt ) const {
        for  ( int  ind = 0; ind < count; ind++ ) 
            if  ( hist[ ind ] == pnt ) 
                return  true;

        return  false;
    }
};

#else   /* __HISTORY__H */
#error  Header file history.h included twice
#endif  /* __HISTORY__H */

/* history.h - End of File ------------------------------------------*/
