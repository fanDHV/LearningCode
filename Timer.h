/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * Timer.h -
 *     
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#ifndef  ___TIMER__H
#define  ___TIMER__H

#include  <time.h>

/* measuring times */
#ifndef  CLOCKS_PER_SEC
#define   CLOCKS_PER_SEC  1000000.0
#endif 

//extern long clock (void);

class  Timer
{
private:
  long int    start_t, end_t;

  static inline float  SECONDS( long int   time ) 
  {
      return  ( (double)time ) / ((double)CLOCKS_PER_SEC);
  }

public:
  void  start( void )
  {
    start_t = clock();
  }

  void  end( void )
  {
    end_t = clock();
  }

  float  seconds( void )
  {
    return  SECONDS( end_t - start_t ); 
  }
};


#else   /* ___TIMER__H */
//#error  Header file timer.h included twice
#endif  /* ___TIMER__H */

/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 *     
 * Timer.h - End of File
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

