/*--------------------------------------------------------------------
 * cutting - compute palanr cuttings
 * Copyright (C) 1998 Sariel Har-Peled
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA  02111-1307, USA.
\*--------------------------------------------------------------------*/

/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * debug.h -
 *    Implement a string logging into a file. It also supply a mechanisim
 *  for entry/exit of functions. Debugging is performed only if
 *  the environment variable MMD_DEBUG is defined.
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#ifndef  __DEBUG__H
#define  __DEBUG__H

//#ifndef  bool
//typedef  int  bool;
//#endif  // bool

#ifndef  TRUE
#define  TRUE  1
#endif 

#ifndef  FALSE
#define  FALSE  0
#endif 


/*======================================================================
 * FunctionExitLog - Reports when a function is being exited.
\*======================================================================*/
struct  FunctionDebugDef {
    char      * fileName;
    int         defLine;
    char      * funcName;
    bool        fLog;
    int         currLine;

    void  DEBUG_lprintf( const char *frmt,... );

/*
    FunctionDebugDef( char  *file, int   line,
                      char   * func, bool fLogExit = FALSE )
    {
        fileName = file;
        funcName = func;
        defLine = line;
        fLog= fLogExit;
    }
*/
};
#define  dprintf  \
      func_pos_def.currLine = __LINE__; \
      func_pos_def.DEBUG_lprintf

#define  FEL_ZERO_SIZE  10
class  FunctionExitLog
{
    FunctionDebugDef        * pDef;
    char                      zero[ FEL_ZERO_SIZE ];

public :
     FunctionExitLog(
        FunctionDebugDef        * pDefinition );
     ~FunctionExitLog();
};


/*======================================================================
 * Defines some basic debug macros
\*======================================================================*/
#ifdef __GNUG__
static char       __file__[];
static char       _func_name[];
#else
extern char       * __file__;
extern char       * _func_name;
#endif  // __GNUG__

//#define  SOURCE__FILE              static char *__file__ = __FILE__
#ifndef  __GNUG__
#define  SOURCE__FILE     static char  __file__[] \
                                         = __FILE__; \
                          static char  _func_name[] \
                                         = __FILE__ 
#else
#define  SOURCE__FILE     static char __file__[] \
                                         = __FILE__; \
                          static char _func_name[]\
                                         = __FILE__ 
#endif  // __GNUG__

#define  __CALLER__                __file__, __LINE__
#ifndef  USE_ARGUMENT
   #define  USE_ARGUMENT( X )         ((void)X)
   #define  USE_VARIABLE( X )         ((void)&(X));
#endif 
#define  __CALLER_PARAMS__         char    * callingFile, int  callingLine
#define  __CALLER_PARAMS_PASS__    callingFile, callingLine

#define __POS__             __CALLER__, _func_name
#define __NULL_POS__        NULL, __LINE__, NULL
#define __REAL_POS__        __FILE__, __LINE__, NULL
#define __POS_PARAMS__      __CALLER_PARAMS__, char  * callingFunc
#define __POS_PASS__        __CALLER_PARAMS_PASS__, callingFunc

#define FUNCTION            FUNC_NAMEX( __PRETTY_FUNCTION__ )

class  DEBUG_dumper {
public:
    char  * str;

    DEBUG_dumper( char  * _str ) {
        str = _str; 
        printf( "[%s] called\n", str );
        fflush( stdout );
    }
    ~DEBUG_dumper() { 
        printf( "[%s] returned\n", str );
        fflush( stdout );
    }
};

#define FUNCTIONPR          DEBUG_dumper  __ddd( __PRETTY_FUNCTION__ )

                                         //#define  NOASSERT 

#ifdef  NOASSERT
#define FUNC_NAMEX( X )      static char *_func_name = __PRETTY_FUNCTION__; \
                            static FunctionDebugDef func_pos_def = \
                                         { __POS__, FALSE, 0 }; 
#define FUNC_NAME( X )      static char _func_name[ sizeof( X) + 1 ] = X; \
                            static FunctionDebugDef func_pos_def = \
                              { __POS__, FALSE, 0 }; 
#else
#define FUNC_NAMEX( X )      static char *_func_name = __PRETTY_FUNCTION__; \
                            static FunctionDebugDef func_pos_def = \
                                   { __POS__, FALSE, 0 }; \
                            FunctionExitLog   exitLogging( &func_pos_def );
#define FUNC_NAME( X )      static char _func_name[ sizeof( X) + 1 ] = X; \
                            static FunctionDebugDef func_pos_def = \
                                   { __POS__, FALSE, 0 }; \
                            FunctionExitLog   exitLogging( &func_pos_def );
#endif
#define FUNC_NAME_LOG( X )  static char _func_name[ sizeof( X) + 1 ] = X; \
                            static FunctionDebugDef func_pos_def = \
                                   { __POS__, TRUE, 0 }; \
                            FunctionExitLog   exitLogging( &func_pos_def );
#define FUNC_NAME_RSC( X )  FUNC_NAME( X ) \
                            ResourceManager         localManager;


/*======================================================================
 * Global Variables.
\*======================================================================*/
extern bool  DEBUG_fDebug;
extern bool  DEBUG_fTrace;
extern bool  DEBUG_fMem;
extern bool  DEBUG_fFiles;
extern bool  DEBUG_fFTE;


/*======================================================================
 * Function declerations
\*======================================================================*/
/* this function are used by the engine */
void DEBUG_init( void );
void DEBUG_exit( void );


/* exported functions */
void  DEBUG_GenLogComment(__CALLER_PARAMS__,char *title,char *comment,bool ishard);
void  DEBUG_lprintf( __POS_PARAMS__,const char *frmt,... );
void  DEBUG_vprintf( __POS_PARAMS__,
                                   const char   * format,
                                   va_list        arglist );
void          DEBUG_hard_lprintf( __POS_PARAMS__,const char *frmt,... );
void         DEBUG_logError( __POS_PARAMS__,int errorCode);
void         DEBUG_assert_fail(__POS_PARAMS__, char *message);
int          DEBUG_trackError( __POS_PARAMS__, int errorCode );
void         DEBUG_printStack( void );
void         DEBUG_position( __POS_PARAMS__ );
void         DEBUG_SystemPanic( __POS_PARAMS__, char * str );
inline void  DUMP_STACK( void ) { DEBUG_printStack(); }
int          stricmp( char     * pA, char      * pB );

/*======================================================================
 * Macros
\*======================================================================*/

#define DEBUG_LogComment(x,y) { if (DEBUG_fDebug) DEBUG_GenLogComment(__CALLER__,x,y,FALSE); }
#define DEBUG_LogCommentHard(x,y) { if (DEBUG_fDebug) DEBUG_GenLogComment(__CALLER__,x,y,TRUE); }

#define   PROPOGATE( X )      DEBUG_trackError( __POS__, X)


#ifndef  NOASSERT
   #define     DASSERT( COND, TEXT ) \
                    ((COND) ? (void)0 : DEBUG_assert_fail( __POS__, TEXT ))
   #ifndef  ASSERT
      #define     ASSERT( COND )  \
                    ((COND) ? (void)0 : DEBUG_assert_fail( __POS__, #COND ))
   #endif // ASSERT
   #define    PANIC( X )  DEBUG_SystemPanic( __POS__, X );
#else
   #define     DASSERT( X, Y )             {}
   #ifndef  ASSERT
       #define     ASSERT( COND)               {}
   #endif 
#endif

#define   LOG_POSITION()   DEBUG_position( __POS__ );

#define   CHECK( X )     { \
                            int   _r; \
                                           \
                            if  ( ( _r = X ) < 0 ) \
                               {                      \
                                  DEBUG_logError( __POS__, _r ); \
                                  return  _r;          \
                               }                      \
                          }
#define   _CHECK( X )     { \
                            int   __r; \
                                           \
                            if  ( ( __r = (X) ) < 0 )  return  __r; }
#ifndef  CHECK_NULL
#define   CHECK_NULL( X ) { \
                              void     * _ptr = (void *)(X);   \
                              \
                              if  ( _ptr == NULL ) {                              \
                                  DEBUG_logError( __POS__, ERR_NOTENOUGHMEMORY ); \
                                  return   ERR_NOTENOUGHMEMORY; \
                              }\
                          }
#endif 
#define   CHECK_MALLOC( X, Y ) { \
              void     * _ptr = (void *)MALLOC( Y );\
              \
              if  ( _ptr == NULL ) {                              \
                  DEBUG_logError( __POS__, ERR_NOTENOUGHMEMORY ); \
                  return   ERR_NOTENOUGHMEMORY; \
              }\
              *((char **)&(X)) = (char *)_ptr; \
          }
#define   CHECK_LOCAL_MALLOC( X, Y ) { \
              void     * _ptr = (void *)localManager.malloc( __POS__, Y );\
              \
              if  ( _ptr == NULL ) {                              \
                  DEBUG_logError( __POS__, ERR_NOTENOUGHMEMORY ); \
                  return   ERR_NOTENOUGHMEMORY; \
              }\
              *((char **)&(X)) = (char *)_ptr; \
          }
#define   LOCAL_MALLOC( X, Y ) { \
              *((char **)&(X)) = (char *) \
                      ((void *)localManager.malloc( __POS__, Y ));\
          }
#define   CHECK_DISPLAY_LOCAL_MALLOC( X, Y ) { \
              void     * _ptr = (void *)localManager.malloc( __POS__, Y );\
              \
              if  ( _ptr == NULL ) {                              \
                  DisplayERROR( ERR_NOTENOUGHMEMORY ); \
                  return; \
              }\
              *((char **)&(X)) = (char *)_ptr; \
          }
#define   CHECK_DISPLAY_MALLOC( X, Y ) { \
              void     * _ptr = (void *)MALLOC( Y );\
              \
              if  ( _ptr == NULL ) {                              \
                  DisplayERROR( ERR_NOTENOUGHMEMORY ); \
                  return; \
              }\
              *((char **)&(X)) = (char *)_ptr; \
          }

char            * strupr( char      * str );
void      myExit( int   code );

#endif   /* __DEBUG__H  */

/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 *
 * debug.h - end of file
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/
