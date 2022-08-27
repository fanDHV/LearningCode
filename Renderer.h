#ifndef _RENDERER_H
#define _RENDERER_H

#include <list>
#include <string>
#include <vector>
#include <deque>

#include  "Functions.h"
//#include "Leda.h"
#include "myCGAL.h"
#include "Path.h"
#include "Face.h"
#include "Vert.h"

// for Mac OSx
//#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
//#include <GLUT/glut.h>
// for linux
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

//#define TRANSLATION_STEP 1
#define  TRANSLATION_STEP 10
#define  ROTATION_STEP 10.0
#define  PI_DIV_180 0.01745329 // for converting from degrees to rads
#define  VEC_MAG(x,y,z) sqrt(x*x+y*y+z*z)

#define  FLY_STEP              0.5 // step distance in fly bys
#define  FLY_DIST_OFF_GROUND  10 // distance we hover above ground in flyby
#define  FLY_LOOK_AHEAD        3 // how many points in the fly by path that we look
                                 // ahead

class SPath;
class SPathNode;


struct FlyByPoint
{
    float pointOnPath [3]; // actual point on path
    float unitVec [3]; // the unit vector made from the origin to the point on
    // the path
    float flyByPoint [3]; // point that we fly through
};

#define  MAX_PATHS  50

class Renderer
{
public:
    Renderer();
    ~Renderer();
    void  MainLoop(std::string filename,double epsilon,
                   std::string infilename, int argc,
                   char ** argv, bool  f_test, bool  _f_print,
                   bool  f_super_test, bool  f_white );

    // for running LiuWong algorithm 2022-8-28
    void  LiuWong(string filename,
                    double epsilon,
                    string infilename, bool _f_LiuWong);

  // has to be here for the wrapper crap
    void  display();      // function called to draw the scene
    void  reshape(int, int);      // called when window resizes
    void  keyboard(unsigned char key, int x, int y); // called when key is pressed
    void  mouse(int button,int state, int x, int y); // called when mouse is
    // clicked
    void  updateDebugDraw();
    void  flyBy(); // called when user wants to to a fly-by of the path

    void  algorithmDisplayMenu(int value);
    void  algorithmAnimateMenu(int value);
    void  otherDisplayMenu(int value);
    void  mainMenu(int value);

    void  useHershSuri(); // tells the renderer to use the hersh-suri method
    void  drawFillFace( const Vert ** vertices, int  face_ind );
    void  drawFrameFace( const Vert ** vertices, int  face_index );
    void  drawNormalFace(  const Vert ** vertices, int  face_index );

    // trang 2013/12: for handling direction key pressing
    void  handleSpecialKeypress(int key, int x, int y);

private:
    void  init();         // initialize the window, lighting etc...

    void  moveViewPoint(int axis, int dir);
    void  rotateViewDir(float angle);
    void  rotateViewPoint(float angle);

    void  drawAxes(float size); // draws 3d axes
    void  drawPolytope(GLenum drawMode); // draw the polytope either in
    // render or select mode
    void  drawPath();
    void  drawGraph();
    void  drawSphere();

    void  drawDebug();
    void  animateDebug();
    void  readSrcDestFile( const char  * infilename );

    void  generateFlyByPoints(); // generates the array of points that we will
                              // travel through in our flyby
    bool  processHits(GLint hits, GLuint buffer[]); // processes the mouse hits

    void  initMainMenu(); // initializes the menu stuff

    ///////////////////////////////////////////////////
    // private members
    ///////////////////////////////////////////////////

    // index for display list of polytope
    int  m_listIndex;

    // data for the faces of the polytope
    Face  ** m_faces;
    int  m_faceNum;
    bool  f_no_draw;
    bool  f_print;

    // the path object that we will use to read the polytope
    // and determine the path
    Path  * m_path;

    SPath  * shortest_path;
    /**********************************
      Hoai - 7/1/2013
    ***********************************/
    std::vector<Point> *multiSP_path;
    std::vector<Point> *multiSP_path_first;
    std::deque< std::deque < Seg > > *slices_;
    /**********************************/

    // trang (2015) --- for comparing with the global sdp
    SPath * globalSGP;

    SPath  * arr_paths[ MAX_PATHS ];
    int  paths_count;


    //GraphNode  * m_traversalPath; // the linked list containing the points on
    // the shortest path

    // the filename of the polytope
    std::string  m_filename;

    float   m_viewPos [4]; // array that holds coords of camera position
    // 4th float is used when we set light position
    // using view position
    float   m_viewDir [3]; // coords of what point camera points at
    float   m_viewUp  [3]; // coords of vector saying which way is up

    float   m_tempVec [3];
    float   m_viewPointDistance; // distance between viewPos and viewDir

    // info concerning the source and destination points of the path
    int   m_srcFace;
    int   m_destFace;
    double   m_srcPoint[3];
    double   m_destPoint[3];

    float   m_length;

    // how many faces are selected
    int   m_facesSelected;

    // the dimensions of the current window
    int   m_windowWidth;
    int   m_windowHeight;

    GLushort * m_viewDepthBuffer; // pointer to the depth buffer

    // info for how much the polytope is currently rotated
    float   m_xRot; // how much rotated around x axis
    float   m_yRot; // how much rotated around y axis

    bool   m_drawNormals; // do we draw the normals?
    bool   m_drawGraph; // do we draw the graph?
    bool   m_drawPath; // do we draw the path?
    bool   m_drawPolytope; // do we draw polytope?
    bool   m_drawPolytopeFilled; // do we draw the polytope filled or lined?
    bool   m_drawAxes; // do we draw axes?

    // fly by stuff
    bool   m_flyBy; // are we in flyBy mode?
    SPathNode * m_currFlyByNode; // the current node in the path that we are

    // flying by
    float   m_currFlyBySrc [3]; // current source point
    float   m_currFlyByDest [3]; // current destination point
    float   m_currFlyByUnitVec [3]; //unit vector made from current source
                                // to current destination
    float   m_flyByDist; // distance between current source and destination
    float   m_currFlyByDist; // how far we've travelled so far from our source
    // to the destination
    float   m_unitVecSrc[3]; // unit vector made from current source to origin
    float   m_unitVecDest[3]; // unit vector made from current dest to origin

    FlyByPoint * m_flyByPoints; // list of points that are on our fly by path
    int   m_numFlyByPoints; // number of flyby points in the array
    int   m_currFlyByPoint; // the current index into the flyByPoints array that
    // holds the point that we are currently on

     // debugging the algorithm stuff
    std::list<AnimateVecInfo *> m_debugPoints; // holds the list of vectors to draw
    bool   m_animateDebug;
    bool   m_drawDebug;
    int   m_debugListIndex;
    bool   m_animatingDebug;

    // more debugging of the algorithm stuff, this is to show what parts of
    // the algothm to display
    bool   m_showLong;
    bool   m_showReg;
    bool   m_showNN;

    // for hersh-suri
    bool   m_useHershSuri;

    /* ------------------- trang 2013/12 ------------------------*/
    float m_transDis[3];
    float scale;
    bool problem_solved;
};

#endif

