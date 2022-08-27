//#include  <vector>
#include <fstream>
#include  <stdio.h>
#include  <memory.h>
#include  <assert.h>
#include  <stdlib.h>

#include  "sarray.h"
#include  "Renderer.h"
#include  "rgeneric.h"
#include  <math.h>
#include  "spath.h"

extern void RendererDisplayWrapper ();
extern void RendererReshapeWrapper (int w, int h);
extern void RendererKeyboardWrapper (unsigned char key, int x, int y);
extern void RendererMouseWrapper (int button, int state, int x, int y);
extern void RendererFlyByWrapper ();
extern void RendererUpdateDebugDrawWrapper ();
extern void RendererHandleSpecialKeypressWrapper(int key, int x, int y);

extern void RendererAlgorithmDisplayMenuWrapper (int value);
extern void RendererAlgorithmAnimateMenuWrapper (int value);
extern void RendererOtherDisplayMenuWrapper (int value);
extern void RendererMainMenuWrapper (int value);
static bool f_solid = false;
static bool f_display_selected = true;
static bool f_draw_sphere = false;
static bool f_all_lights = true;  //false; //true;
static bool f_fly_mode = false;
static bool f_white = false;

using namespace std;

Renderer::Renderer ()
{
    // trang: init translation distances
    m_transDis[0] = -50;
    m_transDis[1] = -10;
    m_transDis[2] = -150;

    scale = 0.85;

    f_no_draw = false;
    m_faces = NULL;
    m_path = NULL;
    m_listIndex = 0;

    memset (arr_paths, 0, sizeof (arr_paths));
    paths_count = 0;

    m_viewPos[0] = 0;
    m_viewPos[1] = 0;
    m_viewPos[2] = 1;

    m_viewDir[0] = 0;
    m_viewDir[1] = 0;
    m_viewDir[2] = 0;

    m_viewUp[0] = 0.0;
    m_viewUp[1] = 1.0;
    m_viewUp[2] = 0.0;

    m_viewPointDistance = sqrt (pow (m_viewDir[0] - m_viewPos[0], 2.0) +
                                pow (m_viewDir[2] - m_viewPos[2], 2.0));

    m_facesSelected = 0;
    m_srcFace = -1;
    m_destFace = -1;

    //m_xRot = 25;
    m_xRot = 0;
    m_yRot = 0;
    //m_yRot = -45;

    m_windowWidth = 1000; // old: 800
    m_windowHeight = 800; // old: 600

    m_viewDepthBuffer = NULL;
    //m_traversalPath = NULL;
    shortest_path = NULL;
    /******** Hoai - 7/1/2013 **************/
    multiSP_path = NULL;
    multiSP_path_first = NULL;
    slices_ = NULL;
    /***************************************/

    m_drawNormals = false;
    m_flyBy = false;
    m_drawGraph = false;
    m_drawPath = true;
    m_drawPolytope = true;
    m_drawPolytopeFilled = true;
    m_drawAxes = false;

    m_currFlyByNode = NULL;
    m_currFlyByDist = 0;
    m_flyByDist = 0;

    m_length = 0.0;

    m_flyByPoints = NULL;
    m_numFlyByPoints = 0;
    m_currFlyByPoint = -1;

    m_drawDebug = false;
    m_debugListIndex = 0;
    m_animatingDebug = false;
    m_animateDebug = false;

    m_showLong = false;
    m_showReg = false;
    m_showNN = false;

    m_useHershSuri = false;


}

Renderer::~Renderer ()
{
  if ( multiSP_path != NULL )
    delete multiSP_path;
  if ( multiSP_path_first != NULL )
    delete multiSP_path_first;
  if ( slices_ != NULL )
    delete slices_;
}


void
Renderer::init ()
{
  if (f_fly_mode)
    {
      glClearColor (0.7, 0.8, 1, 0);
    }
  else
    {
      if (f_white)
        glClearColor (1.0, 1.0, 1.0, 0.0);  // Background color = White
      else if (f_print)
        glClearColor (0.7, 0.7, 1.0, 0.0);
      else
        glClearColor (0.1, 0.1, 0.1, 0.0);
    }
  glShadeModel (GL_SMOOTH);     // Smooth shading

  // lighting defaults enabled below

  GLfloat mat_specular[] = { 0.5, 0.5, 0.5, 1 };
  GLfloat mat_shininess[] = { 0.5 };
  //GLfloat light1_position[] = {2000,0,0,1};

  GLfloat light0_position[] = { 0, 200, 0, 1 };
  GLfloat light1_position[] = { -200, -200, 1000, 1 };
  GLfloat light2_position[] = { 200, 200, 0, 1 };
  GLfloat light3_position[] = { 0, 0, -20000, 1 };

  GLfloat light_ambient[] = { 0.12, 0.12, 0.12, 1.0 };
  //GLfloat light_direction[] = {0,1,0};

  //GLfloat light0_direction[] = {-1,-1,-1};
  //GLfloat light1_direction[] = {1,0,0};
  //GLfloat light2_direction[] = {0,1,0};
  //GLfloat light3_direction[] = {0,0,1};

  GLfloat light0_diffuse[] = { 0.90 / 2, 1.00 / 2, 0.90 / 2, 1.0 };
  GLfloat light1_diffuse[] = { 1.00 / 2, 0.90 / 2, 0.90 / 2, 1.0 };
  GLfloat light2_diffuse[] = { 0.20, 0.20, 0.20, 1.0 };
  GLfloat light3_diffuse[] = { 0.20, 0.20, 0.02, 1.0 };

  //trang: original one
  //double param = 0.1;
  float param = 0.1; // remove some warnings for C++11 compiler ...

  if (f_print && (!f_white))
    param = 0.3;
  else
    param = 0.5;
  GLfloat lmodel_ambient[] = { param, param, param, 1.0 };

  glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
  glMaterialfv (GL_BACK, GL_SPECULAR, mat_specular);
  //glMaterialfv(GL_BACK, GL_SHININESS, mat_shininess);

  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, light_ambient);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, light_ambient);

  glColorMaterial (GL_FRONT_AND_BACK, GL_DIFFUSE);

  //glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light0_direction);
  glLightfv (GL_LIGHT0, GL_POSITION, light0_position);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, light0_diffuse);

  glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv (GL_LIGHT0, GL_SPECULAR, mat_specular);

  glLightfv (GL_LIGHT1, GL_POSITION, light1_position);
  //glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, light1_direction);
  glLightfv (GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv (GL_LIGHT1, GL_AMBIENT, light_ambient);
  glLightfv (GL_LIGHT1, GL_SPECULAR, mat_specular);

  //  glLightf(GL_LIGHT0,GL_SPOT_CUTOFF,150.0);
  glLightfv (GL_LIGHT2, GL_POSITION, light2_position);
  //glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light2_direction);
  glLightfv (GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
  glLightfv (GL_LIGHT2, GL_AMBIENT, light_ambient);

  glLightfv (GL_LIGHT3, GL_POSITION, light3_position);
  //glLightfv(GL_LIGHT3, GL_SPOT_DIRECTION, light3_direction);
  glLightfv (GL_LIGHT3, GL_DIFFUSE, light3_diffuse);
  glLightfv (GL_LIGHT3, GL_AMBIENT, light_ambient);
  glLightfv (GL_LIGHT3, GL_SPECULAR, mat_specular);

  glLightModelfv (GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);

  if (f_all_lights)
    {
      glEnable (GL_LIGHT1);
      //glEnable(GL_LIGHT2);
      //glEnable(GL_LIGHT3);
    }

  glEnable (GL_DEPTH_TEST);
  glDepthMask (GL_TRUE);

  glEnable (GL_COLOR_MATERIAL);
  glEnable (GL_NORMALIZE);

  // glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  glEnable (GL_CULL_FACE);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}


void
Renderer::drawFillFace (const Vert ** vertices, int face_index)
{
  double transparency = 1.0;
  int pt;

  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  glBegin (GL_POLYGON);


  // if we're drawing the path, then set our faces to be transparent
  if (shortest_path && (!f_solid))
    {
      if (!m_drawGraph)
        transparency = 0.7;
    }

  // if this face is either the source or destination face,
  // then set its color to red

  if (f_display_selected
      && ((m_srcFace == face_index) || (m_destFace == face_index)))
    //glColor4f (1.0, 0.0, 0.0, transparency);
    glColor4f (0.0, 0.0, 1.0, transparency);
  else
    glColor4f (0.6, 0.6, 0.7, transparency);

  //first, set the normal vector for the face
  //glNormal3f(m_faces[ face_index ]->UnitNormal().to_vector()[0],
  //           m_faces[ face_index ]->UnitNormal().to_vector()[1],
  //           m_faces[ face_index ]->UnitNormal().to_vector()[2]);
  glNormal3f (m_faces[face_index]->UnitNormal ().x ().to_double (),
              m_faces[face_index]->UnitNormal ().y ().to_double (),
              m_faces[face_index]->UnitNormal ().z ().to_double ());


  for (pt = 0; pt < m_faces[face_index]->VertDegree (); pt++)
    {
      glVertex3f (vertices[pt]->Location ().x ().to_double (),
                  vertices[pt]->Location ().y ().to_double (),
                  vertices[pt]->Location ().z ().to_double ());
    }
  //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnd ();
}


void
Renderer::drawFrameFace (const Vert ** vertices, int face_index)
{
  GLfloat bcolors[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  };
  glColorPointer (3, GL_FLOAT, 0, bcolors);

  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
  glLineWidth (0.1);
  glBegin (GL_POLYGON);
  //glColorPointer(3,GL_FLOAT, 0, bcolors);

  glColor4f (0.0, 0.0, 0.0, 1.0);

  //first, set the normal vector for the face
  glNormal3f (1.0, 1.0, 1.0);


  for (int pt = 0; pt < m_faces[face_index]->VertDegree (); pt++)
    {
      glVertex3f (vertices[pt]->Location ().x ().to_double (),
                  vertices[pt]->Location ().y ().to_double (),
                  vertices[pt]->Location ().z ().to_double ());
    }

  glEnd ();
}


void
Renderer::drawNormalFace (const Vert ** vertices, int face_index)
{
  //if we want to draw the normal for this face
  if (!m_drawNormals)
    return;

  glColor3f (1.0, 1.0, 0.0);

  // now do some work to draw the normal vectors of the face
  // we draw the normals out from the approxiamte centers
  // of the face (it doesn't really matter)
  // easy way to find center of face is to take midpoint
  // of two vertices and then find midpoint of this point
  // and third vertex
  //float center[3];
  //float midpoint[3];
  RT center[3];
  RT midpoint[3];

  // first find midpoint of two vertices
  midpoint[0] = (vertices[1]->Location ().x () +
                 vertices[0]->Location ().x ()) / 2;
  midpoint[1] = (vertices[1]->Location ().y () +
                 vertices[0]->Location ().y ()) / 2;
  midpoint[2] = (vertices[1]->Location ().z () +
                 vertices[0]->Location ().z ()) / 2;

  // now find the center from the midpoint and the third
  // vertex
  center[0] = (vertices[2]->Location ().x () + midpoint[0]) / 2;
  center[1] = (vertices[2]->Location ().y () + midpoint[1]) / 2;
  center[2] = (vertices[2]->Location ().z () + midpoint[2]) / 2;

  // now draw the normal
  glLineWidth (0.1);
  glBegin (GL_LINES);

  glVertex3f (center[0].to_double (), center[1].to_double (),
              center[2].to_double ());
  /*
     glVertex3f(center[0]+
     10*m_faces[ face_index ]->UnitNormal().to_vector()[0],
     center[1]+
     10*m_faces[ face_index ]->UnitNormal().to_vector()[1],
     center[2]+
     10*m_faces[ face_index ]->UnitNormal().to_vector()[2]);
   */
  glVertex3f ((center[0] +
               10 * m_faces[face_index]->UnitNormal ().x ()).to_double (),
              (center[1] +
               10 * m_faces[face_index]->UnitNormal ().y ()).to_double (),
              (center[2] +
               10 * m_faces[face_index]->UnitNormal ().z ()).to_double ());

  glEnd ();
}


// function that generates the polytope display list
void
Renderer::drawPolytope (GLenum drawMode)
{
  int i;                        //, pt;

  glEnable (GL_CULL_FACE);
  glEnable (GL_DEPTH_TEST);

  // go through polytope faces and draw them
  for (i = 0; i < m_faceNum; i++)
    {
      // get the vertices for this face
      const Vert *const *vertices = m_faces[i]->Verts ();

      // if we're in mouse selection mode
      if (drawMode == GL_SELECT)
        {
          // replace name on top of stack with this face index
          glLoadName (i);
        }

      drawFillFace ((const Vert **) vertices, i);
      //drawFrameFace( (const Vert **)vertices, i );
      drawNormalFace ((const Vert **) vertices, i);
    }

  glDisable (GL_DEPTH_TEST);
  glDisable (GL_CULL_FACE);

  // now disable the culling of our faces
  if ((drawMode == GL_SELECT) || (shortest_path != NULL))
    if (!f_solid)
      glDisable (GL_CULL_FACE);
}


void
Renderer::drawGraph ()
{
  //const GraphNode * const * graph = m_path->ReturnGraph();
  const GridGraph & grid (m_path->getGridGraph ());
  int numNodes = grid.getVerticesNum ();
  int i;

  //glEnable(GL_CULL_FACE);

  //glEnable(GL_DEPTH_TEST);


  glLineWidth (1.0);
  for (i = 0; i < numNodes; i++)
    {
      //printf( "drawGraph: %d\n", i ); fflush( stdout );
      //      if(i==0||i==numNodes-1)
      //  glColor3f(1.0,1.0,0.0);
      glColor3f (0.0, 1.0, 0.0);

      glPushMatrix ();
      glTranslatef (grid.getVertex_DoubleXCoord (i),
                    grid.getVertex_DoubleYCoord (i),
                    grid.getVertex_DoubleZCoord (i));
      //printf( " (%g, %g, %g)\n",
      //        grid.getVertex_DoubleXCoord( i ),
      //        grid.getVertex_DoubleYCoord( i ),
      //        grid.getVertex_DoubleZCoord( i ) );


      glutSolidSphere (2, 10, 10);
      glPopMatrix ();

      glColor3f (0.3, 0.9, .5);

      glBegin (GL_LINES);

      //      draw the point (with a sphere or something) graph[i]->m_p
      //NodeSet::iterator si;
      const GraphNode & gn (*(grid.getNode (i)));
      for (int si = 0; si < (int) gn.neighbors.size (); si++)
        {
          //printf( "si: %d\n", si ); fflush( stdout );
          const GraphNode & gn_trg (*(grid.getNode (gn.neighbors[si].m_adj)));

          if (gn.m_p == gn_trg.m_p)
            continue;

          //point_print( gn.m_p );
          //point_print( gn_trg.m_p );

          glVertex3f ((float) gn.m_p.x ().to_double () * 1.01,
                      (float) gn.m_p.y ().to_double () * 1.01,
                      (float) gn.m_p.z ().to_double () * 1.01);

          glVertex3f ((float) gn_trg.m_p.x ().to_double () * 1.01,
                      (float) gn_trg.m_p.y ().to_double () * 1.01,
                      (float) gn_trg.m_p.z ().to_double () * 1.01);

        }

      //printf( "before glEnd!\n" ); fflush( stdout );
      glEnd ();
    }

  glLineWidth (1.0);
  //glDisable(GL_CULL_FACE);

  //glDisable(GL_DEPTH_TEST);
}


void
Renderer::drawSphere ()
{
  //const GraphNode * const * graph = m_path->ReturnGraph();
  const GridGraph & grid (m_path->getGridGraph ());
  int numNodes = grid.getVerticesNum ();
  int i;


  glLineWidth (1.0);
  for (i = 0; i < numNodes; i++)
    {
      glColor3f (0.0, 1.0, 0.0);

      GraphNode *node = grid.getNode (i);
      const Point & p (node->on_sphere);

      glPushMatrix ();
      glTranslatef (p.x ().to_double (),
                    p.y ().to_double (), p.z ().to_double ());
      glutSolidSphere (10, 10, 10);
      glPopMatrix ();
      glColor3f (0.3, 0.9, .5);
    }

  glLineWidth (1.0);
}


static void
draw_path (SPath * path )
{
  //double  length;

  if (path == NULL)
    return;

  SPathNode *ptr = path->first ();

  glBegin (GL_LINE_STRIP);

  while (ptr != NULL)
    {
      glVertex3f ((float) ptr->point ().x ().to_double (),
                  (float) ptr->point ().y ().to_double (),
                  (float) ptr->point ().z ().to_double ());

      ptr = ptr->next ();
    }

  glEnd ();
}

/*********** Hoai - 7/1/2012 ************/
static void draw_multishooting_path( vector<Point> *multiSP_path )
{
  if ( multiSP_path != NULL ){
    glBegin( GL_LINE_STRIP );
    for( int i=0; i<(int)multiSP_path->size(); i++ ){
      glVertex3f(
          (float)(*multiSP_path)[i].x().to_double(),
		  (float)(*multiSP_path)[i].y().to_double(),
		  (float)(*multiSP_path)[i].z().to_double() );
    }
    glEnd();
  }
}



static void draw_slices( deque< deque< Seg > > *slices )
{
  if ( slices != NULL ){
    //glLineWidth ( 3.0 );
    glLineWidth ( 2.0 );
    glLineStipple( 1, 0xFCFC );
    glEnable( GL_LINE_STIPPLE );
    for( int s = 1; s <= (int)slices->size(); s++ ){
      glBegin( GL_LINES );
      for( int i=0; i<(int)(*slices)[s].size(); i++ ){
        glColor3f (0.0, 0.0, 0.0);
        glVertex3f( (float)(*slices)[s][i].source().x().to_double(),
                    (float)(*slices)[s][i].source().y().to_double(),
                    (float)(*slices)[s][i].source().z().to_double() );
        glColor3f (0.0, 0.0, 1.0);
        glVertex3f( (float)(*slices)[s][i].target().x().to_double(),
		    (float)(*slices)[s][i].target().y().to_double(),
                    (float)(*slices)[s][i].target().z().to_double() );
      }
      glEnd();
    }
  }
}
/****************************************/

void
Renderer::drawPath ()
{
  SPathNode *ptr = shortest_path->first ();

  if (ptr == NULL)
    return;

  glEnable (GL_CULL_FACE);

  // draw shortest path by Agarwal
  // trang: disable ...
  /*
  glColor3f (0.0, 0.78, 0.0);
  glLineWidth (3.0);
  draw_path (shortest_path);
  */

  /************ Hoai - 27/8/2013 *************/
  glColor3f (1.0, 0.0, 0.0);
  glLineWidth (3.0);
  draw_multishooting_path( multiSP_path_first );

  /************ Hoai - 7/1/2013 *************/
  glColor3f (0.0, 0.0, 0.0);
  glLineWidth (6.0);
  draw_multishooting_path( multiSP_path );
  draw_slices( slices_ );
  /******************************************/

  // trang (2015/08/30) --- draw global one
  glColor3f (0.0, 1.0, 0.0);
  glLineWidth (10.0);
  draw_path( globalSGP );

  glColor3f (1.0, 1.0, 0.0);
  glLineWidth (3.0);
  for (int ind = 0; ind < paths_count; ind++)
    draw_path (arr_paths[ind]);

#ifdef  BOGI_GOGI
  glBegin (GL_LINE_STRIP);
  // get length of path
  //m_length = m_path->PathLength();
  m_length = shortest_path->getAprxLength ();

  while (ptr != NULL)
    {
      glVertex3f ((float) ptr->point ().x ().to_double (),
                  (float) ptr->point ().y ().to_double (),
                  (float) ptr->point ().z ().to_double ());
      //point_print( ptr->point() );
      //printf( " (%g, %g, %g)\n",
      //                (double)ptr->point().x().to_double(),
      //(double)ptr->point().y().to_double(),
      //    (double)ptr->point().z().to_double() );

      ptr = ptr->next ();
      //printf( "ptr = %p\n", ptr );
      //fflush( stdout );
    }

  //printf( "before gl_end\n" );
  //fflush( stdout );
  glEnd ();
  //printf( "after gl_end\n" );
#endif // BOGI GOGI
  glLineWidth (1.0);
  //printf( "before seting pointer\n" );

  // now draw the vertices of the path
  ptr = shortest_path->first ();
  //printf( "__ptr = %p\n", ptr );
  //fflush( stdout );

  // draw points on shortest path (by Agarwal)
  // trang: disable ...
  /*
  glColor3f (0.0, 1.0, 0.0);
  while (ptr != NULL)
    {
      glPushMatrix ();
      glTranslatef ((float) ptr->point ().x ().to_double (),
                    (float) ptr->point ().y ().to_double (),
                    (float) ptr->point ().z ().to_double ());
      glutSolidCube ( 0.5 );
      glPopMatrix ();

      ptr = ptr->next ();
    }
  */

  // Hoai - 7/1/2013
  // draw points for first multishooting path
  glColor3f (0.0, 0.1, 0.1);
  if ( multiSP_path_first != NULL ){
    for( int i=0; i<(int)multiSP_path_first->size(); i++ ){
      glPushMatrix ();
      glTranslatef( (float)(*multiSP_path_first)[i].x().to_double(),
		    (float)(*multiSP_path_first)[i].y().to_double(),
		    (float)(*multiSP_path_first)[i].z().to_double() );
      glutSolidCube ( 0.8 );
      glPopMatrix ();
    }
  }

  /************ Hoai - 7/1/2013 *************/
  glColor3f( 0.0, 0.0, 0.0 );
  if ( multiSP_path != NULL ){
    for( int i=0; i<(int)multiSP_path->size(); i++ ){
      glPushMatrix ();
      glTranslatef( (float)(*multiSP_path)[i].x().to_double(),
		    (float)(*multiSP_path)[i].y().to_double(),
		    (float)(*multiSP_path)[i].z().to_double() );
      glutSolidSphere (0.8, 5, 3);
      glPopMatrix ();
    }
  }
  /******************************************/

  // trang (2015/08/30) draw points (vertices) for global sdp
  glColor3f( 0.0, 0.0, 1.0 );
  if ( globalSGP != NULL )
  {
      ptr = globalSGP->first();
      //
      glBegin (GL_LINE_STRIP);
      //
      while (ptr != NULL)
        {
          glPushMatrix ();
          glTranslatef( (float) ptr->point().x().to_double(),
                        (float) ptr->point().y().to_double(),
                        (float) ptr->point().z().to_double() );
          glutSolidSphere (1.0, 5, 3);
          glPopMatrix ();

          ptr = ptr->next ();
        }
      //
      glEnd ();
  }
}


void
Renderer::updateDebugDraw ()
{
  //increment our list index, if we are at the beginning
  // our increment will set it to 0
  m_debugListIndex += 1;

  // check to see if we're at the end of the list, then stop
  if (m_debugListIndex >= (int) m_debugPoints.size ())
    {
      glutIdleFunc (NULL);
      return;
    }

  //otherwise, redraw
  glutPostRedisplay ();
}


void
Renderer::animateDebug ()
{
  //list_item li;
  list < AnimateVecInfo * >::iterator li;
  Point p;
  glLineWidth (2.0);

  glBegin (GL_LINES);


  // set our list item to the first line in our debug list
  //li = m_debugPoints.first();
  li = m_debugPoints.begin ();

  // loop through and only draw up to our current count of debug lines
  // we use <= because we are going to draw all lines up to and including
  // that index

  for (int i = 0; i <= m_debugListIndex; i++)
    {

      // set the color for the first vertex in the line depending
      // on what type of debug it is

      // also, if we are not supposed to show this type of line, then
      // just skip to the next line in the list
      //switch ((m_debugPoints.contents(li))->m_vecType) {
      assert( *li != NULL );
      switch ((*li)->m_vecType)
        {
        case (AnimateVecInfo::LONGITUDE):
          glColor3f (1.0, 1.0, 0.0);
          if (!m_showLong)
            {
              //first move to next item in list
              //li = m_debugPoints.succ(li);
              li++;
              continue;
            }
          break;

        case (AnimateVecInfo::REGION):
          glColor3f (0.0, 1.0, 1.0);
          if (!m_showReg)
            {
              //first move to next item in list
              //li = m_debugPoints.succ(li);
              li++;
              continue;
            }
          break;
        case (AnimateVecInfo::NN):
          glColor3f (1.0, 1.0, 0.0);
          if (!m_showNN)
            {
              //first move to next item in list
              //li = m_debugPoints.succ(li);
              li++;
              continue;
            }
          break;
        case (AnimateVecInfo::OTHER):
          glColor3f (0.0, 0.0, 0.0);
          break;
        }

      // get the first point
      //p = (m_debugPoints.contents(li))->m_p1;
      p = (*li)->m_p1;

      // draw the point
      glVertex3f ((float) p.x ().to_double (),
                  (float) p.y ().to_double (), (float) p.z ().to_double ());

      // now set our second color for the second point in our line
      //switch ((m_debugPoints.contents(li))->m_vecType)
      switch ((*li)->m_vecType)
        {
        case (AnimateVecInfo::LONGITUDE):
          glColor3f (1.0, 0.0, 0.0);
          break;
        case (AnimateVecInfo::REGION):
          glColor3f (0.0, 0.0, 1.0);
          break;
        case (AnimateVecInfo::NN):
          glColor3f (0.0, 1.0, 0.0);
          break;
        case (AnimateVecInfo::OTHER):
          glColor3f (1.0, 0.0, 0.0);
          break;
        }

      // now get the second point
      //p = (m_debugPoints.contents(li))->m_p2;
      p = (*li)->m_p2;

      // now draw the second poevôlutinint
      glVertex3f ((float) p.x ().to_double (),
                  (float) p.y ().to_double (), (float) p.z ().to_double ());

      // now move to the next item in the list
      //li = m_debugPoints.succ(li);
      li++;
    }

  glEnd ();
  glLineWidth (1.0);

}

void
Renderer::drawDebug ()
{
  list < AnimateVecInfo * >::iterator li;
  Point p;
  glLineWidth (2.0);

  // first check to make sure that atleast one of the types of debug
  // are to be drawn.  if not, we don't have to draw anything
  if (!m_showLong && !m_showReg && !m_showNN)
    return;

  glBegin (GL_LINES);

  // set our list item to the first line in our debug list
  li = m_debugPoints.begin ();

  // loop through and draw all debug lines
  for (int i = 0; i < (int) m_debugPoints.size (); i++)
    {
      // set our color for the first point in the line, depending on what
      // type of debugging stuff we have

      // also, if we are not supposed to show this type of line, then
      // just skip to the next line in the list
      // switch ((m_debugPoints.contents(li))->m_vecType) {
      switch ((*li)->m_vecType)
        {
        case (AnimateVecInfo::LONGITUDE):
          glColor3f (1.0, 1.0, 0.0);
          if (!m_showLong)
            {
              //first move to next item in list
              //li = m_debugPoints.succ(li);
              li++;
              continue;
            }
          break;
        case (AnimateVecInfo::REGION):
          glColor3f (0.0, 1.0, 1.0);
          if (!m_showReg)
            {
              //first move to next item in list
              //li = m_debugPoints.succ(li);
              li++;
              continue;
            }
          break;
        case (AnimateVecInfo::NN):
          glColor3f (1.0, 1.0, 0.0);
          if (!m_showNN)
            {
              //first move to next item in list
              // li = m_debugPoints.succ(li);
              li++;
              continue;
            }
          break;
        case (AnimateVecInfo::OTHER):
          glColor3f (0.0, 0.0, 0.0);
          break;
        }

      // get the first point
      //p = (m_debugPoints.contents(li))->m_p1;
      p = (*li)->m_p1;

      // draw the point
      glVertex3f ((float) p.x ().to_double (),
                  (float) p.y ().to_double (), (float) p.z ().to_double ());


      // now set our second color for the second point in our line
      // switch ((m_debugPoints.contents(li))->m_vecType)
      switch ((*li)->m_vecType)
        {
        case (AnimateVecInfo::LONGITUDE):
          glColor3f (1.0, 0.0, 0.0);
          break;
        case (AnimateVecInfo::REGION):
          glColor3f (0.0, 0.0, 1.0);
          break;
        case (AnimateVecInfo::NN):
          glColor3f (0.0, 1.0, 0.0);
          break;
        case (AnimateVecInfo::OTHER):
          glColor3f (1.0, 0.0, 0.0);
          break;
        }

      // now get the second point
      //p = (m_debugPoints.contents(li))->m_p2;
      p = (*li)->m_p2;

      // now draw the second point
      glVertex3f ((float) p.x ().to_double (),
                  (float) p.y ().to_double (), (float) p.z ().to_double ());

      /*
       */
      /*
         glPushMatrix();
         glTranslatef( (float)p.x().to_double(),
         (float)p.y().to_double(),
         (float)p.z().to_double());
         //glutSolidSphere( 2, 10, 10 );
         glPopMatrix();

       */

      // now move to the next item in the list
      //li = m_debugPoints.succ(li);
      li++;
    }
  glEnd ();

  glLineWidth (1.0);
}

void
Renderer::flyBy ()
{
  // if we haven't made our path of fly by points, we shouldn't be here
  if (!m_flyByPoints)
    return;

  // first check to see if this is our first loop through the
  // fly by.  if so, m_flyBy will still be false
  if (!m_flyBy)
    {
      m_flyBy = true;
      m_currFlyByPoint = 0;
    }

  // first check to see if we are to the end of the path, which
  // means that where we are looking from is the end of the path
  if (m_currFlyByPoint + FLY_LOOK_AHEAD >= m_numFlyByPoints)
    {
      m_flyBy = false;
      glutIdleFunc (NULL);

      // also, reset our viewpoint,viewdir, and up vector
      m_viewPos[0] = 0;
      m_viewPos[1] = 0;
      m_viewPos[2] = 200;

      m_viewDir[0] = 0;
      m_viewDir[1] = 0;
      m_viewDir[2] = 0;

      m_viewUp[0] = 0.0;
      m_viewUp[1] = 1.0;
      m_viewUp[2] = 0.0;

      return;
    }

  // otherwise, set our new flyby point and direction and up vector

  // set view position to fly by point
  m_viewPos[0] = m_flyByPoints[m_currFlyByPoint].flyByPoint[0];
  m_viewPos[1] = m_flyByPoints[m_currFlyByPoint].flyByPoint[1];
  m_viewPos[2] = m_flyByPoints[m_currFlyByPoint].flyByPoint[2];

  // set view direction the point on the path further down in
  // the path by FLY_LOOK_AHEAD number of points
  m_viewDir[0] =
    m_flyByPoints[m_currFlyByPoint + FLY_LOOK_AHEAD].pointOnPath[0];
  m_viewDir[1] =
    m_flyByPoints[m_currFlyByPoint + FLY_LOOK_AHEAD].pointOnPath[1];
  m_viewDir[2] =
    m_flyByPoints[m_currFlyByPoint + FLY_LOOK_AHEAD].pointOnPath[2];

  // bump up to next point on fly by path
  m_currFlyByPoint++;

  //now, redraw us
  glutPostRedisplay ();

}

void
Renderer::generateFlyByPoints ()
{
  // if we haven't made a path, we shouldn't even be here!
  if (shortest_path == NULL)
    return;

  m_length = shortest_path->getAprxLength ();
  SPathNode *currFlyByNode = NULL;
  float vecSrcDest[3];
  float flyByDist = 0.0;        // distance between source and dest
  float unitVecSrcDest[3];      //unit vector from source to dest
  float currFlyByDist = 0.0;    // our current distance travelled from source to
  // dest
  float remainder = 0.0;        //the extra distance travelled past the destination node
  // this will be added on to the next part of the path so
  // that throughout the path, each step is FLY_STEPs apart
  int i = 0;

  // trang: avoiding uninitialize variable
  unitVecSrcDest[0] = 0; unitVecSrcDest[1] = 0; unitVecSrcDest[2] = 0;

  // we generate the fly by points in our array, such that we have
  // a point for every FLY_STEP on the path, plus a couple additional points
  // at the end (FLY_LOOK_AHEAD amount), so that we reach all the way to the
  // end point in the path with our view position (and thus our view
  // direction point would be FLY_LOOK_AHEAD points past the end of the path

  // see how many points are going to be in the path( we add 1 so that
  // there is a point at the beginning of the path as well as one at the
  // end of the real path)
  m_numFlyByPoints = (int) (m_length / FLY_STEP) + 1 + FLY_LOOK_AHEAD;

  // now allocate flyby points
  if (m_flyByPoints)
    delete[]m_flyByPoints;
  m_flyByPoints = new FlyByPoint[m_numFlyByPoints];

  // now, there's a two step process in filling up the array
  // first, loop through the points in the path and fill up the first
  // m_numFlyByPoints-FLY_LOOK_AHEAD points, which is basically
  // all of them except for those last points that we add at the
  // end of the path
  for (i = 0; i < m_numFlyByPoints - FLY_LOOK_AHEAD; i++)
    {
      // if we've travelled as far as we're going to from one node to the next,
      // then we need to move on to the next node in the path
      if (currFlyByDist >= flyByDist)
        {

          // first bump up pointer to next node (or initialize it if this is
          // the first time through the loop
          if (!currFlyByNode)
            currFlyByNode = shortest_path->first ();
          else
            currFlyByNode = currFlyByNode->next ();

          //then get the remainder so we know how far we need to jump
          //along when starting on the new part of the path
          remainder = currFlyByDist - flyByDist;


          // find the new vector going from our current node to the next node
          vecSrcDest[0] =
            (float) currFlyByNode->next ()->point ().x ().to_double () -
            (float) currFlyByNode->point ().x ().to_double ();
          vecSrcDest[1] =
            (float) currFlyByNode->next ()->point ().y ().to_double () -
            (float) currFlyByNode->point ().y ().to_double ();
          vecSrcDest[2] =
            (float) currFlyByNode->next ()->point ().z ().to_double () -
            (float) currFlyByNode->point ().z ().to_double ();

          // now find the distance between these two points
          flyByDist = VEC_MAG (vecSrcDest[0], vecSrcDest[1], vecSrcDest[2]);

          // now find the unit vector going from source to dest
          unitVecSrcDest[0] = (float) vecSrcDest[0] / flyByDist;
          unitVecSrcDest[1] = (float) vecSrcDest[1] / flyByDist;
          unitVecSrcDest[2] = (float) vecSrcDest[2] / flyByDist;

          // set our currently travelled distance from source to dest as the
          // remainder left over from the last part of the path
          currFlyByDist = remainder;
        }

      // set the point on the path in our current entry.  it's going to be
      // measured from the current source point plus some constant (which is
      // updated by FLY_STEP amount each time) times the unit vector made
      // from the current source node to the current destination node

      m_flyByPoints[i].pointOnPath[0] =
        (float) currFlyByNode->point ().x ().to_double () +
        currFlyByDist * unitVecSrcDest[0];
      m_flyByPoints[i].pointOnPath[1] =
        (float) currFlyByNode->point ().y ().to_double () +
        currFlyByDist * unitVecSrcDest[1];
      m_flyByPoints[i].pointOnPath[2] =
        (float) currFlyByNode->point ().z ().to_double () +
        currFlyByDist * unitVecSrcDest[2];

      // now set the unit vector in our current entry, so first find magnitude
      // of vector (point on path to origin)
      float mag = VEC_MAG (m_flyByPoints[i].pointOnPath[0],
                           m_flyByPoints[i].pointOnPath[1],
                           m_flyByPoints[i].pointOnPath[2]);

      // now find unit vector
      m_flyByPoints[i].unitVec[0] = m_flyByPoints[i].pointOnPath[0] / mag;
      m_flyByPoints[i].unitVec[1] = m_flyByPoints[i].pointOnPath[1] / mag;
      m_flyByPoints[i].unitVec[2] = m_flyByPoints[i].pointOnPath[2] / mag;

      // now set the hovering point above the point on the path that we will
      // be flying at.  this point is found by adding a constant times the unit
      // vector made from the point on the path

      m_flyByPoints[i].flyByPoint[0] = m_flyByPoints[i].pointOnPath[0] +
        FLY_DIST_OFF_GROUND * m_flyByPoints[i].unitVec[0];
      m_flyByPoints[i].flyByPoint[1] = m_flyByPoints[i].pointOnPath[1] +
        FLY_DIST_OFF_GROUND * m_flyByPoints[i].unitVec[1];
      m_flyByPoints[i].flyByPoint[2] = m_flyByPoints[i].pointOnPath[2] +
        FLY_DIST_OFF_GROUND * m_flyByPoints[i].unitVec[2];

      // now, bump up our current distance travelled by our translation
      // constant
      currFlyByDist += FLY_STEP;

    }

  /////////////////////////////////////////////////////////////////////////
  // now, the second and final step is to add those couple of points at the
  // end of the path
  //
  // we do this by working with the vector from the second to last node
  // and the last node (which should be the current vector) and we basically
  // just keep extending point along this vector
  //
  // we need to only set the point on path point for each of these
  // entries since they will only end up being used for view direction
  // (our view position/point stops at the last point on the actual path)
  /////////////////////////////////////////////////////////////////////////

  // now loop through and set point on path
  for (i = m_numFlyByPoints - FLY_LOOK_AHEAD; i < m_numFlyByPoints; i++)
    {
      // set the point on path
      m_flyByPoints[i].pointOnPath[0] =
        (float) currFlyByNode->point ().x ().to_double () +
        currFlyByDist * unitVecSrcDest[0];
      m_flyByPoints[i].pointOnPath[1] =
        (float) currFlyByNode->point ().y ().to_double () +
        currFlyByDist * unitVecSrcDest[1];
      m_flyByPoints[i].pointOnPath[2] =
        (float) currFlyByNode->point ().z ().to_double () +
        currFlyByDist * unitVecSrcDest[2];

      // bump up our distance
      currFlyByDist += FLY_STEP;
    }

  // **** wait to set the up vector for the view point until we are
  // actually going through the loop that will do the fly by since
  // we want to be able to change how far ahead we are looking, and
  // the up vector depends on this
}



// axis: 0 - move in Z direction (towards view point)
//       1 - move in X direction (perpendicular to view point)
// dir: >= 0 in positive direction (towards viewpoint)
//      < 0 in negative direction (away from viewpoint)

void
Renderer::moveViewPoint (int axis, int dir)
{
  float phi = 0.0;
  float deltaX, deltaZ, val;

  // val = (m_viewDir[2] - m_viewPos[2]) / m_viewPointDistance;
  // trang: for x and y moving:
  val = (m_viewDir[1] - m_viewPos[1]) / m_viewPointDistance;

  if (m_viewDir[0] > m_viewPos[0])
    {
      phi = acos (val);
    }
  else
    phi = -acos (val);

  switch (axis)
    {
    case 0:

      deltaZ = TRANSLATION_STEP * cos (phi);
      deltaX = TRANSLATION_STEP * sin (phi);

      if (dir < 0)
        {
          deltaZ = -deltaZ;
          deltaX = -deltaX;
        }
      m_viewDir[0] += deltaX;
      //m_viewDir[2] += deltaZ;
      m_viewDir[1] += deltaZ;                   // trang: for x and y moving

      m_viewPos[0] += deltaX;
      //m_viewPos[2] += deltaZ;
      m_viewPos[2] += deltaZ;                   // trang: for x and y moving

      break;

    case 1:                    // x axis
      if (dir < 0)              // strafing left
        {
          deltaX = TRANSLATION_STEP * cos (phi);
          deltaZ = -TRANSLATION_STEP * sin (phi);
        }
      else                      // strafing right
        {
          deltaX = -TRANSLATION_STEP * cos (phi);
          deltaZ = TRANSLATION_STEP * sin (phi);
        }

      m_viewDir[0] += deltaX;
      //m_viewDir[2] += deltaZ;
      m_viewDir[1] += deltaZ;                   // trang: for x and y moving

      m_viewPos[0] += deltaX;
      //m_viewPos[2] += deltaZ;
      m_viewPos[1] += deltaZ;                   // trang: for x and y moving
      break;
    }
}


void
Renderer::rotateViewDir (float angle)
{
  float phi = 0.0;
  float deltaX, deltaZ, val;

  val = (m_viewDir[2] - m_viewPos[2]) / m_viewPointDistance;
  if (m_viewDir[0] > m_viewPos[0])
    {
      phi = acos (val);
    }
  else
    phi = -acos (val);

  deltaZ = m_viewPointDistance * cos (angle + phi);
  deltaX = m_viewPointDistance * sin (angle + phi);

  m_viewDir[0] = m_viewPos[0] + deltaX;
  m_viewDir[2] = m_viewPos[2] + deltaZ;
}

void
Renderer::rotateViewPoint (float angle)
{
  float phi = 0.0;
  float deltaX, deltaZ, val;

  val = (m_viewPos[2] - m_viewDir[2]) / m_viewPointDistance;
  if (m_viewPos[0] > m_viewDir[0])
    {
      phi = acos (val);
    }
  else
    phi = -acos (val);

  deltaZ = m_viewPointDistance * cos (angle + phi);
  deltaX = m_viewPointDistance * sin (angle + phi);
  m_viewPos[0] = m_viewDir[0] + deltaX;
  m_viewPos[2] = m_viewDir[2] + deltaZ;
}

void
Renderer::keyboard (unsigned char key, int x, int y)
    // checks keypresses: pressing Enter key begin animations, pressing Q quits
{
    // trang 2013/12: avoid unused variables
    x = x + 0; y = y + 0;

    // key pressed ...
    switch (key)
    {
    case 13: // trang: solve the problem (press enter)
        if ( m_facesSelected == 2  )
        {   
            m_path->solve_finish = true;
            m_path->multispCalcPath ( m_srcPoint, m_destPoint,
                                      m_srcFace, m_destFace, m_useHershSuri,
                                      shortest_path, multiSP_path, multiSP_path_first,
                                      slices_ );

            // trang: mark for something
            problem_solved = true;
            m_path->solve_finish = false;
        }

        // trang: set globalSGP to be null for computing the global solution
        globalSGP = NULL;
        //
        glutPostRedisplay ();
        break;

    case 'o': // (2015/08/28) trang: find the global SDP from p to q - 'o' for optimal
        if (problem_solved )
        {
            cout << "%% Computing the global shortest descending path ..." << endl;

            // declare and assign source/destination info
            SourceDest srcdest;
            srcdest.m_srcPoint = m_path->m_srcdest->m_srcPoint;
            srcdest.m_destPoint = m_path->m_srcdest->m_destPoint;
            srcdest.m_srcFace = m_path->m_srcdest->m_srcFace;
            srcdest.m_destFace = m_path->m_srcdest->m_destFace;

            // trang's debug for source/dest points ... OK
            cout << "   - determining source point (p): "
                 << srcdest.m_srcPoint.x().to_double () << ", "
                 << srcdest.m_srcPoint.y().to_double () << ", "
                 << srcdest.m_srcPoint.z().to_double () << endl;
            cout << "   - determining destination point (q): "
                 << srcdest.m_destPoint.x().to_double () << ", "
                 << srcdest.m_destPoint.y().to_double () << ", "
                 << srcdest.m_destPoint.z().to_double () << endl;
            //
            //exit(0);

            // get the cutting slice through the source point
            srcdest.m_srcSlice = m_path->m_srcdest->m_srcSlice;
            srcdest.m_srcSliceFaces = m_path->m_srcdest->m_srcSliceFaces;

            //
            m_path->calcSGP( &srcdest, false, globalSGP );

            //
            // final information
            cout << "====================================== Global solution ===================================" << endl;
            cout << "%% Global SGP(p,q) obtained ..." << endl;
            cout << "%% Source and destination: " << endl;
            cout << "\t "; point_print( srcdest.m_srcPoint  ); cout << endl;
            cout << "\t "; point_print( srcdest.m_destPoint ); cout << endl;
            //
            cout << "%% The final path: " << endl;
            globalSGP->printNodeList();
            //cout << "\t - initial path:" << initiallen << endl;
            cout << "\t - approximate length:" << globalSGP->getAprxLength() << endl;

            problem_solved = false;
        }
        else
        {
            cout << "-- solution of multiple shooting NEEDED to solved first ..." << endl;
        }

        glutPostRedisplay ();
        break;

    case 9: // trang: solve step by step (press tab)
        if ( m_facesSelected == 2  )
        {
            if (!m_path->show_cut_slices && !m_path->finish)
            {
                m_path->show_initial_path = true;
            }
            m_path->multispCalcPath ( m_srcPoint, m_destPoint,
                                      m_srcFace, m_destFace, m_useHershSuri,
                                      shortest_path, multiSP_path, multiSP_path_first,
                                      slices_ );
        }
        glutPostRedisplay ();     
        break;

    case 'z': // trang: zoom in
        scale = scale * 1.05; glutPostRedisplay (); break;
    case 'Z': // trang: zoom out
        scale = scale / 1.05; glutPostRedisplay (); break;

    case 27:
      //case 'q':
      //case 'Q':
      exit (0);

    case 'I':                  // for I, in positive y direction (strafe up)
      m_viewPos[1] += TRANSLATION_STEP;
      m_viewDir[1] += TRANSLATION_STEP;
      glutPostRedisplay ();
      break;
    case 'K':                  // for K, in negative y direction (strafe down)
      m_viewPos[1] -= TRANSLATION_STEP;
      m_viewDir[1] -= TRANSLATION_STEP;
      glutPostRedisplay ();
      break;
    case 'J':                  // for strafe left
      moveViewPoint (1, -1);
      glutPostRedisplay ();
      break;
    case 'L':                  // for strafe right
      moveViewPoint (1, 1);
      glutPostRedisplay ();
      break;
    case 'k':                  // for k, in negative z direction
      moveViewPoint (0, -1);
      glutPostRedisplay ();
      break;
    case 'i':                  // for i, in positive z direction
      moveViewPoint (0, 1);
      glutPostRedisplay ();
      break;
    case 'G':                  // do we draw the whole graph
    case 'g':
      m_drawGraph = !m_drawGraph;
      glutPostRedisplay ();
      break;

    case 'F':                  // do fly by only if we've already created a path
    case 'f':
      m_drawGraph = false;
      if (shortest_path != NULL)
        {
          f_fly_mode = true;
          //m_pathLength
          // first generate the points on the fly by path
          generateFlyByPoints ();

          // now set our idle function to the fly by updater
          glutIdleFunc (RendererFlyByWrapper);
          f_fly_mode = false;
        }
      break;

    case 'B':
    case 'b':                  // animate debug stuff

      if (m_animateDebug)       // if it's on and we want to either
        //stop it momentarily (until we hit 'animate' again, which restarts
        // it), or turn it off
        {
          // if we're done drawing, then set the flag off
          if (m_debugListIndex >= (int) m_debugPoints.size ())
            {
              m_animateDebug = false;
              m_animatingDebug = false;
            }
          else if (m_animatingDebug)  // if we are animating, then
            // just 'freeze' the animation
            {
              glutIdleFunc (NULL);
              m_animatingDebug = false;
            }
          else                  // otherwise, we are currently 'frozen' and we
            // want to continue animating
            {
              glutIdleFunc (RendererUpdateDebugDrawWrapper);
              m_animatingDebug = true;
            }

          glutPostRedisplay ();
        }
      else                      // other wise, we want to animate the debug
        {
          m_animatingDebug = true;

          glutIdleFunc (RendererUpdateDebugDrawWrapper);
        }

      break;
    case 'V':                  // simply draw debug stuff
    case 'v':
      m_drawDebug = !m_drawDebug;

      // also, if we we're animating the debug and they want it to go away,
      // then let's turn it off too
      if (m_animateDebug)
        {
          m_animateDebug = false;
          m_animatingDebug = false;
          glutIdleFunc (NULL);
        }

      glutPostRedisplay ();
      break;

    case 'A':                  // rotate viewpoint CW around y-axis
    case 'a':
      m_yRot += 5.0;
      m_yRot = fmod (m_yRot, (float) 360.0);
      glutPostRedisplay ();
      break;

    case 'D':                  // rotate viewpoint CCW around y-axis
    case 'd':
      m_yRot -= 5.0;
      m_yRot = fmod (m_yRot, (float) 360.0);
      glutPostRedisplay ();
      break;
    case 'W':                  // rotate viewpoint CW around x-axis
    case 'w':
      m_xRot += 5.0;
      m_xRot = fmod (m_xRot, (float) 360.0);
      glutPostRedisplay ();
      break;

    case 's':
      m_xRot -= 5.0;
      m_xRot = fmod (m_xRot, (float) 360.0);
      glutPostRedisplay ();
      break;

    case '1':
      if (shortest_path == NULL)
        break;
      f_no_draw = true;
      printf ("1-Shortcutting. Current length: %g\n",
              shortest_path->getAprxLength ());
      shortest_path->shortcut_on_faces ();
      printf ("1-                  New length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      glutPostRedisplay ();
      f_no_draw = false;
      break;

    case '2':
      if (shortest_path == NULL)
        break;
      f_no_draw = true;
      printf ("2- Shortcutting. on edge Current length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      shortest_path->shortcut_on_edge ();
      printf ("2-                           New length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      glutPostRedisplay ();
      f_no_draw = false;
      break;

    case '3':
      if (shortest_path == NULL)
        break;
      f_no_draw = true;
      printf ("3 - Shortcutting. on vertex Current length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      shortest_path->shortcut_on_vertex ();
      printf ("3 -                           New length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      glutPostRedisplay ();
      f_no_draw = false;
      break;
      /*
         case  '4' :
         if  ( shortest_path == NULL )
         break;
         f_no_draw = true;
         printf( "4 - Shortcutting.by relifting. Current length: %g\n",
         shortest_path->getAprxLength() );
         fflush( stdout );

         shortest_path->shortcut_by_relifting();
         printf( "4 -                               New length: %g\n",
         shortest_path->getAprxLength() );
         fflush( stdout );

         f_no_draw = false;

         glutPostRedisplay();
         break;
       */

    case '5':
      if (shortest_path == NULL)
        break;
      f_no_draw = true;
      printf ("5 - Shortcutting.by relifting (jump). Current length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      shortest_path->shortcut_by_relifting_jump ();
      printf ("5 -                               New length: %g\n",
              shortest_path->getAprxLength ());
      fflush (stdout);

      f_no_draw = false;

      glutPostRedisplay ();
      break;

    case '0':
      if (shortest_path == NULL)
        break;
      shortest_path->dump ();
      break;

    case ' ':
      f_display_selected = !f_display_selected;
      glutPostRedisplay ();
      break;

    case '#':
      f_draw_sphere = !f_draw_sphere;
      glutPostRedisplay ();
      break;


    case '@':
      if ((paths_count >= MAX_PATHS) || (shortest_path == NULL))
        break;
      arr_paths[paths_count++] = new SPath (shortest_path);
      break;

    case '+':
    case '-':
      break;

    case '!':
      f_solid = !f_solid;
      break;

    case 'j':                  // rotate view direction (turn)
      rotateViewDir (.1);
      glutPostRedisplay ();
      break;
    case 'l':                  // rotate view direction (turn)
      rotateViewDir (-.1);
      glutPostRedisplay ();
      break;
    case 'n':                  //toggle normals
      m_drawNormals = !m_drawNormals;
      glutPostRedisplay ();
      break;

    case '?':
        fprintf (stderr, "Using: \n"
               "\t q - quit, \n "
               "\t I - strafe up, K - strafe down \n"
               "\t J - strafe left, L - strafe right \n"
               "\t k - strafe forward, l - strafe backword\n"
               "\t a,d - rotate around y axis\n"
               "\t g - graph on/off  f - fly\n"
               "\t b,v - animate,draw debug stuff\n"
               "\t 1,2,3 - shortcutting methods\n"
               "\t 0 - dump shortest path\n"
               "\t s,j,l,n -???\n");
        break;
    }
}

/*---------------- trang: handling direction keys -----------------*/
void Renderer::handleSpecialKeypress(int key, int x, int y)
{
    // avoiding unused variables
    x = x + 0; y = y + 0;

    // reset view position of eyes
    //m_viewPos[0] = 0; m_viewPos[1] = 0; m_viewPos[2] = 100;

    switch (key)
    {
    case GLUT_KEY_LEFT:
        m_transDis[0] -= 5;    break;
    case GLUT_KEY_RIGHT:
        m_transDis[0] += 5;    break;
    case GLUT_KEY_DOWN:
        m_transDis[1] -= 5;    break;
    case GLUT_KEY_UP:
        m_transDis[1] += 5;    break;
    }

    glutPostRedisplay ();
}
/*-----------------------------------------------------------------*/

void
Renderer::mouse (int button, int state, int x, int y)
{
  GLuint selectBuf[512];
  GLint viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];
  GLint realy;                  //openGL y coordinate position
  GLint hits;
  GLdouble wx, wy, wz;          //world coordinate of the point that we've clicked

  if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN)
    return;

  // trang: avoid moving of object which is caused by direction keys
  //m_transDis[0] = 0; m_transDis[1] = 0; m_transDis[2] = 0;

  ////////////////////////////////////////////////////////////
  // the first thing to do is find the world coordinates of the
  // point that we've clicked....
  // this is a pain in the butt

  // first, we must get the view Depth buffer and store in m_viewDepthBuffer

  if (m_viewDepthBuffer)
    free (m_viewDepthBuffer);
  m_viewDepthBuffer = (GLushort *) malloc (4 * m_windowHeight * m_windowWidth
                                           * sizeof (GLushort));
  glReadPixels (0, 0, m_windowWidth, m_windowHeight, GL_DEPTH_COMPONENT,
                GL_UNSIGNED_SHORT, m_viewDepthBuffer);

  // now get all the matrix info that we need to unproject the 2D
  // point to get our 3D world coord
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  // now get our real y coordinate
  realy = viewport[3] - (GLint) y - 1;

  // now we must find what our depth buffer says about the point we clicked
  // at.  the result in the zViewPlane will be between 0 and 1, where 0
  // represents the near view plane, and 1 represents the far plane

  debug ("depth buffer is: " << m_viewDepthBuffer[x + realy * m_windowWidth]);
  cout << "depth buffer is: " << m_viewDepthBuffer[x + realy * m_windowWidth] << endl;
  float zViewPlane =
    ((float) m_viewDepthBuffer[x + realy * m_windowWidth]) / ((float) 65536);


  debug ("zviewplane is: " << zViewPlane);
  cout << "zviewplane is: " << zViewPlane << endl;

  // finally, we get the world coords for view point, storing in wx,wy
  // and wz
  gluUnProject ((GLdouble) x, (GLdouble) realy, zViewPlane, mvmatrix,
                projmatrix, viewport, &wx, &wy, &wz);

  debug ("world coords are: " << (double) wx << " " << (double) wy << " " <<
         (double) wz);
  fprintf (stderr, "%.10g %.10g %.10g\n", (double) wx, (double) wy,
           (double) wz);

  /////////////////////////////////////////////////////////////
  // now deal with the selection business

  glSelectBuffer (512, selectBuf);
  (void) glRenderMode (GL_SELECT);

  glInitNames ();
  glPushName (0);

  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity ();

  // make a very very small picking window to make sure that we can only
  // select one face
  gluPickMatrix ((GLdouble) x, (GLdouble) (viewport[3] - y), 0.000001,
                 0.000001, viewport);

  gluPerspective (90.0, (GLfloat) m_windowWidth / (GLfloat) m_windowHeight, 1,
                  6000.0);
  // now draw our polytope in select mode
  drawPolytope (GL_SELECT);

  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glFlush ();

  hits = glRenderMode (GL_RENDER);

  // here we check to see what we've actually selected...
  // only proceed if we've actually selected a point on a face
  // and the view plane isn't equal to one
  // we must do this because there are border cases where you
  // can generate a hit on the face but the point clicked is
  // on the far plane
  if (!processHits (hits, selectBuf) ||
      m_viewDepthBuffer[x + realy * m_windowWidth] == 0xffff)
    {
      glutPostRedisplay ();     // redraw anyway, if we wanted to unselect
      // a face
      return;
    }

  // now that we know that our selected point is valid,
  // store these coords in the appropriate array (dest or src)
  if (m_facesSelected == 1)
    {
      m_srcPoint[0] = (double) wx;
      m_srcPoint[1] = (double) wy;
      m_srcPoint[2] = (double) wz;
    }
  else
    {
      m_destPoint[0] = (double) wx;
      m_destPoint[1] = (double) wy;
      m_destPoint[2] = (double) wz;
    }

  // now, if we've clicked on two faces, then we have all the info
  // we need to calculate the path
  if (m_facesSelected == 2)
    {
      // first, force a redraw so we can highlight the second face
      // that we selected before we have to wait to calculate the
      // path
      display ();

      debug ("Source face = " << dec << m_srcFace);
      debug ("Dest face = " << dec << m_destFace);

      // now calculate the path
      printf ("Calling calcpth - A\n");
      fflush (stdout);

      FILE *fp = fopen( "srcdest.txt", "w" );
      fprintf( fp, "%g\t%g\t%g\n", m_srcPoint[0],  m_srcPoint[1],  m_srcPoint[2] );
      fprintf( fp, "%g\t%g\t%g\n", m_destPoint[0], m_destPoint[1], m_destPoint[2] );
      fprintf( fp, "%d\t%d\n", m_srcFace, m_destFace );
      fclose( fp );

      //m_path->CalcPath ( m_srcPoint, m_destPoint,
      // 		 m_srcFace, m_destFace, m_useHershSuri,
      //		 shortest_path );

      // trang: make comment to test step by step visualization ...
      /*
      m_path->multispCalcPath ( m_srcPoint, m_destPoint,
                                m_srcFace, m_destFace, m_useHershSuri,
                                shortest_path, multiSP_path, multiSP_path_first,
                                slices_ );

      m_facesSelected = 0;
      */
    }

  glutPostRedisplay ();
}

bool Renderer::processHits (GLint hits, GLuint buffer[])
{
  int i, j;
  GLuint
  names, *
    ptr;

  ptr = (GLuint *) buffer;

  // if hits is zero, we clicked somewhere on the
  // screen that wasn't a face on the polytope, such as debugging
  // lines, or black space
  // so, let's unselected the other faces that may be selected
  if (hits == 0)
    {
      m_facesSelected = 0;      // now, no faces are selected
      shortest_path = NULL;
      m_srcFace = -1;
      m_destFace = -1;
      return false;
    }

  // else, we have a face that's been selected
  m_facesSelected++;

  // if we already have 2 points selected and thus have calculated
  // a path, then reset everything first
  if (m_facesSelected > 2)
    {
      m_facesSelected = 1;      // now just current point is selected
      //m_traversalPath = NULL;
      m_srcFace = -1;
      m_destFace = -1;
    }

  for (i = 0; i < hits; i++)
    {
      names = *ptr;
      ptr++;
      ptr++;
      ptr++;
      for (j = 0; j < (int) names; j++)
        {
          // set our face... either the source or destination face
          // depending on what our count says
          if (m_facesSelected == 1)
            m_srcFace = (int) *ptr;
          else if (m_facesSelected == 2)
            m_destFace = (int) *ptr;

          ptr++;
        }
      return true;
    }

  return false;
}

void
Renderer::drawAxes (float size)
{
  glBegin (GL_LINES);
  glVertex3f (size, 0, 0);
  glVertex3f (-size, 0, 0);
  glVertex3f (0, size, 0);
  glVertex3f (0, -size, 0);
  glVertex3f (0, 0, size);
  glVertex3f (0, 0, -size);
  glEnd ();

}


// function that draws the polytope
void
Renderer::display ()
{

    //printf( "Rendered::dipslay: %d\n", (int)f_no_draw ); fflush( stdout );

    // clear the buffers for redrawing
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    // say how we want to draw our polys
    if (m_drawPolytopeFilled)
      {
        glPolygonMode (GL_FRONT, GL_FILL);
        glPolygonMode (GL_BACK, GL_FILL);
      }
    else
      {
        glPolygonMode (GL_FRONT, GL_LINE);
        glPolygonMode (GL_BACK, GL_LINE);
      }

    // let us write into the depth buffers
    glEnable (GL_DEPTH_TEST);
    glDepthMask (GL_TRUE);

    // reset our modelview matrix and viewpoint
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
    gluLookAt (m_viewPos[0], m_viewPos[1], m_viewPos[2],
               m_viewDir[0], m_viewDir[1], m_viewDir[2],
               m_viewUp[0], m_viewUp[1], m_viewUp[2]);

    /* ---------- 2013/12, trang: some modifications ------------- */
    // move object in vertical and horizontal directions
    glTranslatef( m_transDis[0], m_transDis[1], m_transDis[2]);

    // zoom in/out
    glScalef(scale, scale, scale);
    /* ----------------------end of trang's----------------------- */

    // now do the lighting updating here, so that our lighting conditions
    // aren't affected when we move the polytope
    // all we really have to update is the position
    /*
       GLfloat light1_position[] = {0000,10000,20000,1};
       GLfloat light2_position[] = {-20000,0,0,1};
       GLfloat light3_position[] = {0,0,20000,1};
       GLfloat light4_position[] = {0,0,-20000,1};

       glLightfv(GL_LIGHT0, GL_POSITION, light1_position);
       glLightfv(GL_LIGHT1, GL_POSITION, light2_position);
       glLightfv(GL_LIGHT2, GL_POSITION, light3_position);
       glLightfv(GL_LIGHT3, GL_POSITION, light4_position);
     */
    // for the rest of here, we're drawing the polytope (and possibly
    // the path, too)

    glColor3f (1, 1, 1);

    if (m_drawAxes)
      drawAxes (500);

    // before we draw the polytope or the path, we apply the rotation
    // transformations that the user may have created from
    // keystrokes. this involves some unfortunate vector
    // manipulation for the second rotation that basically
    // removes the y rotation effect so that when we rotate
    // with respect to x, it is in the original coordinate
    // system instead of the new coordinate system made
    // by the y rotation

    if (!m_flyBy)
      {
        glRotatef (m_yRot, 0.0, 1.0, 0.0);  // do our y rotate
        glRotatef (m_xRot, sin (m_yRot * PI_DIV_180), 0, cos (m_yRot * PI_DIV_180));  // do our x rotate
        /**************** Hoai - 19/6/2013 ******************/
        glRotatef (-60.0, 1.0, 0.0, 0.0 );
        glRotatef (-50.0, 0.0, 0.0, 1.0 );
        /****************************************************/
      }

    if ((shortest_path != NULL) && m_drawPath)
      {
        drawPath ();
      }

    if (m_drawGraph)
      {
        drawGraph ();
      }

    if (f_draw_sphere)
      {
        drawSphere ();
      }


    //printf( "3Rendered::dipslay: - done %d\n", (int)f_no_draw );
    //fflush( stdout );

    if (m_animateDebug)
      {
        animateDebug ();
      }
    else                          // if we're not animating, then just draw our current debug stuff
      // may not draw anything... if none of the debug type flags
      // are set then function just returns
      drawDebug ();

    //printf( "2Rendered::dipslay: - done %d\n", (int)f_no_draw );
    //fflush( stdout );

    if (m_drawPolytope)
      {
        drawPolytope (GL_RENDER);
      }


    //printf( "1Rendered::dipslay: - done %d\n", (int)f_no_draw );
    //fflush( stdout );

    glutSwapBuffers ();

    //printf( "Rendered::dipslay: - done %d\n", (int)f_no_draw );
    //fflush( stdout );
}

// called when window is resized
void
Renderer::reshape (int w, int h)
{
  m_windowWidth = w;
  m_windowHeight = h;

  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  //  glOrtho(-200,200,-200,200,5,600);
  gluPerspective (90.0, (GLfloat) m_windowWidth / (GLfloat) m_windowHeight, 1,
                  6000.0);
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
  gluLookAt (m_viewPos[0], m_viewPos[1], m_viewPos[2],
             m_viewDir[0], m_viewDir[1], m_viewDir[2],
             m_viewUp[0], m_viewUp[1], m_viewUp[2]);
}


void
Renderer::algorithmDisplayMenu (int value)
{
  switch (value)
    {
    case (1):                  // we're showing all debug stuff
      m_showLong = true;
      m_showReg = true;
      m_showNN = true;
      glutPostRedisplay ();
      break;
    case (2):                  // hiding all debug stuff
      m_showLong = false;
      m_showReg = false;
      m_showNN = false;
      glutPostRedisplay ();
      break;
    case (3):                  // show or hide locations
      m_showLong = !m_showLong;
      glutPostRedisplay ();
      break;
    case (4):                  // show or hide Regions
      m_showReg = !m_showReg;
      glutPostRedisplay ();
      break;
    case (5):                  // show or hide NNs
      m_showNN = !m_showNN;
      glutPostRedisplay ();
      break;
    }

}

void
Renderer::algorithmAnimateMenu (int value)
{
  switch (value)
    {
    case (1):                  // start or pause or unpause animation

      if (m_animateDebug)       // if animation is on and we want to either
        // stop it momentarily (until we hit 'animate' again, which restarts
        // it), or turn it off
        {
          // if we're done drawing, then set the flag off
          if (m_debugListIndex >= (int) m_debugPoints.size ())
            {
              m_animateDebug = false;
              m_animatingDebug = false;
            }
          else if (m_animatingDebug)  // if we are animating, then
            // just 'freeze' the animation
            {
              glutIdleFunc (NULL);
              m_animatingDebug = false;
            }
          else                  // otherwise, we are currently 'frozen' and we
            // want to continue animating
            {
              glutIdleFunc (RendererUpdateDebugDrawWrapper);
              m_animatingDebug = true;
            }

          glutPostRedisplay ();
        }
      else                      // other wise, we want to animate the debug
        {
          // if none of the debug drawing flags are set, then set them all
          if (!m_showLong && !m_showReg && !m_showNN)
            {
              m_showLong = true;
              m_showReg = true;
              m_showNN = true;
            }

          if (m_flyBy)          // if fly by is already running, then stop it
            {
              glutIdleFunc (NULL);
              m_flyBy = false;
            }

          m_animateDebug = true;
          m_animatingDebug = true;
          m_debugListIndex = -1;

          glutIdleFunc (RendererUpdateDebugDrawWrapper);
        }
      break;
    case (2):                  // stop animation

      if (m_animateDebug)       // first, make sure that the animation was
        // running.... otherwise, do nothing
        {
          m_animateDebug = false;
          m_animatingDebug = false;
          glutIdleFunc (NULL);
          glutPostRedisplay ();
          break;
        }
    }
}

void
Renderer::otherDisplayMenu (int value)
{
  switch (value)
    {
    case (1):                  // show, hide normals
      m_drawNormals = !m_drawNormals;
      glutPostRedisplay ();
      break;
    case (2):                  // show, hide axes
      m_drawAxes = !m_drawAxes;
      glutPostRedisplay ();
      break;
    case (3):                  // show, hide graph
      m_drawGraph = !m_drawGraph;
      glutPostRedisplay ();
      break;
    case (4):                  // show, hide path
      m_drawPath = !m_drawPath;
      glutPostRedisplay ();
      break;
    case (5):                  // show, hide polytope
      m_drawPolytope = !m_drawPolytope;
      glutPostRedisplay ();
      break;
    case (6):                  // fill, lined of polytope
      m_drawPolytopeFilled = !m_drawPolytopeFilled;
      glutPostRedisplay ();
      break;
    }
}


void
Renderer::mainMenu (int value)
{
  switch (value)
    {
    case (2):                  // user wants to quit
      exit (0);
      break;
    case (1):                  // user wants to do a flyby
      //first check to see if we want to turn it off
      if (m_flyBy)
        {                       //fly by is already running
          glutIdleFunc (NULL);
          m_flyBy = false;
        }
      else if (shortest_path != NULL)
        {
          // make sure that path exists
          // before we do anything, turn off algorithm animation
          // if it's running
          if (m_animateDebug)
            {
              m_animateDebug = false;
              m_animatingDebug = false;
              glutIdleFunc (NULL);
            }

          // first generate the points on the fly by path
          generateFlyByPoints ();

          // now set our idle function to the fly by updater
          glutIdleFunc (RendererFlyByWrapper);
        }
    }
}

// creates a menu attached to the right mouse button that handles
// standard options

void
Renderer::initMainMenu ()
{
  int algorithmDisplayMenuID;
  glutCreateMenu (RendererAlgorithmDisplayMenuWrapper);
  glutAddMenuEntry ("Show All", 1);
  glutAddMenuEntry ("Hide All", 2);
  glutAddMenuEntry ("Show/Hide Longitudes", 3);
  glutAddMenuEntry ("Show/Hide Regions", 4);
  glutAddMenuEntry ("Show/Hide NNs", 5);
  algorithmDisplayMenuID = glutGetMenu ();

  int algorithmAnimateMenuID;
  glutCreateMenu (RendererAlgorithmAnimateMenuWrapper);
  glutAddMenuEntry ("Start/Pause/Unpause Animation", 1);
  glutAddMenuEntry ("Stop Animation", 2);
  algorithmAnimateMenuID = glutGetMenu ();

  int otherDisplayMenuID;
  glutCreateMenu (RendererOtherDisplayMenuWrapper);
  glutAddMenuEntry ("Show/Hide Normals", 1);
  glutAddMenuEntry ("Show/Hide Axes", 2);
  glutAddMenuEntry ("Show/Hide Graph", 3);
  glutAddMenuEntry ("Show/Hide Path", 4);
  glutAddMenuEntry ("Show/Hide Polytope", 5);
  glutAddMenuEntry ("Filled/Line Polytope", 6);
  otherDisplayMenuID = glutGetMenu ();

  glutCreateMenu (RendererMainMenuWrapper);
  glutAddSubMenu ("Algorithm Display", algorithmDisplayMenuID);
  glutAddSubMenu ("Algorithm Animate", algorithmAnimateMenuID);
  glutAddMenuEntry ("Run/Stop FlyBy", 1);
  glutAddSubMenu ("Other Display Options", otherDisplayMenuID);
  glutAddMenuEntry ("Quit", 2);
  glutAttachMenu (GLUT_RIGHT_BUTTON);
}

void
Renderer::readSrcDestFile (const char *infilename)
{
  FILE *fl;
  int params;

  fl = fopen (infilename, "r");
  if (fl == NULL)
    {
      fprintf (stderr, "Unable to read file: [%s]\n", infilename);
      exit (-1);
    }

  params = fscanf (fl, "%lg %lg %lg\n%lg %lg %lg\n",
                   &(m_srcPoint[0]),
                   &(m_srcPoint[1]),
                   &(m_srcPoint[2]),
                   &(m_destPoint[0]),
                   &(m_destPoint[1]),
                   &(m_destPoint[2]));
  if (params != 6)
    {
      fprintf (stderr, "Unable to read coordinates from file [%s]\n",
               infilename);
      exit (-1);
    }
  params = fscanf (fl, "%d %d", &m_srcFace, &m_destFace);
  if (params != 2)
    {
      fprintf (stderr, "Unable to read faces. from file [%s]\n", infilename);
      exit (-1);
    }
  // Le: add to debug
  
  cout <<"here ("<< m_srcPoint[0] <<"," <<  m_srcPoint[1]<< "," << m_srcPoint[2] <<")";
  cout <<"("<< m_destPoint[0] <<"," <<  m_destPoint[1] <<"," <<  m_destPoint[2] <<")";
  cout << m_srcFace <<"," <<  m_destFace << endl;

  // Le: end

  fclose (fl);
}


void
Renderer::MainLoop (string filename,
                    double epsilon,
                    string infilename,
                    int argc, char **argv, bool f_test,
                    bool _f_print, bool f_super_test, bool _f_white)
{
  //printf( "Rendered::MainLoop\n" );
  m_filename = filename;
  f_print = _f_print;
  f_white = _f_white;

  double new_epsilon;

  new_epsilon = 3.1 / (double ((int)(3.14159265 / epsilon)));
  printf ("Fixed epsilon: %g\n", new_epsilon);
  epsilon = new_epsilon;

  //DDD( "a" );
  //Path m_path;
  m_path = new Path (epsilon);

  //DDD( "b" );
  ifstream inFile (infilename.c_str ());
  bool skipdisplay;
  inFile.good ();

  skipdisplay = false;

  // read the polytope data
  //DDD( "c" );
  m_faceNum = m_path->ReadPoly (m_filename, (const Face * const *&) m_faces);

  // do preprocess path stuff
  //DDD ("d");
  if (!m_useHershSuri)
    m_path->ConstructPathGraph (m_debugPoints);
  //DDD ("a");

  // if we read no faces, then there's a problem
  if (m_faceNum <= 0)
    {
      cerr << "read poly failed\n";
      exit (0);
    }

  if (f_super_test)
    skipdisplay = true;

  if (!skipdisplay) // trang's: run in graphical mode ...
    {
      //initialize glut stuff
      glutInit (&argc, argv);
      glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
      glutInitWindowSize (m_windowWidth, m_windowHeight);
      glutInitWindowPosition (0, 0);
      glutCreateWindow ("Multiple shooting approach for shortest gentle paths ... ");

      init ();
      initMainMenu ();

      glutReshapeFunc (RendererReshapeWrapper);
      glutKeyboardFunc (RendererKeyboardWrapper);
      glutDisplayFunc (RendererDisplayWrapper);
      glutMouseFunc (RendererMouseWrapper);

      // trang: handling direction keys
      glutSpecialFunc(RendererHandleSpecialKeypressWrapper);

      if (f_test)
        {
          readSrcDestFile (infilename.c_str ());

          printf ("My Calling calcpth - B\n");
          fflush (stdout);
          cout << "my source from file: (" <<  m_srcPoint[0] << "," << m_srcPoint[1] << "," << m_srcPoint[2] 
              << ") at:"<< m_srcFace << endl;
          cout << "my destination from file: (" << m_destPoint[0] << "," << m_destPoint[1] << "," << m_destPoint[2]
              << ") at:"<< m_destFace <<  endl;

          // Le: 6-8-2022
          m_path->solve_finish = true;
          // Le:    

          m_path->multispCalcPath ( m_srcPoint, m_destPoint,
                                    m_srcFace, m_destFace, m_useHershSuri,
                                    shortest_path, multiSP_path, multiSP_path_first,
                                    slices_ );

          // Le: 6-8-2022
          // run LiuWong algorithm
          m_path->liuWongCalcPath ( m_srcPoint, m_destPoint,
                            m_srcFace, m_destFace);
          
          m_path->solve_finish = false;
          // Le:
        }

      // Le: 6-8-2022
      globalSGP = NULL;
      //
      glutPostRedisplay ();
      // Le: 6-8-2022

      glutMainLoop ();
    }
  else // trang's: run in console mode ...
    {
      readSrcDestFile (infilename.c_str ());
      printf ("Calling calcpth - B\n");
      fflush (stdout);

      cout << " source dest are choosen by click ... my source: (" <<  m_srcPoint[0] << "," << m_srcPoint[1] << "," << m_srcPoint[2] 
            << ") at:"<< m_srcFace << endl;
      cout << "my destination: (" << m_destPoint[0] << "," << m_destPoint[1] << "," << m_destPoint[2]
           << ") at:"<< m_destFace <<  endl;

      // Le: 6-8-2022
      m_path->solve_finish = true;
      // Le:

      m_path->multispCalcPath ( m_srcPoint, m_destPoint,
                                m_srcFace, m_destFace, m_useHershSuri,
                                shortest_path, multiSP_path, multiSP_path_first,
                                slices_ );

      // Le: 6-8-2022
      m_path->solve_finish = false;
      // Le:

    }
    // Le: 6-8-2022
    globalSGP = NULL;
    // Le: 6-8-2022

    delete m_path;

}

void
Renderer::useHershSuri ()
{
  m_useHershSuri = true;
}



void
Renderer::LiuWong (string filename,
                    double epsilon,
                    string infilename, bool _f_LiuWong)
{
  cout << " \t\t\t running LiuWong algorithm.." << endl;
  m_filename = filename;
  //bool f_LiuWong = _f_LiuWong;


  //DDD( "a" );
  //Path m_path;
  m_path = new Path (epsilon);

  //DDD( "b" );
  ifstream inFile (infilename.c_str ());
  inFile.good ();

  //skipdisplay = false;

  // read the polytope data
  //DDD( "c" );
  m_faceNum = m_path->ReadPoly (m_filename, (const Face * const *&) m_faces);

  // if we read no faces, then there's a problem
  if (m_faceNum <= 0)
    {
      cerr << "read poly failed\n";
      exit (0);
    }

  readSrcDestFile (infilename.c_str ());

  fflush (stdout);
  cout << "LiuWong my source from file: (" <<  m_srcPoint[0] << "," << m_srcPoint[1] << "," << m_srcPoint[2] 
      << ") at:"<< m_srcFace << endl;
  cout << "my destination from file: (" << m_destPoint[0] << "," << m_destPoint[1] << "," << m_destPoint[2]
      << ") at:"<< m_destFace <<  endl;

  if (_f_LiuWong)
      m_path->liuWongCalcPath ( m_srcPoint, m_destPoint,
                            m_srcFace, m_destFace);

    delete m_path;

}

/* Renderer.cc - end of file */
