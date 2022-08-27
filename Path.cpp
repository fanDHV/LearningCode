/*****************************************************************************
 * Path.cc
 *
 * Path is the main interface with the outside world.  Holds the path graph
 * and does everything except for construct it.
 *****************************************************************************/

#include  <stdio.h>
#include  <memory.h>
#include  <assert.h>
#include  <stdlib.h>
#include <fstream>

#include  "sarray.h"

#include  "Functions.h"
#include  "Path.h"
#include  "Face.h"
#include  "Edge.h"
#include  "Vert.h"
#include  "Poly.h"
#include  <math.h>
#include  "history.h"

#include  "rgeneric.h"

#include "slice-mngt.h"
#include "global-settings.h"

// for number decimals
#include <iomanip>
// for calculate running time
#include <time.h>


SourceDest::SourceDest (double *src, double *dest,  int srcFaceIdx,
                         int destFaceIdx, Face * const *faces, Poly * poly,
      bool isprojected )
{
    // trang: avoiding unused variable
isprojected = isprojected;

assert (poly != NULL);
printf( "srcdest from input file .in: ");
printf( "src:(%g, %g, %g) at %d\n", src[ 0 ], src[ 1 ], src[ 2 ], srcFaceIdx );
printf( "dst:(%g, %g, %g) at %d\n", dest[ 0 ], dest[ 1 ], dest[ 2 ],  destFaceIdx );
  


  
  cout << "srcface: " << srcFaceIdx << ", dstface: " << destFaceIdx << endl;
  m_srcFace = faces[srcFaceIdx];
  m_destFace = faces[destFaceIdx];

  if ( isprojected ){
    //cout << "Before projected: " << Point( RT (src[0]),
    //                                       RT (src[1]),
    //                                       RT (src[2]) ) << endl;
    //cout << "Before projected: " << Point( RT (dest[0]),
    //                                       RT (dest[1]),
    //                                       RT (dest[2]) ) << endl;

    m_srcPoint = m_srcFace->project (Point (RT (src[0]),
              RT (src[1]),
              RT (src[2])));
    m_destPoint = m_destFace->project (Point (RT (dest[0]),
                RT (dest[1]),
                RT (dest[2])));
    } 
    else {
      m_srcPoint = Point( RT(src[0]), RT(src[1]), RT(src[2]) );
      m_destPoint = Point( RT(dest[0]), RT(dest[1]), RT(dest[2]) );

      assert (m_srcFace->isOn (m_srcPoint));
      assert (m_destFace->isOn (m_destPoint));
  }

  printf ("src/dest Points in input format!\n");
  point_print_input_frmt (m_srcPoint);
  printf ("\n");
  point_print_input_frmt (m_destPoint);
  printf ("\n");

  printf ("Face Idx %d %d\n",
          poly->getFaceIndex (m_srcFace), poly->getFaceIndex (m_destFace));

  printf ("\n\n");
  printf ("src :");
  point_print (m_srcPoint);
  printf (" ");
  printf (" %d \n", m_srcFace->getSideSign (m_srcPoint));

  printf ("dst :");
  point_print (m_destPoint);
  printf (" ");
  printf (" %d \n", m_destFace->getSideSign (m_destPoint));

}


Path::Path (double epsilon):
m_epsilon (epsilon),
m_pickNum (0),
m_extraGraphNodes (0),
m_currExtraNode (0),
m_pathLength (0.0)
{
  int iiBound, jjBound;
  m_poly = new Poly;

  // calculate the number of latitude and longitude lines
  if (fmod (CGAL_PI, m_epsilon) == 0.0)
    {
      iiBound = int (CGAL_PI / m_epsilon - 1);
      jjBound = int (2 * CGAL_PI / m_epsilon);
    }
  else
    {
      iiBound = int (CGAL_PI / m_epsilon);
      jjBound = int (2 * CGAL_PI / m_epsilon + 1);
    }

  grid.init (iiBound, jjBound);

  //m_numNodes = 0;

  // reserve 8 extra nodes for more picks, if we need more, we will increase
  // the size of the array
  // the extra four nodes are for the north/south poles and the source
  // and destination points
  //m_sizeGraph = m_iiBound * m_jjBound + 12;

  // create the graph
  //m_graph = new GraphNode*[m_sizeGraph];

  /*----------- trang: options for solving in step by step ----------*/
  // solve to finish ...
  solve_finish = false;
  // solve in step by step ...
  show_cut_slices = false;
  show_initial_path = false;
  finish = false;
}


Path::~Path ()
{
  delete m_poly;
  //for  ( int  i = 0; i < m_numNodes; i++)
  //  delete [] m_graph[i];
  //delete m_graph;
  grid.term ();

  GraphNode *curr = m_extraGraphNodes, *curr2;
  while (curr)
    {
      curr2 = curr->m_next;
      delete curr;
      curr = curr2;
    }
}


int
Path::ReadPoly (const std::string & filename, const Face * const *&faces)
{
  //DDD( "1" );
  ifstream ifstrIn (filename.c_str ());

  //DDD( "2" );
  if (!ifstrIn.good ())
    {
      cerr << "Cannot open " << filename << "!" << endl;
      exit (1);
    }

  //DDD( "3" );
  m_poly->Read (ifstrIn);

  //DDD( "4" );
  faces = m_poly->Faces ();

  m_poly->doPolyDump (NULL);

  //DDD( "5" );
  return m_poly->NumFaces ();
}


// adds a node to the path graph at the specified index.  if the index is
// larger than the graph, then increase the size of the path graph by 10
void
Path::AddNodeToGraph (GraphNode * gn, int idx)
{
  grid.addNodeToGraph (gn, idx);
#ifdef  OLD_JANK
  if (idx >= m_sizeGraph)
    {
      // this only reserves space for 10 more now because we probably
      // aren't going to pick a bunch, but if so, another method can
      // be used
      GraphNode **temp = new GraphNode *[m_sizeGraph + 10];
      for (int i = 0; i < m_sizeGraph; i++)
        temp[i] = m_graph[i];
      m_sizeGraph += 10;

      delete m_graph;
      m_graph = temp;
    }

  m_graph[idx] = gn;
  m_numNodes++;
#endif // OLD_JUNK
}


// called to preprocess the polytope: generates the path graph and calculates
// distances between connected nodes in the graph
void
Path::ConstructPathGraph (std::list < AnimateVecInfo * >&debugPoints)
{
  m_poly->ConstructPathGraph (this, m_epsilon, grid, debugPoints);
  CalcDists ();
}



// void Path::algHershSuri (double *src, double *dest,
//                     int srcFaceIdx, int destFaceIdx,
//                     bool bUseHershSuri, SPath * &p_path)
void Path::algHershSuri ( bool bUseHershSuri, SPath * &p_path)
{

    // trang: avoiding unused variable
    bUseHershSuri = bUseHershSuri;

  // run the Hershberger-Suri algorithm
  GraphNode *curr, *path;

  m_Hs = Plane (m_srcdest->m_srcPoint, -(m_srcdest->m_srcFace->NormalDir ()));
  m_Ht = Plane (m_srcdest->m_destPoint,
                -(m_srcdest->m_destFace->NormalDir ()));

  // if the angle is less than 2pi/3, then we can use a 2-plane
  // wedge
  //double angle = (m_srcdest->m_srcFace->NormalDir().to_vector()).
  //    angle(m_srcdest->m_destFace->NormalDir().to_vector());
  double angle = angleBetween( m_srcdest->m_srcFace->NormalDir (),
                               m_srcdest->m_destFace->NormalDir ());

  NNLocation snnloc;
  snnloc.m_loc = m_srcdest->m_srcFace;
  snnloc.m_locType = NNLocation::FACE;
  curr = new GraphNode (m_srcdest->m_srcPoint, snnloc, m_srcdest->m_srcPoint);

  if ( angle < 2 * CGAL_PI / 3 )
    {
      // we can use the 2-plane wedge
      AddExtraGraphNode (curr,
                         RidgePoint (m_srcdest->m_srcPoint,
                                     m_srcdest->m_destPoint, m_Hs, m_Ht));

      NNLocation dnnloc;
      dnnloc.m_loc = m_srcdest->m_destFace;
      dnnloc.m_locType = NNLocation::FACE;
      curr->m_next->m_next = new
        GraphNode (m_srcdest->m_destPoint, dnnloc, m_srcdest->m_destPoint);
      curr->m_next->m_next->m_next = 0;

      path = curr;

      ProjectPath (path);

      p_path = new SPath (path, m_poly);
      m_pathLength = p_path->getAprxLength ();
      return;
    }

  // we have to use the 3-plane wedge
  // this is Neill's code

  // make the "i" vector
  Vec n1 = m_srcdest->m_srcFace->NormalDir ();
  Vec n2 = m_srcdest->m_destFace->NormalDir ();

  Vec n1u = Unit (n1);
  Vec n2u = Unit (n2);
  //Vec i = CrossProductDir(CrossProductDir(n1u, n2u,31415),
  //                        n1u+n2u,7888833);
  Vec i = CrossProduct (CrossProduct (n1u, n2u), n1u + n2u); //tich co huong

  Simplify (i);

  list < Edge * >L = m_poly->HorizonEdges (i);

  Edge *e;
  Point *currp, *best;
  RT currLength, bestLength;
  Plane *mid, *tempplane;
  bool first = true;

  currp = new Point[4];
  best = new Point[4];

  //forall (e, L) {
  for (list < Edge * >::iterator ei = L.begin (); ei != L.end (); ei++)
    {
      e = *ei;
      mid = MiddlePlane (e, i);

      // make sure normal points towards polytope
      if (sign (DotProduct (normalVector (*mid),
                            e->V1 ()->Location () - m_poly->Center ())) < 0)
        {
          tempplane = new Plane (mid->point (), -(normalVector (*mid)));
          delete mid;
          mid = tempplane;
        }

      currLength = ThreeWedge (currp, mid);

      if (first || (currLength < bestLength))
        {
          first = false;

          for (int j = 0; j < 4; j++)
            best[j] = currp[j];

          bestLength = currLength;
        }
      delete mid;
    }

  AddExtraGraphNode (curr, best[1]);
  AddExtraGraphNode (curr->m_next, best[2]);

  NNLocation dnnloc;
  dnnloc.m_loc = m_srcdest->m_destFace;
  dnnloc.m_locType = NNLocation::FACE;
  curr->m_next->m_next->m_next = new
    GraphNode (m_srcdest->m_destPoint, dnnloc, m_srcdest->m_destPoint);
  curr->m_next->m_next->m_next->m_next = 0;

  path = curr;

  delete currp;
  delete best;

  //ProjectPath(path, true);
  p_path = new SPath (path, m_poly);
  m_pathLength = p_path->getAprxLength ();
}


// void Path::algHershSuriExt (double *src, double *dest,
//          int srcFaceIdx, int destFaceIdx,
//          bool bUseHershSuri, SPath * &p_path)
void Path::algHershSuriExt ( bool bUseHershSuri, SPath * &p_path)
{
    // trang: avoiding unused variable
    bUseHershSuri = bUseHershSuri;

  //printf( "Path::algHershSurExt called!\n" );
  //fflush( stdout );

  // run the Hershberger-Suri algorithm
  //GraphNode * curr, * path;

  m_Hs = Plane (m_srcdest->m_srcPoint, -(m_srcdest->m_srcFace->NormalDir ()));
  m_Ht = Plane (m_srcdest->m_destPoint,
                -(m_srcdest->m_destFace->NormalDir ()));

  //assert( m_Hs.side_of( m_srcdest->m_srcPoint ) == 0 );
  //assert( m_Ht.side_of( m_srcdest->m_destPoint ) == 0 );
  assert (m_Hs.has_on (m_srcdest->m_srcPoint));
  assert (m_Ht.has_on (m_srcdest->m_destPoint));

  // if the angle is less than 2pi/3, then we can use a 2-plane
  // wedge
  //double angle = (m_srcdest->m_srcFace->NormalDir().to_vector()).
  //    angle(m_srcdest->m_destFace->NormalDir().to_vector());
  double angle = angleBetween( m_srcdest->m_srcFace->NormalDir (),
                               m_srcdest->m_destFace->NormalDir ());

  SPathNode *sp_start, *sp_end;

  sp_start = new SPathNode (m_srcdest->m_srcFace, m_srcdest->m_srcPoint);
  sp_end = new SPathNode (m_srcdest->m_destFace, m_srcdest->m_destPoint);

  p_path = new SPath (sp_start, m_poly);
  p_path->dump ();

  // p_path now contains the starting point.

  if (angle < 2 * CGAL_PI / 3)
    {
      p_path->pushTop (sp_end);
      p_path = p_path->shortcut ();
      printf ("small angle!\n");
      return;
    }

  printf ("big angle!\n");
  fflush (stdout);

  // we have to use the 3-plane wedge
  // this is Neill's code

  // make the "i" vector
  Vec n1 = m_srcdest->m_srcFace->NormalDir ();
  Vec n2 = m_srcdest->m_destFace->NormalDir ();

  Vec n1u = Unit (n1);
  Vec n2u = Unit (n2);
  Vec i = CrossProduct (CrossProduct (n1u, n2u), n1u + n2u);

  Simplify (i);

  list < Edge * >L = m_poly->HorizonEdges (i);

  Edge *e;
  Point *currp, *best;
  RT currLength, bestLength;
  Plane *mid, *tempplane;
  bool first = true;

  currp = new Point[4];
  best = new Point[4];

  //forall (e, L) {
  for (list < Edge * >::iterator ei = L.begin (); ei != L.end (); ei++)
    {
      e = *ei;
      mid = MiddlePlane (e, i);

      // make sure normal points towards polytope
      if (sign (DotProduct (normalVector (*mid),
                            e->V1 ()->Location () - m_poly->Center ())) < 0)
        {
          tempplane = new Plane (mid->point (), -(normalVector (*mid)));
          delete mid;
          mid = tempplane;
        }

      currLength = ThreeWedge (currp, mid);

      if (first || (currLength < bestLength))
        {
          first = false;

          for (int j = 0; j < 4; j++)
            best[j] = approx_pnt (currp[j]);

          bestLength = currLength;
        }
      delete mid;
    }

  NNLocation loc_out;
  Point tmp;

  memset (&loc_out, 0, sizeof (NNLocation));


  p_path->pushTop (loc_out, best[1]);
  p_path->pushTop (loc_out, best[2]);

  p_path->dump ();



  //printf( "BRUTE 1!\n" ); fflush( stdout );

  p_path->pushTop (sp_end);
  printf ("BRUTE! 2\n");
  fflush (stdout);
  p_path->shortcut_brute_force ();
  //p_path->shortcut_brute_force();
  printf ("BRUTE!\n");
  fflush (stdout);
}


// the fun stuff: calculate the path
// don't use the function below? Use spCalcPath(..)
void Path::CalcPath ( double *src, double *dest,
          int srcFaceIdx, int destFaceIdx,
          bool bUseHershSuri,
          SPath * &s_path )
{
  // get source and destination point and face information
  m_srcdest = new SourceDest( src, dest, srcFaceIdx, destFaceIdx,
                              m_poly->Faces (), m_poly);

  printf( "Path::calcPath called\n" );
  if (bUseHershSuri)
    {
      // algHershSuriExt (src, dest, srcFaceIdx, destFaceIdx,
      //                  bUseHershSuri, s_path);
      algHershSuriExt ( bUseHershSuri, s_path);
    }
  else
    {
      // Construct the path graph
      //  if (!m_pickNum)
      //    m_poly->ConstructPathGraph(this, m_epsilon, m_iiBound, m_jjBound);

      // add the source and destination points to the graph
      debug ("Adding Source and Destination Points ...");
      printf( "before adding source and destination points\n" );
      //fflush( stdout );
      m_poly->AddSourceDest (this, m_srcdest, m_epsilon, grid, m_pickNum++);
      printf( "after addition source and destination points\n" );
      //fflush( stdout );

      // calculate the distances to and from every node to all other
      // connected nodes
      printf("Calculating Distances ...\n");
      CalcDists ();
      printf("After calculating distances...\n");
      // run Djikstra's to find the shortest path in the graph
      GraphNode *gn = FindShortestPathInGraph (grid.getVerticesNum () - 2,
                                               grid.getVerticesNum () - 1);
      printf( "After shortest path ...\n" );
      //m_numNodes-2, m_numNodes-1);

      // Shortcut the shortest path and return it
      //s_path = new  SPath( ShortcutPath(gn), m_poly );

      /**********************************/
      /* Hoai - 4/12/2012 - no shortcut */
      /**********************************/
      printf( "Before calculating shortcut...\n" );
      s_path = ShortcutPathSP (gn);
      printf( "After calculating shortcut...\n" );
      s_path->dump();
      //printf( "source - dest -?\n" );

      /*printf( "src :" );
         point_print( m_srcdest->m_srcPoint );
         printf( " " );
         printf( " %d \n", m_srcdest->m_srcFace->getSideSign(
         m_srcdest->m_srcPoint ) );

         printf( "dst :" );
         point_print( m_srcdest->m_destPoint );
         printf( " " );
         printf( " %d \n", m_srcdest->m_destFace->getSideSign(
         m_srcdest->m_destPoint ) );
       */
    }

  m_pathLength = s_path->getAprxLength ();

  cerr << "Path length = " << m_pathLength << endl;
}

/****************** Hoai - 22/1/2013 ********************/
/* calculate the path by multiple shooting method       */
/********************************************************/
void Path::multispCalcPath ( double *src, double *dest,
                             int srcFaceIdx, int destFaceIdx,
                             bool bUseHershSuri,
                             SPath * &apath, vector<Point> *&multisp_path,
                             vector<Point> *&multisp_path_first,
                             deque< deque< Seg > > *&myslices )
{
    // trang: avoiding unused variable
    bUseHershSuri = bUseHershSuri;

    // Le 6-8-2022 ----------------------------------------------------------------------
    // get set of faces of the terrain
    Face * const * faces = m_poly->Faces();
    SourceDest mL_srcdest;
    //mL_srcdest.m_srcPoint = (Point (RT (src[0]), RT (src[1]), RT (src[2])));
    //mL_srcdest.m_destPoint = (Point (RT (dest[0]), RT (dest[1]), RT (dest[2])));
    mL_srcdest.m_srcFace = faces[srcFaceIdx];
    mL_srcdest.m_destFace = faces[destFaceIdx];
    mL_srcdest.m_srcPoint = mL_srcdest.m_srcFace->project (Point (RT (src[0]),
              RT (src[1]),
              RT (src[2])));
    mL_srcdest.m_destPoint = mL_srcdest.m_destFace->project (Point (RT (dest[0]),
                RT (dest[1]),
                RT (dest[2])));

    cout << "   - srcFaceIdx and  destFaceIdx -> " << srcFaceIdx << ";" << destFaceIdx << endl;

    // m_srcdest = new SourceDest ( src, dest, srcFaceIdx, destFaceIdx, m_poly->Faces (), m_poly);
    // avoid error generating by shalow copy Le 5-8-2022
    // SourceDest mL_srcdest( src, dest, srcFaceIdx, destFaceIdx, m_poly->Faces (), m_poly);
    m_srcdest = new SourceDest(mL_srcdest);
    // avoid error generating by shalow copy Le 5-8-2022 --------------------------------

    // for calculate running time
    clock_t start, end;
    double used_time;

    //printf( "Path::multiple shooting calcPath called\n" );
    start = clock();

    printf( "before adding source and destination points\n" );
    m_poly->AddSourceDest (this, m_srcdest, m_epsilon, grid, m_pickNum++);

    // calculate the distances to and from every node to all other
    // connected nodes
    //printf("Calculating Distances ...\n");
    CalcDists ();

    // run Djikstra's to find the shortest path in the graph
    //printf("Run Dijkstra to find shortest path ...\n");
    GraphNode *gn = FindShortestPathInGraph ( grid.getVerticesNum () - 2,
                                              grid.getVerticesNum () - 1);

    // shortcut shortest path
    //printf( "Before calculating shortcut...\n" );

    // trang: obtaining all slices through vertices of polytope???
    //sliceMngt slices( m_poly );
    //Trang's in the order to construct a palne throught source point:
     sliceMngt slices( m_poly, m_srcdest );

    // solves iteratively ...
    if (solve_finish) // trang: solve to finish
    {
        // trang's debug:
        cout << " solve to finish the probem..." << endl; 

        apath = ShortcutPathSP (gn); // trang's comment
        apath->dump();
        cerr << "Agarwal's length: " << apath->getAprxLength() << endl;

        // Hoai - 6/1/2013 - no shortcut
        slices.detectEffectivePlanes( m_srcdest );
        slices.getSlices( myslices );

        //
        slices.detectSlicePaths( apath );

        // hoai's old initializing
        //slices.vibrateSlicePaths( m_srcdest, srcFaceIdx, destFaceIdx, this,
        //                         VIBRATE_RATIO, VIBRATE_NEDGES );

        // trang: initialize a path through slices
        slices.initializePath( m_srcdest, srcFaceIdx, destFaceIdx, this,
                                  VIBRATE_RATIO, VIBRATE_NEDGES );

        // compute approximate shortest path using a given shortest path method
        if ( multisp_path != NULL )
          delete multisp_path;
        if ( multisp_path_first != NULL )
          delete multisp_path_first;
        //
        multisp_path_first = slices.getShootingPath();
        multisp_path = slices.multiShooting( m_srcdest, srcFaceIdx, destFaceIdx, this );

        //
        finish = true;
    }
    else // trang: solve in step by step
    {
      if ( show_initial_path )
      { 
          cout << " solve in step by step... and show initial path..." << endl; 

          apath = ShortcutPathSP (gn);
          apath->dump();
          cerr << "Agarwal's length: " << apath->getAprxLength() << endl;

          //
          show_initial_path = false;
          show_cut_slices = true;
      }
      else // trang: solve in step by step
          if ( show_cut_slices )
          {
              cout << "solve in step by step... DO NOT show initial path... and show cut slices..." << endl;
              apath = ShortcutPathSP (gn);
              apath->dump();
              slices.detectEffectivePlanes( m_srcdest );
              slices.getSlices( myslices );

              //hoai: slices.detectPoorSlicePaths( apath );
              slices.detectSlicePaths( apath );
              slices.vibrateSlicePaths( m_srcdest, srcFaceIdx, destFaceIdx, this,
                      VIBRATE_RATIO, VIBRATE_NEDGES );

              //
              show_cut_slices = false;
              finish = true;
          }
        else // trang: solve in step by step
          {
              //
              cout << "solve in step by step... DO NOT show initial path...and DO NOT show cut slices..." << endl;

              apath = ShortcutPathSP (gn);
              apath->dump();
              cerr << "Agarwal's length: " << apath->getAprxLength() << endl;

              slices.detectEffectivePlanes( m_srcdest );
              slices.getSlices( myslices );

              // hoai: slices.detectPoorSlicePaths( apath );
              slices.detectSlicePaths( apath );
              slices.vibrateSlicePaths( m_srcdest, srcFaceIdx, destFaceIdx, this,
                      VIBRATE_RATIO, VIBRATE_NEDGES );

              // trang: initialize a path through slices
              slices.initializePath( m_srcdest, srcFaceIdx, destFaceIdx, this,
                                        VIBRATE_RATIO, VIBRATE_NEDGES );

              // compute approximate shortest path using a given shortest path method
              if ( multisp_path != NULL )
                delete multisp_path;
              if ( multisp_path_first != NULL )
                delete multisp_path_first;

              //
              multisp_path_first = slices.getShootingPath();
              multisp_path = slices.multiShooting( m_srcdest, srcFaceIdx, destFaceIdx, this );

              //
              finish = false;
        }
    }
    end = clock();
    used_time = (double)(end - start) / CLOCKS_PER_SEC; 

    //%%%%%%%%% Write into file 6-8-2022
    ofstream ofs("result.txt", ios::app);
    if(!ofs)
    {
      cerr << "Error: file not opened." << endl;
      exit(1);
    }
  
    ofs << "\n Running time is....: " << setprecision(5) << fixed << used_time;
      
    //Close file
    ofs.close();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // ----------------------- 2015/08/30: computing the global solution ------------------------
    if (globalComputing_)
    {
      cout << endl <<endl;
      cout << "%% Computing the global shortest gentle path ..." << endl;

    
      SPath * global_SGP = NULL;
      //
      this->calcSGP( m_srcdest, false, global_SGP );

      //
      // final information
      cout << "====================================== Global solution ===================================" << endl;
      cout << "%% Global SGP(p,q) obtained ..." << endl;
      cout << "%% Source and destination: " << endl;
      cout << "\t "; point_print( m_srcdest->m_srcPoint  ); cout << endl;
      cout << "\t "; point_print( m_srcdest->m_destPoint ); cout << endl;
      //
      cout << "%% The final path: " << endl;
      global_SGP->printNodeList();
      //cout << "\t - initial path:" << initiallen << endl;
      cout << "\t - approximate length:" << global_SGP->getAprxLength() << endl;
      // ------------------ 2015/08/30: end of computing the global solution ------------------------
    
      //%%%%%%%%% Write into file 2022-4-6
      ofstream ofs("result.txt", ios::app);
      if(!ofs)
      {
        cerr << "Error: file not opened." << endl;
        exit(1);
      }
        
      ofs << "\n Slope is...." << SIN_THETA ;
        
      ofs << "\n  === Global Solution === SGP(p,q) obtained has length..."
      << setprecision(5) << fixed << global_SGP->getAprxLength() ;
        
      //Close file
      ofs.close();
      //%%%%%%%%%%%%%%%%
    }

    if (globalClassicalComputing_)
    {
        cout << endl <<endl;
        cout << "%% Computing the global classical adjusted shortest path ..." << endl;

        SPath * global_adjusted_SP = NULL;
        //
        this->calcSP( m_srcdest, false, global_adjusted_SP );

        //
        // final information
        cout << "====================================== Global solution ===================================" << endl;
        cout << "%% global classical adjusted shortest path SP(p,q) obtained ..." << endl;
        cout << "%% Source and destination: " << endl;
        cout << "\t "; point_print( m_srcdest->m_srcPoint  ); cout << endl;
        cout << "\t "; point_print( m_srcdest->m_destPoint ); cout << endl;
        //
        cout << "%% The final path: " << endl;
        global_adjusted_SP->printNodeList();
        //cout << "\t - initial path:" << initiallen << endl;
        cout << "\t - approximate length:" << global_adjusted_SP->getAprxLength() << endl;
        // ------------------ 2015/08/30: end of computing the global solution ------------------------
    
      //%%%%%%%%% Ghi vao file 2022-4-6
      //Mở file bằng ofstream 

      ofstream ofs("result.txt", ios::app);
      //Kiểm tra file đã mở thành công hay chưa
      if(!ofs){
        cerr << "Error: file not opened." << endl;
        exit(1);
      }
        
      //Ghi vào file 
        
      ofs << "\n  Adjusted SP(p,q) has length..."
      << setprecision(5) << fixed << global_adjusted_SP->getAprxLength() ;
        
      //Đóng file
      ofs.close();
      //%%%%%%%%%%%%%%%%
    }
    // Le 6-8-2022
    delete m_srcdest;

}

/*--------------------------------------------------------------------------*/
// the fun stuff: calculate the path
void Path::spCalcPath ( SourceDest *srcdest, bool bUseHershSuri, SPath * &s_path )
{
  // printf( "Path::spCalcPath called\n" );
  if (bUseHershSuri)
  {

      algHershSuriExt( bUseHershSuri, s_path);
  }
  else
  {

      m_poly->AddSourceDest (this, srcdest, m_epsilon, grid, m_pickNum++);

      CalcDists ();

      GraphNode *gn = FindShortestPathInGraph ( grid.getVerticesNum () - 2,
                                                grid.getVerticesNum () - 1 );

      s_path = ShortcutPathSP( gn );
  }

  //
  m_pathLength = s_path->getAprxLength ();
}

/*--------------------------------------------------------------------------*/
// for SGP distance
void Path::spCalcSGPPath ( SourceDest *srcdest, bool bUseHershSuri, SPath * &s_path )
{
  // printf( "Path::spCalcSGPPath called\n" );
  if (bUseHershSuri)
  {
      // algHershSuriExt( src, dest, srcFaceIdx, destFaceIdx,
      //         bUseHershSuri, s_path);
      algHershSuriExt( bUseHershSuri, s_path);
  }
  else
  {
      // add the source and destination points to the graph
      // printf( "before adding source and destination points\n" );
      // m_poly->AddSourceDest (this, m_srcdest, m_epsilon, grid, m_pickNum++);
      m_poly->AddSourceDest (this, srcdest, m_epsilon, grid, m_pickNum++);
      // printf( "after addition source and destination points\n" );

      // calculate the distances to and from every node to all other
      // connected nodes
      // printf("Calculating Distances ...\n");
      CalcSGPDists ();
      // printf("After calculating distances...\n");

      // run Djikstra's to find the shortest path in the graph
      // old's one:
      GraphNode *gn = FindShortestGentlePathInGraph ( grid.getVerticesNum () - 2,
                                                grid.getVerticesNum () - 1 );

      // trang: find non increasing shortest path from p' to q. This is not natural!!!
      //GraphNode * gn = FindNonIncreasingShortestPathInGraph( grid.getVerticesNum () - 2,
      //                                                       grid.getVerticesNum () - 1 );

      // printf( "After shortest path ...\n" );

      // shortcut shortest path
      // printf( "Before calculating shortcut...\n" );
      s_path = ShortcutPathSP( gn );
      // printf( "After calculating shortcut...\n" );
      // s_path->dump();
  }

  //
  m_pathLength = s_path->getAprxLength ();
}

/****** Le: Function to find shoting points in case of lateral and vertex (5/10/2019)*******/
Point Path::spCalcSGPPath(SourceDest *srcdest, deque <Seg> slices_s, Point atVertex_s)
{
  // add the source and destination points to the graph
  // printf( "before adding source and destination points\n" );
  // m_poly->AddSourceDest (this, m_srcdest, m_epsilon, grid, m_pickNum++);
  m_poly->AddSourceDest (this, srcdest, m_epsilon, grid, m_pickNum++);
  // printf( "after addition source and destination points\n" );

  // calculate the distances to and from every node to all other
  // connected nodes
  // printf("Calculating Distances ...\n");
  CalcSGPDists ();
  // printf("After calculating distances...\n");

  // run Djikstra's to find the shortest path in the graph
  // old's one:
  GraphNode *gn = FindShortestGentlePathInGraph ( grid.getVerticesNum () - 2,
                                            grid.getVerticesNum () - 1 );

  // trang: find non increasing shortest path from p' to q. This is not natural!!!
  //GraphNode * gn = FindNonIncreasingShortestPathInGraph( grid.getVerticesNum () - 2,
  //                                                       grid.getVerticesNum () - 1 );

  // printf( "After shortest path ...\n" );

  // shortcut shortest path
  // printf( "Before calculating shortcut...\n" );
  SPath path_ = ShortcutPathSP( gn );
  SPathNode *ptr = path_.first();
  deque<Point> *edit_path = new deque<Point>;
  edit_path->push_back( ptr->point() );
  assert ( ptr->point().z() > slices_s.front().source().z() );

  while (ptr != NULL)
  {
    if ( fabs( ptr->point().z().to_double() - slices_s.front().source().z().to_double() )  < 1E-3 )
    {
      return ptr->point();
    }
    else if ( ptr->point().z() > slices_s.front().source().z() )
    {
      // save point to slice path
      edit_path->push_back( ptr->point() );
      // increase pointer
      ptr = ptr->next();
      continue;
    }

    else // if ( ptr->point().z() < slices_s.front().source().z() )
    {
        // find intersection of the slice and path (only current segment)
        Seg pseg( edit_path->back(), ptr->point() );
        int k = slices_s.size();
       //
        for( int i = 0; i < k; i++ )
        {
          CGAL::Object obj = CGAL::intersection( pseg, slices_s[ i ] );
          if ( const Point *pobj = CGAL::object_cast< Point >( &obj ) )
            return *pobj;
        }
        cout << "... have no intersection..." << endl;
        PressEnterToContinue();
        return  atVertex_s;
    }
  }
//return *pobj;
}
/****** Le: Function to find shoting points in case of lateral and vertex (5/10/2019) end *******/

/************ Le modifies SDP->SGP to find shortest gentle paths (5/10/2019)**************/
// old path: s_path -> newpath
//use to find SGP between consesutive slices, srcdest are newpoints are updated
void Path::calcSGP(SourceDest *srcdest, bool bUseHershSuri, SPath *&s_path)
{
    // Point p_prime(CGAL::ORIGIN);
    if (bUseHershSuri)
    {
        algHershSuriExt( bUseHershSuri, s_path);
    }
    else
    {
        // add the source and destination points to the graph
        cout << "AddSourceDest ... " << endl;
        m_poly->AddSourceDest (this, srcdest, m_epsilon, grid, m_pickNum++);

         // if ( isDebug_2) PressEnterToContinue();

        // calculate the distances to and from every node to all other connected nodes
         cout << "calculate the distances ..." << endl;
        CalcSGPDists ();

         // if ( isDebug_2) PressEnterToContinue();

        // run Djikstra's to find the shortest path in the graph
        cout << "FindShortestPathInGraph ..." << endl;
        GraphNode *gn = FindShortestGentlePathInGraph ( grid.getVerticesNum () - 2,
                              grid.getVerticesNum () - 1 );

      //  if ( isDebug_2) PressEnterToContinue();

        // shortcut shortest path
         cout << "shortcut shortest path ..." << endl;
        s_path = ShortcutPathSP( gn );
        //
       //  cout << "shortcut shortest path ->dump..." << endl;
          s_path->dump();
    }
      m_pathLength = s_path->getAprxLength ();

}

void Path::calcSP(SourceDest *srcdest, bool bUseHershSuri, SPath *&s_path)
{
    // Point p_prime(CGAL::ORIGIN);
    if (bUseHershSuri)
    {
        algHershSuriExt( bUseHershSuri, s_path);
    }
    else
    {
        // add the source and destination points to the graph
        cout << "AddSourceDest ... " << endl;
        m_poly->AddSourceDest (this, srcdest, m_epsilon, grid, m_pickNum++);

         // if ( isDebug_2) PressEnterToContinue();

        // calculate the distances to and from every node to all other connected nodes
         cout << "calculate the distances ..." << endl;
        CalcDists ();

         // if ( isDebug_2) PressEnterToContinue();

        // run Djikstra's to find the shortest path in the graph
        cout << "FindShortestPathInGraph ..." << endl;
        GraphNode *gn = FindShortestPathInGraph ( grid.getVerticesNum () - 2,
                              grid.getVerticesNum () - 1 );

      //  if ( isDebug_2) PressEnterToContinue();

        // shortcut shortest path
         cout << "shortcut shortest path ..." << endl;
        s_path = ShortcutPathSP( gn );
        //
       //  cout << "shortcut shortest path ->dump..." << endl;
          s_path->dump();
    }
      m_pathLength = s_path->getAprxLength ();

}



//---------- trang: check a shortest path if it is descending -------------
bool Path::isDescending( SPath s_path )
{
    if ( s_path.first()->point().z() < s_path.first()->next()->point().z() )
    {
        return false;
    }
    //
    return true;
}
//-------------------------------------------------------------------------------


/*--------------------------------------------------------------------------*/
SPath * Path::ShortcutPathSP (GraphNode * path)
{
    SPath *new_path;

    new_path = convertToSPath (path);
    if (new_path == NULL)
      return NULL;

    //------------------------------------------------------------------------
    // check for z-monotonicity: needed to our method ???
    for( SPathNode *ptr = new_path->first(); ptr != NULL; ptr = ptr->next() )
    {
        if ( ptr != new_path->first() )
        {
          if ( ptr->point().z() > ptr->prev()->point().z() )
          {
            point_print( ptr->point() );
            cerr << "-bbb-" << flush;
             cerr << "stop to check... ShortcutPathSP()";
            point_print( ptr->prev()->point() );
            cerr << endl;
          }

          // trang - 22/12/2013 - ALLOW z increasing
          // assert( ptr->point().z() <= ptr->prev()->point().z() );
        }
    }
    //------------------------------------------------------------------------

    SPath *shortcuted_path = new_path->shortcut ();
    delete new_path;
    return shortcuted_path;
    // return new_path->shortcut ();
}

/*--------------------------------------------------------------------------*/
SPath * Path::convertToSPath (GraphNode * path)
{
    SPath *sp;
    GraphNode *curr;

    if (path == NULL)
      return NULL;
    sp = new SPath ();

    curr = path;
    while (curr->m_prev != NULL)
    {
        curr->m_prev->m_next = curr;
        curr = curr->m_prev;

        //------------------------------------------------------------------------
        // check for z-monotonicity: needed to our method ???
        if ( curr->m_prev != NULL ){
          if ( curr->m_p.z() >  curr->m_prev->m_p.z() ){
            point_print( curr->m_p );
            cerr << "-bbb-" << flush;
            cerr << "stop to check 1... convertToSPath()";
            point_print( curr->m_prev->m_p );
            cerr << endl;
          }

          // trang - 22/12/2013 - ALLOW z increasing
          // assert( curr->m_p.z() <= curr->m_prev->m_p.z() );
        }
        //------------------------------------------------------------------------
    }

    sp->init_first (curr, m_poly);

    // shortcut the path
    while (curr->m_next != NULL)
    {
        sp->pushTop (curr);
        curr = curr->m_next;
        //------------------------------------------------------------------------
        // check for z-monotonicity: needed to our method ???
        if ( curr->m_next != NULL ){
          if ( curr->m_p.z() <  curr->m_next->m_p.z() ){
            point_print( curr->m_p );
            cerr << "-bbb-" << flush;
             cerr << "stop to check 2... convertToSPath()";
            point_print( curr->m_prev->m_p );
            cerr << endl;
          }

          // trang - 22/12/2013 - ALLOW z increasing
          // assert( curr->m_p.z() >= curr->m_next->m_p.z() );
        }
        //------------------------------------------------------------------------
    }
    sp->pushTop (curr);

    return sp;
}
/*--------------------------------------------------------------------------*/

GraphNode * Path::ReturnNode (int i) const
{
  return grid.getNode (i);
  //return m_graph[i];
}

void Path::PrintGraph ()
{
  grid.PrintGraph ();
#ifdef  XXX
  //NodeSet::iterator si;
  for (int j = 0; j <  m_numNodes; j++)
    {
      cerr << "[" << j << "]: ";
      //      for (int i=0; i<m_graph[j]->m_deg; i++)
      for (int si = 0; si < (int) m_graph[j]->neighbors.size (); si++)
        {

          cerr << m_graph[j]->neighbors[si].m_adj << " ";
        }
      cerr << endl;
    }
#endif
}

/*--------------------------------------------------------------------------*/
// updates all distances in the path graph that are not already computed
void
Path::CalcDists ()
{
  grid.CalcDists ();
}

void
Path::CalcSGPDists ()
{
  grid.CalcSGPDists ();
}

// runs Djikstra's algorithm to find the shortest path in the path graph
// from s to t
GraphNode * Path::FindShortestPathInGraph (const int src, const int trg)
{
  /* Trang
   return grid.FindShortestPathInGraph (src, trg); */
   // Le return to Hoai's code in order to z-coordinate is descreasing
  return grid.FindShortestPathInGraph (src, trg);
}

GraphNode * Path::FindShortestGentlePathInGraph (const int src, const int trg)
{
  /* Trang
   return grid.FindShortestPathInGraph (src, trg); */
   // Le return to Hoai's code in order to z-coordinate is descreasing
  return grid.FindShortestGentlePathInGraph (src, trg);
}



/*---------------------------------------------------------------------------*/
#ifdef  OLD_CODE
// This shortcuts the path, finds the two-plane wedges, and projects the path
// onto the polytope.  Also removes duplicates
SPath *
Path::ShortcutPathSP (GraphNode * path)
{
  GraphNode *curr;
  GraphNode *curr2;
  Vec n1, n2;
  Point rp;
  SPath *sp;

  if (curr == NULL)
    return NULL;

  sp = new SPath ();

  curr = path;
  while (curr->m_prev != NULL)
    {
      curr->m_prev->m_next = curr;
      curr = curr->m_prev;
    }

  /*printf( "start point: " );
     point_print( curr->m_p );
     printf( "\n" );
     fflush( stdout );
   */
  sp->init_first (curr, m_poly);
  //printf( "X" ); fflush( stdout );
  // shortcut the path
  curr2 = curr;
  while (curr->m_next != NULL)
    {
      // printf( "A\n" ); fflush( stdout );
      // use Hersherberger-Suri to find the 2-plane wedge path
      n1 = curr->m_nnloc.m_loc->exactNormal ();
      n2 = curr->m_next->m_nnloc.m_loc->exactNormal ();

      sp->pushTop (curr);
      debug ("n1 = " << n1.to_vector ());
      debug ("n2 = " << n2.to_vector ());
      if (n1 != n2)
        {
          Plane Hs (curr->m_p, -n1);
          Plane Ht (curr->m_next->m_p, -n2);

          rp = RidgePoint (curr->m_p, curr->m_next->m_p, Hs, Ht);
          /*printf( "\n\nRidge point: " );
             point_print( rp );
             printf( "\n" );

             printf( "curr->m_p: " );
             point_print( curr->m_p );
             printf( "\n" );

             printf( "curr->m_next->m_p: " );
             point_print( curr->m_next->m_p );
             printf( "\n" );
             printf( "-----------\n" );
             sp->dump();
             printf( "-----------\n" );

             fflush( stdout );
           */
          ProjectPathSP (curr, rp, *sp);
        }

      curr = curr->m_next;
    }
  sp->pushTop (curr);

  return sp;
}
#endif /* GOGI */


// This shortcuts the path, finds the two-plane wedges, and projects the path
// onto the polytope.  Also removes duplicates
GraphNode *
Path::ShortcutPath (GraphNode * path)
{
  GraphNode *curr = path;
  GraphNode *curr2;
  GraphNode *gn;
  Vec n1, n2;
  Point rp;

  if (curr == NULL)
    return NULL;

  // add in the forward direction pointers for the path list.  also save
  // information about which nodes are the last ones in the path that lie
  // on any given face, edge, or vertex to be used in shortcutting the path
  while (curr->m_prev != NULL)
    {
      curr->m_prev->m_next = curr;
      curr->m_nnloc.m_loc->SetLastInPath (curr);
      curr = curr->m_prev;
    }

  // shortcut the path
  curr2 = curr;
  while (curr->m_next)
    {
      gn = curr->m_nnloc.m_loc->ReturnLastInPath ();
      if ((gn != NULL) && (gn != curr))
        {
          curr->m_next = gn;
          gn->m_prev = curr;
        }

      // use Hersherberger-Suri to find the 2-plane wedge path
      n1 = curr->m_nnloc.m_loc->NormalDir ();
      n2 = curr->m_next->m_nnloc.m_loc->NormalDir ();

      debug ("n1 = (" << n1.x() << "," << n1.y() << "," << n1.z() << ")" );
      debug ("n2 = (" << n2.x() << "," << n2.y() << "," << n2.z() << ")" );
      if (n1 != n2)
        {
          Plane Hs (curr->m_p, -n1);
          Plane Ht (curr->m_next->m_p, -n2);

          rp = RidgePoint (curr->m_p, curr->m_next->m_p, Hs, Ht);

          AddExtraGraphNode (curr, rp);

          printf ("Ridge point: ");
          point_print (rp);
          printf ("\n");

          curr = ProjectPath (curr);
        }

      curr = curr->m_next;
    }

  // remove duplicates from the path and calculate the length of the
  // path
  m_pathLength = 0.0;
  curr = curr2;
  while (curr->m_next)
    {
      while (PointEq (curr->m_p, curr->m_next->m_p))
        curr->m_next = curr->m_next->m_next;

      if (curr->m_next)
        //m_pathLength += sqrt(curr->m_p.
        //                     sqr_dist(curr->m_next->m_p).to_double());
        m_pathLength +=
          SGP_distance( curr->m_p, curr->m_next->m_p );

        // m_pathLength +=
        //  sqrt (CGAL::squared_distance (curr->m_p, curr->m_next->m_p).to_double ());

      curr = curr->m_next;
    }

  return curr2;
}


// This projects a two-plane wedge onto the polytope
// returns a pointer to one node before the destination
GraphNode *
Path::ProjectPath (GraphNode * path, bool b3Wedge)
{
  GraphNode *curr, *dest;
  int i, j;
  Point p1, p2, q, pprev;
  bool contains;
  Edge *e, *e1, *e2;
  Face *f;
  Vec v;

  curr = path;
  dest = (b3Wedge ? curr->m_next->m_next->m_next : curr->m_next->m_next);

  Point v0 (curr->m_p);
  Point v1 (curr->m_next->m_p);
  Point v2 (curr->m_next->m_next->m_p);
  Point v3 (b3Wedge ? curr->m_next->m_next->m_next->m_p : v2);

  debug ("v0 = " << v0);
  debug ("v1 = " << v1);
  debug ("v2 = " << v2);
  debug (CGAL::squared_distance(v0,v1).to_double ());

  // if any of these points are equal then we have no work to do
  if (PointEq (v0, v1) || PointEq (v1, v2) || PointEq (v2, v0))
    return curr->m_next;

  bool setp1 = true;
  Plane pathPlane (v0, CrossProduct (v1 - v0, v2 - v0));
  //bool tempbool=false;

  e1 = e2 = NULL;
  f = NULL;

  // find the two points at which the plane containing the path
  // intersects the source face.  the function Inter returns the
  // intersection of a line and a plane thus the dot products are used
  // to verify that the point of intersection lies on the edge
  // since we don't know whether we are dealing with a face, edge, or
  // vertex, we must check all adjacent faces of the location for the
  // one that is in the direction of the path
  for (j = 0; j < curr->m_nnloc.m_loc->FaceDegree (); j++)
    {
      debug ("j=" << j);
      f = curr->m_nnloc.m_loc->Faces ()[j];

      if (dest->m_nnloc.m_loc->ContainsFace (f))
        {
          curr->m_next = curr->m_next->m_next;
          return curr;
        }
      for (i = 0; i < f->EdgeDegree (); i++)
        {
          e = f->Edges ()[i];

          if (!pathPlane.has_on (e->V1 ()->Location ())
              && (pathPlane.oriented_side (e->V1 ()->Location ()) ==
                  pathPlane.oriented_side (e->V2 ()->Location ())))
            continue;
          q = IntersectLinePlane (e->V1 ()->Location (),
                                  e->V2 ()->Location (), pathPlane, contains);

          v = e->V1 ()->Location () - e->V2 ()->Location ();
          debug ("q = " << q);
          debug ("v = " << v);
          debug ("eV1 = " << e->V1 ()->Location ());
          debug ("eV2 = " << e->V2 ()->Location ());
          if (sign (DotProduct (e->V1 ()->Location () - q, v)) !=
              sign (DotProduct (e->V2 ()->Location () - q, v)))
            {
              // q is on the middle of the edge e...
              printf ("found q in the midle of the edge!\n");
              if (setp1)
                {
                  p1 = Point (q);
                  e1 = e;
                  setp1 = false;
                }
              else if (!PointEq (p1, q))
                {
                  p2 = Point (q);
                  e2 = e;
                  setp1 = true;
                  break;
                }
            }
        }

      if (!setp1)
        continue;

      if (DotProduct (v0 - p2, v0 - v2) > DotProduct (v0 - p1, v0 - v2))
        {
          p1 = p2;
          e1 = e2;
        }

      if (!PointEq (p1, v0))
        break;
    }

  // now it is necessary to determine which of the points of intersection
  // is in the direction of the path.  this is done by checking to see
  // if <s-p1>*<s-t> or <s-p2>*<s-t> > 0.  the correct point will be
  // stored in p1

  debug ("v0 = " << v0);
  debug ("p1 = " << p1);
  debug ("p2 = " << p2);
  debug ("dot = " << DotProduct (v0 - p2, v0 - v2).to_double ());
  debug ("sign(dot) = " << sign (DotProduct (v0 - p2, v0 - v2)));

  // replace the middle point in the path with this one.  other point
  // will be inserted below
  curr = curr->m_next;
  curr->m_p = p1;

  // remove the next node if we're projecting a 3-wedge
  if (b3Wedge)
    curr->m_next = curr->m_next->m_next;

  // continue the intersections until the destination face is reached
  // the algorithm is the same as above except that you already know
  // which way the path is supposed to go
  pprev = v0;
  while (1)
    {
      printf ("Face: %p\n", f);
      // select the next face
      if (e1->F1 () == f)
        f = e1->F2 ();
      else
        f = e1->F1 ();

      // break out if we've reached the destination face
      if (dest->m_nnloc.m_loc->ContainsEdge (e1))
        break;

      printf ("-----------before loop\n");
      // find the next intersection point
      for (i = 0; i < f->EdgeDegree (); i++)
        {
          e = f->Edges ()[i];
          if (e == e1)
            continue;
          q = IntersectLinePlane (e->V1 ()->Location (),
                                  e->V2 ()->Location (), pathPlane, contains);

          // if the edge is contained in the plane then we don't want
          // to use it because it could cause us to go in a circle about
          // a vertex looking for the next face
          if (contains)
            {
              printf ("comtains!\n");
              continue;
            }

          v = e->V1 ()->Location () - e->V2 ()->Location ();
          if (!pathPlane.has_on (e->V1 ()->Location ())
              && (pathPlane.oriented_side (e->V1 ()->Location ()) ==
                  pathPlane.oriented_side (e->V2 ()->Location ())))
            {
              printf ("Does not intersect!\n");
              continue;
            }
          printf ("probably intersect: %d\n", i);
          if ((sign (DotProduct (e->V1 ()->Location () - q, v)) !=
               sign (DotProduct (e->V2 ()->Location () - q, v))) &&
              (sign (DotProduct (v0 - q, v0 - v2)) >= 0))
            //                && !PointEq(q,curr->m_p))
            {
              printf ("intersect!\n");
              p2 = Point (q);
              e1 = e;
              break;
            }
        }

      // if the path is a 3-plane wedge, determine if we've projected enough
      // with the first 3 points.  this is done by determining if the
      // current point is tangent with respect to the path.  the reasoning
      // is discussed in Hershberger and Suri's paper
      if (b3Wedge)
        {
          if (sign (DotProduct (CrossProduct (v2 - p1, p1 - p2),
                                CrossProduct (v2 - p1, p1 - pprev))) > 0)
            {
              pathPlane = Plane (p1, CrossProduct (v2 - p1, v3 - p1));
              b3Wedge = false;
              continue;
            }
          pprev = p1;
          p1 = p2;
        }

      // add the point to the path
      debug (p2);
      if (i != f->EdgeDegree ())
        {
          AddExtraGraphNodeExt (curr, f->project (p2), e1);
          curr = curr->m_next;
        }

      //if (tempbool)
      //     break;
    }

  return curr;
}



Edge *
Path::findExitEdge (NNLocation & nnloc,
                    Plane & pathPlane, Plane & hPositive,
                    history_t & hist, Point & exit_point, GraphNode * dest)
{
  //bool  f_defined;
  RT val;
  Edge *out_edge = NULL;
  Face *f;
  Edge *e;
  Point q;
  bool contains;

  //f_defined = false;
  for (int j = 0; j < nnloc.m_loc->FaceDegree (); j++)
    {
      f = nnloc.m_loc->Faces ()[j];

      // we dont need to do any smart thing...
      if (dest->m_nnloc.m_loc->ContainsFace (f))
        return NULL;
      for (int i = 0; i < f->EdgeDegree (); i++)
        {
          e = f->Edges ()[i];

          if (!pathPlane.has_on (e->V1 ()->Location ())
              && (pathPlane.oriented_side (e->V1 ()->Location ()) ==
                  pathPlane.oriented_side (e->V2 ()->Location ())))
            continue;
          q = IntersectLinePlane (e->V1 ()->Location (),
                                  e->V2 ()->Location (), pathPlane, contains);
          if (hist.isIn (q))
            continue;
          //if  ( hPositive.side_of( q ) < 0 )
          if (hPositive.has_on_negative_side (q))
            continue;

          if (e->isInSlab (q))
            {
              q = e->project_dbl (q);
              exit_point = q;
              hist.register_pt (q);

              out_edge = e;
            }
        }
    }

  return out_edge;
}



// This projects a two-plane wedge onto the polytope
// returns a pointer to one node before the destination
GraphNode *
Path::ProjectPathSP (GraphNode * path, Point & ridge_point, SPath & sp)
{
    GraphNode *curr, *dest;
    Edge *e1, *e2;
    //Face  * f;
    NNLocation nnloc;

    printf ("Path::ProjectPathSP called\n");
    curr = path;
    dest = curr->m_next;

    sp.pushTop (curr);

    Point v0 (curr->m_p);
    Point & v1 (ridge_point);
    Point v2 (curr->m_next->m_p);
    //Point vec_dir( v2 - v0 );
    Point exit_point (v0);

    // if any of these points are equal then we have no work to do
    if (PointEq (v0, v1) || PointEq (v1, v2) || PointEq (v2, v0))
    {
        //printf( "GOGI\n" );
        sp.pushTop (curr);
        return curr->m_next;
    }

    history_t hist;

    hist.init (v0);
    // all the interesting things happen on the side where
    // hPositive is, well, positive
    Plane hPositive (getPositivePlane (v0, v1, v2));

    Plane pathPlane (v0, CrossProduct (v1 - v0, v2 - v0));

    e1 = e2 = NULL;
    //f = NULL;
    nnloc = curr->m_nnloc;
    //printf( "Current point: " );
    //point_print( curr->m_p );
    do
    {
        //printf( "-!  " );
        //nnloc.dump();
        //printf( "-!\n" );

        e1 = findExitEdge (nnloc, pathPlane, hPositive, hist, exit_point, dest);
        if (e1 == NULL)
          break;
        //sp.pushTop( e1->intersect( pathPlane ), e1 );

        nnloc.m_loc = e1;
        nnloc.m_locType = NNLocation::EDGE;

        /* Hoai 07/08/2012 - RIDICULOUS */
        if (pathPlane.has_on (e1->V1 ()->Location ()) && pathPlane.has_on (e1->V1 ()->Location ()))
        {
            assert (false);
        }
        if (pathPlane.has_on (e1->V1 ()->Location ()))
        {
            nnloc.m_loc = e1->V1 ();
            nnloc.m_locType = NNLocation::VERT;
        }
        if (pathPlane.has_on (e1->V2 ()->Location ()))
        {
            nnloc.m_loc = e1->V2 ();
            nnloc.m_locType = NNLocation::VERT;
        }

        //printf( "exit_point: \n" );
        //point_print( exit_point );
        //fflush( stdout );

        sp.pushTop (nnloc, exit_point);
    }
    while (e1 != NULL);

    sp.pushTop (dest);

    return curr->m_next;
}


// adds another node to the path graph.  extra nodes are those that are
// only used to store information about the path (ie, they do not represent
// a nearest neighbor).  so, these are also stored on another list and
// can be reused for subsequent queries
void
Path::AddExtraGraphNode (GraphNode * curr, const Point & p)
{
  GraphNode *gn;

  if (m_currExtraNode == NULL)
    {
      // all info in this GraphNode is bogus (including
      // NNLocation) except for point
      // solely for use in displaying path
      m_extraGraphNodes = new GraphNode (p, curr->m_nnloc, p);
      m_currExtraNode = m_extraGraphNodes;
      // for later use
      m_currExtraNode->m_prev = new GraphNode (p, curr->m_nnloc, p);
    }
  else
    {
      m_currExtraNode->m_p = p;
      if (!m_currExtraNode->m_prev)
        m_currExtraNode->m_prev = new GraphNode (p, curr->m_nnloc, p);
    }

  // note that we cannot traverse the path list backwards anymore because
  // the prev pointer is used for the extra node list
  gn = curr->m_next;
  curr->m_next = m_currExtraNode;
  m_currExtraNode->m_next = gn;
  m_currExtraNode = m_currExtraNode->m_prev;
}


void
Path::AddExtraGraphNodeExt (GraphNode * curr, const Point & p, Edge * e)
{
  GraphNode *gn;
  NNLocation nnloc;

  printf ("addExtraGraphNodeExt  ");
  point_print (p);
  printf ("\n");
  printf ("e->F1()->sidesign, F2: %d, %d\n",
          e->F1 ()->getSideSign (p), e->F2 ()->getSideSign (p));


  nnloc.m_loc = e;
  nnloc.m_locType = NNLocation::EDGE;
  if (m_currExtraNode == NULL)
    {
      // all info in this GraphNode is bogus (including
      // NNLocation) except for point
      // solely for use in displaying path
      m_extraGraphNodes = new GraphNode (p, nnloc, p);
      m_currExtraNode = m_extraGraphNodes;
      // for later use
      m_currExtraNode->m_prev = new GraphNode (p, nnloc, p);
    }
  else
    {
      m_currExtraNode->m_p = p;
      if (!m_currExtraNode->m_prev)
        m_currExtraNode->m_prev = new GraphNode (p, nnloc, p);
    }

  // note that we cannot traverse the path list backwards anymore because
  // the prev pointer is used for the extra node list
  gn = curr->m_next;
  curr->m_next = m_currExtraNode;
  m_currExtraNode->m_next = gn;
  m_currExtraNode = m_currExtraNode->m_prev;
}


double
     Path::PathLength () const
     {
       return m_pathLength;
     }

/*const GraphNode * const * Path::ReturnGraph() const
{
    return m_graph;
}*/

/*int Path::NumNodes() const
{
    return m_numNodes;
}*/

     Plane *Path::MiddlePlane (Edge * e, const Vec & i)
{
  Vec n = CrossProduct ((e->V1 ()->Location () - e->V2 ()->Location ()), i);

  return new Plane (e->V1 ()->Location (), n);
}

RT Path::ThreeWedge (Point * v, Plane * P)
{
  v[0] = m_srcdest->m_srcPoint;
  v[3] = m_srcdest->m_destPoint;

  Point
    s1 = ExtUnfold (m_poly->Center (), *P, m_Hs, v[0]);
  Point
    t1 = ExtUnfold (m_poly->Center (), *P, m_Ht, v[3]);

  Point
    us1 = Unfold (v[0], m_Hs, *P);
  Point
    ut1 = Unfold (v[3], m_Ht, *P);
  /*
     printf( "Source!\n" );
     point_print( s1 );
     printf( "\n" );
     point_print( us1 );
     printf( "\n" );

     printf( "Target!\n" );
     printf( "\n" );
     point_print( t1 );
     printf( "\n" );
     point_print( ut1 );
     printf( "\n" );
   */

  bool
    temp;
  v[1] = IntersectLinePlane (s1, t1, m_Hs, temp);
  v[2] = IntersectLinePlane (s1, t1, m_Ht, temp);

  double
    len = 0.0;
  for (int i = 0; i < 3; i++)
    //len += SGP_distance(v[i], v[i + 1]);
   // Le test them SP SGP
     if(SGP)
    {
        len += SGP_distance(v[i], v[i + 1]);
    }
    else
    {
        len += sqrt (CGAL::squared_distance (v[i], v[i + 1]).to_double ());
    }
  return len;
}





/* calculate the path by LiuWong method       */
/********************************************************/
void Path::liuWongCalcPath ( double *src, double *dest,
                             int srcFaceIdx, int destFaceIdx)
{   // get set of faces of the terrain
    Face * const * faces = m_poly->Faces();
    SourceDest mL_srcdest;
    //mL_srcdest.m_srcPoint = (Point (RT (src[0]), RT (src[1]), RT (src[2])));
    //mL_srcdest.m_destPoint = (Point (RT (dest[0]), RT (dest[1]), RT (dest[2])));
    mL_srcdest.m_srcFace = faces[srcFaceIdx];
    mL_srcdest.m_destFace = faces[destFaceIdx];
    mL_srcdest.m_srcPoint = mL_srcdest.m_srcFace->project (Point (RT (src[0]),
              RT (src[1]),
              RT (src[2])));
    mL_srcdest.m_destPoint = mL_srcdest.m_destFace->project (Point (RT (dest[0]),
                RT (dest[1]),
                RT (dest[2])));

    cout << "   - srcFaceIdx and  destFaceIdx -> " << srcFaceIdx << ";" << destFaceIdx << endl;

    // m_srcdest = new SourceDest ( src, dest, srcFaceIdx, destFaceIdx, m_poly->Faces (), m_poly);
    // avoid error generating by shalow copy Le 5-8-2022
    // SourceDest mL_srcdest( src, dest, srcFaceIdx, destFaceIdx, m_poly->Faces (), m_poly);
    m_srcdest = new SourceDest(mL_srcdest);
    // avoid error generating by shalow copy Le 5-8-2022 --------------------------------

    // for calculate running time
    clock_t start, end;
    double used_time;

    start = clock();

    //##########################################
    // calculate time of LiuWong algorithm at here...
    //##########################################
    end = clock();
    used_time = (double)(end - start) / CLOCKS_PER_SEC; 

    //%%%%%%%%% Write into file 6-8-2022
    ofstream ofs("result.txt", ios::app);
    if(!ofs)
    {
      cerr << "Error: file not opened." << endl;
      exit(1);
    }
  
    ofs << "\n Running time of LiuWong is....: " << setprecision(5) << fixed << used_time;
      
    //Close file
    ofs.close();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete m_srcdest;

}


/* path.cc - end of file */