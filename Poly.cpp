/*****************************************************************************
 * Poly.cc
 *
 * Represents the polytope.  Does the meat of the work in that it constructs
 * the path graph
 *****************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <stdio.h>
#include  <memory.h>
#include  <assert.h>

#include  "sarray.h"
#include  <vector>

#include "Functions.h"
#include "Poly.h"
#include "Face.h"
#include "Vert.h"
#include "Edge.h"
#include "Path.h"
#include <limits.h>
#include  "rgeneric.h"
#include  "Timer.h"

// used to determine if a double is "equal" to zero

Poly::Poly ()
{
}

Poly::~Poly ()
{
}

void
Poly::Read (istream & istrIn)
{
  Poly *p_poly;

  p_poly = new Poly;

  printf ("Reading polytope.\n");
  p_poly->ReadInner (istrIn);
  printf ("\tDone\n");
  fflush (stdout);


  //printf( "before poly dump\n" );
  //fflush( stdout );

  //p_poly->doPolyDump( NULL );
  // fflush( stdout );

  //printf( "Before Convexification\n" );
  //fflush( stdout );
  printf ("Computing convex hull.\n");
  convexify (p_poly->m_verts, p_poly->m_numVerts);
  printf ("\tDone\n");
  fflush (stdout);
  //printf( "After convexification\n" );
  doPolyDump (NULL);
  doPolyTest (NULL);
  //printf( "After second dump\n" );
}



void
Poly::ReadInner (istream & istrIn)
{
  int i, j;
  Point p;
  string sTemp;
  Edge *e;
  int deg, iCurrVert, iPrevVert, iFirstVert;
  bool readDoubles = false;
  double x, y, z;

  // change format for in stream to ascii
  CGAL::set_ascii_mode (istrIn);

  printf ("Entering Poly::Read\n");
  fflush (stdout);

  // read in the "4OFF" on the first line and discard
  istrIn >> ws >> sTemp;

  // if sTemp == "4OFF", then we are reading the coordinates in
  // x,y,z,w format ... otherwise, read doubles
  if (sTemp != "4OFF")
    readDoubles = true;

  // number of edges is initially incorrect (information in file is wrong)
  istrIn >> m_numVerts >> m_numFaces >> m_numEdges;
  m_numEdges = 0;               // make sure to increment this when adding edges
  cout << "Poly has " << m_numVerts << " verts and " << m_numFaces <<
    " faces.";
  fflush (stdout);
  // set up the bounding box with min and max values
  // 0,2,4: min x,y,z values
  // 1,3,5: max x,y,z values
  m_bBox.init ();

  //m_boundingBox[ 0 ] = m_boundingBox[ 2 ] = m_boundingBox[ 4 ] =
  //    rational(integer(INT_MAX));
  //m_boundingBox[1] = m_boundingBox[3] = m_boundingBox[5] =
  //    rational(integer(INT_MIN));

  // create the arrays to hold the verts and faces
  m_verts = new Vert *[m_numVerts];
  m_faces = new Face *[m_numFaces];

  // read in the vertices
  if (readDoubles)
    {
      for (i = 0; i < m_numVerts; i++)
        {
          printf ("i :%d\n", i);
          istrIn >> x >> y >> z;

          // this is an optimization to avoid an exorbitant number of
          // multiplications that slow down the program immensely
          // it assumes that all input only has 6 decimal places
          // and it multiplies by 100
          //    p = Point(int(x*1000000), int(y*1000000),
          //       int(z*1000000), 10000);

          // let's try more precision
          /*
             p = Point((double)(long long)(x*100000000),
             (double)(long long)(y*100000000),
             (double)(long long)(z*100000000), 1000000);
           */
          RT new_x (x * 100.0);
          RT new_y (y * 100.0);
          RT new_z (z * 100.0);

          //cout << "VV(( " << new_x << ", " << new_y << ", " << new_z << " )\n";
          //new_x.normalize();
          //new_y.normalize();
          //new_z.normalize();

          p = Point (new_x, new_y, new_z);
          //cout << "__" << p << "\n";
          Vec v (CGAL::ORIGIN, p);
          Simplify (v);
          p = Point (v.x (), v.y (), v.z ());
          //cout << "_xx_" << v << "\n";
          //assert( v == p.to_vector() );

          cout << "(( " << new_x << ", " << new_y << ", " << new_z << " )\n";

          //p = Point( x *100000000.0, y * 100000000.0, z*100000000.0,
          //          1000000.0 );
          m_verts[i] = new Vert (this, p);
          m_bBox.bound (p);
        }
    }
  else
    {
      for (i = 0; i < m_numVerts; i++)
        {
          int tt;
          istrIn >> p;
          istrIn >> tt;

          cout << p << endl;
          m_verts[i] = new Vert (this, p);
          m_bBox.bound (p);
        }
    }

  // read in the faces, create the edges, and generate adjacency information
  // trang: calculate the angles of vertices of each face of the terrain ...
  double total_face_angles;
  //
  for (i = 0; i < m_numFaces; i++)
  {
      //printf( "faces: %d\n", i );
      // create a face of the correct degree
      istrIn >> deg;

      Face * ptr_face;
      ptr_face = new Face (this, deg);

      m_faces[i] = ptr_face;

      // read in the index of the first vert
      istrIn >> iFirstVert;
      iPrevVert = iFirstVert;

      m_faces[i]->AddVert (m_verts[iPrevVert]);
      m_verts[iPrevVert]->AddFace (m_faces[i]);

      // read in the rest of the verts
      for (j = 0; j < deg; j++)
      {
          // only read in the extras (we already read in the first one)
          if (j != deg - 1)
            {
              istrIn >> iCurrVert;

              m_faces[i]->AddVert (m_verts[iCurrVert]);
              m_verts[iCurrVert]->AddFace (m_faces[i]);
            }

          // check to see if there is already an edge connecting these two
          // verts, and if not, create one
          e = m_verts[iPrevVert]->Connected (m_verts[iCurrVert]);

          // now create it if it hasn't been created
          if (e == NULL)
            {
              e = new Edge (this, m_verts[iPrevVert], m_verts[iCurrVert]);
              m_verts[iPrevVert]->AddEdge (e);
              m_verts[iCurrVert]->AddEdge (e);
              m_edges.push_back (e);
              m_numEdges++;
            }

          m_faces[i]->AddEdge (e);
          e->addFace (m_faces[i]);

          // make the current vert the previous one
          iPrevVert = iCurrVert;

          // if we're on the last one, then make sure we connect and edge
          // to the first one on the next time around the loop
          if (j == deg - 2)
            iCurrVert = iFirstVert;
      }

      // trang: for holding angle values of the face ---------------------------------
      // trang's debug:
      cout << "\t - Poly::Poly() - vertex angles of face " << i << ": " << endl;
      double face_angles = 0;

      total_face_angles = 0;
      //
      for ( j = 0; j < deg ; j++ )
      {
          if ( j == 0 )
              face_angles = angleBetween
                               (
                          m_faces[ i ]->Verts()[j+1]->Location() - m_faces[ i ]->Verts()[j]->Location(),
                          m_faces[ i ]->Verts()[deg - 1]->Location() - m_faces[ i ]->Verts()[j]->Location()
                               );
          else
              if ( j == deg -1 )
              {
                  face_angles = angleBetween
                                   (
                              m_faces[ i ]->Verts()[0]->Location() - m_faces[ i ]->Verts()[j]->Location(),
                              m_faces[ i ]->Verts()[j-1]->Location() - m_faces[ i ]->Verts()[j]->Location()
                                   );
              }
          else
                  face_angles = angleBetween
                               (
                          m_faces[ i ]->Verts()[j+1]->Location() - m_faces[ i ]->Verts()[j]->Location(),
                          m_faces[ i ]->Verts()[j-1]->Location() - m_faces[ i ]->Verts()[j]->Location()
                               );
          //
          m_faces[ i ]->AddVertAngle( face_angles );

          // trang's debug:
          cout << "\t\t - angle of vertex " << j << ": " << m_faces[ i ]->Angles()[ j ] << endl;
          //
          total_face_angles += face_angles;
      }
      // trang's debug:
      cout << "\t\t - total of angles " << total_face_angles << endl;
  }
  // trang's debug:

  Point p1 (m_bBox.minx (), m_bBox.miny (), m_bBox.minz ());
  Point p2 (m_bBox.maxx (), m_bBox.maxy (), m_bBox.maxz ());

  m_r = RT (sqrt (CGAL::squared_distance (p1, p2).to_double ()));
  //m_C = (p1.to_vector() + p2.to_vector())/2;
  m_C = CGAL::midpoint (p1, p2);
}


// checks the directions of the normals of all of the faces
void
Poly::Check ()
{
    for (int i = 0; i < m_numFaces; i++)
      debug ("Direction of normal is: "
             << sign (DotProduct (m_faces[i]->NormalDir (),
                                  Point (0, 0, 0) -
                                  m_faces[i]->Verts ()[0]->Location ())));
}


Face *const * Poly::Faces () const
{
  return m_faces;
}

int Poly::NumFaces () const
{
  return m_numFaces;
}

/*-------------------------
  Hoai - 6/1/2012
*/
Vert * const * Poly::Verts() const
{
  return m_verts;
}

int Poly::NumVerts() const
{
  return m_numVerts;
}
/*-------------------------*/

Point
Poly::Center () const
{
  return m_C;
}



// returns true if curr_pnt is in the face region of f.  also returns an
// array of inward pointing normals of the walls of the region
//
// Translation: Check if the point lies inside the half prism induced
// by the face and its normal. If so, then it is the closest point.
bool
Poly::CheckNNFace (const Face * f, const Point & curr_pnt,
                   vector < Vec > &vInwardPerps)
{
  int signSum = 0;
  int i;

  assert (f != NULL);
  if (f->EdgeDegree () != 3)
    {
      printf ("nearest neighbor failure!\n");
      fflush (stdout);
      assert (f->EdgeDegree () == 3);
    }

  for (i = 0; i < 3; i++)
    {
      // these vectors should be in their respective places since
      // they are read in the order that they are used here
      //vInwardPerps[i] = CrossProductDir( f->NormalDir(),
      vInwardPerps[i] = CrossProduct (f->NormalDir (),
                                      f->Verts ()[(i + 1) % 3]->Location () -
                                      f->Verts ()[i]->Location ());

      signSum += sign (DotProduct (vInwardPerps[i], curr_pnt -
                                   f->Verts ()[i]->Location ()));
    }

  // the first clause makes sure that the current point and the face
  // are on the same side of the polytope and the second clause
  // makes sure that the current point is in the region of the face
  if ((sign (DotProduct (f->NormalDir (),
                         Vec (CGAL::ORIGIN, curr_pnt))) > 0) &&
      (abs (signSum) == f->EdgeDegree ()))
    return true;
  else
    return false;
}


// returns true if curr_pnt is in the region for this edge.
// the four vectors that CheckNNEdge returns are the inward pointing normals
// of the faces of the wedge created by the edge
bool
  Poly::CheckNNEdge (const Edge * e, const Point & curr_pnt, Vec & vN1,
                     Vec & vN2, Vec & vV1, Vec & vV2)
{
  bool bReturnVal = true;
  Vec vP (curr_pnt - e->V1 ()->Location ());

  // check to see if the current point is within the "skinny" faces
  // of the wedge created by the edge
  vV1 = (e->V2 ()->Location () - e->V1 ()->Location ());
  vV2 = -vV1;

  if (sign (DotProduct (vV1, vP)) <= 0)
    bReturnVal = false;
  if (sign (DotProduct (vV2, curr_pnt - e->V2 ()->Location ())) <= 0)
    bReturnVal = false;

  // check to see if the current point is within the other two faces
  //vN1 = CrossProductDir( e->F1()->NormalDir(), vV1, 1999 );
  //vN2 = CrossProductDir( e->F2()->NormalDir(), vV2, 2000 );
  vN1 = CrossProduct (e->F1 ()->NormalDir (), vV1);
  vN2 = CrossProduct (e->F2 ()->NormalDir (), vV2);

  int s1 = sign (DotProduct (vN1, vP));
  int s2 = sign (DotProduct (vN2, vP));


  //printf( "real s1: %d    s2: %d\n", s1, s2 ); fflush(stdout );

  if (s1 == 0)
    s2 = s1;
  if (s2 == 0)
    s1 = s2;

  if (s1 != s2)
    {
      //printf( "s1: %d    s2: %d\n", s1, s2 ); fflush(stdout );
      bReturnVal = false;

      // make sure that the normals of the faces of the wedge are pointing
      // inwards
      if (sign (DotProduct (e->F1 ()->NormalDir (), vN2)) < 0)
        {
          vN1 = -vN1;
          vN2 = -vN2;
        }
    }
  else
    {
      // make sure that the normals of the faces of the wedge are pointing
      // inwards
      if (s1 < 0)
        {
          vN1 = -vN1;
          vN2 = -vN2;
        }
    }

  return bReturnVal;
}


// returns true if curr_pnt is in the region for this vertex
// returns array of vectors which are the normals of the face of the "cone"
// generated by the vertex
bool
  Poly::CheckNNVertex (const Vert * v, const Point & curr_pnt,
                       vector < Vec > &vV1)
{
  int i;
  const Edge *e;
  Vec vP (curr_pnt - v->Location ());
  bool bReturnVal = true;
  /*
     printf( "CheckNNVertex called( " );
     point_print( v->Location() );
     printf( "): %g", v->Location().sqr_dist( curr_pnt ).to_double() );
   */
  for (i = 0, e = v->Edges ()[i]; i < v->EdgeDegree (); i++)
    {
      e = v->Edges ()[i];
      //printf( " i:: %d, %p\n", i, e );
      vV1[i] = (v->Location () - (v == e->V1 ()? e->V2 ()->Location () :
                                  e->V1 ()->Location ()));

      RT distx;
      //        NNLocation  nnTmp;
      //printf( "UU " );
      //e->getSqrDistance( curr_pnt, distx, nnTmp );

      if (sign (DotProduct (vV1[i], vP)) < 0)
        {
          /*
             printf( "Point is not in slab of edge : %p\n", e );
             printf( "                          F1 : %p %d\n", e->F1(),
             e->F1()->getSideSign( curr_pnt ) );
             printf( "                          F2 : %p %d\n", e->F2(),
             e->F2()->getSideSign( curr_pnt ) );
           */
          bReturnVal = false;
        }
    }

  //printf( "before return!\n" );
  return bReturnVal;
}


// checks all of the faces, edges, and vertices to find the nearest neighbor
// of curr_pnt
void
Poly::NNCheckAll (const Point & curr_pnt, NNLocation & nnloc, Point & ptNN)
{
  int i;
  Face *f;
  Edge *e, *closestEdge = 0;
  RT q, closestQ;
  RT d, closestDist = -1;
  vector < Vec > vInwardPerps (3);

  nnloc.m_loc = 0;

  printf ("Poly::NNCheckAll!\n");
  fflush (stdout);

  debug ("Entered Poly::NNCheckAll");
  debug ("Checking faces ...");

  // check all of the faces for a nearest neighbor for the current point
  for (i = 0; i < m_numFaces; i++)
    {
      f = m_faces[i];
      assert (m_faces[i] != NULL);
      if (CheckNNFace (f, curr_pnt, vInwardPerps))
        {
          printf ("aaa\n");
          debug ("NNLocation is FACE");
          nnloc.m_loc = f;
          nnloc.m_locType = NNLocation::FACE;


          //q = DotProduct(curr_pnt - f->Verts()[0]->Location(),
          //         f->Normal());
          //Point(curr_pnt - q*f->Normal());
          ptNN = f->project (curr_pnt);
        }
      //printf( "bbbbaaa\n" );
    }

  if (nnloc.m_loc != NULL)
    return;

  // if we haven't found the nearest neighbor yet, check all of the edges
  // and vertices
  debug ("Checking edges and vertices ...");

  // calculate the shortest distance from curr_pnt to each edge and use the
  // smallest one as the nearest neighbor
  //forall (e, m_edges) {
  for (list < Edge * >::iterator ei = m_edges.begin (); ei != m_edges.end ();
       ei++)
    {
      e = *ei;
      d = e->SqrDistance (curr_pnt);

      if ((closestDist == -1) || (d < closestDist))
        {
          closestDist = d;
          closestQ = q;
          closestEdge = e;
        }
    }

  q = closestQ;
  e = closestEdge;

  RT _dist;

  e->getSqrDistance (curr_pnt, _dist, nnloc, -1);

  //nnloc.m_loc = e;
  //nnloc.m_locType = NNLocation::EDGE;

  debug ("NNLocation is EDGE");

  ptNN = e->project (curr_pnt);
}


struct walk_nn_t
{
  //walk_nn_t() {
  //    memset( this, 0, sizeof( walk_nn_t ) );
  //}

  NNLocation nnloc;
  Point curr_pnt;
  GridGraph *pGrid;
  bool bNNFound;
  Point nnpt, oldCurr_pnt;
  int phi_ind, theta_ind;       //, phi_bound, theta_bound;
  double theta;
  Path *pthP;
  bool bFirst;
  void *oldNNLoc;
  RT ratFoundDot;
  bool bFound;
  int unique_walk_id;

  int angle_ind () const
  {
    return phi_ind + theta_ind * pGrid->getPhiBound (); //phi_bound;
  }

  bool fIsNNDefined;
  RT faceNNSqrDist;
  Face *pNNFace;
  NNLocation face_walk_nnloc;

  Face *getNNFace () const
  {
    return pNNFace;
  }

  void resetCandidNNValue ()
  {
    fIsNNDefined = false;
    faceNNSqrDist = 0;
    pNNFace = NULL;
    //memset (&face_walk_nnloc, 0, sizeof (NNLocation));
  }

  bool registerFaceValue (Face * f, RT & candid_dist_sq,
                          NNLocation & test_nnloc)
  {
    if ((!fIsNNDefined) || (candid_dist_sq < faceNNSqrDist))
      {
        pNNFace = f;
        faceNNSqrDist = candid_dist_sq;
        fIsNNDefined = true;
        face_walk_nnloc = test_nnloc;

        return true;
      }

    return false;
  }

  void setFaceWalkID (int id)
  {
    unique_walk_id = id;
  }


  bool tryFace (Face * pBaseFace, Face * pNewFace)
  {
    if (pBaseFace == pNewFace)
      return false;

    if ((pNewFace->getLastFaceWalkID () > 0)
        && (unique_walk_id == pNewFace->getLastFaceWalkID ()))
      {
        //printf( "breaking: %d!\n", unique_walk_id );
        return false;
      }

    pNewFace->setLastFaceWalkID (unique_walk_id);

    //printf( "          tryFace: %p\n", pNewFace );
    RT res;

    NNLocation test_nnloc;

    if ((pNewFace->getSideSign (curr_pnt) > 0)
        && (pNewFace->getNearestSqrDistance (curr_pnt, res,
                                             test_nnloc, unique_walk_id)))
      {
        //printf( "       Value: %g\n", res.to_double() );
        if (registerFaceValue (pNewFace, res, test_nnloc))
          return true;
      }


    return false;
  }
};

void
point_print (const Point & pnt)
{
  printf ("(%g, %g, %g)",
          pnt.x ().to_double (),
          pnt.y ().to_double (),
          pnt.z ().to_double ());
}

void err_point_print( const Point &pnt )
{
  fprintf( stderr, "(%g, %g, %g)",
           pnt.x ().to_double (),
           pnt.y ().to_double (),
           pnt.z ().to_double ());
}


void
vec_print (const Vec & v)
{
  Point p (v.x (), v.y (), v.z ());
  point_print (p);
}


void
point_print_input_frmt (const Point & pnt)
{
  printf ("%f %f %f",
          pnt.x ().to_double (),
          pnt.y ().to_double (), pnt.z ().to_double ());
}

void
Poly::doPolyDump (const Point * p_pnt)
{
  //static int  count = 1;

  //if  ( count <= 0 )
  //    return;
  //count--;

  printf ("\n\n-----------------------------------------\n");

  for (int i = 0; i < m_numVerts; i++)
    {
      printf ("\t%3d: ", i);
      point_print (m_verts[i]->Location ());
      if (p_pnt != NULL)
        printf ("  %g  ",
                CGAL::squared_distance (m_verts[i]->Location (),
                                        *p_pnt).to_double ());
      printf ("\n");
    }

  if (p_pnt != NULL)
    {
      for (int i = 0; i < m_numFaces; i++)
        {
          printf ("FACE\t%3d: ", i);

          if (p_pnt != NULL)
            printf ("[pnt:%3d]   ", m_faces[i]->getSideSign (*p_pnt));
          printf (" [m_C:%3d]\n", m_faces[i]->getSideSign (m_C));
          //m_faces[ i ]->checkIfConvex();
        }
    }
  printf ("-----------------------------------------\n\n");
}


void
Poly::doPolyTest (const Point * p_pnt)
{
  printf ("Do polytope test!\n");

  (void) p_pnt;
  {
    for (int i = 0; i < m_numFaces; i++)
      {
        assert (m_faces[i]->EdgeDegree () == m_faces[i]->VertDegree ());
        assert (m_faces[i]->EdgeDegree () == 3);

        m_faces[i]->checkIfConvex ();
      }
  }

  for (int i = 0; i < m_numVerts; i++)
    {
      Vert & v (*(m_verts[i]));

      for (int jnd = 0; jnd < v.FaceDegree (); jnd++)
        {
          Face & f (*(v.Faces ()[jnd]));
          assert (f.isVertexAdj (&v));

          int res = f.getSideSign (v.Location ());
          if (res != 0)
            {
              printf ("jnd = %d\n", jnd);
              assert (res == 0);
            }
        }
    }

}


int
getUniqueID ()
{
  static int id = 100;

  return id++;
}


int
     Poly::getFaceIndex (const Face * face) const
     {
       assert (face != NULL);

       for (int ind = 0; ind < NumFaces (); ind++)
         if (m_faces[ind] == face)
           return ind;

       assert (false);

       return -1;
     }


// Starte\ fro ma visible point to info.curr_pnt and finds the face
// which contains the closest point - by a brute force algorithm that
// simply check all the adjacent faces that are visible. It is easy to
// provew that this is correct. I hope.
void Poly::doFaceWalk (Face * f, walk_nn_t & info)
{
  bool f_try_again = true;
  RT res;
  //    Edge  * e;

  //printf( "doFaceWalk!\n" );
  //doPolyTest( &(info.curr_pnt) );

  //printf( "doFaceWalk( %p, ? )\n", f );

  info.resetCandidNNValue ();
  while (f_try_again)
    {
      f_try_again = false;

      if (info.tryFace (NULL, f))
        f_try_again = true;

      int i;
      //Edge * e;
      /*
         for  ( i = 0, e = f->Edges()[i]; i < f->EdgeDegree();
         e=f->Edges()[++i]) {
         if  ( info.tryFace( f, e->F1() ) )
         f_try_again = true;
         if  ( info.tryFace( f, e->F2() ) )
         f_try_again = true;
         }
       */

      // Now we should try all the faces that are adjacent to vertices
      // of f
      Vert *v;

      assert (f->EdgeDegree () == 3);
      // for (i = 0, v = f->Verts ()[i]; i < 3; v = f->Verts ()[++i])
      //   {
      //     //printf( "___ i : %d ", i );
      //     for (int jnd = 0; jnd < v->FaceDegree (); jnd++)
      //       if (info.tryFace (f, v->Faces ()[jnd]))
      //         f_try_again = true;
      //   }
      for ( i = 0; i < 3; ++i )
        {
          v = f->Verts()[i];
          //printf( "___ i : %d ", i );
          for (int jnd = 0; jnd < v->FaceDegree (); jnd++)
            if (info.tryFace (f, v->Faces ()[jnd]))
              f_try_again = true;
        }

      if (f_try_again)
        f = info.getNNFace ();
      //printf( "\tinfo.faceNNSqrDist: %g\n",
      //       info.faceNNSqrDist.to_double() );
    }
  info.nnloc = info.face_walk_nnloc;
  assert (info.face_walk_nnloc.m_loc != NULL);
  //printf( "nnloc.m_locType: %d\n", info.nnloc.m_locType );
  //printf( "nnloc.m_loc    : %p\n", info.nnloc.m_loc );
  //printf( "nnloc.DISTANCE: %g\n",
  //        info.faceNNSqrDist.to_double() );

  //if  ( info.nnloc.m_locType == NNLocation::VERT ) {
  //Vert  * v = (Vert *)info.nnloc.m_loc;
  //printf( "AAA " );
  //fflush( stdout );
  ///        point_print( v->Location() );
  //printf( ": %g\n", v->Location().sqr_dist( info.curr_pnt ).to_double() );
  //fflush( stdout );
  //}

  //info.nnloc.m_loc = f;
  //info.nnloc.m_locType = NNLocation::FACE;
}


void
Poly::WalkNNFace (walk_nn_t & info, list < AnimateVecInfo * >&debugPoints)
{
  GraphNode *gn;
  Face *f = static_cast < Face * >(info.nnloc.m_loc);
  vector < Vec > vInwardPerps (3);

  //printf( "walknnface( %p )\n", (void *)f );

  debug ("f=0x" << hex << (int) f);

  assert (f != NULL);
  if (CheckNNFace (f, info.curr_pnt, vInwardPerps))
    {
      // the current point is still closest to the face
      // or we just entered the face "region"
      // and the current point is in here (ie the nearest
      // neighbor lies on this face
      info.bNNFound = true;

      //printf( "face is nearest neighbor!\n" );

      debug ("NNloc = FACE");

      // calculate where the nearest neighbor is using
      // a normal projection
      //RT q(DotProduct(info.curr_pnt -
      //                              f->Verts()[0]->Location(),
      //                      f->Normal()));
      info.nnpt = f->project (info.curr_pnt);

      info.nnloc.m_loc = f;
      info.nnloc.m_locType = NNLocation::FACE;
      //printf( "before creating graph node !\n " );
      //printf( "Point sign: %p: %d\n",
      //      f, f->getSideSign( info.nnpt ) );

      // add the nearest neighbor to the graph
      gn = CreateGraphNode (*info.pGrid, info.nnpt, info.nnloc,
                            info.phi_ind, info.theta_ind, info.curr_pnt);
      info.pthP->AddNodeToGraph (gn, info.angle_ind ());

      info.nnloc.m_loc->AddNN (info.angle_ind ());
      debugPoints.push_back (new AnimateVecInfo (info.curr_pnt, info.nnpt,
                                                 AnimateVecInfo::NN));
      return;
    }
  if (f->isFaceOrNeigborFaceVisible (info.curr_pnt))
    {
      doFaceWalk (f, info);
      return;
    }

  NNCheckAll (info.curr_pnt, info.nnloc, info.nnpt);
}



void
Poly::WalkNNEdge (walk_nn_t & info, list < AnimateVecInfo * >&debugPoints)
{
  //printf( "walknnEdge\n" );
  GraphNode *gn;

  Edge *e = static_cast < Edge * >(info.nnloc.m_loc);

  // the following vectors are the normals of the
  // side faces of the wedge
  Vec vN01, vN23, vN02, vN13;

  debug ("e=0x" << hex << (int) e);

  assert (e != NULL);
  if (CheckNNEdge (e, info.curr_pnt, vN01, vN23, vN02, vN13))
    {
      //printf( "nearest neighbor found!\n" );
      // the current point is still closest to the edge
      // or we just entered the edge region and the
      // current point is in here
      info.bNNFound = true;

      debug ("NNLoc = EDGE");

      info.nnloc.m_loc = e;
      info.nnloc.m_locType = NNLocation::EDGE;

      // calculate where the nearest neighbor is
      // RT q(DotProduct( info.curr_pnt - e->V1()->Location(),
      //          Unit(e->V2()->Location() -
      //         e->V1()->Location())));
      info.nnpt = e->project (info.curr_pnt);
      //printf( "Point projected on edge " );
      //point_print( info.nnpt );
      //printf( "\n" );

      // add the nearest neighbor to the path graph
      gn = CreateGraphNode (*(info.pGrid),
                            info.nnpt, info.nnloc,
                            info.phi_ind, info.theta_ind, info.curr_pnt);
      info.pthP->AddNodeToGraph (gn, info.angle_ind ());

      info.nnloc.m_loc->AddNN (info.angle_ind ());
      debugPoints.push_back (new AnimateVecInfo (info.curr_pnt,
                                                 info.nnpt,
                                                 AnimateVecInfo::NN));
      return;
    }

  if (e->F1 ()->isFaceOrNeigborFaceVisible (info.curr_pnt))
    {
      doFaceWalk (e->F1 (), info);
      return;
    }
  if (e->F2 ()->isFaceOrNeigborFaceVisible (info.curr_pnt))
    {
      doFaceWalk (e->F1 (), info);
      return;
    }

  NNCheckAll (info.curr_pnt, info.nnloc, info.nnpt);
}


void
Poly::WalkNNVertex (walk_nn_t & info, list < AnimateVecInfo * >&debugPoints)
{
  GraphNode *gn;
  Vert *v = static_cast < Vert * >(info.nnloc.m_loc);
  vector < Vec > vInwardPerps (v->EdgeDegree ());

  //printf( "walknnVertex( ?, %p)\n", v );

  debug ("v=0x" << hex << (int) v);
  //printf( "before if: %p \n", v );

  assert (v != NULL);
  if (CheckNNVertex (v, info.curr_pnt, vInwardPerps))
    {
      // the current point is still closest to the
      // vertex or we just entered the vertex region
      // and the current point is in here
      info.bNNFound = true;

      debug ("NNLoc = VERT");

      info.nnloc.m_loc = v;
      info.nnloc.m_locType = NNLocation::VERT;

      //printf(" 12345: %p\n", v );
      info.nnpt = v->Location ();

      //printf( "(%d, %d, %d %d)\n",
      //        info.phi_ind, info.theta_ind,
      //        info.pGrid->getPhiBound(), info.pGrid->getThetaBound() );
      gn = CreateGraphNode (*(info.pGrid),
                            info.nnpt, info.nnloc,
                            info.phi_ind, info.theta_ind, info.curr_pnt);
      info.pthP->AddNodeToGraph (gn, info.angle_ind ());
      //printf( "xxx 12345\n" );

      info.nnloc.m_loc->AddNN (info.angle_ind ());
      debugPoints.push_back (new AnimateVecInfo (info.curr_pnt, info.nnpt,
                                                 AnimateVecInfo::NN));
      //printf( "rturning!\n" );
      return;
    }

  int i;
  Edge *e;

  for (i = 0; i < v->EdgeDegree (); i++)
    {
      e = v->Edges ()[i];
      if (e->F1 ()->isFaceOrNeigborFaceVisible (info.curr_pnt))
        {
          doFaceWalk (e->F1 (), info);
          return;
        }
      if (e->F2 ()->isFaceOrNeigborFaceVisible (info.curr_pnt))
        {
          doFaceWalk (e->F1 (), info);
          return;
        }
    }

  NNCheckAll (info.curr_pnt, info.nnloc, info.nnpt);
}


void
Poly::handlePoint (walk_nn_t & info,
                   double phi, double theta,
                   list < AnimateVecInfo * >&debugPoints)
{
    // trang: avoiding unused variable
    theta = theta;

    //Vec vGC;

    // verify that we do not visit each face too much.
    info.setFaceWalkID (getUniqueID ());

    //printf( "\tphi = %g\n", phi );
    debug ("theta = " << info.theta);
    debug ("phi = " << phi);

    //info.curr_pnt = Point( m_C + m_r * Point(
    //                     RT( sin(phi)*cos(info.theta) ),
    //                     RT(sin(phi)*sin(info.theta)),
    //                     RT(cos(phi))).to_vector());
    info.curr_pnt = Point (m_C + m_r * Vec (RT (sin (phi) * cos (info.theta)),
                                            RT (sin (phi) * sin (info.theta)),
                                            RT (cos (phi))));

    debug ("curr_pnt = " << info.curr_pnt);

    // add points to display the walking of the longitude
    debugPoints.push_back (new AnimateVecInfo (info.oldCurr_pnt,
                                               info.curr_pnt,
                                               AnimateVecInfo::LONGITUDE));

    //vGC = CrossProduct( info.oldCurr_pnt.to_vector(),
    // info.curr_pnt.to_vector() );
    info.bNNFound = false;

    // keep walking through regions until we find the nearest neighbor
    while (!info.bNNFound)
      {
        info.bFound = false;
        info.ratFoundDot = 0;

        printf ("iteration...\n");
        fflush (stdout);
        // check for the nearest neighbor and where to move to next
        // based on what type of region we are currently in
        switch (info.nnloc.m_locType)
          {
          case NNLocation::FACE:
            WalkNNFace (info, debugPoints);
            break;

          case NNLocation::EDGE:
            WalkNNEdge (info, debugPoints);
            break;

          case NNLocation::VERT:
            WalkNNVertex (info, debugPoints);
            break;

          default:
            assert (false);
          }                       /* switch */
      }                           /* while */
    info.oldCurr_pnt = info.curr_pnt;
}



// the meat of the program.  construct the path graph by walking the
// longitudes and finding the nearest neighbors.
// phi_bound and theta_bound are bounds for the indices of the nodes in the
// lat/lon grid (ie, they are the number of latitude and longitude lines
// respectively
void
Poly::ConstructPathGraph (Path * pthP, const double epsilon,
                          GridGraph & _grid,
                          //const int _phi_bound, const int _theta_bound,
                          list < AnimateVecInfo * >&debugPoints)
{
  walk_nn_t info;
  Timer timer;

  timer.start ();

  //Point p (m_C + m_r * Vec (0, 0, 1, 1));
  info.curr_pnt = Point (m_C + m_r * Vec (0, 0, 1, 1));

  //fflush( stdout );
  printf ("Poly::ConstructPathGraph\n");
  fflush (stdout);

  Point northPole = info.curr_pnt;
  info.oldCurr_pnt = info.curr_pnt;

  info.pthP = pthP;
  info.pGrid = &_grid;
  GraphNode *gn;


  // check all the faces, edges, and vertices to find the nearest neighbor
  // for the north pole
  NNCheckAll (info.curr_pnt, info.nnloc, info.nnpt);

  // make a graph node out of the point and add it to the graph
  gn = CreateGraphNode (*(info.pGrid),
                        info.nnpt, info.nnloc, -1, -1, info.curr_pnt);
  //    printf( "B\n" ); fflush( stdout );
  info.pGrid->addNorthPole (gn);

  //pthP->AddNodeToGraph( gn, info.pGrid->getPhiBound()
  //               * info.pGrid->getThetaBound() + 1 );
  //printf( "C\n" ); fflush( stdout );

  info.nnloc.m_loc->AddNN (info.pGrid->getNorthPoleIndex ());
  //printf( "D\n" ); fflush( stdout );
  debugPoints.push_back (new AnimateVecInfo (info.curr_pnt, info.nnpt,
                                             AnimateVecInfo::NN));
  //printf( "E\n" ); fflush( stdout );

  double phi;
  NNLocation northPoleNNLoc = info.nnloc;

  info.bFirst = true;
  info.phi_ind = info.theta_ind = 0;

  // walk the longitude lines and determine the nearest neighbors of all
  // of the lat/lon intersection points

  debug ("Walking Longitudes");

  for (info.theta = 0.0;
       info.theta < 2 * CGAL_PI; info.theta += epsilon, info.theta_ind++)
    {
      info.nnloc = northPoleNNLoc;
      info.oldCurr_pnt = northPole;
      info.oldNNLoc = 0;
      info.bFirst = true;
      info.phi_ind = 0;

      //printf( "theta = %g\n", info.theta );
      // walk down the longitude line and if it's the last one before we
      // would wrap back around, take one extra step and calculate the
      // nearest neighbor for the south pole
      for (phi = epsilon; (phi < CGAL_PI); phi += epsilon, info.phi_ind++)
        {
          printf ("theta: %10g, phi: %10g\n", info.theta, phi);
          fflush (stdout);
          handlePoint (info, phi, info.theta, debugPoints);
        }
    }

  //printf( "\n\n\nSOUTH POLE!\n" );
  // south pole
  //fflush( stdout );
  info.phi_ind = 0;
  //info.theta_ind = 0;
  handlePoint (info, CGAL_PI, 0.0, debugPoints);

  timer.end ();

  printf ("Path Graph construction completed!\n");
  printf ("Path graph construction time: %g\n", timer.seconds ());
  fflush (stdout);
}


// return the index of the new node in the grid-graph
int
Poly::addPoint (GridGraph & grid,
                const Point & pnt,
                const Point & pnt_on_sphere, Face * pFace, double epsilon)
{
  double phi, theta;
  int phi_ind, theta_ind, ind_new_vrtx;
  GraphNode *gn;
  NNLocation nnloc;

  phi = CalcPhi (pnt_on_sphere, m_r);
  theta = CalcTheta (pnt_on_sphere, m_r);

  theta_ind = int (theta / epsilon);
  phi_ind = int (phi / epsilon) - 1;

  nnloc.m_loc = pFace;
  nnloc.m_locType = NNLocation::FACE;

  gn = new GraphNode (pnt, nnloc, pnt_on_sphere);
  ind_new_vrtx = grid.AppendVertex (gn);

  //grid.connect( phi_ind, theta_ind, ind_new_vrtx );
  grid.connect (phi_ind, theta_ind, ind_new_vrtx);

  for (int di = -1; di < 2; di++)
    for (int dj = -1; dj < 2; dj++)
      grid.connect (phi_ind + di, theta_ind + dj, ind_new_vrtx);
//     grid.connect( phi_ind + 1, theta_ind    , ind_new_vrtx );
//     grid.connect( phi_ind + 1, theta_ind + 1, ind_new_vrtx );
//     grid.connect( phi_ind    , theta_ind    , ind_new_vrtx );
//     grid.connect( phi_ind    , theta_ind + 1, ind_new_vrtx );
//     grid.connect( phi_ind + 1, theta_ind    , ind_new_vrtx );
//     grid.connect( phi_ind + 1, theta_ind + 1, ind_new_vrtx );

  //set<int>::iterator sj;

  for (int lnd = 0; lnd < pFace->ReturnNNs ().size (); lnd++)
    {
      grid.connectByIndices (ind_new_vrtx, pFace->ReturnNNs ()[lnd]);
    }
  /*    for (sj = pFace->ReturnNNs().begin();
     sj != pFace->ReturnNNs().end(); sj++)
     grid.connectByIndices( ind_new_vrtx, *sj );
   */
  return ind_new_vrtx;
}


// adds the source and destination points to the path graph and adds the
// appropriate edges too
void
Poly::AddSourceDest (Path * pthP, SourceDest * srcdest,
                     const double epsilon,
                     GridGraph & grid, const int pickNum)
{
    cout << " \t start adding SourceDest... " << endl;
    // trang: avoiding unused variable
    pthP->used = true; int using_pickNum = 0; using_pickNum = using_pickNum + pickNum;

    // original one is bellowed
    int source_ind, dest_ind;

    Vec sv = srcdest->m_srcPoint - m_C; //m_C is a center?
    Vec dv = srcdest->m_destPoint - m_C;
    Point src = ProjectPointOntoSphere (Point (sv.x (), sv.y (), sv.z ()),
                                        srcdest->m_srcFace->NormalDir (), m_r);
    Point dst = ProjectPointOntoSphere (Point (dv.x (), dv.y (), dv.z ()),
                                        srcdest->m_destFace->NormalDir (), m_r);

    source_ind = addPoint (grid, srcdest->m_srcPoint,
                           src, srcdest->m_srcFace, epsilon);
    grid.setSourceVertex (source_ind);

    dest_ind = addPoint (grid, srcdest->m_destPoint,
                         dst, srcdest->m_destFace, epsilon);
    grid.setDestVertex (dest_ind);
    cout << " \t after adding SourceDest... " << endl;
}

/************* LE: do not need the first argument Path * pthP ************/
void
Poly::AddSourceDest (SourceDest * srcdest,
                     const double epsilon,
                     GridGraph & grid, const int pickNum)
{
    cout << " \t start adding SourceDest (Path is unnecessary)... " << endl;
    // trang: avoiding unused variable
    int using_pickNum = 0; using_pickNum = using_pickNum + pickNum;

    // original one is bellowed
    int source_ind, dest_ind;

    Vec sv = srcdest->m_srcPoint - m_C; //m_C is a center?
    Vec dv = srcdest->m_destPoint - m_C;
    Point src = ProjectPointOntoSphere (Point (sv.x (), sv.y (), sv.z ()),
                                        srcdest->m_srcFace->NormalDir (), m_r);
    Point dst = ProjectPointOntoSphere (Point (dv.x (), dv.y (), dv.z ()),
                                        srcdest->m_destFace->NormalDir (), m_r);

    source_ind = addPoint (grid, srcdest->m_srcPoint,
                           src, srcdest->m_srcFace, epsilon);
    grid.setSourceVertex (source_ind);

    dest_ind = addPoint (grid, srcdest->m_destPoint,
                         dst, srcdest->m_destFace, epsilon);
    grid.setDestVertex (dest_ind);
}
/************* LE: do not need the first argument Path * pthP -end ************/

     list < Edge * >Poly::HorizonEdges (const Vec & u) const
     {
       Edge *
         e;
       list <
       Edge * >
         hEdges;
       //forall (e, m_edges)
       for (list < Edge * >::const_iterator ei = m_edges.begin ();
            ei != m_edges.end (); ei++)
         {
           e = *ei;
           if (hEdge (e, u))
             hEdges.
             push_back (e);
         }

       return
         hEdges;
     }

     bool Poly::hEdge (const Edge * e, const Vec & u) const
     {
       int
         sg1,
         sg2;

       Vec
         v1 = e->F1 ()->NormalDir ();
       Vec
         v2 = e->F2 ()->NormalDir ();

       sg1 = sign (DotProduct (u, v1));
       sg2 = sign (DotProduct (u, v2));

       return ((sg1 != 0) && (sg1 == -sg2));
     }

Point
     Poly::computeRealCenter () const
     {
       Point p;
       //leda_rat_vector  v( m_verts[ 0 ]->Location().to_vector() );
       Vec v (CGAL::ORIGIN, m_verts[0]->Location ());

       for (int ind = 1; ind < m_numVerts; ind++)
           v = v + Vec (CGAL::ORIGIN, m_verts[ind]->Location ());

         return Point (v.x () / (double) m_numVerts,
                       v.y () / (double) m_numVerts,
                       v.z () / (double) m_numVerts);
     }



     void Poly::sortAdjInfo ()
{
  for (int ind = 1; ind < m_numVerts; ind++)
    m_verts[ind]->sortAdjInfo ();
}


/* Poly.C - End of File ------------------------------------------*/
