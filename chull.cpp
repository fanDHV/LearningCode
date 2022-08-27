#include "myCGAL.h"
#include "sarray.h"
#include "Vert.h"
#include "Edge.h"
#include "Face.h"
#include "Poly.h"

typedef
  KPoly::Facet_iterator
  Facet_iterator;
typedef
  KPoly::Halfedge_around_facet_circulator
  Halfedge_facet_circulator;

void
Poly::convexify (Vert ** arr, int num)
{
  //int i;
  int faces_count;
  Face *ptr_face;
  Edge *ptr_edge;

  m_bBox.init ();

  // save point list for computing convex hull
  vector < Point > ptlist;
  for (int i = 0; i < num; i++)
    ptlist.push_back (arr[i]->Location ());

  // compute 3D convex hull
  // assume that the convex hull is triangulated already
  KPoly poly;
  CGAL::convex_hull_3 (ptlist.begin (), ptlist.end (), poly);

  // let start the real work ...
  m_numVerts = poly.size_of_vertices ();
  m_numFaces = poly.size_of_facets ();

  std::cout << "Facets: " << m_numFaces << endl;

  m_numEdges = 0;
  m_verts = new Vert *[m_numVerts];
  m_faces = new Face *[m_numFaces];
  memset (m_verts, 0, sizeof (Vert *) * m_numVerts);
  memset (m_faces, 0, sizeof (Face *) * m_numFaces);

  // let's create the vertices
  //i = 0;
  for (KPoly::Vertex_iterator vi = poly.vertices_begin ();
       vi != poly.vertices_end (); vi++)
    {
      //m_verts[ i++ ] = new Vert( this, vi->point() );
      m_verts[std::distance (poly.vertices_begin (), vi)] =
        new Vert (this, vi->point ());
      m_bBox.bound (vi->point ());
    }

  // iterate over all facets
  faces_count = 0;
  for (Facet_iterator fi = poly.facets_begin (); fi != poly.facets_end ();
       ++fi)
    {
      Halfedge_facet_circulator j = fi->facet_begin ();
      // Facets in polyhedral surfaces must be triangles.
      CGAL_assertion (CGAL::circulator_size (j) == 3);
      ptr_face = new Face (this, CGAL::circulator_size (j));
      m_faces[faces_count++] = ptr_face;

      // add vertices of a facet and associate that facet to the vertices
      // Hoai - 17/08/2012 - XXX (poor performance)
      do
        {
          int src_id = std::distance (poly.vertices_begin (), j->vertex ());
          ptr_face->AddVert (m_verts[src_id]);
        }
      while (++j != fi->facet_begin ());
      j = fi->facet_begin ();
      do
        {
          int src_id = std::distance (poly.vertices_begin (), j->vertex ());
          m_verts[src_id]->AddFace (ptr_face);
        }
      while (++j != fi->facet_begin ());

      // iterate for creating edges
      j = fi->facet_begin ();
      do
        {
          int src_id =
            std::distance (poly.vertices_begin (), j->opposite ()->vertex ());
          int trg_id = std::distance (poly.vertices_begin (), j->vertex ());
          ptr_edge = m_verts[src_id]->Connected (m_verts[trg_id]);
          if (ptr_edge == NULL)
            {
              ptr_edge = new Edge (this, m_verts[src_id], m_verts[trg_id]);
              m_verts[src_id]->AddEdge (ptr_edge);
              m_verts[trg_id]->AddEdge (ptr_edge);
              m_edges.push_back (ptr_edge);
              m_numEdges++;
            }
          ptr_face->AddEdge (ptr_edge);
          ptr_edge->addFace (ptr_face);
        }
      while (++j != fi->facet_begin ());
    }

  // init m_r, m_C fields
  Point p1 (m_bBox.minx (), m_bBox.miny (), m_bBox.minz ());
  Point p2 (m_bBox.maxx (), m_bBox.maxy (), m_bBox.maxz ());
  m_r = RT (sqrt (CGAL::squared_distance (p1, p2).to_double ()));
  m_C = CGAL::midpoint (p1, p2);

  sortAdjInfo ();

  doPolyTest (NULL);
}
