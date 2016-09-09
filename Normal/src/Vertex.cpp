#include <vector>
using namespace std;

#include "Vertex.h"
#include "Mesh.h"
#include "HalfEdge.h"

namespace DDG
{
   double Vertex::area( void ) const
   // returns the dual area associated with this vertex
   {
      double A = 0.;

      HalfEdgeCIter h = he;
      do
      {
         if (not h->onBoundary) A += h->face->area();
         h = h->flip->next;
      }
      while( h != he );

      return A / 3.;
   }

   int Vertex::normState = 0;

   Vector Vertex::normal( void ) const
   // returns the vertex normal
   {
       switch(normState)
       {
           case EQUAL: return this->normalEquallyWeighted();
           case AREA:  return this->normalAreaWeighted();
           case ANGLE: return this->normalAngleWeighted();
           default: return this->normalEquallyWeighted();
       }
   }
   
   Vector Vertex::normalEquallyWeighted() const
   // returns unite vertex normal using uniform weights
   {
      Vector N;

      HalfEdgeCIter h = he;
      do
      {
         if (not h->onBoundary) N += h->face->normal();
         h = h->flip->next;
      }
      while( h != he );

      return N.unit();
   }
      
   Vector Vertex::normalAreaWeighted() const
   // return unit vertex normal using face area weights
   {
      Vector N;
      double A = 0;

      HalfEdgeCIter h = he;
      do
      {
         if (not h->onBoundary) 
         { 
             N += h->face->normal() * h->face->area(); 
             A += h->face->area();
         }
         h = h->flip->next;
      }
      while( h != he );
    
      N /= A;

      return N.unit();
   }
 
   Vector Vertex::normalAngleWeighted() const
   // return uni vertex normal using tip angle weights
   {
      Vector N;
      double A = 0;

      HalfEdgeCIter h = he;
      do
      {
         HalfEdgeCIter nexth = h->flip->next;
         if (not nexth->onBoundary) 
         {
             double angle = nexth->next->angle();
             N += nexth->face->normal() * angle;
             A += angle;
         }
         h = nexth;
      }
      while( h != he );

      N /= A;
      
      return N.unit();
   }
 
   Vector Vertex::normalMeanCurvature() const
   // return unit mean curvature normal \lapla f
   {
       return this->normal();
   }

   Vector Vertex::normalSphereInscribed() const
   // returns unit sphere-inscribed normal
   {
       return this->normal();
   }

   vector<HalfEdge> isolated; // all isolated vertices point to isolated.begin()

   bool Vertex::isIsolated( void ) const
   // returns true if the vertex is not contained in any face or edge; false otherwise
   {
      return he == isolated.begin();
   }

   int Vertex :: valence( void ) const
   // returns the number of incident faces
   {
      int n = 0;

      HalfEdgeCIter h = he;
      do
      {
         n++;
         h = h->flip->next;
      }
      while( h != he );

      return n;
   }
   
   void Vertex :: toggleTag()
   {
      tag = !tag;
   }
}

