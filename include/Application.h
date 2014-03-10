/*
 * Spectral Conformal Parameterization
 * Patrick Mullen, Yiying Tong, Pierre Alliez, Mathieu Desbrun.
 * Symposium of Geometry Processing, 2008.
 */

#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.h"
#include "Complex.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"
#include "Vector.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <iterator>

//Vector defining Constant vector field (in 2D)
#define CONSTANT_VECTOR (12.1, 74.1, 0.)

#define SEPARATOR "SEPARATOR"

//Constant Flux Factor
#define CONSTANT_FLUX_SCALE 100.

#define PI 3.1415926

namespace DDG
{
   class Application
   {
   public:

      void run(Mesh& mesh)
      {
         if (mesh.boundaries.empty())
         {
            std::cout << "Mesh has no boundary" << std::endl;
            return;
         }
         double initial_area = mesh.area();

         // create matrix for Lc
         int nV = mesh.vertices.size();
         SparseMatrix<Complex> Lc(nV, nV);
         
         // build energy matrix
         buildEnergy( mesh, Lc );

         Lc = Lc + SparseMatrix<Complex>::identity(mesh.vertices.size()) * Complex(0.00000001, 0.);

         mesh.Lc = Lc;
         
         // compute the solution; 
         DenseMatrix<Complex> x(nV, 1);
         SparseMatrix<Complex> star0;
         HodgeStar0Form<Complex>::build( mesh, star0 );

         // using one of the following functions
         //smallestEig(Lc, x);
         //smallestEig(Lc, star0, x);
         smallestEigPositiveDefinite(Lc, x);
         //smallestEigPositiveDefinite(Lc, star0, x);

         // for (int i = 0; i < nV; i++)
         // {
         //   printComplex(x(i, 0));
         //   cout << " ";
         // }

         // cout << "separation here" << endl;

         cout << "min energy: " << dot(x, Lc * x) << endl;
                  
         // then assign the solution
         assignSolution(x, mesh); // see below
         
         // rescale mesh for display convenience
         double scale = std::sqrt( initial_area / mesh.area() );
         normalizeMesh(scale, mesh);
      }

      void run2(Mesh& mesh)
      {
        //After we already have a parametrization, try to apply a flux through our boundary using a constant vector field

         double initial_area = mesh.area();

         // create matrix for Lc
         int nV = mesh.vertices.size();
         SparseMatrix<Complex> Lc(nV, nV);
         
         // build energy matrix
         buildEnergy( mesh, Lc );

         Lc = Lc + SparseMatrix<Complex>::identity(mesh.vertices.size()) * Complex(0.00000001, 0.);

         setBoundaryFluxes(mesh);

         DenseMatrix<Complex> g (nV, 1);
         buildFluxVector(mesh, g);

         // for (int i = 0; i < nV; i++)
         // {
         //  cout << g(i, 0) << endl;
         // }

         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
          printComplex(v->texture.x, v->texture.y);
         }
         cout << SEPARATOR << endl;

         DenseMatrix<Complex> xB (nV, 1);
         solvePositiveDefinite(Lc, xB, g);

         for (int i = 0; i < nV; i++)
         {
           cout << xB(i, 0) << endl;
         }
         cout << SEPARATOR << endl;

         assignSolution(xB, mesh);

         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
          printComplex(v->texture.x, v->texture.y);
         }
         cout << SEPARATOR << endl;

         // rescale mesh for display convenience
         double scale = std::sqrt( initial_area / mesh.area() );
         normalizeMesh(scale, mesh);

         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
          printComplex(v->texture.x, v->texture.y);
         }
         cout << SEPARATOR << endl;
      }

      void run3(Mesh& mesh)
      {
        double initial_area = mesh.area();

        // create matrix for Lc
         int nV = mesh.vertices.size();
         // SparseMatrix<Complex> Lc(nV, nV);
         // // build energy matrix
         // buildEnergy( mesh, Lc );
         // Lc = Lc + SparseMatrix<Complex>::identity(mesh.vertices.size()) * Complex(0.00001, 0.);

         //set the angles (represented as complex numbers) of each boundary edge
         setAbsolutePhaseOnBoundary(mesh);

         //build up the matrix M which converts unknown lengths to potentials for boundary vertices
         SparseMatrix<Complex> M (nV, nV);
         buildPositionCurvatureMatrix(mesh, M);

         printM(M, nV);

         // DenseMatrix<Complex> g (nV, 1);
         // buildSolutionVector(mesh, g);

         // SparseMatrix<Complex> EcM (nV, nV);
         // EcM = Lc * M;

         DenseMatrix<Complex> x (nV, 1);
         smallestEigPositiveDefinite(mesh.Lc, x, false);
         cout << "successfully solved for regular x, without M" << endl;


         SparseMatrix<Complex> MEcM (nV, nV);
         MEcM = M.transpose() * mesh.Lc * M;
         MEcM = MEcM + SparseMatrix<Complex>::identity(mesh.vertices.size()) * Complex(0.00001, 0.);

         //hold the solution
         DenseMatrix<Complex> xB (nV, 1);
         //solve(MEcM, xB, g);
         smallestEigPositiveDefinite(MEcM, xB, false);

         // for (int i = 0; i < nV; i++)
         // {
         //   cout << xB(i, 0) << endl;
         // }
         // cout << SEPARATOR << endl;

        DenseMatrix<Complex> MxB (nV, 1);
        MxB = M * xB;
        //printMxB(MxB, nV);

        // for (int i = 0; i < nV; i++)
        //  {
        //    cout << MxB(i, 0) << endl;
        //  }
        //  cout << SEPARATOR << endl;

        //assign solution
        assignSolution(MxB, mesh);

        // for( VertexIter v = mesh.vertices.begin();
        //      v != mesh.vertices.end();
        //      v ++ )
        //  {
        //   printComplex(v->texture.x, v->texture.y);
        //  }
        //  cout << SEPARATOR << endl;

        //printCurvatures(mesh);

        // rescale mesh for display convenience
         double scale = std::sqrt( initial_area / mesh.area() );
         normalizeMesh(scale, mesh);

         // for( VertexIter v = mesh.vertices.begin();
         //     v != mesh.vertices.end();
         //     v ++ )
         // {
         //  printComplex(v->texture.x, v->texture.y);
         // }
      }

      double sign(double a) {
        if (a < 0) return -1.;
        if (a > 0) return 1.;
        return 0.;
      }

      void printM(const SparseMatrix<Complex>& M, int nV)
      {
        cout << "printing M: " << endl;
        for (int i = 0; i < nV; i++)
        {
          for (int j = 0; j < nV; j++)
          {
            cout << M(i, j) << " ";
          }
          cout << endl;
        }
      }

      void printMxB(const DenseMatrix<Complex>& M, int nV)
      {
        cout << "printing MxB: " << endl;
        for (int i = 0; i < nV; i++)
        {
          cout << M(i, 0) << endl;
        }
        cout << endl;
      }
      
   protected:
      void buildEnergy(const Mesh& mesh, SparseMatrix<Complex>& A) const
      {
         // L 
         SparseMatrix<Complex> d0, star0, star1, Delta;
         ExteriorDerivative0Form<Complex>::build( mesh, d0 );
         HodgeStar0Form<Complex>::build( mesh, star0 );
         HodgeStar1Form<Complex>::build( mesh, star1 );
         SparseMatrix<Complex> L = d0.transpose() * star1 * d0;
            
         // Lc = L - A 
         A = L;
         for( FaceCIter f  = mesh.boundaries.begin();
                        f != mesh.boundaries.end();
                        f ++ )
         {
             HalfEdgeCIter he = f->he;
             do
             {
                int i = he->vertex->index;
                int j = he->next->vertex->index;
                A(i, j) = L(i, j) - Complex(0, 0.5);
                A(j, i) = L(j, i) + Complex(0, 0.5);
                
                he = he->next;
             }
             while( he != f->he );
         }
      }

      void buildPositionCurvatureMatrix(const Mesh& mesh, SparseMatrix<Complex>& M)
      {
        std::vector<int> bi; //boundary indices
        std::vector<Complex> ba; //boundary angles
        bool cont = false;
        for ( FaceCIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          do
          {
            if (!cont) {
              cont = true;
              he = he->next;
              continue;
            }
            bi.push_back(he->vertex->index);
            ba.push_back(he->phase);
            //cout << he->vertex->index << endl;
            for (unsigned int i = 0; i < bi.size(); i++)
            {
              M(he->vertex->index, bi.at(i)) = ba.at(i);
              //cout << "entered " << ba.at(i) << " at " << "(" << he->vertex->index << ", " << bi.at(i) << ")" << endl;
            }
            he = he->next;
          }
          while( he != f->he );
        }

        for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
        {
          if ( find(bi.begin(), bi.end(), v->index)==bi.end() )
          {
            M(v->index, v->index) = Complex(1., 0.);
          }
          // M(v->index, v->index) = Complex(1., 0.);
        }


        copy(bi.begin(), bi.end(), ostream_iterator<int>(cout, " "));

        // cout << "bi size: " << bi.size() << endl;
        // cout << "ba size: " << ba.size() << endl;
      }

      // void buildSolutionVector(const Mesh& mesh, DenseMatrix<Complex>& g)
      // {
      //   for ( FaceCIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
      //   {
      //     HalfEdgeIter he = f->he;
      //     do
      //     {
      //       g(he->vertex->index, 0) = he->length;   
             
      //       he = he->next;
      //     }
      //     while( he != f->he );
      //   }
      // }
      
      void assignSolution(const DenseMatrix<Complex>& x, Mesh& mesh)
      {
         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
            int i = v->index;
            const Complex& z = x(i,0);
            //cout << z << endl;
            v->texture.x = z.re;
            v->texture.y = z.im;
            v->texture.z = 0.0;
            // v->tag = true;

            // if (v->tag)
            // {
            //   cout << "tagged vertex position: " << "(" << v->texture.x << ", " << v->texture.y << ")" << endl;
            // }
            //cout << "vertex position: " << "(" << v->texture.x << ", " << v->texture.y << ")" << endl;
            /*v->position.x = z.re;
            v->position.y = z.im;
            v->position.z = 0;*/
         }
      }

      // void assignCurvatureSolution(const DenseMatrix<Complex>& x, Mesh& mesh)
      // {
      //   for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
      //   {
      //     HalfEdgeIter he = f->he;
      //     Complex z (0.0, 0.0);
      //     do
      //     {
      //       z += x(he->vertex->index, 0) * he->phase;

      //       printComplex(z);

      //       he->vertex->texture.x = z.re;
      //       he->vertex->texture.y = z.im;
      //       he->vertex->texture.z = 0.0;

      //       he = he->next;
      //     }
      //     while( he != f->he );
      //   }

      //   for( VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v ++ )
      //   {
      //     if (!v->he->onBoundary)
      //     {
      //       int i = v->index;
      //       const Complex& z = x(i,0);
      //       //printComplex(z);
      //       v->texture.x = z.re;
      //       v->texture.y = z.im;
      //       v->texture.z = 0.0;
      //       /*v->position.x = z.re;
      //       v->position.y = z.im;
      //       v->position.z = 0;*/
      //     }
      //   }
      // }
      
      void normalizeMesh(const double scale, Mesh& mesh)
      {
         cout << "scale: " << scale << endl;
         Vector avg;
         cout << "textures: " << endl;
         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
            avg += v->texture;
            //cout << v->texture << endl;
         }
         avg /= mesh.vertices.size();
         cout << "avg: " << avg << endl;
         
         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
            //cout << "before: " << v->texture << endl;
            //cout << "diff: " << v->texture - avg << endl;
            v->texture = scale*(v->texture - avg);
            //cout << "after: " << v->texture << endl;
         }
      }

      void printComplex(double a, double b)
      {
        cout << a << " + " << b << "i" << endl;
      }

      void printComplex(Complex a)
      {
        printComplex(a.re, a.im);
      }

      //Please provide the parametrized mesh
      void setBoundaryFluxes(Mesh& mesh) 
      {
        int count = 0;
        //cout << "Number of boundaries: " << mesh.boundaries.size() << endl;
        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          //Define some constant complex vector (will act as our constant vector field)
          Vector C CONSTANT_VECTOR;
          do
          {
             VertexIter v1 = he->vertex;
             VertexIter v2 = he->next->vertex;

             Vector texturedEdge = v2->texture - v1->texture;
             double L = texturedEdge.norm();
             Vector N (-1. * texturedEdge.y, texturedEdge.x, texturedEdge.z);

             //Find perpendicular component to each boundary edge (in both real and complex setting)
             //   and integrate over edge (multiply by edge length)
             //Extra sign term there to negate outward fluxes
             Vector F = (dot(C, N) / N.norm2()) * N * L * sign(dot(C, N));
             he->flux.re = F.x;
             he->flux.im = F.y;
             
             he = he->next;
             count++;
          }
          while( he != f->he );
        }
      }

      //Please provide parametrized mesh, and two tagged adjacent boundary vertices (defining a boundary edge)
      void setBoundaryFluxes2(Mesh& mesh)
      {
        cout << "Number of boundaries: " << mesh.boundaries.size() << endl;
        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;

          do
          {
            //Only set boundary flux through this edge, and the one right after
            if (he->vertex->tag && he->next->vertex->tag)
            {
              VertexIter v1 = he->vertex;
              VertexIter v2 = he->next->vertex;
              VertexIter v3 = he->next->next->vertex;
              Vector edge1 = v2->texture - v1->texture;
              Vector edge2 = v3->texture - v2->texture;
              // double L1 = edge1.norm();
              // double L2 = edge2.norm();
              Vector N1 (-1. * edge1.y, edge1.x, edge1.z);
              Vector N2 (-1. * edge2.y, edge2.x, edge1.z);

              he->flux.re = CONSTANT_FLUX_SCALE * N1.x / (N1.norm());
              he->flux.im = CONSTANT_FLUX_SCALE * N1.y / (N1.norm());

              he->next->flux.re = -1. * CONSTANT_FLUX_SCALE * N2.x / (N2.norm());
              he->next->flux.im = -1. * CONSTANT_FLUX_SCALE * N2.y / (N2.norm());

            } 
             
            he = he->next;
          }
          while( he != f->he );
        }
      }

      //Please provide parametrized mesh, and two sets of non-intersecting boundary edges (defined by 2 tagged boundary vertices)
      void setBoundaryFluxes3(Mesh& mesh)
      {
        cout << "Number of boundaries: " << mesh.boundaries.size() << endl;
        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          bool first = true;
          do
          {
            //Only set boundary flux through this edge, and the one right after
            if (he->vertex->tag && he->next->vertex->tag)
            {
              VertexIter v1 = he->vertex;
              VertexIter v2 = he->next->vertex;
              Vector e = v2->texture - v1->texture;
              // double L1 = edge1.norm();
              // double L2 = edge2.norm();
              Vector N (-1. * e.y, e.x, e.z);

              if (first)
              {
                he->flux.re = CONSTANT_FLUX_SCALE * N.x / (N.norm());
                he->flux.im = CONSTANT_FLUX_SCALE * N.y / (N.norm());
                first = false;
              }
              else
              {
                he->flux.re = -1. * CONSTANT_FLUX_SCALE * N.x / (N.norm());
                he->flux.im = -1. * CONSTANT_FLUX_SCALE * N.y / (N.norm());
              } 

            } 
             
            he = he->next;
          }
          while( he != f->he );
        }
      }

      void buildFluxVector(Mesh& mesh, DenseMatrix<Complex>& g)
      {
        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          do
          {
            g(he->next->vertex->index, 0) = 0.5 * (he->flux.conj() + he->next->flux.conj());   
             
            he = he->next;
          }
          while( he != f->he );
        }
      }

      //for unit complex vector a
      double angle(Vector a)
      {

        if (sign(a.y) == 0)
        {
          return 0.;
        }
        else if (sign(a.y) == -1)
        {
          return 2. * PI - acos(a.x);
        }
        else
        {
          return acos(a.x);
        }
      }

      double angleDiff(Vector a, Vector b)
      {
        return angle(a) - angle(b);
      }

      void printCurvatures(Mesh& mesh)
      {
        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          do
          {
             VertexIter v1 = he->vertex;
             VertexIter v2 = he->next->vertex;
             VertexIter v3 = he->next->next->vertex;

             Vector e1 = (v2->texture - v1->texture);
             Vector e2 = (v3->texture - v2->texture);
             //he->length = e1.norm();
             e1.normalize();
             e2.normalize();

             double angle = angleDiff(e1, e2);
             //printComplex(cos(angle), sin(angle));
             //cout << angle/(2*PI) << " pi" << endl;

             he->next->phase.re = cos(angle);
             he->next->phase.im = sin(angle);
             
             he = he->next;
          }
          while( he != f->he );
        }
      }

      //For now, just assumes a previously parametrized mesh- so that we can set curvatures at the boundary
      void setAbsolutePhaseOnBoundary(Mesh& mesh)
      {

        // bool first = true;
        // for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        // {
        //   HalfEdgeIter he = f->he;
        //   he->phase.re = 0.;
        //   he->phase.im = 0.;
        //   double angle = 0.;
        //   do
        //   {
        //     if (first)
        //     {
        //       first = false;
        //       he->vertex->tag = true;
        //     }
        //      VertexIter v1 = he->vertex;
        //      VertexIter v2 = he->next->vertex;
        //      VertexIter v3 = he->next->next->vertex;

        //      Vector e1 = (v2->texture - v1->texture);
        //      Vector e2 = (v3->texture - v2->texture);
        //      //he->length = e1.norm();
        //      e1.normalize();
        //      e2.normalize();

        //      angle += angleDiff(e1, e2);
        //      //cout << angle/(2*PI) << " pi" << endl;

        //      he->next->phase.re = cos(angle);
        //      he->next->phase.im = sin(angle);
             
        //      he = he->next;
        //   }
        //   while( he != f->he );
        // }

        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          // he->phase.re = 0.;
          // he->phase.im = 0.;
          double angle = 0.;
          do
          {
             VertexIter v = he->vertex;
             if (v->tag)
             {
               angle += PI / 2.;
             }
             //else don't add anything to angle

             //cout << angle/(PI) << " pi" << endl;

             he->next->phase.re = cos(angle);
             he->next->phase.im = sin(angle);

             //printComplex(he->next->phase);
             
             he = he->next;
          }
          while( he != f->he );
        }

        // int nV = 83;
        // bool first = true;

        // for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        // {
        //   HalfEdgeIter he = f->he;
        //   he->phase.re = 0.;
        //   he->phase.im = 0.;
        //   double angle = 0.;
        //   do
        //   {
        //     if (first)
        //     {
        //       first = false;
        //       he->vertex->tag = true;
        //     }
        //      angle += (2.*PI) / (double) nV;

        //      he->next->phase.re = cos(angle);
        //      he->next->phase.im = sin(angle);
             
        //      he = he->next;
        //   }
        //   while( he != f->he );
        // }

        for ( FaceIter f = mesh.boundaries.begin(); f != mesh.boundaries.end(); f++)
        {
          HalfEdgeIter he = f->he;
          do
          {
             printComplex(he->next->phase);
             he = he->next;
          }
          while( he != f->he );
        }
      }
   };
}

#endif
