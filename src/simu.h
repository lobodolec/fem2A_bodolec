#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }
        
        double sin_bump( vertex v )
        {
            
            return 
        }

        //#################################
        //  Simulations
        //#################################
        
        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }            
            Mesh mesh;
            mesh.load(mesh_filename);
            // on parcourt le maillage :
            DenseMatrix Ke;
            SparseMatrix K = SparseMatrix(mesh.nb_vertices());
            for ( int triangle = 0; triangle < mesh.nb_triangles(); triangle++ ) { 
            /* assemble_elt_matrix est définie sur des triangles uniquement !*/
                // calcul de Ke : avec k = 1, donc unit_fc
                ElementMapping eltmap = ElementMapping(mesh, false, triangle);
                ShapeFunctions ref_func(2, 1);
                Quadrature quad = Quadrature::get_quadrature(2, false);
                assemble_elementary_matrix( eltmap, ref_func, quad, unit_fct, Ke );
                local_to_global_matrix( mesh, triangle, Ke, K );               
            }
            
            // set_attribute(..., 1, true pour les bords puis false pour les triangles, pas besoin de boucler sur les elements); = pas encore pour cette simu; les bords positifs seront de dirichlet (ex : attr=1 => unit_fct, coefficient permettra de différencier les régions de l'enveloppe) 
            
            // application de la condition Dirichlet :            
            std::vector< double > values (mesh.nb_vertices()); 
            for (int vert = 0; vert < mesh.nb_vertices(); vert++) { 
                values[vert] = xy_fct(mesh.get_vertex(vert)); /* la condition de Dirichlet est u = x + y */
            }
            mesh.set_attribute(unit_fct, 1, true);
            std::vector< bool > attr_is_dir (2, false);
            attr_is_dir[1] = true; /* = [0, 1] 1 si on applique Dirichlet et 0 sinon */
            std::vector< double > F (mesh.nb_vertices(), 0); /* 0 car pas de terme source*/
            
            apply_dirichlet_boundary_conditions( mesh, attr_is_dir, values, K, F );
            
            // resolution du systeme lineaire
            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);
            
            // enregistrement :
            std::string export_name = "pure_dirichlet_square";
            mesh.save(export_name + ".mesh"); /* le maillage avec les nouveaux attributs */
            save_solution(u, export_name + ".bb"); /* la solution du système */
            
        }
        
        void dirichlet_source_term_pb( const std::string& mesh_filename, bool verbose ) 
        {
            std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }            
            Mesh mesh;
            mesh.load(mesh_filename);
            // on parcourt le maillage :
            DenseMatrix Ke;
            SparseMatrix K = SparseMatrix(mesh.nb_vertices());
            std::vector< double > F (mesh.nb_vertices(), 0);
            for ( int triangle = 0; triangle < mesh.nb_triangles(); triangle++ ) {
                // calcul de Ke : avec k = 1, donc unit_fc
                ElementMapping eltmap = ElementMapping(mesh, false, triangle);
                ShapeFunctions ref_func(2, 1);
                Quadrature quad = Quadrature::get_quadrature(2, false);
                assemble_elementary_matrix( eltmap, ref_func, quad, unit_fct, Ke );
                local_to_global_matrix( mesh, triangle, Ke, K );
                // calcul de Fe : avec f = 1, donc unit_fct 
                std::vector< double > Fe;
                assemble_elementary_vector( eltmap, ref_func, quad, unit_fct, Fe );
                local_to_global_vector( mesh, false, triangle, Fe, F );
            } 
            
            // application de la condition Dirichlet :            
            std::vector< double > values (mesh.nb_vertices()); 
            for (int vert = 0; vert < mesh.nb_vertices(); vert++) { 
                /* la condition de Dirichlet est u = 0 */
                values[vert] = zero_fct(mesh.get_vertex(vert));
            }
            mesh.set_attribute(unit_fct, 1, true);
            std::vector< bool > attr_is_dir (2, false);
            attr_is_dir[1] = true; /* = [0, 1] , 1 si on applique Dirichlet et 0 sinon */
                  
            apply_dirichlet_boundary_conditions( mesh, attr_is_dir, values, K, F );
            
            // resolution du systeme lineaire
            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);
            
            // enregistrement :
            std::string export_name = "dirichlet_source_term_square_fine";
            mesh.save(export_name + ".mesh"); /* le maillage avec les nouveaux attributs */
            save_solution(u, export_name + ".bb"); /* la solution du système */
        }
        
        void sinus_bump_pb ( const std::string& mesh_filename, bool verbose ) 
        {
            std::cout << "Solving a sinus bump problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }            
            /* Mesh mesh;
            mesh.load(mesh_filename);
            // on parcourt le maillage :*/
        }

    }

}
