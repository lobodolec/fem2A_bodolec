#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#define M_PIL
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
            return 2 * M_PI * M_PI * sin(M_PI * v.x) * sin(M_PI * v.x);
        }
        
        double sol_ex_sin( vertex v)
        {
            return sin(M_PI * v.x) * sin(M_PI * v.y);
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
            DenseMatrix Ke;
            SparseMatrix K = SparseMatrix(mesh.nb_vertices());
            for ( int triangle = 0; triangle < mesh.nb_triangles(); triangle++ ) { 
                ElementMapping eltmap = ElementMapping(mesh, false, triangle);
                ShapeFunctions ref_func(2, 1);
                Quadrature quad = Quadrature::get_quadrature(2, false);
                assemble_elementary_matrix( eltmap, ref_func, quad, unit_fct, Ke );
                local_to_global_matrix( mesh, triangle, Ke, K );               
            }
            
            std::vector< double > values (mesh.nb_vertices()); 
            for (int vert = 0; vert < mesh.nb_vertices(); vert++) { 
                values[vert] = xy_fct(mesh.get_vertex(vert));
            }
            mesh.set_attribute(unit_fct, 1, true);
            std::vector< bool > attr_is_dir (2, false);
            attr_is_dir[1] = true; 
            std::vector< double > F (mesh.nb_vertices(), 0); 
            
            apply_dirichlet_boundary_conditions( mesh, attr_is_dir, values, K, F );

            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);
            
            std::string export_name = "pure_dirichlet_square";
            mesh.save(export_name + ".mesh"); 
            save_solution(u, export_name + ".bb"); 
            
        }
        
        void dirichlet_source_term_pb( const std::string& mesh_filename, bool verbose ) 
        {
            std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }            
            Mesh mesh;
            mesh.load(mesh_filename);
            DenseMatrix Ke;
            SparseMatrix K = SparseMatrix(mesh.nb_vertices());
            std::vector< double > F (mesh.nb_vertices(), 0);
            for ( int triangle = 0; triangle < mesh.nb_triangles(); triangle++ ) {
                ElementMapping eltmap = ElementMapping(mesh, false, triangle);
                ShapeFunctions ref_func(2, 1);
                Quadrature quad = Quadrature::get_quadrature(2, false);
                assemble_elementary_matrix( eltmap, ref_func, quad, unit_fct, Ke );
                local_to_global_matrix( mesh, triangle, Ke, K );
                std::vector< double > Fe;
                assemble_elementary_vector( eltmap, ref_func, quad, unit_fct, Fe );
                local_to_global_vector( mesh, false, triangle, Fe, F );
            } 
                      
            std::vector< double > values (mesh.nb_vertices()); 
            for (int vert = 0; vert < mesh.nb_vertices(); vert++) { 
                values[vert] = zero_fct(mesh.get_vertex(vert));
            }
            mesh.set_attribute(unit_fct, 1, true);
            std::vector< bool > attr_is_dir (2, false);
            attr_is_dir[1] = true;
                  
            apply_dirichlet_boundary_conditions( mesh, attr_is_dir, values, K, F );

            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);

            std::string export_name = "dirichlet_source_term_square";
            mesh.save(export_name + ".mesh"); 
            save_solution(u, export_name + ".bb");
        }
        
        void sinus_bump_pb ( const std::string& mesh_filename, bool verbose ) 
        {
            std::cout << "Solving a sinus bump problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }
            Mesh mesh;
            mesh.load(mesh_filename);
            DenseMatrix Ke;
            SparseMatrix K = SparseMatrix(mesh.nb_vertices());
            std::vector< double > F (mesh.nb_vertices(), 0);
            for ( int triangle = 0; triangle < mesh.nb_triangles(); triangle++ ) {
                ElementMapping eltmap = ElementMapping(mesh, false, triangle);
                ShapeFunctions ref_func(2, 1);
                Quadrature quad = Quadrature::get_quadrature(2, false);
                assemble_elementary_matrix( eltmap, ref_func, quad, unit_fct, Ke );
                local_to_global_matrix( mesh, triangle, Ke, K );
                std::vector< double > Fe;
                assemble_elementary_vector( eltmap, ref_func, quad, sin_bump, Fe );
                local_to_global_vector( mesh, false, triangle, Fe, F );
            } 
                     
            std::vector< double > values (mesh.nb_vertices()); 
            for (int vert = 0; vert < mesh.nb_vertices(); vert++) { 
                values[vert] = zero_fct(mesh.get_vertex(vert));
            }
            mesh.set_attribute(unit_fct, 1, true);
            std::vector< bool > attr_is_dir (2, false);
            attr_is_dir[1] = true; 
                  
            apply_dirichlet_boundary_conditions( mesh, attr_is_dir, values, K, F );
            
            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);
            
            std::string export_name = "sinus_bump_pb_square_fine";
            mesh.save(export_name + ".mesh"); 
            save_solution(u, export_name + ".bb"); 
        }
        
        void solexact_sin_pb ( const std::string& mesh_filename, bool verbose ) 
        {
            std::cout << "Exact solution for a sinus bump problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }
            Mesh mesh;
            mesh.load(mesh_filename);
            int N = mesh.nb_vertices();
            std::vector< double > u_exact(N, 0);
            for (int i = 0; i < N; i++){
                u_exact[i] = sol_ex_sin(mesh.get_vertex(i));
            }
    
            std::string export_name = "solexact_sin_pb_square_fine";
            mesh.save(export_name + ".mesh"); /* le maillage avec les nouveaux attributs */
            save_solution(u_exact, export_name + ".bb"); /* la solution du système */
        }
        
        void diff_sinexact ( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Difference between analytic and numeric solutions for a sinus bump problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details... " << std::endl;
            }
            Mesh mesh;
            mesh.load(mesh_filename);
            
            // numeric solution :
            DenseMatrix Ke;
            SparseMatrix K = SparseMatrix(mesh.nb_vertices());
            std::vector< double > F (mesh.nb_vertices(), 0);
            for ( int triangle = 0; triangle < mesh.nb_triangles(); triangle++ ) {
                ElementMapping eltmap = ElementMapping(mesh, false, triangle);
                ShapeFunctions ref_func(2, 1);
                Quadrature quad = Quadrature::get_quadrature(2, false);
                assemble_elementary_matrix( eltmap, ref_func, quad, unit_fct, Ke );
                local_to_global_matrix( mesh, triangle, Ke, K );
                std::vector< double > Fe;
                assemble_elementary_vector( eltmap, ref_func, quad, sin_bump, Fe );
                local_to_global_vector( mesh, false, triangle, Fe, F );
            } 
            std::vector< double > values (mesh.nb_vertices()); 
            for (int vert = 0; vert < mesh.nb_vertices(); vert++) { 
                values[vert] = zero_fct(mesh.get_vertex(vert));
            }
            mesh.set_attribute(unit_fct, 1, true);
            std::vector< bool > attr_is_dir (2, false);
            attr_is_dir[1] = true;   
            apply_dirichlet_boundary_conditions( mesh, attr_is_dir, values, K, F );
            
            std::vector< double > u_num(mesh.nb_vertices());
            solve(K, F, u_num);
            
            // difference = numeric - analytic solution
            int N = mesh.nb_vertices();
            std::vector< double > u_diff(N, 0);
            for (int i = 0; i < N; i++){
                u_diff[i] = u_num[i] - sol_ex_sin(mesh.get_vertex(i));
            }

            std::string export_name = "diff_sinexact_square_fine";
            mesh.save(export_name + ".mesh"); 
            save_solution(u_diff, export_name + ".bb"); 
        }

    }

}
