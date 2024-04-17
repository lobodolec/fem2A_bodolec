#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature(int ordre, bool bord)
        {
            Quadrature quad;
            quad = Quadrature::get_quadrature(ordre, bord);
            std::cout << quad.nb_points() << std::endl;
            double weight_sum = 0;
            for (int i = 0; i < quad.nb_points(); i++) {
                std::cout << quad.weight(i) << std::endl;
            	weight_sum =+ quad.weight(i);
            }
            std::cout << "somme des poids : " << weight_sum << std::endl;
            return true;
        }

        bool test_element_mapping(bool border, int i)
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping(mesh, border, i);
            return true;
        }
        
        bool test_transform()
        {	
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping map = ElementMapping(mesh, false, 4);
            vertex point ;
            point.x = 0.2;
            point.y = 0.4;
            std::cout << map.transform(point).x << " " 
            << map.transform(point).y << std::endl;
            return true;
        }
        
        bool test_jacobian_matrix()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping map = ElementMapping(mesh, false, 4);
            vertex point ;
            point.x = 0.2;
            point.y = 0.4;
            DenseMatrix mat = map.jacobian_matrix(point);
            mat.print();
            return true;
        }
        
        bool test_det_jacobian()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping map = ElementMapping(mesh, true, 4);
            vertex point ;
            point.x = 0.2;
            point.y = 0.4;
            std::cout << map.jacobian(point) << std::endl;
            return true;
        }
        
        bool test_shape_functions()
        {
            vertex point ;
            point.x = 0.2;
            point.y = 0.4;
            
            ShapeFunctions shapefunc1 = ShapeFunctions(1, 1);
            std::cout << shapefunc1.nb_functions() << std::endl;
            std::cout << shapefunc1.evaluate(1, point) << std::endl;
            vec2 g1 = shapefunc1.evaluate_grad(0, point);
            std::cout << g1.x << " " << g1.y << std::endl;
            
            ShapeFunctions shapefunc2 = ShapeFunctions(2, 1);
            std::cout << shapefunc2.nb_functions() << std::endl;
            std::cout << shapefunc2.evaluate(2, point) << std::endl;
            vec2 g2 = shapefunc2.evaluate_grad(2, point);
            std::cout << g2.x << " " << g2.y << std::endl;
            
            return true;
        }
        
        double diffusion ( vertex x_r )
        {
            double k = 1;
            return k;
        }
        
        bool test_K()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping map = ElementMapping(mesh, false, 4);
            ShapeFunctions ref_func(2, 1);
            Quadrature quad = Quadrature::get_quadrature(2, false); // tester avec 0 et 2
            
            DenseMatrix Ke;
            assemble_elementary_matrix( map, ref_func, quad, diffusion, Ke );
            Ke.print(); 
            
            SparseMatrix K = SparseMatrix( mesh.nb_vertices() );
            local_to_global_matrix(mesh, 4, Ke, K);
            K.print();
            
            return true;
        }
        
    }
}
