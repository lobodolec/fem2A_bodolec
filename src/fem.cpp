#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        std::cout << "[ElementMapping] constructor for element " << i << " ";
        if ( border ) std::cout << "(border)";
        std::cout << '\n';

        /* on remplit les attributs (constructeur !) avec les indexes locaux+globaux puis les coord*/
        // trouver le nombre de points
        // reshape vertices_ pour les contenir
        // le remplir avec les méthodes de mesh : avoir les indices locaux(calcul) et globaux(i n'est pas l'indice GLOBAL mais l'indice de l'élément entier !)
        
        if ( border ) {
            for (int i_local = 0; i_local<2; ++i_local) {
                vertices_.push_back(M.get_edge_vertex(i, i_local));
            }
            /* for (int i_local = 0; i_local<2; ++i_local) {
                std::cout << vertices_[i_local].x << " " << vertices_[i_local].y << std::endl
                ;} */
        }
        else {
            for (int i_local = 0; i_local<3; ++i_local) {
                vertices_.push_back(M.get_triangle_vertex(i, i_local));
            }
            /* for (int i_local = 0; i_local<3; ++i_local) {
                std::cout << vertices_[i_local].x << " " << vertices_[i_local].y << std::endl
                ;}*/
        } 
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        std::cout << "[ElementMapping] transform reference to world space" << '\n';
        
        /* les coordonnées de références sont dans l'attribut vertices_
        on applique les fonctions phi^ */
        vertex r ;
        if ( border_ ) {
            r.x = (1-x_r.x) * vertices_[0].x + x_r.x * vertices_[1].x;
            r.y = (1-x_r.x) * vertices_[0].y + x_r.x * vertices_[1].y;
        }
        else {
            r.x = (1-x_r.x-x_r.y) * vertices_[0].x + x_r.x 
                  * vertices_[1].x + x_r.y * vertices_[2].x ;
  	    r.y = (1-x_r.x-x_r.y) * vertices_[0].y + x_r.x 
  	          * vertices_[1].y + x_r.y * vertices_[2].y ;
  	}
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        
        DenseMatrix J ;
        if ( border_ ) {
            J.set_size(2, 1);
            J.set(0, 0, - vertices_[0].x + vertices_[1].x);
            J.set(1, 0, - vertices_[0].y + vertices_[1].y);

        }
        else {
       	    J.set_size(2, 2);
            J.set(0, 0, - vertices_[0].x + vertices_[1].x);
            J.set(1, 0, - vertices_[0].y + vertices_[1].y);
            J.set(0, 1, - vertices_[0].x + vertices_[2].x);
            J.set(1, 1, - vertices_[0].y + vertices_[2].y);
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        DenseMatrix J = jacobian_matrix( x_r );
        double det_jac;
        if ( border_ ) {
            double JTJ = J.get(0,0)*J.get(0,0) + J.get(1,0)*J.get(1,0);
            det_jac = sqrt(JTJ);
        }
        else {
            det_jac = J.det_2x2();
        }     
        return det_jac;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        // (elle concerne les fonctions d'interpolation!) 
        bool shape_func_constr = true;
        if ( dim != 1 && dim != 2 ) {
            std::cout << "ShapeFunctions is only implemented in 1D or 2D" << std::endl;
            shape_func_constr = false;
        }
        if ( order != 1 ) {
            std::cout << "Only order-1 functoons are implemented" << std::endl;
            shape_func_constr = false;
        assert (shape_func_constr);
        }
    }

    int ShapeFunctions::nb_functions() const
    {
        std::cout << "[ShapeFunctions] number of functions" << '\n';
        return dim_ + 1 ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        // qu'apporte un switch case par rapport à un if ici ? 
        // juste plus simple à lire et à écrire
        if ( dim_ == 1 ) {
            switch (i) {
            	case 0 :
            	    return 1 - x_r.x ;
  
                case 1 :
            	    return x_r.x ;
            }
        }
        else {
            switch (i) {
                case 0 :
                    return 1 - x_r.x - x_r.y ;
                case 1 :
                    return x_r.x ;
                case 2 :
                    return x_r.y ;
            }
        }
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        // quel intérêt d'utiliser le type vec2 plutot que la structure vertex, dans DenseMatrix ?
        // le gradient de phi dans l'espace réel (g) est le gradient de phi^ dans l'espace de  	référence (g_ref), multiplié par la matrice jacobienne            
        vec2 g;
        vec2 g_ref;
        DenseMatrix J = jacobian_matrix( x_r );
        // PB : jacobian_matrix est une méthode de ElementMaillage, donc on doit introduire un mapping !, mais aucun n'est passé en argument : comment savoir quels paramètres choisir ?

        if ( dim_ == 1 ) {
            switch (i) {
            	case 0 :
            	    g_ref.x = -1 ;
            	    g_ref.y = 0 ;
            	    g = J.mult_2x2_2(g_ref);
                case 1 :
                    g_ref.x = 1 ;
            	    g_ref.y = 0 ;
            	    g = J.mult_2x2_2(g_ref);
            }
        }
        else {
            switch (i) {
                case 0 :
                    g_ref.x = -1 ;
            	    g_ref.y = -1 ;
            	    g = J.mult_2x2_2(g_ref);
                case 1 :
                    g_ref.x = 1 ;
            	    g_ref.y = 0 ;
            	    g = J.mult_2x2_2(g_ref);
                case 2 :
                    g_ref.x = 0 ;
            	    g_ref.y = 1 ;
            	    g = J.mult_2x2_2(g_ref);
            }
        }
        
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        std::cout << "compute elementary matrix" << '\n';
        // TODO (utiliser les méthodes de la classe DenseMatrix; 
        // double (*coefficient)(vertex) est un pointeur de fonction : help.md dans la doc
        //pour un produit scalaire avec des vecteurs : dot (.,.)
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        std::cout << "Ke -> K" << '\n';
        // TODO
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        // TODO
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        // TODO
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        // TODO
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
