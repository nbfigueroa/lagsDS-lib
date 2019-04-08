/*
 * Copyright (C) 2018 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author:  Nadia Figueroa
 * email:   nadia.figueroafernandez@epfl.ch
 * website: lasa.epfl.ch
 *
 * This work was supported by the EU project Cogimon H2020-ICT-23-2014.
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#include <stdlib.h>
#include <string>
#include "eigen3/Eigen/Dense"
#include "ros/ros.h"
#include "utils.h"

/* To use the class in the ds_motion_generators package */
#include "MathLib.h"

using namespace std;
using namespace Eigen;

const double PI_ = 3.14159265358979323846264338327950288419716939937510;

class lagsDS
{	
private:

        int 		K_;
        int 		M_;
        VectorXd	gamma_;
        double 		*Prior_;
        VectorXd 	*Mu_;
        MatrixXd 	*Sigma_;
        MatrixXd 	*A_g_matrix_;
        VectorXd    att_g_;
        MatrixXd    *A_l_matrix_;
        MatrixXd    *A_d_matrix_;
        VectorXd    *att_l_;
        VectorXd    *w_l_;
        double 		*b_l_;
        double      scale_;

public:             

        /* For Agnostic C++ use */
        lagsDS(const char  *path_dims);
        lagsDS(const char  *path_dims, const char  *path_prior,const char  *path_mu,const char  *path_sigma, const char  *path_Ag, const char  *path_Al, const char  *path_Ad, const char  *path_att_l, const char  *path_w_l, const char  *path_b_l, const char  *path_scale);
        lagsDS(const int K, const int M, const MatrixXd Priors_fMatrix, const MatrixXd Mu_fMatrix, const MatrixXd Sigma_fMatrix, const MatrixXd Ag_fMatrix, const MatrixXd Al_fMatrix, const MatrixXd Ad_fMatrix, const MatrixXd attl_fMatrix, const MatrixXd wl_fMatrix, const MatrixXd bl_fMatrix, const MatrixXd scale_fMatrix );

        /* For YAML/ROS Interface */
        lagsDS(const int K, const int M, const vector<double> Priors_vec, const vector<double> Mu_vec, const vector<double> Sigma_vec, const vector<double> Ag_vec);
        ~lagsDS(void);

        /*Initialization of Global Components */
        void        initialize_Ag(const char  *path_Ag);
        void        initialize_gamma(const char  *path_prior, const char  *path_mu, const char  *path_sigma);

        /*Initialization of Local Components */
        void        initialize_Al(const char  *path_Al);
        void        initialize_Ad(const char  *path_Ad);
        void        initialize_local_params(const char  *path_att_l, const char  *path_w_l, const char  *path_b_l, const char  *path_scale);


        /* For Agnostic C++ use */
        MatrixXd         compute_Ag(VectorXd xi);
        void             set_att_g(VectorXd att_g);
        VectorXd         compute_fg(VectorXd xi, VectorXd att);
        VectorXd         compute_gamma(VectorXd xi);
        VectorXd         compute_fl(VectorXd xi);

        /* For YAML/ROS Interface with DS motion generators */
        MathLib::Vector  compute_fg(MathLib::Vector xi, MathLib::Vector att);


private:

        /* General Stuff */
        void        setup_params();
        fileUtils   fileUtils_;
        void 	 	ERROR();

        /*Initialization of Global Components */
        void        initialize_Ag(const MatrixXd fMatrix);        
        void        initialize_Priors(const MatrixXd fMatrix);
        void        initialize_Mu(const MatrixXd fMatrix);
        void        initialize_Sigma(const MatrixXd fMatrix);

        /*Initialization of Local Components */
        void        initialize_Al(const MatrixXd fMatrix);
        void        initialize_Ad(const MatrixXd fMatrix);
        void        initialize_att_l(const MatrixXd fMatrix);
        void        initialize_w_l(const MatrixXd fMatrix);
        void        initialize_b_l(const MatrixXd fMatrix);

        void        initialize_Priors_vec(const vector<double> Priors_vec);
        void        initialize_Mu_vec(const vector<double> Mu_vec);
        void        initialize_Sigma_vec(const vector<double> Sigma_vec);
        void        initialize_Ag_vec(const vector<double> A_vec);

        /* Private Computations (for gaussian and local components) */
        VectorXd    compute_flk(VectorXd xi, int k);
        double      compute_hk(VectorXd xi, int k);
        VectorXd    compute_grad_hk(VectorXd xi, int k);
        double      compute_lambda_k(VectorXd xi, int k);
        double 		GaussianPDF(VectorXd x,VectorXd Mu,MatrixXd Sigma);
        double 		compute_rbf(double b, VectorXd xi, VectorXd center);


};





