/*
 * Copyright (C) 2019 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
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


#include "ros/ros.h"
#include "lagsDS.h"
#include "utils.h"
#include "eigen3/Eigen/Dense"
#include <vector>


using namespace std;

double         K;
double         M;
vector<double> Priors;
vector<double> Mu;
vector<double> Sigma;
vector<double> A_g;
vector<double> att_g;
vector<double> A_l;
vector<double> A_d;
vector<double> att_l;
vector<double> w_l;
vector<double> b_l;
double         scale;
double         b_g;
string         gpr_path;
string         path_model;

bool parseParams(const ros::NodeHandle& nh) {
    bool ret = true;

    if (!nh.getParam("model_path", path_model)) {
        ROS_ERROR("Couldn't retrieve model_path. ");
        ret = false;
    } else {
        cout << "Model path: "<< path_model << endl;
    }


    if (!nh.getParam("K", K))   {
        ROS_ERROR("Couldn't retrieve the number of guassians. ");
        ret = false;
    } else {
        cout << "Number of Components K: "<< (int)K << endl;
    }

    if (!nh.getParam("M", M))  {
        ROS_ERROR("Couldn't retrieve dimension of state. ");
        ret = false;
    } else {
        cout << "Dimensionality of state M: "<< (int)M << endl;
    }

    if (!nh.getParam("Priors", Priors))   {
        ROS_ERROR("Couldn't retrieve Priors. ");
        ret = false;
    } else {
        cout << "Priors: " << endl;
        for (int k = 0; k< int(K); k++)
            cout << Priors.at(k) << " ";
        cout << endl;
    }

    if (!nh.getParam("Mu", Mu))   {
        ROS_ERROR("Couldn't retrieve Mu. ");
        ret = false;
    } else {
        cout << "Mu: " << endl;
        for (int m = 0; m< int(K)*int(M); m++)
            cout << Mu.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("Sigma", Sigma))  {
        ROS_ERROR("Couldn't retrieve Sigma. ");
        ret = false;
    } else {
        cout << "Sigma [0]: " << endl;
        for (int m = 0; m< int(M)*int(M); m++)
            cout << Sigma.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("A_g", A_g))  {
        ROS_ERROR("Couldn't retrieve A_g. ");
        ret = false;
    } else {
        cout << "A_g [0]: " << endl;
        for (int m = 0; m< int(M)*int(M); m++)
            cout << A_g.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("att_g", att_g))  {
        ROS_ERROR("Couldn't retrieve att_g. ");
        ret = false;
    } else {
        cout << "att_g: " << endl;
        for (int m = 0; m< int(M); m++)
            cout << att_g.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("A_l", A_l))  {
        ROS_ERROR("Couldn't retrieve A_l. ");
        ret = false;
    } else {
        cout << "A_l [0]: " << endl;
        for (int m = 0; m< int(M)*int(M); m++)
            cout << A_l.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("A_d", A_d))  {
        ROS_ERROR("Couldn't retrieve A_d. ");
        ret = false;
    } else {
        cout << "A_d [0]: " << endl;
        for (int m = 0; m< int(M)*int(M); m++)
            cout << A_d.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("att_l", att_l))   {
        ROS_ERROR("Couldn't retrieve att_l. ");
        ret = false;
    } else {
        cout << "att_l: " << endl;
        for (int m = 0; m< int(K)*int(M); m++)
            cout << att_l.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("w_l", w_l))   {
        ROS_ERROR("Couldn't retrieve w_l. ");
        ret = false;
    } else {
        cout << "w_l: " << endl;
        for (int m = 0; m< int(K)*int(M); m++)
            cout << w_l.at(m) << " ";
        cout << endl;
    }

    if (!nh.getParam("b_l", b_l))   {
        ROS_ERROR("Couldn't retrieve b_l. ");
        ret = false;
    } else {
        cout << "b_l: " << endl;
        for (int k = 0; k< int(K); k++)
            cout << b_l.at(k) << " ";
        cout << endl;
    }


    if (!nh.getParam("scale", scale))   {
        ROS_ERROR("Couldn't retrieve the scale. ");
        ret = false;
    } else {
        cout << "scale: "<< scale << endl;
    }


    if (!nh.getParam("b_g", b_g))   {
        ROS_ERROR("Couldn't retrieve b_g. ");
        ret = false;
    } else {
        cout << "b_g: "<< b_g << endl;
    }

    if (!nh.getParam("gpr_path", gpr_path))   {
        ROS_ERROR("Couldn't retrieve gpr_path. ");
        ret = false;
    } else {
        cout << "gpr_path: "<< gpr_path << endl;
    }

    return ret;
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "test_lagsDS_node");
    ros::NodeHandle nh;
    ros::NodeHandle _nh("~");


    if(!parseParams(_nh)) {
        ROS_ERROR("Errors while parsing arguments.");
        return 1;
    }

    /* Instantiate lags-DS Model with parameters read from Yaml file*/
    lagsDS lagsDS_((int)K, (int)M, Priors, Mu, Sigma, A_g, att_g, A_l, A_d, att_l, w_l, b_l, scale, b_g, gpr_path);

    /* Testing the LAGS-DS on training data from MATLAB */
    cout << endl << "******** Testing Accuracy of model... ******** " << endl;

    string path_data        = path_model +  "Data";
    string path_xi_dot_g    = path_model +  "xi_dot_g";
    string path_xi_dot_l    = path_model +  "xi_dot_l";
    string path_xi_alpha    = path_model +  "xi_alpha";
    string path_xi_dot_lags = path_model +  "xi_dot_lags";
    string path_att_g    = path_model +  "att_g";
    MatrixXd attractor, Data, xi_dot_g, xi_dot_l, xi_alpha, xi_dot_lags;

    fileUtils fileUtils_;
    attractor    = fileUtils_.readMatrix(path_att_g.c_str());
    Data         = fileUtils_.readMatrix(path_data.c_str());
    xi_dot_g     = fileUtils_.readMatrix(path_xi_dot_g.c_str());
    xi_dot_l     = fileUtils_.readMatrix(path_xi_dot_l.c_str());
    xi_alpha     = fileUtils_.readMatrix(path_xi_alpha.c_str());
    xi_dot_lags  = fileUtils_.readMatrix(path_xi_dot_lags.c_str());
    int samples  = Data.cols();


    /* Fill in attractor */
    VectorXd att; att.resize(M);
    att = attractor.col(0);

    /* Fill in reference trajectories */
    MatrixXd xi_ref;  xi_ref.resize(M,samples);
    for (int i=0; i<M; i++)
        xi_ref.row(i) = Data.row(i);

    /* Compute estimated velocities from model */

    /* For Global Component Estimates */
    VectorXd xi_ref_test;  xi_ref_test.resize(M);
    VectorXd xi_dot_test;  xi_dot_test.resize(M);
    VectorXd xi_dot_mat;   xi_dot_mat.resize(M);
    MatrixXd A_matrix;     A_matrix.resize(M,M);

    /* For Local Component Estimates */
    VectorXd xi_dot_test_l;  xi_dot_test_l.resize(M);
    VectorXd xi_dot_mat_l; xi_dot_mat_l.resize(M);


    /* For LAGS Estimates */
    VectorXd xi_dot_test_lags;  xi_dot_test_lags.resize(M);
    VectorXd xi_dot_mat_lags; xi_dot_mat_lags.resize(M);

    /* For Error Estimates */
    VectorXd xi_dot_error;  xi_dot_error.resize(M);
    VectorXd est_error_1; est_error_1.resize(samples);
    VectorXd est_error_2; est_error_2.resize(samples);
    VectorXd est_error_3; est_error_3.resize(samples);
    VectorXd est_error_4; est_error_4.resize(samples);
    VectorXd est_error_5; est_error_5.resize(samples);
    double alpha_mat(0.0), alpha_test(0.0);


    lagsDS_.set_att_g(att);
    for (int i=0; i<samples; i++){

        /* Computing desired velocity */
        xi_ref_test  = xi_ref.col(i);
        xi_dot_mat   = xi_dot_g.col(i);

        /* Computing error between this estimate (using A-matrix) and MATLAB */
        A_matrix       =  lagsDS_.compute_Ag(xi_ref_test);
        xi_dot_test    =  A_matrix*(xi_ref_test - att);
        xi_dot_error   =  xi_dot_test-xi_dot_mat;
        est_error_1[i] =  xi_dot_error.norm();

        /* Computing error between this estimats (using f_g) and MATLAB */
        xi_dot_test    =  lagsDS_.compute_fg(xi_ref_test);
        xi_dot_error   =  xi_dot_test-xi_dot_mat;
        est_error_2[i] =  xi_dot_error.norm();


        /* Computing error between this estimats (using f_l) and MATLAB */
        xi_dot_mat_l   =  xi_dot_l.col(i);
        xi_dot_test_l  =  lagsDS_.compute_fl(xi_ref_test);
        xi_dot_error   =  xi_dot_test_l - xi_dot_mat_l;
        est_error_3[i] =  xi_dot_error.norm();

        /* Computing error between this estimate (using alpha) and MATLAB */
        alpha_mat    = xi_alpha.coeff(i,0);
        alpha_test   = lagsDS_.compute_alpha(xi_ref_test);
        est_error_4[i] = fabs(alpha_test-alpha_mat);

        /* Computing error between this estimate (using f-lags) and MATLAB */
        xi_dot_mat_lags   =  xi_dot_lags.col(i);
        xi_dot_test_lags  =  lagsDS_.compute_f(xi_ref_test);
        xi_dot_error      =  xi_dot_test_lags - xi_dot_mat_lags;
        est_error_5[i]    = xi_dot_error.norm();

    }

    /* Stats on Estimation error between MATLAB-C++ model */
    cout << "Average Estimation Error for Global Componnet f_g(x)" << " (Norm of predicted Matlab and C++ velocities): " << est_error_1.mean() << endl;
    cout << "Average Estimation Error for Global Componnet f_g(x)" << " (Norm of predicted Matlab and C++ velocities): " << est_error_2.mean() << endl;
    cout << "Average Estimation Error for Local Componnet f_l(x)"  << " (Norm of predicted Matlab and C++ velocities): " << est_error_3.mean() << endl;
    cout << "Average Estimation Error for Activation Function  alpha(x)"  << " (Absolute error of predicted Matlab and C++ alphas): " << est_error_4.mean() << endl;
    cout << "Average Estimation Error for LAGS-DS  f(x)"  << " (Norm of predicted Matlab and C++ velocities): " << est_error_5.mean() << endl;




    // Stop the node's resources
    ros::shutdown();
    // Exit tranquilly
    return 0;

}
