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


#include <stdio.h>
#include <fstream>
#include <time.h>
#include "eigen3/Eigen/Dense"
#include "lagsDS.h"
#include "utils.h"
#include "gaussian_process_regression/gaussian_process_regression.h"
#include "armadillo"
#include "GPRwrap.h"

using namespace arma;
using namespace std;

struct GPRTestData{
    int D;         // Input dimensions
    int M;         // samples
    mat x;         // input
    vec y;         // output
};

void loadTestData(string& f_GPRTestData, GPRTestData& data)
{

    cout << "\nGPR Testing File: " << f_GPRTestData << endl;
    ifstream fin(f_GPRTestData.c_str());
    unsigned int d , s;

    // Load GPRTest Data
    fin >> data.D
        >> data.M;

    cout << "data.D: " << data.D << endl
         << "data.M: " << data.M << endl;

    data.x        = zeros<mat>(data.D, data.M);
    data.y        = zeros<colvec>(data.M);

    for( d=0; d<data.D; d++ ){
        for( s=0; s<data.M; s++ ){
            fin >> data.x(d,s);
        }
    }
    for( s=0; s<data.M; s++ ){
           fin >> data.y(s);
    }
    fin.close();
}

int
main (int argc, char **argv)
{
    /* These are the parameters I have to load from a file */
    int D_in(2); float length_scale (0.25), sigma_f(1.0), sigma_n (0.01);

    cout << "*********** Starting Random test ***********" << endl;
    // input_dim and output_dim are integers specifying the dimensionality of input and output data respectively
    GaussianProcessRegression<float> myGPR(D_in, 1);
    
    // the hyperparameters are float values specifying length-scale, signal variance and observation noise variance
    myGPR.SetHyperParams(length_scale, sigma_f, sigma_n);
    Eigen::VectorXf train_input(D_in);
    Eigen::VectorXf train_output(1);

    // add some training data
    int n_train = 100;
    for(int k=0; k<n_train; k++){
        train_input.setRandom();
        train_output.setRandom();
        myGPR.AddTrainingData(train_input,train_output);
    }

    // get output at a testing point
    Eigen::VectorXf test_input(D_in);
    Eigen::VectorXf test_output(1);
    test_input.setRandom();
    test_output = myGPR.DoRegression(test_input);
    cout << "**************** finished ******************" << endl << endl;



    cout << "*********** Testing Class Wrapper with Loaded Model/Test-Data ***********" << endl;
    /* Start testing my class wrapper */
    string model_fileName = "/home/nbfigueroa/Dropbox/PhD_papers/LAGS-paper/new-code/lagsDS-opt/models/GPR-Models/iCub-Narrow-Passage_model.txt";
    // string data_fileName = "/home/nbfigueroa/Dropbox/PhD_papers/LAGS-paper/new-code/lagsDS-opt/models/GPR-Models/iCub-Narrow-Passage_data.txt";

    /* Instatiate Class and Load Model Parameters */
    GPRwrap gpr_(model_fileName);


    // GPRTestData data_;
    // loadTestData(data_fileName,data_);


    // string path_dim    = path_model +  "dimensions";
    // string path_Priors = path_model +  "Priors";
    // string path_Mu     = path_model +  "Mu";
    // string path_Sigma  = path_model +  "Sigma";
    // string path_A      = path_model +  "A_k";

    // /* Instantiate an LPV-DS class Option 1 */
    // cout << "Initialization Test 1: " << endl;
    // lagsDS lagsDS_test1(path_dim.c_str());
    // lagsDS_test1.initialize_gamma(path_Priors.c_str(), path_Mu.c_str(), path_Sigma.c_str());
    // lagsDS_test1.initialize_A(path_A.c_str());

    // /* Instantiate an LPV-DS class Option 2 */
    // cout << "Initialization Test 2: " << endl;
    // lagsDS lagsDS_test2 (path_dim.c_str(), path_Priors.c_str(), path_Mu.c_str(), path_Sigma.c_str(), path_A.c_str());

    // /* Instantiate an LPV-DS class Option 3 */
    // cout << "Initialization Test 3: " << endl;
    // fileUtils fileUtils_;
    // MatrixXd dim, Priors, Mu, Sigma, A;
    // dim     = fileUtils_.readMatrix(path_dim.c_str());
    // Priors  = fileUtils_.readMatrix(path_Priors.c_str());
    // Mu      = fileUtils_.readMatrix(path_Mu.c_str());
    // Sigma   = fileUtils_.readMatrix(path_Sigma.c_str());
    // A       = fileUtils_.readMatrix(path_A.c_str());
    // int K = (int)dim(0,0);
    // int M = (int)dim(1,0);
    // lagsDS lagsDS_test3 (K, M, Priors, Mu, Sigma, A);

    // /* Testing the LPV-DS on training data from MATLAB */
    // cout << "Testing Accuracy of model..." << endl;
    // string path_att    = path_model +  "attractor";
    // string path_data   = path_model +  "Data";
    // string path_xi_dot = path_model +  "xi_dot";
    // MatrixXd attractor, Data, xi_dot;
    // attractor = fileUtils_.readMatrix(path_att.c_str());
    // Data      = fileUtils_.readMatrix(path_data.c_str());
    // xi_dot    = fileUtils_.readMatrix(path_xi_dot.c_str());
    // int samples = Data.cols();

    // /* Fill in attractor */
    // VectorXd att; att.resize(M);
    // att = attractor.col(0);

    // /* Fill in reference trajectories */
    // MatrixXd xi_ref;  xi_ref.resize(M,samples);
    // for (int i=0; i<M; i++)
    //     xi_ref.row(i) = Data.row(i);


    // /* Compute estimated velocities from model */
    // VectorXd xi_ref_test;  xi_ref_test.resize(M);
    // VectorXd xi_dot_test;  xi_dot_test.resize(M);
    // VectorXd xi_dot_mat;   xi_dot_mat.resize(M);
    // VectorXd xi_dot_error;  xi_dot_error.resize(M);
    // MatrixXd A_matrix; A_matrix.resize(M,M);
    // VectorXd  est_error; est_error.resize(samples);
    // for (int i=0; i<samples; i++){

    //     /* Computing desired velocity */
    //     xi_ref_test = xi_ref.col(i);
    //     A_matrix = lagsDS_test3.compute_A(xi_ref_test);
    //     xi_dot_test = A_matrix*(xi_ref_test - att);

    //     /* Computing error between this estimate and MATLAB */
    //     xi_dot_mat = xi_dot.col(i);
    //     xi_dot_error =  xi_dot_test-xi_dot_mat;
    //     est_error[i] = xi_dot_error.norm();

    // }

    /* Stats on Estimation error between MATLAB-C++ model */
    // cout << "Average Estimation Error" << " (Norm of predicted Matlab and C++ velocities): " << est_error.mean() << endl;
    return 0;
}
