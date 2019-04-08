/*
 * Copyright (C) 2018 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author:  Sina Mirrazavi and Nadia Figueroa
 * email:   {sina.mirrazavi,nadia.figueroafernandez}@epfl.ch
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


#include "lagsDS.h"

lagsDS::lagsDS(const char  *path_dims):K_(0),M_(0) {

    /* Declare the number of the components
     * and the dimension  of the state  */
    MatrixXd fMatrix(1,1);fMatrix.setZero();

    if (fileUtils_.is_file_exist(path_dims))
        fMatrix=fileUtils_.readMatrix(path_dims);
    else{
        cout<<"The provided path for Dimensions does not exist:"<< path_dims << endl;
        ERROR();
    }
    if ((fMatrix.rows()!=2)){
        cout<<"Initialization of the Dimensions is wrong."<<endl;
        ERROR();
    }

    K_ = (int)fMatrix.coeff(0,0);
    M_ = (int)fMatrix.coeff(1,0);
    setup_params();
}


lagsDS::lagsDS(const char  *path_dims, const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_, const char  *path_Ag_, const char *path_Al_, const char *path_Ad_, const char *path_att_l_, const char *path_w_l_, const char *path_b_l_, const char *path_scale_):K_(0),M_(0) {

    /* Declare the number of the components
     * and the dimension  of the state  */
    MatrixXd fMatrix(1,1);fMatrix.setZero();

    if ( fileUtils_.is_file_exist(path_dims))
        fMatrix= fileUtils_.readMatrix(path_dims);
    else{
        cout<<"The provided path for Dimensions does not exist:"<< path_dims << endl;
        ERROR();
    }
    if ((fMatrix.rows()!=2)){
        cout<<"Initialization of the Dimensions is wrong."<<endl;
        ERROR();
    }

    K_ = (int)fMatrix.coeff(0,0);
    M_ = (int)fMatrix.coeff(1,0);
    setup_params();
    initialize_gamma(path_prior_, path_mu_, path_sigma_);
    initialize_Ag(path_Ag_);
    initialize_Al(path_Al_);
    initialize_Ad(path_Ad_);
    initialize_local_params(path_att_l_, path_w_l_, path_b_l_, path_scale_);

}


lagsDS::lagsDS(int K, int M, const MatrixXd Priors_fMatrix, const MatrixXd Mu_fMatrix, const MatrixXd Sigma_fMatrix, const MatrixXd Ag_fMatrix, const MatrixXd Al_fMatrix, const MatrixXd Ad_fMatrix, const MatrixXd attl_fMatrix, const MatrixXd wl_fMatrix, const MatrixXd bl_fMatrix, const MatrixXd scale_fMatrix):K_(K),M_(M) {

    /* Given the parameters directly as MatrixXd, initialize all matrices*/
    setup_params();
    initialize_Priors(Priors_fMatrix);
    initialize_Mu(Mu_fMatrix);
    initialize_Sigma(Sigma_fMatrix);
    initialize_Ag(Ag_fMatrix);
    initialize_Al(Al_fMatrix);
    initialize_Ad(Ad_fMatrix);
    initialize_att_l(attl_fMatrix);
    initialize_w_l(wl_fMatrix);
    initialize_b_l(bl_fMatrix);

    cout<<"** Initializing scale **"<< endl;
    scale_ = scale_fMatrix.coeff(0,0);
}


lagsDS::lagsDS(const int K, const int M, const vector<double> Priors_vec, const vector<double> Mu_vec, const vector<double> Sigma_vec, const vector<double> Ag_vec):K_(K),M_(M){


    /* Given the parameters directly as vector<double>, initialize all matrices*/
    setup_params();
    initialize_Priors_vec(Priors_vec);
    initialize_Mu_vec(Mu_vec);
    initialize_Sigma_vec(Sigma_vec);
    initialize_Ag_vec(Ag_vec);
}


lagsDS::~lagsDS(){

}

void lagsDS::ERROR()
{
    while(ros::ok())
    {

    }
}

void lagsDS::setup_params()
{

    /* Setup matrices for Global Component */
    A_g_matrix_ = new MatrixXd[K_]; for(int s=0; s<K_; s++ ){A_g_matrix_[s].resize(M_,M_);}
    Prior_      = new double[K_];
    Mu_         = new VectorXd[K_]; for(int s=0; s<K_; s++ ){	Mu_[s].resize(M_);	}
    Sigma_      = new MatrixXd[K_]; for(int s=0; s<K_; s++ ){	Sigma_[s].resize(M_,M_);	}

    att_g_.resize(M_);
    att_g_.setZero();

    gamma_.resize(K_);
    gamma_.setZero();

    /* Setup matrices for Local Component */
    A_l_matrix_ = new MatrixXd[K_]; for(int s=0; s<K_; s++ ){A_l_matrix_[s].resize(M_,M_);}
    A_d_matrix_ = new MatrixXd[K_]; for(int s=0; s<K_; s++ ){A_d_matrix_[s].resize(M_,M_);}
    att_l_      = new VectorXd[K_]; for(int s=0; s<K_; s++ ){att_l_[s].resize(M_);	}
    w_l_        = new VectorXd[K_]; for(int s=0; s<K_; s++ ){w_l_[s].resize(M_);	}
    b_l_        = new double[K_];

    cout << "Initialized an M:" << M_ << " dimensional LAGS-DS with K: " << K_ << " Global Components" << endl;
    cout << "Initialized an M:" << M_ << " dimensional LAGS-DS with K: " << K_ << " Locally-Active Components" << endl;
    cout << "Initialized an M:" << M_ << " dimensional LAGS-DS with K: " << K_ << " Locally-Deflective Components" << endl;
}


/***************************************************/
/* Initialization functions with MatrixXd as input */
/***************************************************/

void lagsDS::initialize_Priors(const MatrixXd fMatrix ){

    cout<<"** Initializing Priors **"<< endl;
    if ((fMatrix.cols()!=K_)||(fMatrix.rows()!=1))
    {
        cout<<"Initialization of Prior is wrong."<<endl;
        cout<<"Number of components is: "<<K_<<endl;
        cout<<"Dimension of states of Prior is: "<<fMatrix.cols()<<endl;
        ERROR();
    }

    for (int i=0; i<K_; i++)
        Prior_[i]=fMatrix(0,i);
}


void lagsDS::initialize_Mu(const MatrixXd fMatrix ){

    cout<<"** Initializing Mu **"<< endl;
    if ((fMatrix.cols()!=K_)||(fMatrix.rows()!=M_)){
        cout<<"Initialization of Mean is wrong."<<endl;
        cout<<"Number of components is: "<<K_<<endl;
        cout<<"Dimension of states is: "<<M_<<endl;
        cout<<"Dimension of states of Mean is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
        ERROR();
    }

    for(int s=0; s<K_; s++ )
        Mu_[s]=fMatrix.col(s);
}

void lagsDS::initialize_Sigma(const MatrixXd fMatrix ){
    cout<<"** Initializing Sigma **"<< endl;
    if ((fMatrix.rows()!=K_*M_)||(fMatrix.cols()!=M_))
    {
        cout<<"Initialization of the covariance matrix is wrong."<<endl;
        cout<<"the covariance matrix : "<<endl;cout<<fMatrix<<endl;
        cout<<"Number of components is: "<< K_ <<endl;
        cout<<"Dimension of states is: "<< M_ <<endl;
        cout<<"Dimension of states of the covariance matrix is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
        ERROR();
    }

    int s=0;

    for(int i=0; i<K_; i++ ){
        for(int j=0; j<M_; j++ ){
            Sigma_[i].row(j)=fMatrix.row(s);
            s++;
        }
    }
}

void lagsDS::initialize_Ag(const MatrixXd fMatrix ){

    cout<<"** Initializing A_g's' **"<< endl;
    if ((fMatrix.rows()!=K_*M_)||(fMatrix.cols()!=M_))
    {
        cout<<"Initialization of the A_g matrices is wrong!!"<<endl;
        cout<<"Ag_k: "<< endl << fMatrix<< endl;
        cout<<"[Proposed Dimensionality] K : "<<K_<<" M:"<<M_<<endl;
        cout<<"[Actual Dimensionality] K : "<< fMatrix.cols()/M_ <<" M:" << fMatrix.rows() << endl;
        ERROR();
    }
    int j = 0;
    for(int s=0; s<K_; s++ ){
        for(int i=0; i<M_; i++ ){
            A_g_matrix_[s].row(i)=fMatrix.row(j);
            j++;
        }
    }
}

void lagsDS::initialize_Al(const MatrixXd fMatrix ){

    cout<<"** Initializing A_l's' **"<< endl;
    if ((fMatrix.rows()!=K_*M_)||(fMatrix.cols()!=M_))
    {
        cout<<"Initialization of the A_l matrices is wrong!!"<<endl;
        cout<<"Al_k: "<< endl << fMatrix<< endl;
        cout<<"[Proposed Dimensionality] K : "<<K_<<" M:"<<M_<<endl;
        cout<<"[Actual Dimensionality] K : "<< fMatrix.cols()/M_ <<" M:" << fMatrix.rows() << endl;
        ERROR();
    }
    int j = 0;
    for(int s=0; s<K_; s++ ){
        for(int i=0; i<M_; i++ ){
            A_l_matrix_[s].row(i)=fMatrix.row(j);
            j++;
        }
    }
}


void lagsDS::initialize_Ad(const MatrixXd fMatrix ){

    cout<<"** Initializing A_d's' **"<< endl;
    if ((fMatrix.rows()!=K_*M_)||(fMatrix.cols()!=M_))
    {
        cout<<"Initialization of the A_d matrices is wrong!!"<<endl;
        cout<<"Ad_k: "<< endl << fMatrix<< endl;
        cout<<"[Proposed Dimensionality] K : "<<K_<<" M:"<<M_<<endl;
        cout<<"[Actual Dimensionality] K : "<< fMatrix.cols()/M_ <<" M:" << fMatrix.rows() << endl;
        ERROR();
    }
    int j = 0;
    for(int s=0; s<K_; s++ ){
        for(int i=0; i<M_; i++ ){
            A_d_matrix_[s].row(i)=fMatrix.row(j);
            j++;
        }
    }
}


void lagsDS::initialize_att_l(const MatrixXd fMatrix ){

    cout<<"** Initializing att_l **"<< endl;
    if ((fMatrix.cols()!=K_)||(fMatrix.rows()!=M_)){
        cout<<"Initialization of Local Virtual Attractors is wrong."<<endl;
        cout<<"Number of components is: "<<K_<<endl;
        cout<<"Dimension of states is: "<<M_<<endl;
        cout<<"Dimension of states of attractors is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
        ERROR();
    }

    for(int s=0; s<K_; s++ )
        att_l_[s]=fMatrix.col(s);
}


void lagsDS::initialize_w_l(const MatrixXd fMatrix ){

    cout<<"** Initializing w_l **"<< endl;
    if ((fMatrix.cols()!=K_)||(fMatrix.rows()!=M_)){
        cout<<"Initialization of Local w vectors is wrong."<<endl;
        cout<<"Number of components is: "<<K_<<endl;
        cout<<"Dimension of states is: "<<M_<<endl;
        cout<<"Dimension of states of attractors is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
        ERROR();
    }

    for(int s=0; s<K_; s++ )
        w_l_[s]=fMatrix.col(s);
}


void lagsDS::initialize_b_l(const MatrixXd fMatrix ){

    cout<<"** Initializing b_l **"<< endl;
    if ((fMatrix.cols()!=K_)||(fMatrix.rows()!=1))
    {
        cout<<"Initialization of b_l is wrong."<<endl;
        cout<<"Number of components is: "<<K_<<endl;
        cout<<"Dimension of states of b_l is: "<<fMatrix.cols()<<endl;
        ERROR();
    }

    for (int i=0; i<K_; i++)
        b_l_[i]=fMatrix(0,i);
}

/*********************************************************/
/* Initialization functions with vector<double> as input */
/*********************************************************/

void lagsDS::initialize_Priors_vec(const vector<double> Priors_vec){
    cout<<"** Initializing Priors **"<< endl;
    if (Priors_vec.size() != K_){
        cout<<"Initialization of Prior is wrong."<<endl;
        cout<<"Number of components is: "<<K_<<endl;
        cout<<"Dimension of states of Prior is: "<<Priors_vec.size()<<endl;
        ERROR();
    }

    for (int k=0; k<K_; k++)
        Prior_[k] = Priors_vec[k];
}


void lagsDS::initialize_Mu_vec(const vector<double> Mu_vec){

    cout<<"** Initializing Mu **"<< endl;
    if (Mu_vec.size() != K_*M_){
        cout<<"Initialization of Mean is wrong."<<endl;
        cout<<"Size of vector should be K ("<<K_ << ")*M("<< M_<< ") =" << K_*M_ <<endl;
        cout<<"Given vector is of size "<< Mu_vec.size() << endl;
        ERROR();
    }


    for(int k=0; k<K_; k++ ){
        VectorXd Mu_k; Mu_k.resize(M_); Mu_k.setZero();
        for (int m = 0; m < M_; m++)
            Mu_k[m] = Mu_vec[k * M_ + m];
        Mu_[k]=Mu_k;
    }

    /* For Debugging */
//    for(int k=0; k<K_; k++ )
//        cout << "Mu["<< k << "]"<< endl << Mu_[k] << endl;
}


void lagsDS::initialize_Sigma_vec(const vector<double> Sigma_vec){

    cout<<"** Initializing Sigma **"<< endl;
    if (Sigma_vec.size() != K_*M_*M_){
        cout<<"Initialization of Sigma is wrong."<<endl;
        cout<<"Size of vector should be K ("<<K_ << ")*M("<< M_<< ")*M(" << M_<< ") =" << K_*M_*M_ <<endl;
        cout<<"Given vector is of size "<< Sigma_vec.size() << endl;
        ERROR();
    }

    for(int k=0; k<K_; k++ ){
        MatrixXd Sigma_k; Sigma_k.resize(M_,M_); Sigma_k.setZero();
        for (int row = 0; row < M_; row++) {
            for (int col = 0; col < M_; col++) {
                int ind = k * M_ * M_ + row * M_ + col;
                Sigma_k(col,row)  = Sigma_vec[ind];
            }
        }
        Sigma_[k] = Sigma_k;
    }

    /* For Debugging */
//    for(int k=0; k<K_; k++ )
//        cout << "Sigma["<< k << "]"<< endl << Sigma_[k] << endl;

}


void lagsDS::initialize_Ag_vec(const vector<double> Ag_vec){
    cout<<"** Initializing A **"<< endl;
    if (Ag_vec.size() != K_*M_*M_){
        cout<<"Initialization of A-matrices is wrong."<<endl;
        cout<<"Size of vector should be K ("<<K_ << ")*M("<< M_<< ")*M(" << M_<< ") =" << K_*M_*M_ <<endl;
        cout<<"Given vector is of size "<< Ag_vec.size() << endl;
        ERROR();
    }

    for(int k=0; k<K_; k++ ){
        MatrixXd A_k; A_k.resize(M_,M_); A_k.setZero();
        for (int row = 0; row < M_; row++) {
            for (int col = 0; col < M_; col++) {
                int ind = k * M_ * M_ + row * M_ + col;
                A_k(col,row)  = Ag_vec[ind];
            }
        }
        A_g_matrix_[k] = A_k;
    }

    /* For Debugging */
//    for(int k=0; k<K_; k++ )
//        cout << "A["<< k << "]"<< endl << A_g_matrix_[k] << endl;
}



/************************************************************/
/* Initialization functions with path to text file as input */
/************************************************************/

void lagsDS::initialize_Ag(const char  *path_Ag_){

    /* Initialize A
     * path_Ag_ is the path of A matrix*/

    MatrixXd fMatrix(1,1);fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_Ag_))
        fMatrix=fileUtils_.readMatrix(path_Ag_);
    else{
        cout<<"The provided path does not exist: "<<path_Ag_<<endl;
        ERROR();
    }

    initialize_Ag(fMatrix);
}

void lagsDS::initialize_gamma(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_){

	/* Initialize scheduling/activation function parameters
	 *  path_prior_ is the path of the prior matrix
	 *	path_mu_ is the path of the mean matrix
	 *	path_sigma_ is the path of the covariance matrix */

    /* Initializing Priors*/
    MatrixXd fMatrix;
    if (fileUtils_.is_file_exist(path_prior_))
        fMatrix=fileUtils_.readMatrix(path_prior_);
	else{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_prior_lagsDS "<<endl;
		cout<<path_prior_<<endl;
		ERROR();
	}
    initialize_Priors(fMatrix);


    /* Initializing Mu*/
	fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_mu_))
        fMatrix=fileUtils_.readMatrix(path_mu_);
	else{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_mu_lagsDS "<<endl;
		cout<<path_mu_<<endl;
		ERROR();
	}
    initialize_Mu(fMatrix);

    /* Initializing Sigma*/
	fMatrix.resize(1,1);fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_sigma_))
        fMatrix=fileUtils_.readMatrix(path_sigma_);
	else{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_sigma_lagsDS "<<endl;
		cout<<path_sigma_<<endl;
		ERROR();
	}
    initialize_Sigma( fMatrix );
}


void lagsDS::initialize_Al(const char  *path_Al_){

    /* Initialize A
     * path_Al_ is the path of A matrix*/

    MatrixXd fMatrix(1,1);fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_Al_))
        fMatrix=fileUtils_.readMatrix(path_Al_);
    else{
        cout<<"The provided path does not exist: "<<path_Al_<<endl;
        ERROR();
    }

    initialize_Al(fMatrix);
}


void lagsDS::initialize_Ad(const char  *path_Ad_){

    /* Initialize A
     * path_Ad_ is the path of A matrix*/

    MatrixXd fMatrix(1,1);fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_Ad_))
        fMatrix=fileUtils_.readMatrix(path_Ad_);
    else{
        cout<<"The provided path does not exist: "<<path_Ad_<<endl;
        ERROR();
    }

    initialize_Ad(fMatrix);
}


void lagsDS::initialize_local_params(const char  *path_att_l_, const char  *path_w_l_, const char  *path_b_l_, const char  *path_scale_){

    /* Initialize local attractors and activation function parameters
     *	path_att_l_ is the path of the local virtual attractors
     *	path_w_l_ is the path of the local w vectors for hyper-plane functions
     *  path_b_l_ is the path of the params for lambda functions
     *  path_scale_ is the path of the scaling parameter used during estimation
    */


    /* Initializing att_l*/
    MatrixXd fMatrix;
    fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_att_l_))
        fMatrix=fileUtils_.readMatrix(path_att_l_);
    else{
        cout<<"The provided path does not exist."<<endl;
        cout<<"path_att_l "<<endl;
        cout<<path_att_l_<<endl;
        ERROR();
    }
    initialize_att_l(fMatrix);


    /* Initializing w_l*/
    fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_w_l_))
        fMatrix=fileUtils_.readMatrix(path_w_l_);
    else{
        cout<<"The provided path does not exist."<<endl;
        cout<<"path_w_l "<<endl;
        cout<<path_w_l_<<endl;
        ERROR();
    }
    initialize_w_l(fMatrix);


    /* Initializing b_l*/
    fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_b_l_))
        fMatrix=fileUtils_.readMatrix(path_b_l_);
    else{
        cout<<"The provided path does not exist."<<endl;
        cout<<"path_b_l "<<endl;
        cout<<path_b_l_<<endl;
        ERROR();
    }
    initialize_b_l(fMatrix);


    /* Initializing scale*/
    fMatrix.setZero();
    if (fileUtils_.is_file_exist(path_scale_))
        fMatrix=fileUtils_.readMatrix(path_scale_);
    else{
        cout<<"The provided path for scale does not exist:"<< path_scale_ << endl;
        ERROR();
    }
    if ((fMatrix.rows()!=1)){
        cout<<"Initialization of scale is wrong."<<endl;
        ERROR();
    }

    cout<<"** Initializing scale **"<< endl;
    scale_ = fMatrix.coeff(0,0);

}



/****************************************/
/*     Actual computation functions     */
/****************************************/

/***************************************************/
/******** Computations for Global Component ********/
/***************************************************/
MatrixXd lagsDS::compute_Ag(VectorXd X){

    /* Calculating the weighted sum of A matrices */

	if ((X.rows()!=M_)){
		cout<<"The dimension of X in compute_A is wrong."<<endl;
		cout<<"Dimension of states is: "<<M_<<endl;
		cout<<"Dimension of X "<<X.rows()<<endl;
		ERROR();
	}

	MatrixXd A; A.resize(M_,M_);A.setZero();
	if (K_>1)
        gamma_= compute_gamma(X);
	else
		gamma_(K_-1)=1;

	for (int i=0;i<K_;i++)
        A = A + A_g_matrix_[i]*gamma_(i);

	return A;
}

VectorXd lagsDS::compute_gamma(VectorXd X){
    VectorXd gamma;
    gamma.resize(K_);
    gamma.setZero();

	for (int i=0;i<K_;i++)
        gamma(i)=Prior_[i]*GaussianPDF(X,Mu_[i],Sigma_[i]);

    double sum = gamma.sum();
	if (sum<1e-100){
		for (int i=0;i<K_;i++)
            gamma(i)=1.0/K_;
	}
	else
        gamma = gamma/sum;

    return gamma;
}

VectorXd   lagsDS::compute_fg(VectorXd xi, VectorXd att){
    MatrixXd A_matrix; A_matrix.resize(M_,M_); A_matrix.setZero();
    VectorXd xi_dot;     xi_dot.resize(M_);    xi_dot.setZero();

    A_matrix = compute_Ag(xi);
    xi_dot = A_matrix*(xi - att);

    return xi_dot;
}


void lagsDS::set_att_g(VectorXd att_g){
    att_g_ = att_g;

    cout << "Global attractor set:" << endl << att_g_ << endl;
}

MathLib::Vector lagsDS::compute_fg(MathLib::Vector xi, MathLib::Vector att){

    /* Check size of input vectors */

    if (xi.Size() != M_){
        cout<<"The dimension of X in compute_f is wrong."<<endl;
        cout<<"Dimension of states is: "<<M_<<endl;
        cout<<"You provided a vector of size "<< xi.Size()<<endl;
        ERROR();
    }
    if (att.Size() != M_){
        cout<<"The dimension of X in compute_f is wrong."<<endl;
        cout<<"Dimension of states is: "<<M_<<endl;
        cout<<"You provided a vector of size "<< att.Size()<<endl;
        ERROR();
    }

    /* Fill in VectorXd versions of xi and att */
    VectorXd xi_;  xi_.resize(M_);   xi_.setZero();
    VectorXd att_; att_.resize(M_);  att_.setZero();
    for (int m=0;m<M_;m++){
        xi_[m]  = xi[m];
        att_[m] = att[m];
    }

    /* Compute Desired Velocity */
    VectorXd xi_dot_;  xi_dot_.resize(M_);    xi_dot_.setZero();
    xi_dot_ = compute_fg(xi_,att_);

    /* Transform Desired Velocity to MathLib form */
    MathLib::Vector xi_dot; xi_dot.Resize(M_);
    for (int m=0;m<M_;m++)
        xi_dot[m] = xi_dot_[m];

    return xi_dot;
}


/***************************************************/
/******** Computations for Local Component *********/
/***************************************************/

double lagsDS::compute_hk(VectorXd xi, int k){

    /* Calculating each hyper-plane function h_k(x) */
    double h, h_tilde, bias;
    bias    = 1 - w_l_[k].transpose()*att_l_[k];
    h       = w_l_[k].transpose()*xi + bias;
    h_tilde = 0.5 * (h + abs(h));

    return h_tilde;
}


VectorXd lagsDS::compute_grad_hk(VectorXd xi, int k){

    /* Calculating the gradient of each hyper-plane function h_k(x) */
    VectorXd grad_h;  grad_h.resize(M_);  grad_h.setZero();

    double h = compute_hk(xi, k);
    if (h > 0.5)
         grad_h = w_l_[k];

    return grad_h;
}


double lagsDS::compute_lambda_k(VectorXd xi, int k){

    /* Calculating each hyper-plane function h_k(x) */
    double lambda(0.0) ;
    lambda = 1 - compute_rbf(b_l_[k], xi, att_l_[k]);
    return lambda;

}


VectorXd lagsDS::compute_flk(VectorXd xi, int k){

    /* Calculating each locally active component f_l^k(x) */
    VectorXd xi_dot;  xi_dot.resize(M_);  xi_dot.setZero();
    VectorXd att_diff;  att_diff.resize(M_);  att_diff.setZero();

    /* Compute value of local hyper-plane function */
    double hk = compute_hk(xi, k);
    double h_set(1.0), corr_scale(0.0), d_att_g(0.0);
    att_diff = att_l_[k] - att_g_;
    d_att_g = att_diff.norm();

    /* Compute values of h_set and corr_scale */

    if (hk > 1.0)
        hk = 1.0;
    else
        hk = hk*h_set;

    /* No local deflective DS necessary */
    if (d_att_g < 0.1)
        xi_dot = hk*(A_l_matrix_[k])*(xi - att_l_[k]);
    else{
        /* Compute full local Dynamics Components */
        xi_dot = (hk*A_l_matrix_[k] + (1-hk)*A_d_matrix_[k])*(xi - att_l_[k]);

        /* Sum of components + modulation/correction */
        xi_dot = xi_dot - corr_scale*compute_lambda_k(xi,k)*compute_grad_hk(xi,k);
    }

    return xi_dot;

}


VectorXd   lagsDS::compute_fl(VectorXd xi){

    /* Calculating the weighted sum of local components f_l(x) = sum f_l^k(x) */
    VectorXd xi_dot;     xi_dot.resize(M_);    xi_dot.setZero();

    if (K_>1)
        gamma_= compute_gamma(xi);
    else
        gamma_(K_-1)=1;

    for (int k=0;k<K_;k++)
        xi_dot = xi_dot + gamma_(k)*compute_flk(xi, k);

    return xi_dot;
}



double lagsDS::GaussianPDF(VectorXd x, VectorXd Mu, MatrixXd Sigma){

	double p;
	MatrixXd gfDiff;gfDiff.resize(1,M_);
	MatrixXd gfDiff_T;gfDiff_T.resize(M_,1);
	MatrixXd SigmaIIInv;SigmaIIInv.resize(M_,M_);
	double detSigmaII=0;
	MatrixXd gfDiffp;gfDiffp.resize(1,1);gfDiffp.setZero();

	detSigmaII=Sigma.determinant();
	SigmaIIInv=Sigma.inverse();

    if (detSigmaII<0)
        detSigmaII=0;

	gfDiff=(x - Mu).transpose();
	gfDiff_T=x - Mu;
	gfDiffp =gfDiff*SigmaIIInv* gfDiff_T;
	gfDiffp(0,0)=fabs(0.5*gfDiffp(0,0));
	p = exp(-gfDiffp(0,0)) / sqrt(pow(2.0*PI_, M_)*( detSigmaII +1e-50));

    return p;
}

double lagsDS::compute_rbf(double b_r, VectorXd xi, VectorXd center){
    double r;
    VectorXd xiDiff(M_);
    xiDiff = xi - center;
//    r = 1 - exp(-b*xiDiff.norm());
    r = 1 - exp(-b_r* xiDiff.squaredNorm());
    return r;
}


