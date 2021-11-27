#include <iostream>
#include "mpc2qp.h"

namespace {
    const long MAX = 100000;
}
LMPC_QP::LMPC_QP(int statedim, int inputdim, int windows):
    dim_x(statedim), dim_u(inputdim),
    windowsize(windows), 
{
    A_.resize(dim_x, dim_x);
    B_.resize(dim_x, dim_u);
    x0_.resize(dim_x, 1);
    Q_.resize(dim_x, dim_x);
    Q_.setIdentity();
    R_.resize(dim_u, dim_u);
    R_.setIdentity();
    xref_.resize(windowsize, dim_x);
    u_lowbound.resize(dim_u, 1);
    x_lowbound.resize(dim_x, 1);
    for(int i = 0; i < dim_u; i++) {
        u_lowbound[i] = -MAX;
    }
    for(int i = 0; i < dim_x; i++) {
        x_lowbound[i] = -MAX;
    }
    x_upbound = -x_lowbound;
    u_upbound = -u_lowbound;
}

LMPC_QP::LMPC_QP(Eigen::MatrixXd &A, Eigen::MatrixXd &B, int windows):
    A_(A), B_(B),
    dim_x(A.cols()), dim_u(B.cols()),
    windowsize(windows),
{
    x0_.resize(dim_x, 1);
    Q_.resize(dim_x, dim_x);
    Q_.setIdentity();
    R_.resize(dim_u, dim_u);
    R_.setIdentity();
    xref_.resize(windowsize, dim_x);
    u_lowbound.resize(dim_u, 1);
    x_lowbound.resize(dim_x, 1);
    for(int i = 0; i < dim_u; i++) {
        u_lowbound[i] = -MAX;
    }
    for(int i = 0; i < dim_x; i++) {
        x_lowbound[i] = -MAX;
    }
    x_upbound = -x_lowbound;
    u_upbound = -u_lowbound;
}

void LMPC_QP::setWeightMatrices(Eigen::MatrixXd &Q, Eigen::MatrixXd &R, Eigen::MatrixXd &Qn)
{
    Q_ = Q;
    R_ = R;
    Qn_ = Qn;
}

void LMPC_QP::setDynamicsMatrices(Eigen::MatrixXd &A, Eigen::MatrixXd &B)
{
    A_ = A;
    B_ = B;
}

void LMPC_QP::setInitState(Eigen::VectorXd& x0)
{
    x0_ = x0;
}

void LMPC_QP::setRefState(Eigen::MatrixXd& xref)
{
    if(xref.cols() != dim_x) {
        std::cout << " Ref State dim Error: " << " expected " << dim_x << " got " << xref.cols() << std::endl;
        return;
    }
    if(xref.rows() < windowsize) {
        std::cout << " Ref State Length Error: " << "expected " << windowsize << " got " << xref.rows() << std::endl;
        for(int i = 0; i < xref.rows(); i++) {
            xref_.row(i) = xref.row(i);
        }
        for(int i = xref.rows(); i < windowsize; i++) {
            xref_.row(i) = xref.row(xref.rows()-1);
        }
    }
    else {
        for(int i = 0; i < xref.rows(); i++) {
            xref_.row(i) = xref.row(i);
        }
    }
}

void LMPC_QP::setStateBound(Eigen::VectorXd& xlowbound, Eigen::VectorXd& xupbound)
{
    x_lowbound = xlowbound;
    x_upbound = xupbound;
}

void LMPC_QP::setInputBound(Eigen::VectorXd& ulowbound, Eigen::VectorXd uupbound)
{
    u_lowbound = ulowbound;
    u_upbound = uupbound;
}

void LMPC_QP::updateHessianAndGradient()
{
    Eigen::MatrixXd hessianDense;
    constexpr int variableNum = dim_x * (windowsize + 1) + dim_u * windowsize;
    hessianDense.resize(variableNum, variableNum);
    gradient_.resize(variableNum, 1);
    hessianDense.setZero();
    gradient_.setZero();

    hessianDense.block<dim_x, dim_x>(0, 0) = Q_;
    gradient_.block<dim_x, 1>(0, 0) = -Q_ * xref_.block<1, x_dim>(0, 0).transpose();
    for(int i = 1; i <= windowsize; i++) {
        hessianDense.block<dim_x, dim_x>(dim_x * i, dim_x * i) = Q;
        hessianDense.block<dim_u, dim_u>(dim_x * (windowsize + 1) + i - 1, dim_x * (windowsize + 1) + i - 1) = R_;
        gradient_.block<dim_x, 1>(dim_x * i , 0) = -Q * xref_.block<dim_x, 1>(i, 0).transpose();
    }
    hessian_ = hessianDense.sparseView();
}

void LMPC_QP::updateConstraintAndBound()
{
    Eigen::MatrixXd constraintDense;
    constexpr int variableNum = dim_x * (windowsize + 1) + dim_u * windowsize;
    
    constraintDense.resize(variableNum + dim_x * (windowsize + 1), variableNum);
    lowerBound_.resize(variableNum, 1);
    upperBound_.resize(variableNum, 1);
    constraintDense.setZero();
    lowerBound_.setZero();
    upperBound_.setZero();
    lowerBound.block<dim_x, 1>(0, 0) = -x0_;
    upperBound.block<dim_x, 1>(0, 0) = -x0_; 

    constaints.block<dim_x * (windowsize + 1), dim_x * (windowsize + 1)>(0, 0) = -Eigen::Matrix<double, dim_x * (windowsize + 1), dim_x * (windowsize + 1)>::Identity();
    constaints.block<variableNum, variableNum>(dim_x * (windowsize + 1), 0) = Eigen::Matrix<double, variableNum, variableNum>::Identity();

    for(int i = 1; i <= windowsize; i++) {
        // lowerBound
    }

}

void LMPC_QP::solve()
{

}