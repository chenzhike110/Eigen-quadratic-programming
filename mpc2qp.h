#include <OsqpEigen/OsqpEigen.h>
#include <Eigen/Dense>

class LMPC_QP
{
// convert Model Predictive Control to Quadratic Programming and solve
public:
    LMPC_QP(int statedim, int inputdim, int windows=12);
    LMPC_QP(Eigen::MatrixXd &A, Eigen::MatrixXd &B, int windows=12);
    // set system
    void setDynamicsMatrices(Eigen::MatrixXd &A, Eigen::MatrixXd &B);
    void setInitState(Eigen::VectorXd& x0);
    void setRefState(Eigen::MatrixXd& xref);

    // set qp parameters
    void setWeightMatrices(Eigen::MatrixXd &Q, Eigen::MatrixXd &R, Eigen::MatrixXd &Qn);
    void setStateBound(Eigen::VectorXd& xlowbound, Eigen::VectorXd& xupbound);
    void setInputBound(Eigen::VectorXd& ulowbound, Eigen::VectorXd uupbound);
    // mpc to qp
    void updateHessianAndGradient();
    void updateConstraintAndBound();
    void solve();

private:
    // system parameters
    int dim_x, dim_u; /*state dimension, input dimension*/
    Eigen::MatrixXd A_, B_;
    Eigen::VectorXd x0_;
    // qp parameters
    int windowsize; /*preview window size*/
    Eigen::MatrixXd Q_, Qn_, R_; /*cost matrix, Q for state, Qn for final state, R for input*/
    Eigen::MatrixXd xref_;
    Eigen::VectorXd u_lowbound, u_upbound, x_lowbound, x_upbound;
    // qp transform
    Eigen::SparseMatrix<double> hessian_, constraint_;
    Eigen::VectorXd gradient_, lowerBound_, upperBound_;
}