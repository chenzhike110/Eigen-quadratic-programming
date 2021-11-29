#include <OsqpEigen/OsqpEigen.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <locomotion/ModelPredictiveControl.h>

using namespace lipm_walking;

bool SimpleSystem()
{
    double dt = 0.1;
    double u_min = -1.0;
    double u_max = 1.0;
    Eigen::Vector3d x0{0, 0, 0};
    const int windows = 8;
    
    Eigen::Matrix<double, 3, 9> X_ref;
    X_ref << 0., 0., 0.2, 0.4, 0.6, 0.8, 0.9, 1., 1.,
             0., 0., 2.0, 2.0, 2.0, 2.0, 1.0, 0., 0.,
             0., 20., 0.,  0.,  0.,  0.,-10.,-10., 0.; 
    std::cout << "X_ref: \n" << X_ref << std::endl;
    Eigen::Matrix3d A, Q, R;
    A << 1,   dt, dt*dt,
         0,    1,    dt,
         0,    0,     1;

    Q << 100,   0,  0,
           0,   1,  0,
           0,   0,  0.1;
    R.resize(1, 1);
    R(0, 0) = 1;
    Eigen::Vector3d B;
    B << 0, 0, 1; 

    Eigen::MatrixXd hessian, constaints;
    Eigen::SparseMatrix<double> hessianSparse, constaintsSparse;
    Eigen::VectorXd gradient, lowerBound, upperBound;
    
    // set hessian
    hessian.resize(3 * (windows + 1) + windows, 3 * (windows + 1) + windows);
    gradient.resize(3 * (windows + 1) + windows, 1);
    hessian.setZero();
    gradient.setZero();
    hessian.block<3, 3>(0, 0) = Q;
    gradient.block<3, 1>(0, 0) = Q * X_ref.block<3, 1>(0, 0);
    for(int i = 1; i <= windows; i++)
    {
        hessian.block<3, 3>(3*i, 3*i) = Q;
        hessian(3*(windows + 1) + i - 1, 3*(windows + 1) + i - 1) = R(0, 0);
        gradient.block<3, 1>(3*i , 0) = -Q * X_ref.block<3, 1>(0, i);
    }
    std::cout << "hessian: \n" << hessian << std::endl;
    std::cout << "sparse hessian: \n" << hessian.sparseView() << std::endl;
    std::cout << "gradient: \n" << gradient << std::endl;

    // set containts and bound
    constaints.resize(3 * (windows + 1) + windows, 3 * (windows + 1) + windows);
    lowerBound.resize(3 * (windows + 1) + windows, 1);
    upperBound.resize(3 * (windows + 1) + windows, 1);
    lowerBound.setZero();
    upperBound.setZero();
    lowerBound.block<3, 1>(0, 0) = -x0;
    upperBound.block<3, 1>(0, 0) = -x0; 
    constaints.block<3 * (windows + 1), 3 * (windows + 1)>(0, 0) = -Eigen::Matrix<double, 3 * (windows + 1), 3 * (windows + 1)>::Identity();
    constaints.block<3 * (windows + 1), windows>(0, 3 * (windows + 1)).setZero();
    constaints.block<windows, windows>(3 * (windows + 1), 3 * (windows + 1)) = Eigen::Matrix<double, windows, windows>::Identity();
    for(int i = 1; i <= windows; i++)
    {
        constaints.block<3, 3>(3 * i, 3*(i - 1)) = A;
        constaints.block<3, 1>(3 * i, 3 * (windows + 1) + (i - 1)) = B;
        lowerBound(3 * (windows + 1) + i - 1 , 0) = u_min;
        upperBound(3 * (windows + 1) + i - 1, 0) = u_max;
    }
    
    // debug 
    std::cout << "constaints: \n" << constaints << std::endl;
    std::cout << "uppercontaint: \n" << upperBound << std::endl;
    std::cout << "lowercontaint: \n" << lowerBound << std::endl;

    hessianSparse = hessian.sparseView();
    constaintsSparse = constaints.sparseView();

    OsqpEigen::Solver solver;
    // settings
    // solver.settings()->setWarmStart(true);
    
    // set initial data of the QP solver
    solver.data()->setNumberOfVariables(3 * (windows + 1) + windows);
    solver.data()->setNumberOfConstraints(3 * (windows + 1) + windows);
    if(!solver.data()->setHessianMatrix(hessianSparse))
    {
        std::cout << hessian.cols() << " " << hessian.rows() << std::endl;
        std::cout << "hessian error" << std::endl;
    }
    if(!solver.data()->setGradient(gradient))
    {
        std::cout << "gradient error" << std::endl;
    }
    if(!solver.data()->setLinearConstraintsMatrix(constaintsSparse))
    {
        std::cout << "constaintsSparse error" << std::endl;
    }
    if(!solver.data()->setLowerBound(lowerBound))
    {
        std::cout << "lowerBound error" << std::endl;
    }
    if(!solver.data()->setUpperBound(upperBound))
    {
        std::cout << "upperBound error" << std::endl;
    }

    if(!solver.initSolver())
    {
        std::cout << "solver init error" << std::endl;
    }

    if(!solver.solve())
    {
        std::cout << "solver solve error" << std::endl;
    }

    Eigen::VectorXd QPSolution, ctrl;
    QPSolution = solver.getSolution();
    ctrl = QPSolution.block<windows, 1>(3 * (windows + 1), 0);
    Eigen::Map<Eigen::MatrixXd> states(QPSolution.block<3 * (windows + 1), 1>(0, 0).data(), 3, windows + 1);
    std::cout << "preview states: \n" << states.transpose() << std::endl;
    std::cout << "control: \n" << ctrl << std::endl;

    // calculate loss

    return true;
}

void Locomotion()
{
    Contact init(sva::PTransformd(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero()));
    init.supportState_ = ContactState::DoubleSupport; 
    init.remainSupportTime_ = DoubleSupportTime;
    Contact middle(sva::PTransformd(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0, 0.2, 0)));
    Contact Target(sva::PTransformd(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0, -0.2, 0)));

    Contact end = init;
    end.remainSupportTime_ = 100;
    std::vector<Contact> targetPlan = {init, middle, Target, middle, end};

    ModelPredictiveControl mpc_;
    mpc_.addContactSequence<std::vector<Contact>>(targetPlan);
    // mpc_.phaseDurations(0, 0.4, 0.8, 0.4);
    for (long i = 0; i < 50; i++) {
        mpc_.buildAndSolve();
        Eigen::VectorXd nextState = mpc_.stateTraj_.col(1);
        mpc_.updateEstimation(nextState);
        std::cout << "Traj " << i << " : " <<  nextState.transpose() << std::endl;
    }
    
}

int main(int argc, char** argv)
{
    Locomotion();
    // SimpleSystem();
    return 0;
}