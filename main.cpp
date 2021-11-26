#include <OsqpEigen/OsqpEigen.h>
#include <Eigen/Dense>
#include <iostream>

void SimpleSystem()
{
    double dt = 0.1;
    double u_min = -1.0;
    double u_max = 1.0;
    Eigen::Vector3d x0{0, 0, 0};
    const int windows = 8;
    Eigen::Matrix<double, 8, 3> X_ref;
    X_ref << 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1.,
             0., 0., 2.0, 2.0, 2.0, 2.0, 0., 0.,
             0., 20., 0.,  0.,  0.,  0.,-20.,0.; 
    std::cout << "X_ref: \n" << X_ref << std::endl;
    Eigen::Matrix3d A, Q;
    A << 1,   dt, dt*dt,
         0,    1,    dt,
         0,    0,     1;

    Q = Eigen::Matrix3d::Identity();
    double R = 1.0;
    Eigen::Vector3d B;
    B << 0, 0, 1; 

    Eigen::MatrixXd hessian, gradient, constaints, lowerBound, upperBound;
    
    // set hessian
    hessian.resize(3 * (windows + 1) + windows, 3 * (windows + 1) + windows);
    gradient.resize(3 * (windows + 1) + windows, 3 * (windows + 1) + windows);
    for(int i = 1; i <= windows + 1; i++)
    {
        hessian.block<3, 3>(3*(i - 1), 3*(i - 1)) = Q;
        hessian.block<1, 1>(3*windows + i - 1, 3*windows + i - 1) = R;
        gradient.block<3, 1>(3*(i - 1), 0) = Q * X_ref.block<1, 3>()
    }


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
        lowerBound(3 * (windows + 1) + i, 0) = u_min;
        upperBound(3 * (windows + 1) + i, 0) = u_max;
    }
    
    // debug 
    std::cout << "constaints: \n" << constaints << std::endl;
    std::cout << "uppercontaint: \n" << upperBound << std::endl;
    std::cout << "lowercontaint: \n" << lowerBound << std::endl;
    OsqpEigen::Solver solver;
    // settings
    solver.settings()->setWarmStart(true);
    
    // set initial data of the QP solver
    solver.data()->setNumberOfVariables(3 * (windows + 1) + (windows + 1));
    // solver.data()->


}

int main(int argc, char** argv)
{

    SimpleSystem();
    return 0;
}