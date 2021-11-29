#pragma once
#include <SpaceVecAlg/SpaceVecAlg>

namespace lipm_walking
{

enum class ContactState
{
    DoubleSupport,
    SingleSupport,
};

constexpr int DoubleSupportTime = 20;
constexpr int SingleSupportTime = 30;

struct Contact
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Contact() {}
    Contact(const sva::PTransformd & pose) : pose_(pose) 
    {
        makeConstraint();
    }

    const Eigen::Vector3d & position() const
    {
        return pose_.translation();
    }

    void makeConstraint()
    {
        double yaw = pose_.rotation().eulerAngles(2, 1, 0)[0];
        zmp_constraint = Eigen::Matrix<double, 2, 2>::Zero();
        zmp_constraint(0, 0) = 1.;
        zmp_constraint(0, 1) = tan(yaw);
        zmp_constraint(1, 0) = -tan(yaw);
        zmp_constraint(1, 1) = 1.;

        Eigen::Vector2d center = pose_.translation().head(2);
        upper_bound.resize(2);
        lower_bound.resize(2);
        upper_bound[0] = center[0]+tan(yaw)*center[1]+halfLength/abs(cos(yaw));
        upper_bound[1] = center[1]-tan(yaw)*center[0]+halfWidth/abs(cos(yaw));
        lower_bound[0] = center[0]+tan(yaw)*center[1]-halfLength/abs(cos(yaw));
        lower_bound[1] = center[1]-tan(yaw)*center[0]-halfWidth/abs(cos(yaw));
    
        // std::cout << zmp_constraint << '\n' << upper_bound << "\n" << lower_bound << std::endl;
    }

    void getConstraint(Eigen::MatrixXd& linear_matrix, Eigen::VectorXd& low_bound, Eigen::VectorXd& up_bound)
    {
        linear_matrix = zmp_constraint;
        low_bound = lower_bound;
        up_bound = upper_bound;
    }

    bool containPoint(Eigen::Vector2d& point)
    {
        double error_commit = 0.01;
        point =  zmp_constraint * point;
        if(point[0] < lower_bound[0] - error_commit || point[0] > upper_bound[0] + error_commit) {
            return false;
        }
        if(point[1] < lower_bound[1] - error_commit || point[1] > upper_bound[1] + error_commit) {
            return false;
        }
        return true;
    }

    Eigen::Matrix<double, 4, 2> visual() const
    {
        Eigen::Matrix<double, 4, 2> visualMatrix;
        sva::PTransformd p[4];
        p[0].rotation() = Eigen::Matrix3d::Identity();
        p[0].translation() = Eigen::Vector3d(halfLength, halfWidth, 0);
        p[1] = p[0];
        p[1].translation()[1] = -halfWidth;
        p[2] = p[1];
        p[2].translation()[0] = -halfLength;
        p[3] = p[2];
        p[3].translation()[1] = halfWidth;
        for(int i = 0; i < 4; i++)
        {
            p[i] = pose_ * p[i];
            visualMatrix(i, 0) = p[i].translation()[0];
            visualMatrix(i, 1) = p[i].translation()[1];
        }
        return visualMatrix;
    }

public:

    Eigen::Vector3d refVel = Eigen::Vector3d::Zero();
    sva::PTransformd pose_;
    double halfLength = 0.06;
    double halfWidth = 0.04;
    ContactState supportState_ = ContactState::SingleSupport;
    int remainSupportTime_ = SingleSupportTime;
    Eigen::MatrixXd zmp_constraint;
    Eigen::VectorXd upper_bound, lower_bound;
};



}