#pragma once
#include <SpaceVecAlg/SpaceVecAlg>

namespace lipm_walking
{

enum class ContactState
{
    DoubleSupport,
    SingleSupport,
};

enum SupportTime
{
    DoubleSupport = 4,
    SingleSupport = 8,
};

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
        zmp_constraint = Eigen::Matrix<double, 2, 6>::Zero();
        zmp_constraint(0, 0) = 1.;
        zmp_constraint(0, 1) = tan(yaw);
        zmp_constraint(1, 0) = -tan(yaw);
        zmp_constraint(1, 1) = 1.;

        Eigen::Vector2d center = pose_.translation().head(2);
        upper_bound = Eigen::Vector2d(center[0]+tan(yaw)*center[1]+halfLength/abs(cos(yaw)),
                                      center[1]-tan(yaw)*center[0]+halfWidth/abs(cos(yaw)));
        lower_bound = Eigen::Vector2d(center[0]+tan(yaw)*center[1]-halfLength/abs(cos(yaw))
                                      center[1]-tan(yaw)*center[0]-halfWidth/abs(cos(yaw)));
    }

    void getConstraint(Eigen::MatrixXd& linear_matrix, Eigen::VectorXd& low_bound, Eigen::VectorXd& up_bound)
    {
        linear_matrix = zmp_constraint;
        low_bound = lower_bound;
        up_bound = upper_bound;
    }

public:
    Eigen::Vector3d refVel = Eigen::Vector3d::Zero();
    sva::PTransformd pose_;
    double halfLength = 0.06;
    double halfWidth = 0.04;
    ContactState supportState_ = ContactState::SingleSupport;
    SupportTime remainSupportTime_ = SupportTime::DoubleSupport;
    Eigen::MatrixXd zmp_constraint;
    Eigen::VectorXd upper_bound, lower_bound;
};



}