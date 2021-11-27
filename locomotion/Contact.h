#pragma once
#include <SpaceVecAlg/SpaceVecAlg>

namespace lipm_walking
{

enum class ContactState
{
  DoubleSupport,
  SingleSupport,
};

struct Contact
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Contact() {}
    Contact(const sva::PTransformd & pose) : pose_(pose) {}

    const Eigen::Vector3d & position() const
    {
        return pose_.translation();
    }

public:
    Eigen::Vector3d refVel = Eigen::Vector3d::Zero();
    sva::PTransformd pose_;
    double halfLength = 0.06;
    double halfWidth = 0.04;
    ContactState supportState_ = ContactState::SingleSupport;
    Eigen::MatrixXd zmp_constraint;
    Eigen::VectorXd upper_bound, lower_bound;
};



}