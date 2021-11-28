#pragma once

#include <Eigen/Dense>
#include <OsqpEigen/OsqpEigen.h>

#include <locomotion/Contact.h>

namespace lipm_walking{

constexpr double WALK_HEIGHT = 0.5;
constexpr double GRAVITY = 9.8; /** Gravity **/
constexpr double SAMPLING_PERIOD = 0.1; /** Duration of delta T **/
constexpr unsigned INPUT_SIZE = 2; /** jeck_x, jeck_y **/
constexpr unsigned NB_STEPS = 20; /** window size for model predictive control preview **/
constexpr unsigned STATE_SIZE = 6; /** CoM positions, velocities and accelerations **/

class ModelPredictiveControl
{
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// Model Predictive Control for ZMP tracking
// Ref ZMP is linespace
public:
    ModelPredictiveControl();
    void updateContact(Contact init, Contact middle, Contact final);
    void phaseDurations(double initSupportDuration, double doubleSupportDuration, double targetSupportDuration, double targetDoubleDuration);
    void computeZMP_ref();

    void updateEstimation(Eigen::VectorXd & xest);
    void updateHessian();
    void updateFixedConstraint();

    void updateGradient();
    void updateBound();

    void buildAndSolve();
    void calibrate(Eigen::MatrixXd states);

public:
    Eigen::Vector2d velWeights = {10., 10.}; /**< Weights of CoM velocity tracking cost */
    double jerkWeight = 0.1; /**< Weight of CoM jerk regularization cost */
    double zmpWeight = 1000.; /**< Weight of reference ZMP tracking cost */

private:
    // Foot Plan
    Contact initContact_;
    Contact nextContact_;
    Contact targetContact_;
    // System transition
    Eigen::MatrixXd A_, B_, C_;

    // Constraints
    Eigen::VectorXd initState_;
    unsigned indexToHrep_[NB_STEPS + 1]; /**< Mapping from timestep index to ZMP inequality constraints */

    // mpc to qp
    Eigen::MatrixXd constraintDense;
    Eigen::SparseMatrix<double> hessian_, constraint_;
    Eigen::VectorXd gradient_, lowerBound_, upperBound_;

    // locomotion
    using RefVector = Eigen::Matrix<double, 2 * (NB_STEPS + 1), 1>; /** x, y **/
    using StateMatrix = Eigen::Matrix<double, STATE_SIZE, 2>;
    RefVector zmpRef_ = RefVector::Zero();
    RefVector comRef_ = RefVector::Zero(); /**< Stacked vector of reference CoM state */

    StateMatrix stateTozmp_ = StateMatrix::Zero(); /**< Linear map to compute the ZMP of a CoM state (position, velocity, acceleration) */

    unsigned nbInitSupportSteps_ = 0; /**< Number of sampling steps in the preview spent in the first single-support phase */
    unsigned nbDoubleSupportSteps_ = 0; /**< Number of discretization steps for the double support phase */
    unsigned nbTargetSupportSteps_ = 0; /**< Number of sampling steps in the preview spent on the second single-support phase */
    unsigned nbNextDoubleSupportSteps_ = 0; /**< Number of sampling steps in the preview spent on a potential second double-support phase */
    unsigned nbFinalSupportSteps_ = 0;
    // qp solver
    OsqpEigen::Solver solver_;
};

}
