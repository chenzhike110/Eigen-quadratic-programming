#include <iomanip>
#include <locomotion/ModelPredictiveControl.h>
#include <locomotion/util.h>

namespace lipm_walking
{

ModelPredictiveControl::ModelPredictiveControl()
{
    // constructe system
    constexpr double T = SAMPLING_PERIOD;
    double zeta = WALK_HEIGHT / GRAVITY;
    double S = T * T / 2.;
    double C = T * T * T / 6.;
    A_.resize(STATE_SIZE, STATE_SIZE);
    B_.resize(STATE_SIZE, INPUT_SIZE);

    A_ << 
        1, 0, T, 0, S, 0,
        0, 1, 0, T, 0, S,
        0, 0, 1, 0, T, 0,
        0, 0, 0, 1, 0, T,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1;
    B_ <<
        C, 0,
        0, C,
        S, 0,
        0, S,
        T, 0,
        0, T;
    
    stateTozmp_ << 
        1, 0,
        0, 1,
        0, 0,
        0, 0,
        -zeta, 0,
        0, -zeta;
    
    initState_ = Eigen::VectorXd::Zero(STATE_SIZE);
    updateFixedConstraint();
    updateHessian();
}

void ModelPredictiveControl::phaseDurations(double initSupportDuration,
                                            double doubleSupportDuration,
                                            double targetSupportDuration,
                                            double targetDoubleDuration)
{
    // initSupport -> doubleSupport -> targetSupport -> nextDoubleSupport 
    const double T = SAMPLING_PERIOD;
    unsigned nbStepSoFar = 0;
    nbInitSupportSteps_ = std::min(static_cast<unsigned>(std::round(initSupportDuration / T)), NB_STEPS - nbStepSoFar);
    nbStepSoFar += nbInitSupportSteps_;
    nbDoubleSupportSteps_ = std::min(static_cast<unsigned>(std::round(doubleSupportDuration / T)), NB_STEPS - nbStepSoFar);
    nbStepSoFar += nbDoubleSupportSteps_;
    nbTargetSupportSteps_ = std::min(static_cast<unsigned>(std::round(targetSupportDuration / T)), NB_STEPS - nbStepSoFar);
    nbStepSoFar += nbTargetSupportSteps_;
    nbNextDoubleSupportSteps_ = std::min(static_cast<unsigned>(std::round(targetDoubleDuration / T)), NB_STEPS - nbStepSoFar);
    nbStepSoFar += nbNextDoubleSupportSteps_;
    nbFinalSupportSteps_ = NB_STEPS - nbStepSoFar;

    for (long i = 0; i <= NB_STEPS; i++) {
        if(i < nbInitSupportSteps_ || (i > 0 && i==nbInitSupportSteps_)) {
            indexToHrep_[i] = 0; // single support on first contact
        }
        else if(i - nbInitSupportSteps_ < nbDoubleSupportSteps_) {
            indexToHrep_[i] = 1; // double support between first contact and second contact
        }
        else if(i - nbInitSupportSteps_ - nbDoubleSupportSteps_ < nbTargetSupportSteps_){
            indexToHrep_[i] = 2; // single support on second contact
        }
        else if(i - nbInitSupportSteps_ - nbDoubleSupportSteps_ - nbTargetSupportSteps_ < nbNextDoubleSupportSteps_) {
            indexToHrep_[i] = 3; // double support between second and third contact
        }
        else {
            indexToHrep_[i] = 4;
        }
    }
}

void ModelPredictiveControl::updateContact(Contact init, Contact middle, Contact final)
{
    initContact_ = init;
    nextContact_ = middle;
    targetContact_ = final;
}

void ModelPredictiveControl::computeZMP_ref()
{
    zmpRef_.setZero();
    Eigen::Vector2d p_0 = initContact_.position().head<2>();
    Eigen::Vector2d p_2 = targetContact_.position().head<2>();
    Eigen::Vector2d p_1 = nextContact_.position().head<2>();

    for(long i = 0; i <= NB_STEPS; i++) {
        if( indexToHrep_[i] <= 1) {
            long j = i - nbInitSupportSteps_;
            double x = (nbDoubleSupportSteps_ > 0) ? static_cast<double>(j) / static_cast<double>(nbDoubleSupportSteps_) : 0;
            x = clamp(x, 0., 1.);
            zmpRef_.segment<2>(2 * i) = (1. - x) * p_0 + x * p_1;
        }
        else {
            long j = i - nbInitSupportSteps_ - nbDoubleSupportSteps_ - nbTargetSupportSteps_;
            double x = (nbNextDoubleSupportSteps_ > 0) ? static_cast<double>(j) / static_cast<double>(nbNextDoubleSupportSteps_) : 0;
            x = clamp(x, 0., 1.);
            zmpRef_.segment<2>(2 * i) = (1. - x) * p_1 + x * p_2;
        }
    }
}

void ModelPredictiveControl::updateEstimation(Eigen::VectorXd & xest)
{
    initState_ = xest;
}

void ModelPredictiveControl::updateHessian()
{
    Eigen::MatrixXd hessianDense;
    constexpr int variableNum = STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE * NB_STEPS;
    hessianDense.resize(variableNum, variableNum);
    hessianDense.setZero();
    Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> Q;
    Eigen::Matrix<double, INPUT_SIZE, INPUT_SIZE> R;
    R.setIdentity();
    // R.setZero();
    Eigen::Matrix2d zmpCost = Eigen::Matrix2d::Identity() * zmpWeight;
    Q = stateTozmp_ * zmpCost * stateTozmp_.transpose();
    R = R * jerkWeight;

    // Q(2, 2) = velWeights[0];
    // Q(3, 3) = velWeights[1];

    hessianDense.block<STATE_SIZE, STATE_SIZE>(0, 0) = Q;
    for(int i = 1; i <= NB_STEPS; i++)
    {
        hessianDense.block<STATE_SIZE, STATE_SIZE>(STATE_SIZE * i, STATE_SIZE * i) = Q;
        hessianDense.block<INPUT_SIZE, INPUT_SIZE>(STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE*(i - 1), STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE*(i - 1)) = R;
    }
    hessian_ = hessianDense.sparseView();

    solver_.data()->setNumberOfVariables(variableNum);
    // std::cout << "hessian: \n" << hessian_ << std::endl;

    if(!solver_.data()->setHessianMatrix(hessian_))
    {
        std::cout << hessian_.cols() << " " << hessian_.rows() << std::endl;
        std::cout << "hessian error" << std::endl;
    }
}

void ModelPredictiveControl::updateGradient()
{
    constexpr int variableNum = STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE * NB_STEPS;
    gradient_.resize(variableNum, 1);
    gradient_.setZero();
    Eigen::Matrix<double, STATE_SIZE, 2> Q_;
    Eigen::Matrix2d zmpCost = Eigen::Matrix2d::Identity() * zmpWeight;
    Q_ = stateTozmp_ * zmpCost;

    gradient_.segment<STATE_SIZE>(0) = -Q_ * zmpRef_.segment<2>(0);
    for(int i = 1; i <= NB_STEPS; i++)
    {
        gradient_.segment<STATE_SIZE>(STATE_SIZE * i) = -Q_ * zmpRef_.segment<2>(2 * i);
    }
    
}

void ModelPredictiveControl::updateFixedConstraint()
{
    constexpr int variableNum = STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE * NB_STEPS;
    constexpr int constraintNum = STATE_SIZE * (NB_STEPS + 1) + NB_STEPS * 2;
    constraintDense.resize(constraintNum, variableNum);
    lowerBound_.resize(constraintNum, 1);
    upperBound_.resize(constraintNum, 1);
    constraintDense.setZero();
    lowerBound_.setZero();
    upperBound_.setZero();

    constraintDense.block<STATE_SIZE * (NB_STEPS + 1), STATE_SIZE * (NB_STEPS + 1)>(0, 0).setIdentity();
    constraintDense = -constraintDense;
    lowerBound_.segment<STATE_SIZE>(0) = -initState_;
    upperBound_.segment<STATE_SIZE>(0) = -initState_; 

    for(int i = 1; i <= NB_STEPS; i++){
        constraintDense.block<STATE_SIZE, STATE_SIZE>(STATE_SIZE * i, STATE_SIZE*(i - 1)) = A_;
        constraintDense.block<STATE_SIZE, INPUT_SIZE>(STATE_SIZE * i, STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE*(i - 1)) = B_;
    }

    solver_.data()->setNumberOfConstraints(constraintNum);
}

void ModelPredictiveControl::updateBound()
{
    int constraintIndex = STATE_SIZE * (NB_STEPS + 1);
    Eigen::VectorXd initLow, initUp, nextLow, nextUp, targetLow, targetUp;
    Eigen::MatrixXd initConstraint, nextConstraint, targetConstraint;
    initContact_.getConstraint(initConstraint, initLow, initUp);
    nextContact_.getConstraint(nextConstraint, nextLow, nextUp);
    targetContact_.getConstraint(targetConstraint, targetLow, targetUp);

    constraintDense.block<NB_STEPS * 2, STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE * NB_STEPS>(constraintIndex, 0).setZero();
    lowerBound_.segment<NB_STEPS * 2>(constraintIndex).setZero();
    upperBound_.segment<NB_STEPS * 2>(constraintIndex).setZero();

    initConstraint = initConstraint * stateTozmp_.transpose();
    nextConstraint = nextConstraint * stateTozmp_.transpose();
    targetConstraint = targetConstraint * stateTozmp_.transpose();

    for(int i = 1; i <= NB_STEPS; i++)
    {
        if( indexToHrep_[i] == 0) {
            constraintDense.block<2, STATE_SIZE>(constraintIndex, STATE_SIZE * i) = initConstraint;
            lowerBound_.segment<2>(constraintIndex) = initLow;
            upperBound_.segment<2>(constraintIndex) = initUp;
            constraintIndex += 2;
        }
        if( indexToHrep_[i] == 2) {
            constraintDense.block<2, STATE_SIZE>(constraintIndex, STATE_SIZE * i) = nextConstraint;
            lowerBound_.segment<2>(constraintIndex) = nextLow;
            upperBound_.segment<2>(constraintIndex) = nextUp;
            constraintIndex += 2;
        }
        if( indexToHrep_[i] == 4) {
            constraintDense.block<2, STATE_SIZE>(constraintIndex, STATE_SIZE * i) = targetConstraint;
            lowerBound_.segment<2>(constraintIndex) = targetLow;
            upperBound_.segment<2>(constraintIndex) = targetUp;
            constraintIndex += 2;
        }
    }
}

void ModelPredictiveControl::buildAndSolve()
{
    computeZMP_ref();
    updateGradient();
    updateBound();

    constraint_ = constraintDense.sparseView();
    
    if(!solver_.data()->setGradient(gradient_))
    {
        std::cout << "gradient error" << std::endl;
    }
    if(!solver_.data()->setLinearConstraintsMatrix(constraint_))
    {
        std::cout << "constaintsSparse error" << std::endl;
    }
    if(!solver_.data()->setLowerBound(lowerBound_))
    {
        std::cout << "lowerBound error" << std::endl;
    }
    if(!solver_.data()->setUpperBound(upperBound_))
    {
        std::cout << "upperBound error" << std::endl;
    }

    if(!solver_.initSolver())
    {
        std::cout << "solver init error" << std::endl;
    }
    if(!solver_.solve())
    {
        std::cout << "solver solve error" << std::endl;
    }

    // std::cout << "gradient: \n" << gradient_ << std::endl;
    // std::cout << "constaints: \n" << constraint_ << std::endl;
    // std::cout << "uppercontaint: \n" << upperBound_ << std::endl;
    // std::cout << "lowercontaint: \n" << lowerBound_ << std::endl;

    Eigen::VectorXd QPSolution = solver_.getSolution();
    Eigen::Map<Eigen::MatrixXd> states(QPSolution.block<STATE_SIZE * (NB_STEPS + 1), 1>(0, 0).data(), STATE_SIZE, NB_STEPS + 1);
    std::cout << "preview states: \n" << states.transpose() << std::endl;
    Eigen::MatrixXd statePlan = states;
    calibrate(statePlan);
}

void ModelPredictiveControl::calibrate(Eigen::MatrixXd states)
{
    Eigen::Matrix<double, 2, 1> zmp;
    for(int i = 0; i < states.cols(); i++)
    {
        zmp = stateTozmp_.transpose() * states.col(i);
        std::cout << zmp.transpose() << " | ";
        std::cout << zmpRef_.segment<2>(i * 2).transpose() << std::endl; 
        if(indexToHrep_[i] == 0) {
            std::cout << "zmp in rectangle: " << initContact_.containPoint(zmp) << std::endl;
        }
        else if(indexToHrep_[i] == 2) {
            std::cout << "zmp in rectangle: " << nextContact_.containPoint(zmp) << std::endl;
        }
        else if(indexToHrep_[i] == 4) {
            std::cout << "zmp in rectangle: " << targetContact_.containPoint(zmp) << std::endl;
        }
    }
}

}