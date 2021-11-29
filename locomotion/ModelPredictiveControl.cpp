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
    updateHessian();
    updateFixedConstraint();

    Contact init;
    init.supportState_ = ContactState::DoubleSupport;
    init.remainSupportTime_ = 0;
    footPlan_.push_back(init);

    solver_.settings()->setWarmStart(true);
    // solver_.settings()->setVerbosity(false);

    if(!solver_.initSolver())
    {
        std::cout << "solver init error" << std::endl;
    }
}

void ModelPredictiveControl::updateMeasurement(Eigen::VectorXd & estimate)
{
    for(int i = 0; i < estimate.size() && i < initState_.size(); i++)
    {
        initState_[i] = estimate[i];
    }
}

void ModelPredictiveControl::addContact(Contact & target)
{
    footPlan_.push_back(target);
}

void ModelPredictiveControl::computeZMP_ref()
{   
    totalRef.resize(2 * (totalForward + 1));
    totalRef.setZero();

    Eigen::Vector2d p_0 = footPlan_[0].position().head<2>();
    Eigen::Vector2d p_1 = footPlan_[1].position().head<2>();
    Eigen::Vector2d p_2 = footPlan_[2].position().head<2>();

    for(long i = 0; i <= totalForward; i++) {
        if( totalIndex[i] <= 1) {
            long j = i - nbInitSupportSteps_;
            double x = (nbDoubleSupportSteps_ > 0) ? static_cast<double>(j) / static_cast<double>(nbDoubleSupportSteps_) : 0;
            x = clamp(x, 0., 1.);
            totalRef.segment<2>(2 * i) = (1. - x) * p_0 + x * p_1;
        }
        else {
            long j = i - nbInitSupportSteps_ - nbDoubleSupportSteps_ - nbTargetSupportSteps_;
            double x = (nbNextDoubleSupportSteps_ > 0) ? static_cast<double>(j) / static_cast<double>(nbNextDoubleSupportSteps_) : 0;
            x = clamp(x, 0., 1.);
            totalRef.segment<2>(2 * i) = (1. - x) * p_1 + x * p_2;
        }
    }
}

void ModelPredictiveControl::generateplan()
{
    totalForward = 0;
    int tmpIndex = 0;

    nbInitSupportSteps_ = footPlan_[tmpIndex].supportState_ == ContactState::SingleSupport ? footPlan_[tmpIndex++].remainSupportTime_ : 0;
    totalForward += nbInitSupportSteps_;
    nbDoubleSupportSteps_ = footPlan_[tmpIndex].supportState_ == ContactState::DoubleSupport ? footPlan_[tmpIndex++].remainSupportTime_ : DoubleSupportTime;
    totalForward += nbDoubleSupportSteps_;
    nbTargetSupportSteps_ = footPlan_[tmpIndex].supportState_ == ContactState::SingleSupport ? footPlan_[tmpIndex++].remainSupportTime_ : 0;
    totalForward += nbTargetSupportSteps_;
    nbNextDoubleSupportSteps_ = footPlan_[tmpIndex].supportState_ == ContactState::DoubleSupport ? footPlan_[tmpIndex++].remainSupportTime_ : DoubleSupportTime;
    totalForward += nbNextDoubleSupportSteps_;
    nbFinalSupportSteps_ = footPlan_[tmpIndex].supportState_ == ContactState::SingleSupport ? footPlan_[tmpIndex++].remainSupportTime_ : 0;
    totalForward += nbFinalSupportSteps_;

    totalIndex.resize(totalForward + 1);

    if (totalForward < (NB_STEPS + 1)) {
        std::cout << "【ERROR】ZMP_Ref plan forward size not enough!" << std::endl;
        return;
    }

    for (int i = 0; i <= totalForward; i ++) {
        if (i < nbInitSupportSteps_) {
            totalIndex[i] = 0;
        }
        else if (i < nbInitSupportSteps_ + nbDoubleSupportSteps_)
        {
            totalIndex[i] = 1;
        }
        else if (i < nbInitSupportSteps_ + nbDoubleSupportSteps_ + nbTargetSupportSteps_)
        {
            totalIndex[i] = 2;
        }
        else if (i < nbInitSupportSteps_ + nbDoubleSupportSteps_ + nbTargetSupportSteps_ + nbNextDoubleSupportSteps_)
        {
            totalIndex[i] = 3;
        }
        else {
            totalIndex[i] = 4;
        }
    }
    std::cout << "generate time plan: " << nbInitSupportSteps_ << " " << nbDoubleSupportSteps_
              << " " << nbTargetSupportSteps_ << " " << nbNextDoubleSupportSteps_ << " " << nbFinalSupportSteps_ << std::endl;

    timeCircle_ = 0;
    computeZMP_ref();
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
    gradient_.resize(variableNum, 1);
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
    if(!solver_.data()->setGradient(gradient_)) {
        std::cout << "gradient error" << std::endl;
    }
}

void ModelPredictiveControl::updateGradient()
{ 
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

    for(int i = 1; i <= NB_STEPS; i++){
        constraintDense.block<STATE_SIZE, STATE_SIZE>(STATE_SIZE * i, STATE_SIZE*(i - 1)) = A_;
        constraintDense.block<STATE_SIZE, INPUT_SIZE>(STATE_SIZE * i, STATE_SIZE * (NB_STEPS + 1) + INPUT_SIZE*(i - 1)) = B_;
    }

    constraint_ = constraintDense.sparseView();
    solver_.data()->setNumberOfConstraints(constraintNum);
    if (!solver_.data()->setLinearConstraintsMatrix(constraint_)) {
        std::cout << "set constraint error" << std::endl;
    }
    if (!solver_.data()->setLowerBound(lowerBound_)) {
        std::cout << "set low bound error" << std::endl;
    }
    if (!solver_.data()->setUpperBound(upperBound_)) {
        std::cout << "set up bound error" << std::endl;
    }
}

void ModelPredictiveControl::updateBound()
{
    int constraintIndex = STATE_SIZE * (NB_STEPS + 1);
    Eigen::VectorXd initLow, initUp, nextLow, nextUp, targetLow, targetUp;
    Eigen::MatrixXd initConstraint, nextConstraint, targetConstraint;
    footPlan_[0].getConstraint(initConstraint, initLow, initUp);
    footPlan_[1].getConstraint(nextConstraint, nextLow, nextUp);
    footPlan_[2].getConstraint(targetConstraint, targetLow, targetUp);

    lowerBound_.segment<STATE_SIZE>(0) = -initState_;
    upperBound_.segment<STATE_SIZE>(0) = -initState_; 

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

bool ModelPredictiveControl::checkContact()
{
    if(footPlan_.size() < 3) {
        footPlan_.push_back(footPlan_[footPlan_.size()-1]);
    }
    if(timeCircle_ >= footPlan_[0].remainSupportTime_) {
        // if(footPlan_[0].supportState_ == ContactState::SingleSupport) {
        //     footPlan_[0].supportState_ = ContactState::DoubleSupport;
        //     footPlan_[0].remainSupportTime_ = DoubleSupportTime;
        // }
        // else {
            footPlan_.pop_front();
        // }
        return true;
    }
    return false;
}

void ModelPredictiveControl::buildAndSolve()
{
    if (checkContact()) {
        generateplan();
    }

    zmpRef_.setZero();
    zmpRef_ = totalRef.segment<2 * (NB_STEPS + 1)>(2 * timeCircle_);
    for (long i = 0; i <= NB_STEPS; i++) {
        indexToHrep_[i] = totalIndex[i + timeCircle_];
    }

    updateGradient();
    updateBound();

    constraint_ = constraintDense.sparseView();
    
    if(!solver_.updateGradient(gradient_))
    {
        std::cout << "gradient error" << std::endl;
    }
    if(!solver_.updateLinearConstraintsMatrix(constraint_))
    {
        std::cout << "constaintsSparse error" << std::endl;
    }
    if(!solver_.updateBounds(lowerBound_, upperBound_))
    {
        std::cout << "updateBound error" << std::endl;
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
    stateTraj_ = states;
    calibrate(stateTraj_);
    timeCircle_ ++;
}

void ModelPredictiveControl::calibrate(Eigen::MatrixXd states)
{
    Eigen::Matrix<double, 2, 1> zmp;
    for(int i = 0; i < states.cols(); i++)
    {
        zmp = stateTozmp_.transpose() * states.col(i);
        // std::cout << zmp.transpose() << " | ";
        // std::cout << zmpRef_.segment<2>(i * 2).transpose() << std::endl; 
        if(indexToHrep_[i] % 2 == 0) {
            if(!footPlan_[int(indexToHrep_[i]/2)].containPoint(zmp)) {
                std::cout << "zmp out!" << std::endl;
                break;
            }
        }
    }
}

}