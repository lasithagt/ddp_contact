
#ifndef OPTIMIZERBASE_H
#define OPTIMIZERBASE_H

#include <numeric>
#include <sys/time.h>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>

namespace optimizer {

/* structure for the trajectory */
struct Trajectory
{
    StateTrajectory states;
    ControlTrajectory controls;
    unsigned int iter;
    double final_cost;
    double final_grad;
    double final_lambda;
    Eigen::VectorXd time_forward, time_backward, time_derivative; //computation time?
};

/* structure for the options */
struct OptionSet {
    int n_hor;
    State xInit;

};



template <class DynamicsT>
class OptimizerBase
{

public:
    OptimizerBase(const DynamicsT& model, CostFunctionKukaArm& iiwaCostFunction, bool fullDDP=0,bool QPBox=0);
    virtual ~OptimizerBase() {
        delete TrajOptions_;
    };
    void FirstInitSolver(stateVec_t& iiwaxInit, stateVec_t& iiwaxDes, commandVecTab_t initialTorque, unsigned int& iiwaN,
                    double& iiwadt, unsigned int& iiwamax_iter, double& iiwatolFun, double& iiwatolGrad);
    // Implemetation required
    virtual void Solve() = 0;
    virtual void InitTraj() = 0;
    virtual void DoForwardPass() = 0;
    virtual void DoBackwardPass() = 0;

    // Implementation optional
    virtual void SetOptionsParameters(const TrajectoryOptionSet& options) {
        TrajOptions_ = new TrajectoryOptionSet(options);  
    };

    virtual Trajectory GetLastSolvedTrajectory() {}; // TODO: don't return 
    virtual bool IsPositiveDefinite(const commandMat_t& Quu) {}; 

private:
    Scalar dt_;
    int H_;
    int iter_;
    util::Logger *logger_;
    bool verbose_;

    util::EigenAlignedVector<Scalar, DynamicsT::ControlSize, DynamicsT::StateSize> Lk_;
    util::EigenAlignedVector<Scalar, DynamicsT::StateSize, DynamicsT::StateSize + DynamicsT::ControlSize> df_;
    util::EigenAlignedVector<Scalar, DynamicsT::StateSize + DynamicsT::ControlSize,
                                     DynamicsT::StateSize + DynamicsT::ControlSize> d2L_;
    util::EigenAlignedVector<Scalar, DynamicsT::StateSize, DynamicsT::StateSize> Vxx_;
    util::EigenAlignedVector<Scalar, DynamicsT::StateSize, DynamicsT::StateSize> VxxT_;
    util::BoxQP<Scalar, DynamicsT::ControlSize> boxqp_;
    StateTrajectory x_;
    ControlTrajectory u_;
    Eigen::Matrix<Scalar, DynamicsT::StateSize, DynamicsT::StateSize> Is_;
    Eigen::Matrix<Scalar, DynamicsT::StateSize + DynamicsT::ControlSize, Eigen::Dynamic> dL_;
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> L_;
    Eigen::Matrix<Scalar, DynamicsT::StateSize, Eigen::Dynamic> Vx_;
    Eigen::Matrix<Scalar, 1, Eigen::Dynamic> V_;
    State qx_;
    Control qu_;
    Eigen::Matrix<Scalar, DynamicsT::ControlSize, DynamicsT::StateSize> qux_;
    Eigen::Matrix<Scalar, DynamicsT::StateSize, DynamicsT::StateSize> qxx_;
    Eigen::Matrix<Scalar, DynamicsT::ControlSize, DynamicsT::ControlSize> quu_;
    Eigen::Matrix<Scalar, DynamicsT::ControlSize, Eigen::Dynamic> lk_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> cost_;
    bool warm_start_;

    // Trajectory options
    OptionSet* Op;
    
    // dynamical model
    DynamicsT dynamic_model_; 

    // cost model
    CostT costFunction_;

    int state_n;
    int control_n;

    stateVec_t xInit; //matrix of <statesize, 1> = essentially a vector
    stateVec_t xgoal;

    unsigned int N;
    unsigned int iter_;
    double dt_;
    commandVecTab_t initCommand;


    commandVecTab_t uListFull;
    commandVec_t u_NAN; //matrix of <commandsize, 1> = essentially a vector
    stateVecTab_t updatedxList;
    commandVecTab_t updateduList;
    stateVecTab_t FList;
    costVecTab_t costList;
    costVecTab_t costListNew;

    stateVecTab_t Vx;
    stateMatTab_t Vxx;

    commandMat_t QuuF;
    commandVec_t k;
    commandR_stateC_t K;
    commandVecTab_t kList;
    commandR_stateC_tab_t KList;
    double alpha;

    stateMat_t lambdaEye;


    /* QP variables */
    //QProblemB* qp;
    bool enableQPBox;
    bool enableFullDDP;


    Eigen::Vector2d dV;
    int newDeriv;
    double g_norm_i, g_norm_max, g_norm_sum;
    bool isUNan;


};

}  // namespace

#endif // OPTIMIZERASE_H

