
#ifndef ILQR_H
#define ILQR_H

#include <numeric>
#include <sys/time.h>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>

// Move to config file
// #define ENABLE_QPBOX 0
// #define ENABLE_FULLDDP 0

struct ILQGOptionSet : OptionSet {

    double new_cost, cost, dcost, lambda, dlambda, g_norm, expected;
    double **p;           // parallel line search
    const double* alpha; // backtracking coefficients
    int n_alpha;

    double lambda_max, lambda_min, lambda_init, dlambda_init, lambda_factor; // max lambda value, min labda value, lambda init value.


    unsigned int max_iter;
    double tol_grad, tol_fun, tol_constraint; // gradient exit criterion, reduction exit criterion,  

    double z_min; // minimal accepted reduction ratio
    int reg_type;  // 
    int iterations;
    int *log_linesearch;
    double *log_z;
    double *log_cost;
    
    double w_pen_l, w_pen_f, w_pen_max_l, w_pen_max_f, w_pen_init_l, w_pen_init_f, w_pen_fact1, w_pen_fact2;

    ILQGOptionSet() : n_alpha(11), tol_fun(1e-4), tol_constraint(1e-7), tol_grad(1e-4), max_iter(500), lambda_init(1), dlambda_init(1), lambda_min(1e-6), lambda_max(1e10), reg_type(1), z_min(0.0)
    {

    }

    // Implement copy from other 
    ILQGOptionSet(const ILQGOptionSet& other)
    {

    }

    // Eigen::VectorXd alphaList;

};


#define PRINT(x) do { if (DEBUG_ILQR) printf(x);} while (0)


namespace optimizer {
template <class Dynamics>
class ILQR :: public OptimizerBase
{

public:

    ILQR(const DynamicsT& model, CostFunctionKukaArm& iiwaCostFunction, bool fullDDP=0,bool QPBox=0);
    ~ILQR() = default;
    void InitSolver(stateVec_t& iiwaxInit, stateVec_t& iiwaxDes, commandVecTab_t initialTorque, unsigned int& iiwaN,
                    double& iiwadt, unsigned int& iiwamax_iter, double& iiwatolFun, double& iiwatolGrad);
    void Solve() = 0;
    void InitTrajectory() = 0;
    void DoForwardPass() = 0;
    void DoBackwardPass() = 0;

    // Implementation optional
    // void StandardizeParameters(TrajectoryOptSet* options) {};
    Trajectory GetLastSolvedTrajectory() {}; // TODO: don't return 
    bool IsPositiveDefinite(const commandMat_t& Quu) {}; 

};

}  // namespace

#endif // ILQR_H

