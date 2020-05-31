
#ifndef KUKAARM_H
#define KUKAARM_H

#include "config.h"
#include "cost_function_kuka_arm.h"
#include "SoftContactModel.h"

#include "KukaModel.h"


#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NumericalDiff>

#include <math.h>

#include <memory>
#include <functional>
// #include <mutex>
// std::mutex mtx;

#define pi 3.141592653

#ifndef DEBUG_KUKA_ARM
#define DEBUG_KUKA_ARM 1
#else
    #if PREFIX1(DEBUG_KUKA_ARM)==1
    #define DEBUG_KUKA_ARM 1
    #endif
#endif

#define TRACE_KUKA_ARM(x) do { if (DEBUG_KUKA_ARM) printf(x);} while (0)

using namespace Eigen;
using namespace std;


class KukaArm
{
    using Jacobian = Eigen::Matrix<double, stateSize, stateSize + commandSize>;

    template <class T, int S, int C>
    struct Differentiable
    {
        /*****************************************************************************/
        /*** Replicate Eigen's generic functor implementation to avoid inheritance ***/
        /*** We only use the fixed-size functionality ********************************/
        /*****************************************************************************/
        enum { InputsAtCompileTime = S + C, ValuesAtCompileTime = S };
        using Scalar        = T;
        using InputType     = Eigen::Matrix<T, InputsAtCompileTime, 1>;
        using ValueType     = Eigen::Matrix<T, ValuesAtCompileTime, 1>;
        using JacobianType  = Eigen::Matrix<T, ValuesAtCompileTime, InputsAtCompileTime>;
        int inputs() const { return InputsAtCompileTime; }
        int values() const { return ValuesAtCompileTime; }
        int operator()(const Eigen::Ref<const InputType> &xu, Eigen::Ref<ValueType> dx) const
        {
            dx =  dynamics_(xu.template head<S>(), xu.template tail<C>());
            return 0;
        }
        /*****************************************************************************/

        // using DiffFunc = std::function<Eigen::Matrix<T, S, 1>(const Eigen::Matrix<T, S, 1>&, const Eigen::Matrix<T, C, 1>&)>;
        using DiffFunc = std::function<int(int&)>;
        Differentiable(const DiffFunc &dynamics) : dynamics_(dynamics) {}
        Differentiable() = default;

    private:
        DiffFunc dynamics_;
    };

    /* --------------- to evaluate time profiles ------------------------*/
    struct timeprofile
    {
        double time_period1, time_period2, time_period3, time_period4;
        unsigned int counter0_, counter1_, counter2_;
    };


protected:
    unsigned int stateNb;
    unsigned int commandNb;
    commandVec_t lowerCommandBounds;
    commandVec_t upperCommandBounds;

    stateMat_t fx_;
    stateR_commandC_t fu_;

    stateMatTab_t fxList;
    stateR_commandC_tab_t fuList;

private:
    double dt;
    unsigned int N;
    bool initial_phase_flag_;
    struct timeprofile finalTimeProfile;
    struct timeval tbegin_period, tend_period, tbegin_period4, tend_period4; //tbegin_period2, tend_period2, tbegin_period3, tend_period3, 

public:
    static const double mc, mp, l, g;
    unsigned int globalcnt;
    std::vector<Eigen::Matrix<double,6,1> > fk_ref;

private:
    
    ContactModel::SoftContactModel* contact_model0;
    stateMat_half_t H, C; // inertial, corillois dynamics
    stateVec_half_t G; // gravity? what is this?
    stateR_half_commandC_t Bu; //input mapping
    stateVec_t Xdot_new;
    stateVec_half_t vd;
    stateVecTab_half_t vd_thread;
    stateVecTab_t Xdot_new_thread;

    stateVec_t Xdot1, Xdot2, Xdot3, Xdot4;
    stateMat_t A1, A2, A3, A4, IdentityMat;
    stateR_commandC_t B1, B2, B3, B4;
    stateVec_t Xp, Xp1, Xp2, Xp3, Xp4, Xm, Xm1, Xm2, Xm3, Xm4;
    
    std::unique_ptr<KUKAModelKDL> kukaRobot_;

    Eigen::VectorXd q;
    Eigen::VectorXd qd;
    Eigen::VectorXd qdd;
    Eigen::Vector3d force_current;

    std::vector<Eigen::VectorXd> q_thread, qd_thread;
    bool debugging_print;

    Jacobian j_;
    Differentiable<double, stateSize, commandSize> diff_;
    Eigen::NumericalDiff<Differentiable<double, stateSize, commandSize>, Eigen::Central> num_diff_;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    KukaArm() : diff_([x](int& x) -> int { return x; }),
      num_diff_(diff_) {}

  // diff_([this](const stateVec_t& x, const commandVec_t& u) -> stateVec_t{ return this->kuka_arm_dynamics(x, u); }),
  //     num_diff_(diff_) {}

    ~KukaArm(){};
    KukaArm(double& iiwa_dt, unsigned int& iiwa_N,  std::unique_ptr<KUKAModelKDL>& kukaRobot, ContactModel::SoftContactModel& contact_model, std::vector<Eigen::Matrix<double,6,1> >& iiwa_fk_ref);
    stateVec_t kuka_arm_dynamics(const stateVec_t& X, const commandVec_t& tau);
    void compute_dynamics_jacobian(const stateVecTab_t& xList, const commandVecTab_t& uList);
    void update_fxu(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    struct timeprofile getFinalTimeProfile();

    unsigned int getStateNb();
    unsigned int getCommandNb();
    commandVec_t& getLowerCommandBounds();
    commandVec_t& getUpperCommandBounds();
    stateMatTab_t& getfxList();
    stateR_commandC_tab_t& getfuList();
};



#endif // KUKAARM_H
