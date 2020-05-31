#pragma once

#ifndef KUKAARM_TRK_H
#define KUKAARM_TRK_H

#include "config.h"
#include "cost_function_kuka_arm_track.h"
#include "SoftContactModel.h"
#include "KukaModel.h"

// #include "drake/common/drake_path.h"
// #include "drake/common/drake_assert.h"
// #include "drake/common/find_resource.h"
// #include "drake/common/trajectories/piecewise_polynomial.h"
// #include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
// #include "drake/multibody/joints/floating_base_types.h"
// #include "drake/multibody/parsers/urdf_parser.h"
// #include "drake/multibody/rigid_body_tree.h"
// #include "drake/multibody/rigid_body_plant/rigid_body_plant.h"
// #include "drake/manipulation/util/world_sim_tree_builder.h"
// #include "drake/manipulation/planner/constraint_relaxing_ik.h"

#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

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


class KukaArm_TRK
{
public:
    KukaArm_TRK();
    KukaArm_TRK(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal);
    KukaArm_TRK(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal,ContactModel::SoftContactModel& contact_model,std::vector<Eigen::Matrix<double,6,1> >& iiwa_fk_ref);

    KukaArm_TRK(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal, stateVecTab_t iiwa_xtrack, std::unique_ptr<KUKAModelKDL>& kukaRobot, ContactModel::SoftContactModel& contact_model,
            std::vector<Eigen::Matrix<double,6,1> >& iiwa_fk_ref);


    ~KukaArm_TRK(){};
private:
protected:
    // attributes
    unsigned int stateNb;
    unsigned int commandNb;
    commandVec_t lowerCommandBounds;
    commandVec_t upperCommandBounds;

    stateMat_t fx;
    stateTens_t fxx;
    stateR_commandC_t fu;
    stateR_commandC_commandD_t fuu;
    stateR_stateC_commandD_t fxu;
    stateR_commandC_stateD_t fux;

    stateMatTab_t fxList;
    stateR_commandC_tab_t fuList;
    stateTensTab_t fxxList;
    stateTensTab_t fxuList;
    stateR_commandC_Tens_t fuuList;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    struct timeprofile
    {
        double time_period1, time_period2, time_period3, time_period4;
        unsigned int counter0_, counter1_, counter2_;
    };

private:
    double dt;
    unsigned int N;
    bool initial_phase_flag_;
    struct timeprofile finalTimeProfile;
    struct timeval tbegin_period, tend_period, tbegin_period4, tend_period4; //tbegin_period2, tend_period2, tbegin_period3, tend_period3, 

public:
    static const double mc, mp, l, g;
    
    //######
    unsigned int globalcnt;
    //######
    
    stateVec_t xgoal;
    std::vector<Eigen::Matrix<double,6,1> > fk_ref;
    stateVecTab_t xtrack;

private:
    
    ContactModel::SoftContactModel* contact_model0;
    stateMat_half_t H, C; // inertial, corillois dynamics
    stateVec_half_t G; // gravity? what is this?
    stateR_half_commandC_t Bu; //input mapping
    stateVec_half_t velocity;
    stateVec_half_t accel;
    stateVec_t Xdot_new;
    stateVec_half_t vd;
    stateVecTab_half_t vd_thread;
    stateVecTab_t Xdot_new_thread;

    stateVec_t Xdot1, Xdot2, Xdot3, Xdot4;
    stateMat_t A1, A2, A3, A4, IdentityMat;
    stateR_commandC_t B1, B2, B3, B4;
    stateVec_t Xp, Xp1, Xp2, Xp3, Xp4, Xm, Xm1, Xm2, Xm3, Xm4;
    
    bool debugging_print;
    stateMat_t AA;
    stateR_commandC_t BB;
    stateMatTab_t A_temp;
    stateR_commandC_tab_t B_temp;
    

    std::unique_ptr<KUKAModelKDL> kukaRobot_;

    Eigen::VectorXd q;
    Eigen::VectorXd qd;
    Eigen::Vector3d force_current;
    Eigen::VectorXd qdd;

    std::vector<Eigen::VectorXd> q_thread, qd_thread;
protected:
    // methods
public:

    scalar_t forwardkin_cost(stateVec_t x, commandVec_t u, Eigen::Matrix<double,6,1> fk_goal, 
                                CostFunctionKukaArm_TRK*& costFunction, unsigned int last);
    scalar_t cost_func_expre(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k, 
                                const stateVec_t& xList_bar_k, const commandVec_t& uList_bar_k, CostFunctionKukaArm_TRK*& costFunction);
    stateVec_t finite_diff_cx(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k, const stateVec_t& xList_bar_k, const commandVec_t& uList_bar_k, CostFunctionKukaArm_TRK*& costFunction);
    commandVec_t finite_diff_cu(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k, const stateVec_t& xList_bar_k, const commandVec_t& uList_bar_k, CostFunctionKukaArm_TRK*& costFunction);
    stateVec_t kuka_arm_dynamics(const stateVec_t& X, const commandVec_t& tau);
    
    void kuka_arm_dyn_cst_ilqr(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, const stateVecTab_t& xList_bar, const commandVecTab_t& uList_bar, CostFunctionKukaArm_TRK*& costFunction);
    void kuka_arm_dyn_cst_min_output(const unsigned int& index_k, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr,  const stateVec_t& xList_cur_bar, const commandVec_t& uList_cur_bar, const bool& isUNan, stateVec_t& xList_next, CostFunctionKukaArm_TRK*& costFunction);
    void kuka_arm_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunctionKukaArm_TRK*& costFunction);
    // void kuka_arm_dyn_cst_v3(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList, CostFunctionKukaArm*& costFunction);
    stateVec_t update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    void grad(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    void hessian(const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx_p, stateR_stateC_commandD_t& fxu_p, stateR_commandC_commandD_t& fuu_p);    
    struct timeprofile getFinalTimeProfile();

    unsigned int getStateNb();
    unsigned int getCommandNb();
    commandVec_t& getLowerCommandBounds();
    commandVec_t& getUpperCommandBounds();
    stateMatTab_t& getfxList();
    stateR_commandC_tab_t& getfuList();
private:
protected:
        // accessors //
public:
};


#endif // KUKAARM_H