#ifndef KUKAARM_H
#define KUKAARM_H

#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>


template<int S, int C>
class Plant
{
    enum { StateSize = S, ControlSize = C };
    using State             = Eigen::Matrix<double, StateSize, 1>;
    using Control           = Eigen::Matrix<double, ControlSize, 1>;
    using StateTrajectory   = Eigen::Matrix<double, StateSize, Eigen::Dynamic>;
    using ControlTrajectory = Eigen::Matrix<double, ControlSize, Eigen::Dynamic>;

public:
    Plant() = default;
    Plant(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal);
    ~Plant() = default;


    /**
     * @brief   Pass information to the Plant in order to apply controls and obtain the new state.
     *
     * The user must provide an implementation that takes in the system state and the control calculated
     * by the optimizer, applies one control to the system, performs any desired calculations or updates,
     * and returns the new state of the system after the application of the control.
     *
     * Note that unlike the Dynamics functor, this function must return a STATE rather than a state transition.
     * This is in the interests of compatibility with real robotic systems whose state estimators will usually
     * return a state rather than a transition.
     *
     * @param x The state calculated by the optimizer for the current time window.
     * @param u The control calculated by the optimizer for the current time window.
     * @return  The new state of the system.
     */
    virtual State dynamics(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) = 0;


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
    static const double mc, mp, l, g;    
    unsigned int globalcnt;
    stateVec_t xgoal;

private:
    double dt;
    unsigned int N;
    bool initial_phase_flag_;
    struct timeprofile finalTimeProfile;
    struct timeval tbegin_period, tend_period, tbegin_period4, tend_period4; //tbegin_period2, tend_period2, tbegin_period3, tend_period3, 

    
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
    

    Eigen::VectorXd q;
    Eigen::VectorXd qd;
    std::vector<Eigen::VectorXd> q_thread, qd_thread;


};

}  
#endif // KUKAARM_H
