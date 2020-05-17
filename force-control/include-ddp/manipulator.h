
#ifndef MANIPULATOR_H
#define MANIPULATOR_H

#include "config.h"

#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>


using namespace Eigen;
using namespace std;


class Manipulator
{
public:
    // TODO: set the command and joint limits if used in the planner
    Manipulator() = default;
    Manipulator(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal);
    Manipulator(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal, std::unique_ptr<RigidBodyTree<double>>& totalTree_);
    ~Manipulator() = default;

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
    virtual void GetJointsInfo(JointStates* q, JointStates* qd);
    virtual void GetMassMatrix(const JointStates& q, MassMatrix* M) = 0;
    virtual void GetCoriolisMatrix(const JointStates& q, const JointStates& qd, CoriolisMatrix* C) = 0;
    virtual void GetGravityVec(const JointStates& q, JointStates* G) = 0;

    virtual State Dynamics(const JointStates& q, const JointStates& qd, State* qdd)
    {
        // TODO: to give qdd using the robot dynamics equation
        // qdd = 
    };

    int GetJoint_N() {return state_n_};
    int GetCommand_N() {return command_n_};

    // stateVec_t Update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);

protected:
    /* variables */
    int state_n_;
    int command_n_;


private:
    double dt;
    unsigned int N;
    bool initial_phase_flag_;
    struct timeprofile finalTimeProfile;
    struct timeval tbegin_period, tend_period, tbegin_period4, tend_period4; //tbegin_period2, tend_period2, tbegin_period3, tend_period3, 

    
    stateMat_half_t H, C; // inertial, corillois dynamics
    stateVec_half_t G; // gravity? what is this?

    Eigen::VectorXd q;
    Eigen::VectorXd qd;
    std::vector<Eigen::VectorXd> q_thread, qd_thread;

};

#endif // MANIPULATOR_H
