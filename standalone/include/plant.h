#ifndef PLANT_H
#define PLANT_H

#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "KukaModel.h"
#include "models.h"


template<int S, int C>
class KukaPlant
{
    enum { StateSize = S, ControlSize = C };
    using State             = Eigen::Matrix<double, StateSize, 1>;
    using Scalar            = double;
    using Control           = Eigen::Matrix<double, ControlSize, 1>;
    using StateTrajectory   = Eigen::Matrix<double, StateSize, Eigen::Dynamic>;
    using ControlTrajectory = Eigen::Matrix<double, ControlSize, Eigen::Dynamic>;

public:
    KukaPlant(std::unique_ptr<KUKAModelKDL>& kukaRobot, Scalar dt, Scalar state_var, Scalar control_var)
    : kukaRobot__(kukaRobot), dt_(dt), sdist_(State::Zero(), state_var * StateNoiseVariance::Identity()),
      cdist_(Control::Zero(), control_var * ControlNoiseVariance::Identity()) {}

    KukaPlant() = default;
    KukaPlant(const Plant &other) = default;
    KukaPlant(Plant &&other) = default;
    KukaPlant& operator=(const Plant &other) = default;
    KukaPlant& operator=(Plant &&other) = default;
    ~KukaPlant() = default;


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
    State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u);
private:
    std::unique_ptr<KUKAModelKDL> kukaRobot_;
    Scalar dt_;
    Eigen::EigenMultivariateNormal<double, StateSize> sdist_;
    Eigen::EigenMultivariateNormal<double, ControlSize> cdist_;
};
  
#endif // KUKAARM_H
