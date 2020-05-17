
#ifndef COST_H
#define COST_H

#include "config.h"

#ifdef DEBUG
#include <iostream>
#endif

#include <Eigen/Dense>

class Cost
{
public:
    Cost();

protected:
	 Q;
	stateMat_t Qf;
	commandMat_t R;

	stateVec_t QDiagElementVec;
	stateVec_t QfDiagElementVec;
	commandVec_t RDiagElementVec;
	double pos_scale;
    double vel_scale;
    double pos_f_scale;
    double vel_f_scale;
    double torqoe_scale;
    
	stateVecTab_t cx_new;
	commandVecTab_t cu_new; 
	stateMatTab_t cxx_new; 
	commandR_stateC_tab_t cux_new; 
	commandMatTab_t cuu_new;
	double c_new;

public:
	stateMat_t& getQ();
	stateMat_t& getQf();
	commandMat_t& getR();
	stateVecTab_t& getcx();
	commandVecTab_t& getcu();
	stateMatTab_t& getcxx();
	commandR_stateC_tab_t& getcux();
	commandMatTab_t& getcuu();
	double& getc();

	unsigned int N;

};



template <class Dynamics>
struct CostFunction
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Control   = typename Dynamics::Control;
    using Gradient  = Eigen::Matrix<Scalar, Dynamics::StateSize + Dynamics::ControlSize, 1>;
    using Hessian   = Eigen::Matrix<Scalar, Dynamics::StateSize + Dynamics::ControlSize, Dynamics::StateSize + Dynamics::ControlSize>;

    CostFunction() = default;
    CostFunction(const CostFunction &other) = default;
    CostFunction(CostFunction &&other) = default;
    CostFunction& operator=(const CostFunction &other) = default;
    CostFunction& operator=(CostFunction &&other) = default;
    virtual ~CostFunction() = default;

    CostFunction(const Eigen::Ref<const State> &target)
    : xf(target) {}

    const Eigen::Ref<const State>& target() const { return xf; }
    Eigen::Ref<State> target() { return xf; }
    const Scalar& target(int idx) const { return xf(idx); }
    Scalar& target(int idx) { return xf(idx); }

    virtual Scalar c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) = 0;
    virtual Gradient dc(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) = 0;
    virtual Hessian d2c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) = 0;

    State xf;
};

template <class Dynamics>
struct TerminalCostFunction
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Scalar    = typename Dynamics::Scalar;
    using State     = typename Dynamics::State;
    using Gradient  = Eigen::Matrix<Scalar, Dynamics::StateSize, 1>;
    using Hessian   = Eigen::Matrix<Scalar, Dynamics::StateSize, Dynamics::StateSize>;

    TerminalCostFunction() = default;
    TerminalCostFunction(const TerminalCostFunction &other) = default;
    TerminalCostFunction(TerminalCostFunction &&other) = default;
    TerminalCostFunction& operator=(const TerminalCostFunction &other) = default;
    TerminalCostFunction& operator=(TerminalCostFunction &&other) = default;
    virtual ~TerminalCostFunction() = default;

    TerminalCostFunction(const Eigen::Ref<const State> &target)
    : xf(target) {}

    const Eigen::Ref<const State>& target() const { return xf; }
    Eigen::Ref<State> target() { return xf; }
    const Scalar& target(int idx) const { return xf(idx); }
    Scalar& target(int idx) { return xf(idx); }

    virtual Scalar c(const Eigen::Ref<const State> &x) = 0;
    virtual Gradient dc(const Eigen::Ref<const State> &x) = 0;
    virtual Hessian d2c(const Eigen::Ref<const State> &x) = 0;

    State xf;
};



#endif // COST_H
