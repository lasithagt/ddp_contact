#pragma once

#ifndef COSTFUNCTIONKUKAARM_H
#define COSTFUNCTIONKUKAARM_H

#include "config.h"
#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;


class CostFunctionKukaArm
{
	using Jacobian = Eigen::Matrix<double, 1, stateSize + commandSize>;
	using Hessian = Eigen::Matrix<double, 1, stateSize + commandSize>;

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
            dx =  cost_(xu.template head<S>(), xu.template tail<C>());
            return 0;
        }
        /*****************************************************************************/

        using DiffFunc = std::function<Eigen::Matrix<T, S, 1>(const Eigen::Matrix<T, S, 1>&, const Eigen::Matrix<T, C, 1>&)>;
        Differentiable(const DiffFunc &cost) : cost_(cost) {}
        Differentiable() = default;

    private:
        DiffFunc cost_;
    };

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CostFunctionKukaArm();
    scalar_t forwardkin_cost(stateVec_t x, commandVec_t u, Eigen::Matrix<double,6,1> fkgoal, unsigned int last);
    scalar_t cost_func_expre(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k, const stateVec_t& xList_bar_k);

	stateVec_t finite_diff_cx(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k, const stateVec_t& xList_bar_k);
	commandVec_t finite_diff_cu(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k, const stateVec_t& xList_bar_k);

    void computeDerivatives(const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVecTab_t& xList_bar);

	Eigen::Matrix<double,6,6>& getT();
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

protected:
	stateMat_t Q;
	stateMat_t Qf;
	commandMat_t R;
	Eigen::Matrix<double,6,6> T;

	stateVec_t QDiagElementVec;
	stateVec_t QfDiagElementVec;
	commandVec_t RDiagElementVec;
	Eigen::Matrix<double,6,1> TDiagElementVec;

	double pos_scale;
    double vel_scale;
    double pos_f_scale;
    double vel_f_scale;
    double torqoe_scale;
	double fk_pos_scale;
	double fk_orien_scale;
	double force_scale_x;
	double force_scale_y;
    double force_scale_z;
    double force_f_scale_x;
    double force_f_scale_y;
    double force_f_scale_z;

	stateVecTab_t cx_new;
	commandVecTab_t cu_new; 
	stateMatTab_t cxx_new; 
	commandR_stateC_tab_t cux_new; 
	commandMatTab_t cuu_new;
	double c_new;

};

// class TerminalCostFunction
// {
// public:
// 	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// 	TerminalCostFunction() = default;
// 	~TerminalCostFunction() = default;

// };

// class ADMMCostFunction
// {
// public:
// 	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// 	ADMMCost() = default;
// 	~ADMMCost() = default;

// };


#endif // COSTFUNCTIONKUKAARM_H
