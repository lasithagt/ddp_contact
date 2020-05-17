#pragma once

#ifndef COSTFUNCTIONKUKAARM_H
#define COSTFUNCTIONKUKAARM_H

#include "config.h"
#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;


class CostFunctionKukaArm
{
public:
    CostFunctionKukaArm();
private:
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
    // attributes
public:
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
private:

protected:
    // methods
public:
private:
protected:
    // accessors
public:

};


#endif // COSTFUNCTIONKUKAARM_H
