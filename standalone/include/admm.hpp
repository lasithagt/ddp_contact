#pragma once

/// @file
///
/// kuka_plan_runner is designed to wait for LCM messages contraining
/// a robot_plan_t message, and then execute the plan on an iiwa arm
/// (also communicating via LCM using the
/// lcmt_iiwa_command/lcmt_iiwa_status messages).
///
/// When a plan is received, it will immediately begin executing that
/// plan on the arm (replacing any plan in progress).
///
/// If a stop message is received, it will immediately discard the
/// current plan and wait until a new plan is received.


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <string>
#include <list>

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector2d;
using Eigen::Vector3d;

#include "config.h"
#include "spline.h"
#include "ilqrsolver_admm.h"
#include "kuka_arm.h"
#include "SoftContactModel.h"
#include "KukaModel.h"
#include "models.h"
#include "cost_function_kuka_arm.h"


using namespace std;
using namespace Eigen;

/* DDP trajectory generation */

// static std::list< const char*> gs_fileName;
// static std::list< std::string > gs_fileName_string;


/* -------------------- Soft_contact_state = 17(14+3) ------------------------*/
class ADMM
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ADMM();

  void ADMM(stateVec_t xinit, stateVec_t xgoal, stateVecTab_t xtrack, const stateVec_t& cList_bar, const stateVec_t& xList_bar, const stateVec_t& uList_bar);

private:
  stateVecTab_t joint_state_traj;
  commandVecTab_t torque_traj;
  stateVecTab_t joint_state_traj_interp;
  commandVecTab_t torque_traj_interp;

protected:
  optimizer::ILQRSolver::traj lastTraj;

};





