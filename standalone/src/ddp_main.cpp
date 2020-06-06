#include <iostream>
#include <memory>
#include <Eigen/Dense>

#include "config.h"
#include "spline.h"
#include "ilqrsolver.h"
#include "kuka_arm.h"
#include "SoftContactModel.h"
#include "KukaModel.h"
#include "models.h"

#include "ddp.h"

/* ------------- Eigen print arguments ------------------- */
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
 /* ------------------------------------------------------- */

int main(int argc, char *argv[]) 
{
 
  DDP optimizer;
  stateVec_t xinit, xgoal;
  stateVecTab_t xtrack;
  xtrack.resize(stateSize, NumberofKnotPt + 1);

  xinit << 0, 0.5, 0, 1.0, 0, 0.5, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0;
  xgoal << 1.14, 1.93, -1.48, -1.78, 0.31, 0.13, 1.63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // The trajectory of contact force to be tracked
  xtrack.setZero();
 

  KUKAModelKDLInternalData robotParams;
  robotParams.numJoints = 7;
  robotParams.Kv = Eigen::MatrixXd(7,7);
  robotParams.Kp = Eigen::MatrixXd(7,7);

  /* initialize xinit, xgoal, xtrack - for the hozizon*/
  KDL::Chain robot = KDL::KukaDHKdl();
  std::unique_ptr<KUKAModelKDL> kukaRobot = std::unique_ptr<KUKAModelKDL>(new KUKAModelKDL(robot, robotParams));

  /* initialize xinit, xgoal, xtrack - for the hozizon*/
	optimizer.Run(xinit, xgoal, xtrack);



  // TODO : saving data file
  // for (unsigned int i = 0; i < N; i++)
  // {
  //   saveVector(joint_state_traj[i], "joint_trajectory");
  //   saveVector(torque_traj[i], "joint_torque_command");
  // }

  // saveVector(lastTraj.xList[N], "joint_trajectory");

  // for (unsigned int i=0; i<=N*InterpolationScale;i++){
  //   saveVector(joint_state_traj_interp[i], "joint_trajectory_interpolated");
  // }

	/* TODO : publish to the robot */
  // e.g. ROS, drake, bullet


  return 0;
}
