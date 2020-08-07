#include <iostream>
#include <memory>
#include <Eigen/Dense>


#include "admm.hpp"

/* ------------- Eigen print arguments ------------------- */
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
 /* ------------------------------------------------------- */

int main(int argc, char *argv[]) 
{
 
  int ADMMiterMax = 5;
  double dt = TimeStep;

  ADMM::ADMMopt ADMMopts(dt, 1e-5, 1e-5, 15, ADMMiterMax);

  ADMM optimizerADMM(ADMMopts);
  stateVec_t xinit, xgoal;
  stateVecTab_t xtrack;
  xtrack.resize(stateSize, NumberofKnotPt + 1);

  xinit << 0, 0.5, 0, 1.0, 0, 0.5, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0;
  xgoal << 1.14, 1.93, -1.48, -1.78, 0.31, 0.13, 1.63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;


  KUKAModelKDLInternalData robotParams;
  robotParams.numJoints = 7;
  robotParams.Kv = Eigen::MatrixXd(7,7);
  robotParams.Kp = Eigen::MatrixXd(7,7);

  /* ---------------------------------- Define the robot and contact model ---------------------------------- */
  KDL::Chain robot = KDL::KukaDHKdl();
  std::shared_ptr<KUKAModelKDL> kukaRobot = std::shared_ptr<KUKAModelKDL>(new KUKAModelKDL(robot, robotParams));


  ContactModel::ContactParams cp_;
  cp_.E = 1000;
  cp_.mu = 0.5;
  cp_.nu = 0.4;
  cp_.R  = 0.005;
  cp_.R_path = 1000;
  cp_.Kd = 10;
  ContactModel::SoftContactModel contactModel(cp_);
  kukaRobot->initRobot();

  // Eigen::VectorXd q_pos_init(7);

  // Eigen::VectorXd gravityTorque(7);
  // kukaRobot->getGravityVector(q_pos_init.data(), gravityTorque);

  // dynamic model of the manipulator and the contact model
  unsigned int N = NumberofKnotPt + 1;
  KukaArm KukaArmModel(dt, N, kukaRobot, contactModel);


  /* ------------------------------------------------------------------------------------------------------ */

  /* ---------------------------------- State and Control Limits ---------------------------------- */
  ADMM::Saturation LIMITS;
  Eigen::VectorXd x_limits_lower(stateSize);
  Eigen::VectorXd x_limits_upper(stateSize);
  Eigen::VectorXd u_limits_lower(commandSize);
  Eigen::VectorXd u_limits_upper(commandSize);
  x_limits_lower << -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -10, -10, -10;    
  x_limits_upper << M_PI, M_PI, M_PI, M_PI, M_PI, M_PI, M_PI, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 10, 10, 10;      
  u_limits_lower << -10, -10, -10, -10, -10, -10, -10;
  u_limits_upper << 10, 10, 10, 10, 10, 10, 10;

  LIMITS.stateLimits.row(0) = x_limits_lower;
  LIMITS.stateLimits.row(1) = x_limits_upper;
  LIMITS.controlLimits.row(0) = u_limits_lower; 
  LIMITS.controlLimits.row(1) = u_limits_upper; 

  /* ----------------------------------------------------------------------------------------------- */

  /* State Tracking. Force tracking */
  // Eigen::MatrixXd F(3, NumberofKnotPt + 1);
  // F.row
  // xtrack.block() = F;
 

  /* Cartesian Tracking */
  ADMM::IKopt IK_OPT(7);
  models::KUKA robotIK = models::KUKA();
  Eigen::MatrixXd Slist(6,7);
  Eigen::MatrixXd M(4,4);
  robotIK.getSlist(&Slist); 
  robotIK.getM(&M);


  Eigen::MatrixXd joint_lims(2,7);
  double eomg = 0.00001;
  double ev   = 0.00001;
  // parameters for ADMM, penelty terms
  Eigen::VectorXd rho(5);
  rho << 1, 0.01, 1, 0, 1;
  
  IK_OPT.joint_limits = joint_lims;
  IK_OPT.ev = ev;
  IK_OPT.eomg = eomg;
  IK_OPT.Slist = Slist;
  IK_OPT.M = M;

  IKTrajectory<IK_FIRST_ORDER> IK_traj = IKTrajectory<IK_FIRST_ORDER>(Slist, M, joint_lims, eomg, ev, rho, N);


  rho << 0, 0, 0, 0, 0;
  Eigen::MatrixXd R(3,3);
  R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  double Tf = 2 * M_PI;

  std::vector<Eigen::MatrixXd> cartesianPoses = IK_traj.generateLissajousTrajectories(R, 0.8, 1, 3, 0.08, 0.08, N, Tf);

  /* initialize xinit, xgoal, xtrack - for the hozizon*/
  optimizerADMM.run(kukaRobot, KukaArmModel, xinit, xgoal, xtrack, cartesianPoses, rho, LIMITS, IK_OPT);

  // TODO : saving data file

	/* TODO : publish to the robot */
  // e.g. ROS, drake, bullet


  return 0;
}
