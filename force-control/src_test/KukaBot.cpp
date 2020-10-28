
#include "KukaModel.h"
#include "AdmittanceForceController.h"
#include "SoftContactModel.hpp"
#include "RobotAnalytical.h"

#include <stdio.h>
#include <iostream>
#include <kdl/chain.hpp>

#include <ctime>
#include <ratio>
#include <chrono>
#include <iostream>


using namespace std::chrono;


int main(int argc, char* argv[])
{

	KDL::Chain robot = KDL::KukaDHKdl();

	
	KUKAModelKDLInternalData robotParams;
	robotParams.numJoints = 7;
	robotParams.Kv = Eigen::MatrixXd(7,7);
	robotParams.Kp = Eigen::MatrixXd(7,7);

	// KUKAModelKDL kukaRobot = KUKAModelKDL(robot, robotParams);
	KUKAModelKDL* kukaRobot = new KUKAModelKDL(robot, robotParams);

	RobotAnalyticalInternalData kukaAnalyticalParams;
	kukaAnalyticalParams.numJoints = 7;
	kukaAnalyticalParams.Kv = Eigen::MatrixXd(7,7);
	kukaAnalyticalParams.Kp = Eigen::MatrixXd(7,7);

	RobotAnalytical* kukaAnalytical = new RobotAnalytical(kukaAnalyticalParams);

	kukaAnalytical->initRobot();
	kukaRobot->initRobot();


	Eigen::VectorXd gravityTorque(7);
	Eigen::MatrixXd massMatrix(7, 7);
	Eigen::MatrixXd jacobian(6, 7);

	Eigen::Matrix<double,3,3> poseM;
	Eigen::Vector3d poseP;
	Eigen::Vector3d vel;
	Eigen::Vector3d accel;

	Eigen::VectorXd coriolis(7);
	Eigen::VectorXd force_ext(7);
	Eigen::VectorXd qdd(7);

	double* q = new double[7];
	double* qd = new double[7];
	double* qdd_ = new double[7];
	for (int i = 0; i < 7; i++) {q[i] = 0.5; qd[i] = 0.1;};

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	kukaRobot->getGravityVector(q, gravityTorque);
	kukaRobot->getMassMatrix(q, massMatrix);
	// kukaRobot.getCoriolisMatrix(q, qd, coriolis);
	kukaRobot->getForwardDynamics(q, qd, force_ext, qdd);
	kukaRobot->getSpatialJacobian(q, jacobian);
	// kukaRobot->getSpatialJacobianDot(q, qd, jacobian);
	// kukaRobot->getForwardKinematics(q, qd, qdd_, poseM, poseP, vel, accel, false);
	// kukaAnalytical->getForwardKinematics(q, qd, qdd_, poseM, poseP, vel, accel, false);
	// kukaAnalytical->getSpatialJacobian(q, jacobian);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	// milliseconds time_span = duration_cast<milliseconds>(t2 - t1);
	duration<double, std::micro> time_span = t2 - t1;
	std::cout << "It took me " << time_span.count() << " micros\n";



	/* ------------- Eigen print arguments ------------------- */
	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	/* ------------------------------------------------------- */

	/* Test inverse dynamics */
	std::cout << gravityTorque.format(CleanFmt) << std::endl;
	std::cout << massMatrix.format(CleanFmt) << std::endl;
	// std::cout << coriolis.format(CleanFmt) << std::endl;
	std::cout << qdd.format(CleanFmt) << std::endl;
	std::cout << jacobian.format(CleanFmt) << std::endl;
	std::cout << poseM.format(CleanFmt) << std::endl;
	std::cout << poseP.format(CleanFmt) << std::endl;
	// std::cout << vel.format(CleanFmt) << std::endl;
	// std::cout << accel.format(CleanFmt) << std::endl;


	// Test Admittace Force Controller
	ForceControl::ContactParams cp_FC;
	cp_FC.K = 100;
	cp_FC.D = 10;
	double dt = 0.01;
	bool is_sim = true;

	// double* force_current = new double[3];
	Eigen::Vector3d force_current(0,0,0);
	Eigen::Vector3d force_desired(0,0,5);

	double gains[3] = {1,1,1};
	Eigen::VectorXd q_curr(7);
	q_curr << 0.5,0.5,0.5,0.5,0,0,0;
	Eigen::VectorXd update_q(7);
	Eigen::Vector3d position_ref(0,0,1.1); 
	Eigen::Vector3d orientation_ref(0,0,0);
	Eigen::Vector3d poseQ(0,0,0);


	// std::cout << w << std::endl;
	Eigen::VectorXd q_desired(7);
	q_desired.setZero();

	ForceControl::AdmittanceForceController AdmittanceControl = ForceControl::AdmittanceForceController(cp_FC, dt);
	AdmittanceControl.update(*kukaRobot, q_curr, poseP, poseQ, q_desired, force_current, force_desired, gains, update_q);

	std::cout << update_q.transpose().format(CleanFmt) << std::endl;


	/* Test Soft Contact Model */
	ContactModel::ContactParams cp_;
	cp_.E = 1000;
	cp_.mu = 0.5;
	cp_.nu = 0.4;
	cp_.R  = 0.005;
	cp_.R_path = 1000;
	cp_.Kd = 10;
	ContactModel::SoftContactModel contactModel = ContactModel::SoftContactModel(cp_);

	Eigen::Matrix3d mass_matrix_cart; 
	mass_matrix_cart << 1, 0, 0, 0, 1, 0, 0, 0, 1;

	Eigen::Vector3d position(1, 1, 0); 
	Eigen::Vector3d orientation(1, 1, 0); 
	Eigen::Vector3d velocity(0, 0, 0); 
	Eigen::Vector3d acceleration(1, 1, 0); 
	// Eigen::Vector3d force_current(0, 0, 10); 
	Eigen::Vector3d force_dot;

	contactModel.df(mass_matrix_cart, position, orientation, velocity, acceleration, force_current, force_dot);
	std::cout << force_dot.format(CleanFmt) << std::endl;


	delete[] q;
	delete[] qd;
	delete[] qdd_;

	delete kukaAnalytical;
	delete kukaRobot;

}