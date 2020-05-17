#ifndef ADMITTANCE_FORCE_CONTROLLER_H
#define ADMITTANCE_FORCE_CONTROLLER_H

#include <Eigen/Dense>
#include "KukaModel.h"


/* Admittance Force Controller */
namespace ForceControl {
struct ContactParams {
	double K;
	double D;
};

class AdmittanceForceController
{
	ContactParams cp_;
	double dt_;
	Eigen::Vector3d error_integral;
	
	public:
		AdmittanceForceController() = default;
		AdmittanceForceController(ContactParams cp, double dt);

		~AdmittanceForceController() = default;

		void test();

		// void update(const class RobotAbstract& robot, const Eigen::Vector3d& poseP);
		void update(class KUKAModelKDL& robot, double* q, const Eigen::Vector3d& poseP, const Eigen::Vector4d& poseQ, const Eigen::VectorXd& q_desired, const Eigen::Vector3d& force_current, double* gains, double* update_q);
		
		void EstimateSurfaceNormal(const Eigen::VectorXd& force_current, Eigen::Vector3d& surface_normal);
		
		bool IsContact(const Eigen::Vector3d& force);
		
		void JacobianInvSolve(class KUKAModelKDL& robot, double* q, const Eigen::VectorXd& dp, Eigen::VectorXd& dq); 
};

}

#endif //ADMITTANCE_FORCE_CONTROLLER_H