#include "AdmittanceForceController.h"
#include <iostream>




namespace ForceControl {

AdmittanceForceController::AdmittanceForceController(ContactParams cp, double dt) : cp_(cp), dt_(dt)
{	
	std::cout << "Initilized Admittance Controller" << std::endl;
}

void AdmittanceForceController::test()
{

}

void AdmittanceForceController::update(class KUKAModelKDL& robot, double* q, const Eigen::Vector3d& poseP, const Eigen::Vector4d& poseQ, const Eigen::VectorXd& q_desired, const Eigen::Vector3d& force_current, double* gains, double* update_q)
{
	bool sim = true;
	// Transformation from the EE to the force sensor alignment.
	Eigen::Matrix3d transform_EE_FT;
	if (!sim) 
	{
		// Currently attached force sensor on the real KUKA
    	transform_EE_FT = Eigen::AngleAxisd(0.9325, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX());
    } else 
    {
    	transform_EE_FT = Eigen::Matrix3d::Identity();
    }

    Eigen::Vector3d surface_normal;

    /* transform force to world frame */
	Eigen::Vector3d f_data_transformed = transform_EE_FT.transpose() * force_current;
	Eigen::Vector3d f_test(0,0,1);


	if (IsContact(force_current))
	{
		EstimateSurfaceNormal(force_current, surface_normal);

		// admittance force control law error
		Eigen::Vector3d f_data_normalized = f_data_transformed.normalized();
		Eigen::Vector3d position_err      = f_data_transformed - 5 * f_data_normalized;

		Eigen::Quaterniond orientation_err;
		// orientation error (attempting to keep the orientation normal to the surface)
		// if (f_data_transformed.norm() < 0.1 || f_data_transformed.norm() > 10) 
		// {
		// 	orientation_err = Eigen::Quaterniond::FromTwoVectors(orientation_ref, f_test); 
		// } else 
		// {
		// 	orientation_err = Eigen::Quaterniond::FromTwoVectors(orientation_ref, f_data_transformed); 
		// }

	 	/* --------------------------- integral errror ---------------------------- */
		// integrating the error, make it optional
		error_integral += position_err * dt_;

		// to mitigate integral windup
		if (error_integral.norm() > 1.3 ) 
		{
			error_integral.setZero();
		}
	 	/* ------------------------------------------------------------------------ */

		// concatanate positional and orientation error
		double dt = 0.01;
		Eigen::VectorXd dp(6);
		dp << position_err  * dt, orientation_err.vec() * dt;

		// compute the jacobian inverse (generalized) and make use of the null space
		Eigen::VectorXd dq(7);
		JacobianInvSolve(robot, q, dp, dq);
		

		/* ----------------------------  safety handelling. In case -----------------*/
		for (int i = 0; i < robot.robotParams_.numJoints; i++) 
		{
			if (std::abs(dq[i]) > 0.5) {
				dq.setZero();
			}
		}
		/* ------------------------------------------------------------------------- */

		/* --------------------------update the next pose -------------------------- */
		// admittance force control law
		for (int i = 0; i < robot.robotParams_.numJoints; i++) 
		{
			update_q[i]  = q[i] + dq(i) + dt_ * (q_desired(i) - q[i]);
		}

	} else {
		memcpy(update_q, q, robot.robotParams_.numJoints * sizeof(double));
	}
}

/* estimate surface normal and transform */
void AdmittanceForceController::EstimateSurfaceNormal(const Eigen::VectorXd& force_current, Eigen::Vector3d& surface_normal)
{
	surface_normal = force_current.normalized();
}

bool AdmittanceForceController::IsContact(const Eigen::Vector3d& force)
{
	if (force.norm() < 0.1 || force.norm() > 10) 
	{
		return true;
	}

	return false;
}

void AdmittanceForceController::JacobianInvSolve(class KUKAModelKDL& robot, double* q, const Eigen::VectorXd& dp, Eigen::VectorXd& dq) 
{
	Eigen::MatrixXd Jac;
	robot.getSpatialJacobian(q, Jac); 
	dq = Jac.completeOrthogonalDecomposition().solve(dp);
}

}
