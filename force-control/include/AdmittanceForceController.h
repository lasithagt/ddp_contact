#ifndef ADMITTANCE_FORCE_CONTROLLER_H
#define ADMITTANCE_FORCE_CONTROLLER_H

#include <Eigen/Dense>
#include <iostream>
#include <Eigen/QR>

/* Admittance Force Controller */
namespace ForceControl {

struct ContactParams {
	double K;
	double D;
};

template<class Robot>
class AdmittanceForceController
{
	ContactParams cp_;
	double dt_;
	Eigen::Vector3d error_integral;
	public:
		AdmittanceForceController() = default;
		AdmittanceForceController(Robot &robot_, ContactParams cp, double dt);
		~AdmittanceForceController() = default;
		void update(const Eigen::VectorXd& q, const Eigen::Vector3d& poseP, const Eigen::Vector3d& poseQ, const Eigen::VectorXd& q_desired, const Eigen::Vector3d& force_current, const Eigen::Vector3d& force_desired, double* gains, Eigen::VectorXd& update_q);
		void EstimateSurfaceNormal(const Eigen::VectorXd& force_current, Eigen::Vector3d& surface_normal);
		bool IsContact(const Eigen::Vector3d& force);
		inline Eigen::VectorXd JacobianInvSolve(const Eigen::VectorXd& q, const Eigen::VectorXd& dp); 
	private:
		Robot &robot;
};


/* Methods */
template<class Robot>
AdmittanceForceController<Robot>::AdmittanceForceController(Robot &robot_, ContactParams cp, double dt) : robot(robot_), cp_(cp), dt_(dt)
{	
	std::cout << "Initilized Admittance Controller" << std::endl;
}


template<class Robot>
void AdmittanceForceController<Robot>::update(const Eigen::VectorXd& q, const Eigen::Vector3d& poseP, const Eigen::Vector3d& poseQ, const Eigen::VectorXd& q_desired, const Eigen::Vector3d& force_current, const Eigen::Vector3d& force_desired, double* gains, Eigen::VectorXd& update_q)
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
    /* transform force to world frame 
		Note: Bullet gives it in the world frame. Real KUKA has a force sensor and given in local frame.
    */
	Eigen::Vector3d f_data_transformed = transform_EE_FT.transpose() * force_current;
	Eigen::Vector3d f_test(0,0,1);


	if (IsContact(force_current))
	{
		Eigen::Vector3d surface_normal;
		EstimateSurfaceNormal(force_current, surface_normal);

		// admittance force control law error
		Eigen::Vector3d f_data_normalized = f_data_transformed.normalized();
		Eigen::Vector3d position_err      = force_desired.norm() * f_data_normalized - f_data_transformed;

		Eigen::Quaterniond orientation_err;

		// orientation error (attempting to keep the orientation normal to the surface)
		if (f_data_transformed.norm() < 0.1 || f_data_transformed.norm() > 10) 
		{
			Eigen::Vector3d norm_surface(0,0,1);
			orientation_err = Eigen::Quaterniond::FromTwoVectors(poseQ, norm_surface); 
		} else 
		{
			orientation_err = Eigen::Quaterniond::FromTwoVectors(poseQ, f_data_transformed); 
		}

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
		Eigen::VectorXd dp(6);
		auto euler_err = orientation_err.toRotationMatrix().eulerAngles(0, 1, 2);

		dp << position_err  * dt_, euler_err * dt_;

		// compute the jacobian inverse (generalized) and make use of the null space

		Eigen::VectorXd dq(7);
		dq = JacobianInvSolve(q, dp);


		/* ----------------------------  safety handelling. In case -----------------*/
		// for (int i = 0; i < robot.robotParams_.numJoints; i++) 
		// {
		// 	if (std::abs(dq[i]) > 0.5) {
		// 		dq.setZero();
		// 	}
		// }
		/* ------------------------------------------------------------------------- */

		/* --------------------------update the next pose -------------------------- */
		// admittance force control law

		update_q  = q + dq + dt_ * (q_desired - q);

	} else {
		update_q = q;
	}
}

/* estimate surface normal and transform */
template<class Robot>
void AdmittanceForceController<Robot>::EstimateSurfaceNormal(const Eigen::VectorXd& force_current, Eigen::Vector3d& surface_normal)
{
	surface_normal = force_current.normalized();
}

template<class Robot>
bool AdmittanceForceController<Robot>::IsContact(const Eigen::Vector3d& force)
{
	if (force.norm() < 0.1 || force.norm() > 10) 
	{
		return false;
	}

	return true;
}

template<class Robot>
inline Eigen::VectorXd AdmittanceForceController<Robot>::JacobianInvSolve(const Eigen::VectorXd& q, const Eigen::VectorXd& dp) 
{
	Eigen::MatrixXd Jac(6, 7);
	Eigen::VectorXd dq(7);
	robot.getSpatialJacobian(const_cast<double*> (q.data()), Jac); 
	// Jac.completeOrthogonalDecomposition().setThreshold(0);
	dq = Jac.completeOrthogonalDecomposition().solve(dp);

	return dq;
}


}

#endif //ADMITTANCE_FORCE_CONTROLLER_H