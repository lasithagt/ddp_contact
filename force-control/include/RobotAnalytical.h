#ifndef ROBOT_ANALYTICAL_H
#define ROBOT_ANALYTICAL_H

#include "RobotAbstract.h"
#include "KUKAAnalyticalSolutions.h"

#include <Eigen/Dense>


struct RobotParams 
{	
	double params[84] = {0,0,0,0,0,0,0,0,0,0,0,0, \
							0,0,0,0,0,0,0,0,0,0,0,0, \
								0,0,0,0,0,0,0,0,0,0,0,0, \
									0,0,0,0,0,0,0,0,0,0,0,0, \
										0,0,0,0,0,0,0,0,0,0,0,0, \
											0,0,0,0,0,0,0,0,0,0,0,0, \
												0,0,0,0,0,0,0,0,0,0,0,0};

	// L_1xx, L_1xy, L_1xz, L_1yy, L_1yz, L_1zz, l_1x, l_1y, l_1z, m_1, fv_1, fc_1, 
	// L_2xx, L_2xy, L_2xz, L_2yy, L_2yz, L_2zz, l_2x, l_2y, l_2z, m_2, fv_2, fc_2, 
	// L_3xx, L_3xy, L_3xz, L_3yy, L_3yz, L_3zz, l_3x, l_3y, l_3z, m_3, fv_3, fc_3, 
	// L_4xx, L_4xy, L_4xz, L_4yy, L_4yz, L_4zz, l_4x, l_4y, l_4z, m_4, fv_4, fc_4, 
	// L_5xx, L_5xy, L_5xz, L_5yy, L_5yz, L_5zz, l_5x, l_5y, l_5z, m_5, fv_5, fc_5, 
	// L_6xx, L_6xy, L_6xz, L_6yy, L_6yz, L_6zz, l_6x, l_6y, l_6z, m_6, fv_6, fc_6, 
	// L_7xx, L_7xy, L_7xz, L_7yy, L_7yz, L_7zz, l_7x, l_7y, l_7z, m_7, fv_7, fc_7]

	
};

/* structure for the active robot */
struct RobotAnalyticalInternalData : RobotAbstractInternalData
{
    int numJoints;
    Eigen::MatrixXd Kv; // joint dynamic coefficient
    Eigen::MatrixXd Kp;


};

/* abstract class for a serial robot */
class RobotAnalytical : public RobotAbstract
{	
public:

	using JointState = Eigen::Matrix<double, 7, 1>;
	using MassMatrix = Eigen::Matrix<double, 7, 7>;
	using Jacobian   = Eigen::Matrix<double, 6, 7>;
	using Frame      = Eigen::Matrix<double, 4, 4>;

	struct RobotAnalyticalInternalData* m_data;

	RobotAnalytical() = default;

	RobotAnalytical(const RobotAnalyticalInternalData& robotParams);

	~RobotAnalytical();

	int initRobot(); 

	void getForwardKinematics(double* q, double* qd, double *qdd, Eigen::Matrix<double,3,3>& poseM, Eigen::Vector3d& poseP, Eigen::Vector3d& vel, Eigen::Vector3d& accel, bool computeOther);


	/* given q, qdot, qddot, outputs torque output*/
	void getInverseDynamics(double* q, double* qd, double* qdd, Eigen::VectorXd& torque);

	void getForwardDynamics(double* q, double* qd, const Eigen::VectorXd& force_ext, Eigen::VectorXd& qdd);

	virtual void getForceTorque() {};

	void getMassMatrix(double* jointPositions, Eigen::MatrixXd& massMatrix);

	void getCoriolisMatrix(double* q, double* qd, Eigen::VectorXd& coriolis);

	void getGravityVector(double* q, Eigen::VectorXd& gravityTorque);

	void getSpatialJacobian(double* q, Eigen::MatrixXd& jacobian);

	void getSpatialJacobianDot(double* q, double* qd, Eigen::MatrixXd& jacobianDot);

	void ik();
private:
    // patch variables for speed
    JointState q_;   
    JointState qd_;   
    JointState qdd_;
    MassMatrix inertia_mat_;     // Interia Matrix
    JointState coriolis_;        // CoriolisVector
    JointState gravity_;         // GravityVector
    Jacobian jacobian_;
    Frame frame_;
    Frame frame_vel_;

    Eigen::Matrix<double, 7, 7> Kv_;
    RobotAnalyticalInternalData robotParams_;
    RobotParams intertialParams_;
    KUKAAnalyticalSolutions* KukaAnalytical_;
    double* FK;;

};
#endif  //ROBOT_ANALYTICAL_H
