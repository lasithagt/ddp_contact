#include "kuka.h"

KUKA::KUKA(){

}

KUKA::KUKA(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal)
{
    //#####
    globalcnt = 0;
    //#####
    stateNb = 14;
    commandNb = 7;


    H.setZero();
    C.setZero();
    G.setZero();
    Bu.setZero();
    velocity.setZero();
    accel.setZero();
    Xdot_new.setZero();
    // Xdot_new_thread.resize(NUMBER_OF_THREAD);
    // vd_thread.resize(NUMBER_OF_THREAD);
}


State KUKA::Dynamics(const stateVec_t& X, const commandVec_t& tau)
{


    if(WHOLE_BODY){
        Eigen::Matrix<double,stateSize/2+2,1> q_full;
        Eigen::Matrix<double,stateSize/2+2,1> qd_full;
        q_full.setZero();
        qd_full.setZero();
        q_full.topRows(stateSize/2)=q;
        qd_full.topRows(stateSize/2)=qd;
        // VectorX<double> q_0 = VectorX<double>::Zero(9);
        // VectorX<double> qd_0 = VectorX<double>::Zero(9);
        // VectorX<double> qd_0 = VectorX<double>::Zero(9);

        //    KinematicsCache<double> cache_ = robot_thread_->doKinematics(q, qd);

        //Yuki
        KinematicsCache<double> cache_ = robot_thread_->CreateKinematicsCache();
        cache_.initialize(q_full,qd_full);
        robot_thread_->doKinematics(cache_, true);

        //const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;
        //gettimeofday(&tbegin_period,NULL);
        MatrixX<double> M_ = robot_thread_->massMatrix(cache_); // Inertial matrix
        //gettimeofday(&tend_period,NULL);
        //finalTimeProfile.time_period2 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

        //gettimeofday(&tbegin_period,NULL);

        drake::WrenchVector<double> tw;
        tw  <<  0,0,0,0,0,100;        
        RigidBody<double> const & rb = (*robot_thread_).get_body(10); 
        drake::eigen_aligned_std_unordered_map< RigidBody<double> const*, drake::TwistVector<double> > f_ext;
        f_ext[&rb] = tw;

        // drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext;
        //gettimeofday(&tend_period,NULL);
        //finalTimeProfile.time_period3 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;
        
        //gettimeofday(&tbegin_period,NULL);
        
        // VectorX<double> bias_term_ = robot_thread_->dynamicsBiasTerm(cache_, f_ext);  // Bias term: M * vd + h = tau + J^T * lambda
        VectorX<double> bias_term_ = robot_thread_->dynamicsBiasTerm(cache_, f_ext);  // Bias term: M * vd + h = tau + J^T * lambda

        vd = (M_.inverse()*(tau - bias_term_)).head(stateSize/2);

        Xdot_new << qd, vd;
        
        
        if (globalcnt < 40)
            globalcnt += 1;

        // vdot is newly calculated using q, qdot, u
        // (qdot, vdot) = f((q, qdot), u) ??? Makes sense??
    }
    else{
        KinematicsCache<double> cache_ = robot_thread_->CreateKinematicsCache();
        cache_.initialize(q,qd);
        robot_thread_->doKinematics(cache_, true);

        //const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;
        //gettimeofday(&tbegin_period,NULL);
        MatrixX<double> M_ = robot_thread_->massMatrix(cache_); // Inertial matrix
        //gettimeofday(&tend_period,NULL);
        //finalTimeProfile.time_period2 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

        //Set false for doing only gravity comp
        //     VectorX<double> gtau = robot_thread_->inverseDynamics(cache_, f_ext, qd_0, false);
        //=============================================
        vd = (M_.inverse()*(tau - bias_term_));
        Xdot_new << qd, vd;

        // vdot is newly calculated using q, qdot, u
        // (qdot, vdot) = f((q, qdot), u) ??? Makes sense??
    }
    return Xdot_new;
}
