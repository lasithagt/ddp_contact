#include "kuka_arm.h"
#include <Eigen/Geometry>


KukaArm::KukaArm(double& iiwa_dt, unsigned int& iiwa_N, std::unique_ptr<KUKAModelKDL>& kukaRobot, ContactModel::SoftContactModel& contact_model) 
: diff_([this](const stateVec_t& x, const commandVec_t& u) -> stateVec_t{ return this->kuka_arm_dynamics(x, u); }),
      num_diff_(diff_) 
{
    //#####
    globalcnt = 0;
    //#####

    if (SOFT_CONTACT)
    {
        stateNb = 20;
        q.resize((stateSize-3)/2);
        qd.resize((stateSize-3)/2);
    }
    else
    {
        stateNb = 14;
        q.resize(stateSize/2);
        qd.resize(stateSize/2);
    }
    commandNb = 7;
    dt = iiwa_dt;
    N = iiwa_N;
    fxList.resize(N);
    fuList.resize(N);

    
    contact_model0 = &contact_model;

    H.setZero();
    C.setZero();
    G.setZero();
    Bu.setZero();
    Xdot_new.setZero();

    
    debugging_print = 0;
    finalTimeProfile.counter0_ = 0;
    finalTimeProfile.counter1_ = 0;
    finalTimeProfile.counter2_ = 0;

    initial_phase_flag_ = 1;
    q.resize(stateSize/2);
    qd.resize(stateSize/2);


    finalTimeProfile.time_period1 = 0;
    finalTimeProfile.time_period2 = 0;
    finalTimeProfile.time_period3 = 0;
    finalTimeProfile.time_period4 = 0;

    if (initial_phase_flag_ == 1)
    {
        kukaRobot_    = std::move(kukaRobot);      
        initial_phase_flag_ = 0;
    }
}

stateVec_t KukaArm::kuka_arm_dynamics(const stateVec_t& X, const commandVec_t& tau)
{
    // struct timeval tbegin_dynamics, tend_dynamics;
    // gettimeofday(&tbegin_dynamics,NULL);
    finalTimeProfile.counter0_ += 1;

    if(finalTimeProfile.counter0_ == 10)
        gettimeofday(&tbegin_period,NULL);

    if (WHOLE_BODY)
    {
        if (SOFT_CONTACT)
        {
            q  = X.head((stateSize-3)/2);
            qd = X.segment((stateSize-3)/2, (stateSize-3)/2);
            force_current = X.tail(3);

            //LW-Test
            Eigen::Matrix<double, (stateSize-3)/2+2,1> q_full;
            Eigen::Matrix<double, (stateSize-3)/2+2,1> qd_full;
            Eigen::Matrix<double, (stateSize-3)/2+2,1> vd_full;
            q_full.setZero();
            qd_full.setZero();
            vd_full.setZero();
            q_full.topRows((stateSize-3)/2) = q;
            qd_full.topRows((stateSize-3)/2) = qd;

            Eigen::Vector3d force_dot;

            // dynamics vector
            vd.setZero();
            force_dot.setZero();

            //LW---------------

            force_current.setZero();

            kukaRobot_->getForwardDynamics(q.data(), qd.data(), tau, qdd);

            Xdot_new << qd, qdd, force_dot;
            // std::cout << Xdot_new << std::endl;
            // Xdot_new.setZero();

            if (finalTimeProfile.counter0_ == 10)
            {
                gettimeofday(&tend_period,NULL);
                finalTimeProfile.time_period1 += (static_cast<double>(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;
            }
            

            if (globalcnt < 40) {
                globalcnt += 1;
            }

        }
       
    }

    return Xdot_new;
}


KukaArm::timeprofile KukaArm::getFinalTimeProfile()
{    
    return finalTimeProfile;
}

void KukaArm::compute_dynamics_jacobian(const stateVecTab_t& xList, const commandVecTab_t& uList)
{

    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    if(debugging_print) TRACE_KUKA_ARM("initialize dimensions\n");
    unsigned int Nl = xList.cols();

    if(debugging_print) TRACE_KUKA_ARM("compute cost function\n");

    for (unsigned int k=0; k < Nl-1; k++) 
    {
        /* Numdiff Eigen */
        num_diff_.df((typename Differentiable<double, stateSize, commandSize>::InputType() << xList.col(k), uList.col(k)).finished(), j_);

        fxList[k] = j_.leftCols(stateSize) * dt + Eigen::Matrix<double, stateSize, stateSize>::Identity();
        fuList[k] = j_.rightCols(commandSize) * dt;
    }

    
    if(debugging_print) TRACE_KUKA_ARM("finish kuka_arm_dyn_cst\n");
}



unsigned int KukaArm::getStateNb()
{
    return stateNb;
}

unsigned int KukaArm::getCommandNb()
{
    return commandNb;
}

commandVec_t& KukaArm::getLowerCommandBounds()
{
    return lowerCommandBounds;
}

commandVec_t& KukaArm::getUpperCommandBounds()
{
    return upperCommandBounds;
}

stateMatTab_t& KukaArm::getfxList()
{
    return fxList;
}

stateR_commandC_tab_t& KukaArm::getfuList()
{
    return fuList;
}

