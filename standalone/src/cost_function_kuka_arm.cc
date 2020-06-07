#include "cost_function_kuka_arm.h"


// TODO : incoperate x_track into 
CostFunctionKukaArm::CostFunctionKukaArm(const stateVec_t &x_goal, const stateVecTab_t &x_track) : x_track_(x_track)
{    
    if (SOFT_CONTACT)
    {
        pos_scale = 10;
        vel_scale = 10;
        pos_f_scale = 100;//0.001;
        vel_f_scale = 100;//10;
        torqoe_scale = 1;//100;

        //contact force
        force_scale_x = 0;
        force_scale_y = 0;
        force_scale_z = 0;
        force_f_scale_x = 0;
        force_f_scale_y = 0;
        force_f_scale_z = 0;
        fk_pos_scale = 0;
        fk_orien_scale = 0;

        // torque cost
        // l = sigma(xQx+uRu) + xfQfxf
        QDiagElementVec << pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100,  
                            vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, force_scale_x, force_scale_y, force_scale_z;
        QfDiagElementVec << pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0,
                            vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, force_f_scale_x, force_f_scale_y, force_f_scale_z;
        RDiagElementVec << torqoe_scale*0.005, torqoe_scale*0.005, torqoe_scale*0.007, torqoe_scale*0.007, torqoe_scale*0.02, torqoe_scale*0.02, torqoe_scale*0.05;
        TDiagElementVec << fk_pos_scale*1000, fk_pos_scale*1000, fk_pos_scale*1000, fk_orien_scale*1000, fk_orien_scale*1000, fk_orien_scale*1000;
        
        Q = QDiagElementVec.asDiagonal();
        Qf = QfDiagElementVec.asDiagonal();
        R = RDiagElementVec.asDiagonal();
        T = TDiagElementVec.asDiagonal();
    }
    else
    {
        pos_scale = 10;
        vel_scale = 10;
        pos_f_scale = 100;//0.001;
        vel_f_scale = 100;//10;
        torqoe_scale = 1;//100;
        
        fk_pos_scale = 1;
        fk_orien_scale = 1;
        // initial, final costs (pos ,vel)
        // torque cost
        // l = sigma(xQx+uRu) + xfQfxf

        QDiagElementVec << pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100,  
                            vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10;
        QfDiagElementVec << pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0,
                            vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0;
        RDiagElementVec << torqoe_scale*0.005, torqoe_scale*0.005, torqoe_scale*0.007, torqoe_scale*0.007, torqoe_scale*0.02, torqoe_scale*0.02, torqoe_scale*0.05;
        TDiagElementVec << fk_pos_scale*1000, fk_pos_scale*1000, fk_pos_scale*1000, fk_orien_scale*1000, fk_orien_scale*1000, fk_orien_scale*1000;
        
        Q = QDiagElementVec.asDiagonal();
        Qf = QfDiagElementVec.asDiagonal();
        R = RDiagElementVec.asDiagonal();
        T = TDiagElementVec.asDiagonal();
    }
    
    // TimeHorizon = total time 
    // TimeStep = time between two timesteps
    // N = number of knot
    N = TimeHorizon / TimeStep;
    cx_new.resize(stateSize, N+1);
    cu_new.resize(commandSize, N+1);
    cxx_new.resize(N+1);
    cux_new.resize(N+1);
    cuu_new.resize(N+1);
}


scalar_t CostFunctionKukaArm::forwardkin_cost(stateVec_t x, commandVec_t u, Eigen::Matrix<double,6,1> fkgoal, unsigned int last)
{
    // forward dynamics here xList_curr is a vector of states
    if (SOFT_CONTACT)
    {
        // ---------------------------------------------------------------------------
        // Eigen::VectorXd qq = Eigen::VectorXd((stateSize-3)/2+2);
        // Eigen::VectorXd qqd = Eigen::VectorXd((stateSize-3)/2+2);
        // qq.setZero();
        // qqd.setZero();
        // qq.topRows((stateSize-3)/2) = x.head((stateSize-3)/2);
        // qqd.topRows((stateSize-3)/2) = x.segment((stateSize-3)/2, (stateSize-3)/2);    

        // /*---------------------------------------*/
        // Eigen::Matrix<double,3,3> poseM;
        // Eigen::Vector3d poseP;
        // Eigen::Vector3d vel;
        // Eigen::Vector3d accel;

        // kukaRobot_->getForwardKinematics(qq.data(), qqd.data(), qdd.data(), poseM, poseP, vel, accel, false);

        // // ----------------------------------------------------

        // Eigen::Matrix3d m;
        // m = Eigen::AngleAxisd(poseM);
        // Eigen::Vector3d ea = m.eulerAngles(2, 1, 0); 

        // Eigen::Matrix<double, 6, 1> fk_current;
        // fk_current << poseP(0), poseP(1), poseP(2), ea(0), ea(1), ea(2);

        // // cout << "fk_cur " << fk_cur.transpose() << " fkgoal " << fkgoal.transpose() << endl;
        // // if last element, only add state cost

        // scalar_t c_mat_to_scalar;

        // if (last == 1)
        // {
        //     c_mat_to_scalar = 0.5 * (fk_current.transpose() - fkgoal.transpose()) * costFunction->getT() * (fk_current - fkgoal);
        // }
        // else 
        // {
        //     c_mat_to_scalar = 0.5 * u.transpose() * costFunction->getR() * u;        
        // }

        // std::cout << "in this function" << std::endl;

        scalar_t c_mat_to_scalar;
        // c_mat_to_scalar.setZero();

        return c_mat_to_scalar;
    }

}

scalar_t CostFunctionKukaArm::cost_func_expre(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k)
{
    scalar_t c_mat_to_scalar;
    unsigned int Nl = NumberofKnotPt;

    if (index_k == Nl)
    {
        c_mat_to_scalar = 0.5 * (xList_k.transpose() - x_track_.col(index_k).transpose()) * Qf * (xList_k - x_track_.col(index_k));

    }
    else
    {
        c_mat_to_scalar = 0.5 * (xList_k.transpose() - x_track_.col(index_k).transpose()) * Q * (xList_k - x_track_.col(index_k));
        c_mat_to_scalar += 0.5 * uList_k.transpose() * R * uList_k;
    }
    return c_mat_to_scalar;
}

stateVec_t CostFunctionKukaArm::finite_diff_cx(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k)
{
    stateVec_t cx_fd_k;
    unsigned int n = stateSize;
    // unsigned int m = uList_k.size();
    scalar_t cp1;
    scalar_t cm1;
    double delta = 1e-5;
    stateMat_t Dx;
    Dx.setIdentity();
    Dx = delta * Dx;
    // State perturbation for cost

    for (unsigned int i = 0; i < n; i++)
    {
        cp1 = cost_func_expre(index_k, xList_k+Dx.col(i), uList_k);
        cm1 = cost_func_expre(index_k, xList_k-Dx.col(i), uList_k);
        cx_fd_k(i) = (cp1 - cm1) / (2 * delta);
    }
    return cx_fd_k;
}

commandVec_t CostFunctionKukaArm::finite_diff_cu(const unsigned int& index_k, const stateVec_t& xList_k, const commandVec_t& uList_k)
{
    commandVec_t cu_fd_k;
    // unsigned int n = xList_k.size();
    unsigned int m = commandSize;
    scalar_t cp1;
    scalar_t cm1;
    double delta = 1e-5;
    commandMat_t Du;
    Du.setIdentity();
    Du = delta * Du;
    // State perturbation for cost

    for (unsigned int i=0; i < m; i++)
    {
        cp1 = cost_func_expre(index_k, xList_k, uList_k+Du.col(i));
        cm1 = cost_func_expre(index_k, xList_k, uList_k-Du.col(i));
        cu_fd_k(i) = (cp1 - cm1) / (2 * delta);
    }
    return cu_fd_k;
}

void CostFunctionKukaArm::computeDerivatives(const stateVecTab_t& xList, const commandVecTab_t& uList)
{
    unsigned int Nl = xList.cols();

    stateVec_t cx1;
    stateVec_t cxm1;
    commandVec_t cu1;
    commandVec_t cum1;

    double delta = 1e-5;
    commandMat_t Du;
    Du.setIdentity();
    Du = delta * Du;

    stateMat_t Dx;
    Dx.setIdentity();
    Dx = delta * Dx;

    unsigned int n = stateSize;
    unsigned int m = commandSize;

    for (unsigned int k = 0; k < Nl; k++)
    {
        // if (k < Nl-1)


        // scalar_t c_mat_to_scalar;
        // if (k == Nl - 1) 
        // {
        //     c_mat_to_scalar = cost_func_expre(k, xList.col(k), uList.col(k));
        //     c_new += c_mat_to_scalar(0,0);
        // }
        // else 
        // {
        //     c_mat_to_scalar = cost_func_expre(k, xList.col(k), uList.col(k));
        //     c_new += c_mat_to_scalar(0,0); // TODO: to be checked
        // }

        // cx_new[k] = finite_diff_cx(k, xList[k], uList[k], xList_bar[k]);
        // cu_new[k] = finite_diff_cu(k, xList[k], uList[k], xList_bar[k]);

        // // State perturbation for cxx
        // for (unsigned int i = 0; i < n; i++)
        // {
        //     cx1 = finite_diff_cx(k, xList[k]+Dx.col(i), uList[k], xList_bar[k]);
        //     cxm1 = finite_diff_cx(k, xList[k]-Dx.col(i), uList[k], xList_bar[k]);
        //     cxx_new[k].col(i) = (cx1 - cxm1)/(2*delta);
        // }

        // // Control perturbation for cuu
        // for (unsigned int i = 0; i < m; i++)
        // {
        //     cu1 = finite_diff_cu(k, xList[k], uList[k]+Du.col(i), xList_bar[k]);
        //     cum1 = finite_diff_cu(k, xList[k], uList[k]-Du.col(i), xList_bar[k]);
        //     cuu_new[k].col(i) = (cu1 - cum1)/(2*delta);
        // }

        // Analytical derivatives given quadratic cost
        cx_new.col(k) = Q * xList.col(k);
        cu_new.col(k) = R * uList.col(k);
        cxx_new[k] = Q;

        // costFunction->getcux()[k].setZero();
        cuu_new[k] = R; 

        //Note that cu , cux and cuu at the final time step will never be used (see ilqrsolver::doBackwardPass)
        cux_new[k].setZero();
    } 

    c_new = 0;
    
    // if(debugging_print) TRACE_KUKA_ARM("finish kuka_arm_dyn_cst\n");
}


