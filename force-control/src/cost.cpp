#include "cost.h"

	
Cost::Cost()
{    
    pos_scale = 10;
    vel_scale = 10;
    pos_f_scale = 100;//0.001;
    vel_f_scale = 100;//10;
    torqoe_scale = 1;//100;

    // Move all these to a final main
    QDiagElementVec << pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100, pos_scale*100,  
                        vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10, vel_scale*10;
    QfDiagElementVec << pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0, pos_f_scale*1000.0,
                        vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0, vel_f_scale*100.0;
    RDiagElementVec << torqoe_scale*0.005, torqoe_scale*0.005, torqoe_scale*0.007, torqoe_scale*0.007, torqoe_scale*0.02, torqoe_scale*0.02, torqoe_scale*0.05;

    Q = QDiagElementVec.asDiagonal();
    Qf = QfDiagElementVec.asDiagonal();
    R = RDiagElementVec.asDiagonal();

    // TimeHorizon = total time 
    // TimeStep = time between two timesteps
    // N = number of knot
    N = TimeHorizon/TimeStep;
    cx_new.resize(N+1);
    cu_new.resize(N+1);
    cxx_new.resize(N+1);
    cux_new.resize(N+1);
    cuu_new.resize(N+1);
}

stateMat_t& CostFunctionKukaArm::getQ()
{
    return Q;
}

stateMat_t& CostFunctionKukaArm::getQf()
{
    return Qf;
}

commandMat_t& CostFunctionKukaArm::getR()
{
    return R;
}

stateVecTab_t& CostFunctionKukaArm::getcx()
{
    return cx_new;
}

commandVecTab_t& CostFunctionKukaArm::getcu()
{
    return cu_new;
}

stateMatTab_t& CostFunctionKukaArm::getcxx()
{
    return cxx_new;
}

commandR_stateC_tab_t& CostFunctionKukaArm::getcux()
{
    return cux_new;
}

commandMatTab_t& CostFunctionKukaArm::getcuu()
{
    return cuu_new;
}

double& CostFunctionKukaArm::getc()
{
    return c_new;
}

