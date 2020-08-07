
#include <iostream>
#include <memory>

/* ADMM trajectory generation */
#include "admm.hpp"

using namespace Eigen;

ADMM::ADMM(const ADMMopt& ADMM_opt) : ADMM_OPTS(ADMM_opt)
{
    /* Initalize Primal and Dual variables */
    N = NumberofKnotPt;

    // primal parameters
    xnew.resize(stateSize, N + 1);
    unew.resize(commandSize, N);

    xbar.resize(stateSize, N + 1);
    cbar.resize(stateSize, N + 1);
    ubar.resize(commandSize, N);

    xbar_old.resize(stateSize, N + 1); // "old" for last ADMM iteration 
    cbar_old.resize(stateSize, N + 1);
    ubar_old.resize(commandSize, N); 
    
    // dual parameters
    x_lambda.resize(stateSize, N + 1);
    c_lambda.resize(stateSize, N + 1);
    u_lambda.resize(commandSize, N);

    x_temp.resize(stateSize, N + 1);
    c_temp.resize(stateSize, N + 1);
    u_temp.resize(commandSize, N);

    x_temp2.resize(stateSize, N + 1);
    c_temp2.resize(stateSize, N + 1);
    u_temp2.resize(commandSize, N);

    u_0.resize(commandSize, N);

    xubar.resize(stateSize + commandSize, N); // for projection

    // primal residual
    res_x.resize(ADMM_opt.ADMMiterMax, 0);
    res_u.resize(ADMM_opt.ADMMiterMax, 0);
    res_c.resize(ADMM_opt.ADMMiterMax, 0);

    // dual residual
    res_xlambda.resize(ADMM_opt.ADMMiterMax, 0);
    res_ulambda.resize(ADMM_opt.ADMMiterMax, 0);
    res_clambda.resize(ADMM_opt.ADMMiterMax, 0);

    final_cost.resize(ADMM_opt.ADMMiterMax + 1, 0);

    // joint_positions_IK
    joint_positions_IK.resize(7, N + 1);
    thetalistd0.resize(7);
}
// template<class T>
void ADMM::run(std::shared_ptr<KUKAModelKDL>& kukaRobot, KukaArm& KukaArmModel, const stateVec_t& xinit, const stateVec_t& xgoal,
  const stateVecTab_t& xtrack, const std::vector<Eigen::MatrixXd>& cartesianTrack,
   const Eigen::VectorXd& rho, const Saturation& L, const IKopt& IK_opt) {


    struct timeval tbegin,tend;
    double texec = 0.0;
    
    unsigned int iterMax = 15; //DDP iteration max
 
    thetalistd0.setZero();

    xbar.setZero();
    cbar.setZero();
    ubar.setZero();

    x_temp.setZero();
    c_temp.setZero();
    u_temp.setZero();

    x_temp2.setZero();
    c_temp2.setZero();
    u_temp2.setZero();

    /* -------------------- orocos kdl robot initialization-------------------------*/
    KUKAModelKDLInternalData robotParams;
    robotParams.numJoints = 7;
    robotParams.Kv = Eigen::MatrixXd(7,7);
    robotParams.Kp = Eigen::MatrixXd(7,7);

    // ILQRSolver::traj lastTraj;

    Eigen::VectorXd q_pos_init( (stateSize-3) / 2);
    Eigen::VectorXd q_vel_init( (stateSize-3) / 2);
    Eigen::VectorXd q_pos_goal( (stateSize-3) / 2);
    Eigen::VectorXd q_vel_goal( (stateSize-3) / 2);

    q_pos_init.setZero();
    q_pos_goal.setZero();
    q_vel_goal.setZero();
    q_vel_init.setZero();

    q_pos_init = xinit.head((stateSize-3)/2);
    q_vel_init = xinit.segment((stateSize-3)/2, (stateSize-3)/2);

    // std::cout << q_pos_init.transpose().format(CleanFmt) << std::endl;

    q_pos_goal = xgoal.head(stateSize/2);
    q_vel_goal = xgoal.segment((stateSize-3)/2, (stateSize-3)/2); 


    /*----------------------warm-start-control-----------------------------*/
    // Get the gravity compensation for warm-start
    VectorXd q(9);
    VectorXd qd(9);
    q.setZero();

    qd.setZero();
    Eigen::VectorXd gravityTorque(7);
    kukaRobot->getGravityVector(q_pos_init.data(), gravityTorque);

    /*------------------initialize control input-----------------------*/
    // commandVecTab_t u_0;

    for (unsigned i=0; i < N; i++)
    {
      u_0.col(i) = gravityTorque;
    }


    /* Initialize Solver */
    optimizer::ILQRSolverADMM::traj lastTraj;

    // cost function. TODO: make this updatable
    CostFunctionADMM costFunction_admm(xgoal, xtrack, rho);


    /* -------------------- Optimizer Params ------------------------ */
    optimizer::ILQRSolverADMM::OptSet solverOptions;
    solverOptions.n_hor    = N;
    solverOptions.tolFun   = ADMM_OPTS.tolFun;
    solverOptions.tolGrad  = ADMM_OPTS.tolGrad;
    solverOptions.max_iter = iterMax;

    // TODO: make this updatable
    optimizer::ILQRSolverADMM solverDDP(KukaArmModel, costFunction_admm, solverOptions, N, ADMM_OPTS.dt, ENABLE_FULLDDP, ENABLE_QPBOX);


    // // Initialize Trajectory to get xnew with u_0 
    solverDDP.solve(xinit, u_0, cbar, xbar, ubar);


    /* ---------------------------------------- Initialize IK solver ---------------------------------------- */
    
    IKTrajectory<IK_FIRST_ORDER> IK_solve = IKTrajectory<IK_FIRST_ORDER>(IK_opt.Slist, IK_opt.M, IK_opt.joint_limits, IK_opt.eomg, IK_opt.ev, rho, N);
    IK_solve.getTrajectory(cartesianTrack, xnew.col(0).head(7), thetalistd0, xbar.block(0, 0, 7, N + 1), xbar.block(0, 0, 7, N + 1),  &joint_positions_IK);

    /* ------------------------------------------------------------------------------------------------------ */


    lastTraj = solverDDP.getLastSolvedTrajectory();
    xnew = lastTraj.xList;
    unew = lastTraj.uList;
    final_cost[0] = lastTraj.finalCost;

    for (unsigned int k = 0; k < N; k++) {
      x_lambda.col(k) = xnew.col(k) - xbar.col(k);
      u_lambda.col(k) = unew.col(k) - ubar.col(k);
    }

    x_lambda.col(N) = xnew.col(N) - xbar.col(N);
    double cost = 0.0;
   
    /* ------------------------------------------------ Run ADMM ---------------------------------------------- */
    std::cout << "\n ================================= begin ADMM =================================" << std::endl;
    gettimeofday(&tbegin,NULL);

    for (unsigned int i = 0; i < ADMM_OPTS.ADMMiterMax; i++) {
        // TODO: Stopping criterion is needed

        for (unsigned int k = 0;k < N; k++) {
            x_temp.col(k) = xbar.col(k) - x_lambda.col(k);
            u_temp.col(k) = ubar.col(k) - u_lambda.col(k);
        }
        x_temp.col(N) = xbar.col(N) - x_lambda.col(N);

        std::cout << "\n ================================= ADMM iteration " << i + 1 << " ================================= \n";

        // iLQR solver block
        solverDDP.solve(xinit, u_0, cbar, xbar, ubar);
        lastTraj = solverDDP.getLastSolvedTrajectory();


        xnew = lastTraj.xList;
        unew = lastTraj.uList;

        // IK block
        IK_solve.getTrajectory(cartesianTrack, xnew.col(0).head(7), thetalistd0, xbar.block(0, 0, 7, N + 1), xbar.block(0, 0, 7, N + 1),  &joint_positions_IK);


        /* ADMM update */
        // Projection block to feasible sets (state and control contraints)
        xbar_old = xbar;
        ubar_old = ubar;

        for (unsigned int k = 0;k < N; k++) {
            x_temp2.col(k) = xnew.col(k) + x_lambda.col(k);
            u_temp2.col(k) = unew.col(k) + u_lambda.col(k);
        }

        x_temp2.col(N) = xnew.col(N) + x_lambda.col(N);

        /* Projection */
        xubar = projection(x_temp2, u_temp2, L);


        /* Dual variables update */
        for (unsigned int j = 0;j < N; j++) {

            cbar.col(j) = xubar.col(j).head(stateSize);
            xbar.col(j) = xubar.col(j).head(stateSize);
            ubar.col(j) = xubar.col(j).tail(commandSize);
            // cout << "u_bar[" << j << "]:" << ubar[j].transpose() << endl;

            c_lambda.col(j) += xnew.col(j) - xbar.col(j);
            x_lambda.col(j) += xnew.col(j) - xbar.col(j);
            u_lambda.col(j) += unew.col(j) - ubar.col(j);
            // cout << "u_lambda[" << j << "]:" << u_lambda[j].transpose() << endl;

            // Save residuals for all iterations
            res_c[i] += (xnew.col(j) - xbar.col(j)).norm();
            res_x[i] += (xnew.col(j) - xbar.col(j)).norm();
            res_u[i] += (unew.col(j) - ubar.col(j)).norm();

            // res_xlambda[i] += vel_weight*(xbar.col(j) - xbar_old.col(j)).norm();
            res_clambda[i] += (cbar.col(j) - cbar_old.col(j)).norm();
            res_xlambda[i] += (xbar.col(j) - xbar_old.col(j)).norm();
            res_ulambda[i] += 0*(ubar.col(j) - ubar_old.col(j)).norm();
        }

        xbar.col(N) = xubar.col(N).head(stateSize);
        x_lambda.col(N) += xnew.col(N) - xbar.col(N);

        res_x[i] += (xnew.col(N) - xbar.col(N)).norm();
        res_xlambda[i] += 0*(xbar.col(N) - xbar_old.col(N)).norm();

        // get the cost without augmented Lagrangian terms
        cost = 0;
        for (int i = 0;i < N;i++) {
            cost = cost + costFunction_admm.cost_func_expre(i, xnew.col(i), unew.col(i));
        }

        final_cost[i + 1] = cost;
    }

    gettimeofday(&tend,NULL);    

    // testSolverKukaArm.firstInitSolver(xinit, xgoal, xbar, ubar, unew, N, dt, iterMax, tolFun, tolGrad);
    // testSolverKukaArm.initializeTraj();
    // xnew = testSolverKukaArm.updatedxList;
    

    // joint_state_traj.resize(N+1);
    // joint_state_traj_interp.resize(N * InterpolationScale + 1);

    // for(unsigned int i = 0; i <= N; i++) {
    //   joint_state_traj[i] = xnew[i];
    // }

    // torque_traj = unew;

    // //linear interpolation to 1ms
    // for(unsigned int i = 0; i < stateSize; i++) {

    //   for(unsigned int j = 0; j < N * InterpolationScale; j++) {
    //    unsigned int index = j / 10;
    //    joint_state_traj_interp[j](i,0) =  joint_state_traj[index](i,0) + (static_cast<double>(j)-static_cast<double>(index*10.0))*(joint_state_traj[index+1](i,0) - joint_state_traj[index](i,0))/10.0;
    //   }

    //   joint_state_traj_interp[N*InterpolationScale](i,0) = joint_state_traj[N](i,0);
    // }

    // texec=(static_cast<double>(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
    // // texec /= Num_run;

    cout << endl;
    // cout << "Number of iterations: " << lastTraj.iter + 1 << endl;
    cout << "Final cost: " << lastTraj.finalCost << endl;
    // cout << "Final gradient: " << lastTraj.finalGrad << endl;
    // cout << "Final lambda: " << lastTraj.finalLambda << endl;
    // cout << "Execution time by time step (second): " << texec/N << endl;
    // cout << "Execution time per iteration (second): " << texec/lastTraj.iter << endl;
    cout << "Total execution time of the solver (second): " << texec << endl;
    // cout << "\tTime of derivative (second): " << lastTraj.time_derivative.sum() << " (" << 100.0*lastTraj.time_derivative.sum()/texec << "%)" << endl;
    // cout << "\tTime of backward pass (second): " << lastTraj.time_backward.sum() << " (" << 100.0*lastTraj.time_backward.sum()/texec << "%)" << endl;




    cout << "lastTraj.xList[" << N << "]:" << xnew.col(N).transpose() << endl;
    cout << "lastTraj.uList[" << N-1 << "]:" << unew.col(N - 1).transpose() << endl;

    cout << "lastTraj.xList[0]:" << xnew.col(0).transpose() << endl;
    cout << "lastTraj.uList[0]:" << unew.col(0).transpose() << endl;



    std::cout << "================================= ADMM Trajectory Generation Finished! =================================" << std::endl;


    for(unsigned int i = 0; i < ADMM_OPTS.ADMMiterMax; i++) {
      cout << "res_x[" << i << "]:" << res_x[i] << endl;
      // cout << "res_xlambda[" << i << "]:" << res_xlambda[i] << " ";
      // cout << "res_u[" << i << "]:" << res_u[i] << endl;
      // cout << "res_xlambda[" << i << "]:" << res_xlambda[i] << " ";
      cout << "final_cost[" << i << "]:" << final_cost[i] << endl;
    }

}


/* Projection Block 

Projects the states and commands to be within bounds

*/
projStateAndCommandTab_t ADMM::projection(const stateVecTab_t& xnew, const commandVecTab_t& unew, const ADMM::Saturation& L) {
    projStateAndCommandTab_t xubar;

    for(int i = 0;i < NumberofKnotPt + 1; i++) {

        for (int j = 0;j < stateSize + commandSize; j++) {

            if(j < stateSize) { //postion + velocity + force constraints
                if (xnew(j,i) > L.stateLimits(1, j)) {
                    xubar(j,i) = L.stateLimits(1, j);
                }
                else if(xnew(j,i) < L.stateLimits(0, j)) {
                    xubar(j,i) = L.stateLimits(0, j);
                }
                else {
                    xubar(j,i) = xnew(j,i);
                }
            }

            else
            { //torque constraints
                if(unew(j-stateSize, i) > L.controlLimits(1, j)) {
                    xubar(j,i) = L.controlLimits(1, j);
                }
                else if(unew(j-stateSize,i) < L.controlLimits(0, j)) {
                    xubar(j,i) = L.controlLimits(1, j);
                }
                else {
                    xubar(j,i) = unew(j-stateSize, i);
                }
            }

        }
    }
    return xubar;
}



