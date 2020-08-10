
#include <iostream>
#include <memory>

/* ADMM trajectory generation */
#include "admm.hpp"

using namespace Eigen;

ADMM::ADMM(const ADMMopt& ADMM_opt, const IKTrajectory<IK_FIRST_ORDER>::IKopt& IK_opt) : ADMM_OPTS(ADMM_opt), IK_OPT(IK_opt)
{
    /* Initalize Primal and Dual variables */
    N = NumberofKnotPt;

    // primal parameters
    xnew.resize(stateSize, N + 1);
    cnew.resize(2, N + 1);
    unew.resize(commandSize, N);

    xbar.resize(stateSize, N + 1);
    cbar.resize(2, N + 1);
    ubar.resize(commandSize, N);

    x_avg.resize(stateSize, N + 1);
    x_lambda_avg.resize(stateSize, N + 1);
    q_lambda.resize(7, N + 1);

    // dual parameters
    x_lambda.resize(stateSize, N + 1);
    c_lambda.resize(2, N + 1);
    u_lambda.resize(commandSize, N);

    x_temp.resize(stateSize, N + 1);
    c_temp.resize(2, N + 1);
    u_temp.resize(commandSize, N);


    u_0.resize(commandSize, N);

    xubar.resize(stateSize + commandSize + 2, N); // for projection

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

    Eigen::VectorXd rho_init(5);
    rho_init << 0, 0, 0, 0, 0;
    IK_solve = IKTrajectory<IK_FIRST_ORDER>(IK_opt.Slist, IK_opt.M, IK_opt.joint_limits, IK_opt.eomg, IK_opt.ev, rho_init, N);

    // initialize curvature object
    curve = Curvature();

    L.resize(N + 1);
    R_c.resize(N + 1);
    k.resize(N + 1 ,3);

    X_curve.resize(3, N + 1);

    robotIK = models::KUKA();


}

/* optimizer execution */
void ADMM::run(std::shared_ptr<KUKAModelKDL>& kukaRobot, KukaArm& KukaArmModel, const stateVec_t& xinit,
  const stateVecTab_t& xtrack, const std::vector<Eigen::MatrixXd>& cartesianTrack,
   const Eigen::VectorXd& rho, const Saturation& L) {


    struct timeval tbegin,tend;
    double texec = 0.0;
    
    unsigned int iterMax = 15; // DDP iteration max


    /* -------------------- orocos kdl robot initialization-------------------------*/
    KUKAModelKDLInternalData robotParams;
    robotParams.numJoints = 7;
    robotParams.Kv = Eigen::MatrixXd(7,7);
    robotParams.Kp = Eigen::MatrixXd(7,7);


    /*------------------initialize control input-----------------------*/

    // cost function. TODO: make this updatable
    CostFunctionADMM costFunction_admm(xtrack.col(N), xtrack, rho);


    /* -------------------- Optimizer Params ------------------------ */
    optimizer::ILQRSolverADMM::OptSet solverOptions;
    solverOptions.n_hor    = N;
    solverOptions.tolFun   = ADMM_OPTS.tolFun;
    solverOptions.tolGrad  = ADMM_OPTS.tolGrad;
    solverOptions.max_iter = iterMax;

    // TODO: make this updatable
    optimizer::ILQRSolverADMM solverDDP(KukaArmModel, costFunction_admm, solverOptions, N, ADMM_OPTS.dt, ENABLE_FULLDDP, ENABLE_QPBOX);

    /* ---------------------------------------- Initial Trajectory ---------------------------------------- */
    // Initialize Trajectory to get xnew with u_0 
    optimizer::ILQRSolverADMM::traj lastTraj;
    solverDDP.initializeTraj(xinit, u_0);

    lastTraj = solverDDP.getLastSolvedTrajectory();
    xnew = lastTraj.xList;
    unew = lastTraj.uList;
    // std::cout << xnew << std::endl;

    final_cost[0] = lastTraj.finalCost;

    /* ---------------------------------------- Initialize IK solver ---------------------------------------- */
    
    IK_solve.getTrajectory(cartesianTrack, xnew.col(0).head(7), xnew.col(0).segment(6,7), xbar.block(0, 0, 7, N + 1), xbar.block(0, 0, 7, N + 1), rho, &joint_positions_IK);
    X_curve = IK_solve.getFKCurrentPos(); 


    /* ------------------------------------------------------------------------------------------------------- */
    double error_fk = 0.0;

    /* ----------------------------------------------- TESTING ----------------------------------------------- */
    for (int i = 0;i < cartesianTrack.size()-1;i++) {
        error_fk = error_fk + (cartesianTrack.at(i) - mr::FKinSpace(IK_OPT.M, IK_OPT.Slist, joint_positions_IK.col(i))).norm();
    }
    std::cout << error_fk << std::endl; 
    /* ----------------------------------------- END TESTING ----------------------------------------- */

    /* ---------------------------------------- Initialize xbar, cbar,ubar ---------------------------------------- */
    xbar.block(0,0,7,xbar.cols()) = joint_positions_IK;
    /* calculates contact terms */
    contact_update(kukaRobot, xnew, &cnew);
    cbar = cnew;
    ubar.setZero();

    x_lambda.setZero();
    c_lambda.setZero();
    u_lambda.setZero();
    q_lambda.setZero();


    double cost = 0.0;

   
    /* ------------------------------------------------ Run ADMM ---------------------------------------------- */
    std::cout << "\n ================================= begin ADMM =================================" << std::endl;

    gettimeofday(&tbegin, NULL);

    for (unsigned int i = 0; i < ADMM_OPTS.ADMMiterMax; i++) {
        // TODO: Stopping criterion is needed
        std::cout << "\n ================================= ADMM iteration " << i + 1 << " ================================= \n";

       /* ---------------------------------------- iLQRADMM solver block ----------------------------------------   */
        solverDDP.solve(xinit, u_0, cbar, xbar, ubar);

        lastTraj = solverDDP.getLastSolvedTrajectory();
        xnew = lastTraj.xList;
        unew = lastTraj.uList;

        /* ----------------------------- update cnew. TODO: variable path curves -----------------------------  */
        contact_update(kukaRobot, xnew, &cnew);

        /* ---------------------------------------- IK block update ----------------------------------------   */
        IK_solve.getTrajectory(cartesianTrack, xnew.col(0).head(7), xnew.col(0).segment(6,7), xbar.block(0, 0, 7, N + 1) - q_lambda, 
            xbar.block(0, 0, 7, N + 1), rho,  &joint_positions_IK);


        /* ------------------------------------- Average States ------------------------------------   */

        x_avg.block(0, 0, 7, xnew.cols()) = (xnew.block(0, 0, 7, xnew.cols())  + joint_positions_IK) / 2;
        x_lambda_avg.block(0, 0, 7, x_lambda_avg.cols()) = (q_lambda + x_lambda.block(0, 0, 7, x_lambda.cols())) / 2;

        /* ---------------------------------------- Projection ----------------------------------------   */
        // Projection block to feasible sets (state and control contraints)

        x_temp.block(0, 0, 7, xnew.cols()) = x_avg.block(0, 0, 7, xnew.cols())  + x_lambda_avg.block(0, 0, 7, x_lambda_avg.cols());//
        c_temp = cnew + c_lambda;
        u_temp = unew + u_lambda;


        xubar = projection(x_temp, c_temp, u_temp, L);


        // for (unsigned int k = 0;k < N; k++) {
        //     x_temp.col(k) = xbar.col(k) - x_lambda.col(k);
        //     u_temp.col(k) = ubar.col(k) - u_lambda.col(k);
        // }
        // x_temp.col(N) = xbar.col(N) - x_lambda.col(N);


        /* Dual variables update */
        for (unsigned int j = 0;j < N; j++) {

            cbar.col(j) = xubar.col(j).segment(stateSize, 2);
            xbar.col(j) = xubar.col(j).head(stateSize);
            ubar.col(j) = xubar.col(j).tail(commandSize);
            // cout << "u_bar[" << j << "]:" << ubar[j].transpose() << endl;


            c_lambda.col(j) += cnew.col(j) - cbar.col(j);
            x_lambda.col(j) += xnew.col(j) - xbar.col(j);
            u_lambda.col(j) += unew.col(j) - ubar.col(j);
            q_lambda.col(j) += joint_positions_IK.col(j) - xbar.col(i).head(7);
            // cout << "u_lambda[" << j << "]:" << u_lambda[j].transpose() << endl;

            // Save residuals for all iterations
            res_c[i] += (cnew.col(j) - cbar.col(j)).norm();
            res_x[i] += (xnew.col(j) - xbar.col(j)).norm();
            res_u[i] += (unew.col(j) - ubar.col(j)).norm();

            // res_xlambda[i] += vel_weight*(xbar.col(j) - xbar_old.col(j)).norm();
            // res_clambda[i] += 0*(cbar.col(j) - cbar_old.col(j)).norm();
            // res_xlambda[i] += 0*(xbar.col(j) - xbar_old.col(j)).norm();
            // res_ulambda[i] += 0*(ubar.col(j) - ubar_old.col(j)).norm();
        }

        

        xbar.col(N - 1) = xubar.col(N - 1).head(stateSize);
        x_lambda.col(N) += xnew.col(N) - xbar.col(N);
        

        res_x[i] += (xnew.col(N) - xbar.col(N)).norm();
        // res_xlambda[i] += 0*(xbar.col(N) - xbar_old.col(N)).norm();

        /* ------------------------------- get the cost without augmented Lagrangian terms ------------------------------- */
        cost = 0;
        for (int i = 0;i < N;i++) {
            cost = cost + costFunction_admm.cost_func_expre(i, xnew.col(i), unew.col(i));
        }

        final_cost[i + 1] = cost;
    }

    gettimeofday(&tend,NULL);    

    solverDDP.initializeTraj(xinit, unew);

    lastTraj = solverDDP.getLastSolvedTrajectory();
    xnew = lastTraj.xList;
    unew = lastTraj.uList;


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
Eigen::MatrixXd ADMM::projection(const stateVecTab_t& xnew, const Eigen::MatrixXd& cnew, const commandVecTab_t& unew, const ADMM::Saturation& L) {

    for(int i = 0;i < NumberofKnotPt ; i++) {

        for (int j = 0;j < stateSize + commandSize + 2; j++) {
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
            } else if((j >= stateSize) && (j < (stateSize + 2))) {
                // std::cout << j << std::endl;
                // if((cnew(0, j)) > 0.3 * std::abs(cnew(1, j))) {
                //     xubar(j,i) = 0.3 * std::abs(cnew(1, i));
                // }

            } else { //torque constraints

                if(unew(j - stateSize - 2, i) > L.controlLimits(1, j - stateSize - 2)) {
                    xubar(j, i) = L.controlLimits(1, j - stateSize - 2);
                }
                else if(unew(j - stateSize - 2, i) < L.controlLimits(0, j - stateSize - 2)) {
                    xubar(j, i) = L.controlLimits(0, j - stateSize - 2);
                }
                else {
                    xubar(j, i) = unew(j-stateSize-2, i);
                }
            }

        }
    }
    // std::cout << "projection end" << std::endl;
    return xubar;
}

/* Computes contact terms
*/
void ADMM::contact_update(std::shared_ptr<KUKAModelKDL>& kukaRobot, const stateVecTab_t& xnew, Eigen::MatrixXd* cnew) {
    double vel = 0.0;
    double m = 0.3; 
    double R = 0.4;

    Eigen::MatrixXd jacobian(6, 7);

    curve.curvature(X_curve.transpose(), L, R_c, k);

    for (int i = 0; i < xnew.cols(); i++) {
        // TODO: spatialJacobian won't work
        kukaRobot->getSpatialJacobian(const_cast<double*>(xnew.col(i).head(7).data()), jacobian);
        vel = (jacobian * xnew.col(i).segment(6, 7)).norm();
        (*cnew)(0,i) = m * vel * vel / R_c(i);
        (*cnew)(1,i) = *(const_cast<double*>(xnew.col(i).tail(1).data()));
    }

}



