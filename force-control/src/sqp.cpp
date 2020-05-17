#include "ilqr.h"

#ifdef DEBUG
#include <iostream>
#endif


namespace optimizer {


void ILQR::ILQR(const DynamicsT& DynamicModel, const CostT& CostFunction, bool fullDDP, bool QPBox)
{

    PRINT("initialize dynamic model and cost function\n");

    dynamicModel = &DynamicModel;
    costFunction = &CostFunction;
    stateNb = iiwaDynamicModel.getStateNb();
    commandNb = iiwaDynamicModel.getCommandNb();
    enableQPBox = QPBox;
    enableFullDDP = fullDDP;

    (enableQPBox) ? PRINT("Box QP is enabled\n") : PRINT("Box QP is disabled\n");
    (enableFullDDP) ? PRINT("Full DDP is enabled\n") :  PRINT("Full DDP is disabled\n");

    //tOptSet Op = INIT_OPTSET;
}

iLQR::~ILQR()
{


}

void ILQR::InitSolver(State& iiwaxInit, State& iiwaxgoal, Control initialTorque, unsigned int& iiwaN,
                       double& iiwadt, unsigned int& iiwamax_iter, double& iiwatolFun, double& iiwatolGrad)
{
    // TODO: double check opt params
    xInit = iiwaxInit; // removed iiwaxgoal. Double check whether this makes sense.
    xgoal = iiwaxgoal;
    N = iiwaN;
    dt = iiwadt;
    initCommand = initialTorque;
    
    PRINT("Initializing the trajectory solver options : ILQG\n");

    Op = new ILQGOptionSet(OpSet);

    Op.xInit = iiwaxInit;
    Op.n_hor = N;
    Op.tolFun = iiwatolFun;
    Op.tolGrad = iiwatolGrad;
    Op.max_iter = iiwamax_iter;

    xList.resize(N+1);
    uList.resize(N);
    uListFull.resize(N+1);
    updatedxList.resize(N+1);
    updateduList.resize(N);
    costList.resize(N+1);
    costListNew.resize(N+1);
    kList.resize(N);
    KList.resize(N);
    FList.resize(N+1);
    Vx.resize(N+1);
    Vxx.resize(N+1);
    
    for (int i = 0; i < N; i++) {
        xList[i].setZero();
        uList[i].setZero();
        uListFull[i].setZero();
        updatedxList[i].setZero();
        updateduList[i].setZero();
        costList[i] = 0;
        costListNew[i] = 0;
        kList[i].setZero();
        KList[i].setZero();
        FList[i].setZero();    
        Vx[i].setZero();
        Vxx[i].setZero();
    }

    // xList[N].setZero();
    // uListFull[N].setZero();
    // updatedxList[N].setZero();
    // costList[N] = 0;
    // costListNew[N] = 0;
    // FList[N].setZero();
    // Vx[N].setZero();
    // Vxx[N].setZero();
    
    // k.setZero();
    // K.setZero();
    // dV.setZero();

    /* parameters for line search */
    // Op.alphaList = new double[11];
    // Op.alphaList << 1.0, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;


}

void ILQR::Solve()
{
    #ifdef DEBUG
    PRINT("Initializing the trajectory...");
    #endif

    // Initialize the trajectory
    InitializeTrajectory();

    /* ------------------------------------------------------------------------------ */


    Op.lambda = Op.lambdaInit;
    Op.dlambda = Op.dlambdaInit;
    
    for (int iter = 0; iter < Op.max_iter; iter++)
    {
        #ifdef
        PRINT("======================================== STEP 1: differentiate dynamics and cost along new trajectory ============================================\n");
        #endif

        if(newDeriv){
            int nargout = 7;//fx,fu,cx,cu,cxx,cxu,cuu
            for(int i = 0;i < u_NAN.size(); i++)
                u_NAN(i,0) = sqrt(-1.0); // control vector = Nan for last time step
            
            for (int i = 0; i < uList.size(); i++) {
                uListFull[i] = uList[i];

            }
            uListFull[uList.size()] = u_NAN;

            // FList is empty here on the first interation
            // initial x, u, cost is passed in
            dynamicModel->kuka_arm_dyn_cst_ilqr(nargout, xList, uListFull, FList, costFunction);
            //dynamicModel->kuka_arm_dyn_cst(nargout, dt, xList, uListFull, xgoal, FList, costFunction->getcx(), costFunction->getcu(), costFunction->getcxx(), costFunction->getcux(), costFunction->getcuu(), costFunction->getc());

            newDeriv = 0;
        }


        #ifdef DEBUG
        PRINT("======================================= STEP 2: backward pass, compute optimal control law and cost-to-go ==========================================\n");
        #endif
        doBackwardPass();

        // if te path diverges, do this
        if (diverge) {
            if (Op.debug_level > 1) printf("Cholesky failed at timestep %d.\n",diverge);
            Op.dlambda   = max(Op.dlambda * Op.lambdaFactor, Op.lambdaFactor);
            Op.lambda    = max(Op.lambda * Op.dlambda, Op.lambdaMin);
            if (Op.lambda > Op.lambdaMax) break;

                continue;
        }

        const int backPassDone = 1;
        

        // TODO: add constraint tolerance check
        if (Op.g_norm < Op.tolGrad && Op.lambda < 1e-5){
            Op.dlambda= min(Op.dlambda / Op.lambdaFactor, 1.0 / Op.lambdaFactor);
            Op.lambda= Op.lambda * Op.dlambda * (Op.lambda > Op.lambdaMin);
            if (Op.debug_level >= 1){
                PRINT(("\nSUCCESS: gradient norm < tolGrad\n"));
            }
            break;
        }

        #ifdef DEBUG
        PRINT("======================================= STEP 3: line-search to find new control sequence, trajectory, cost =========================================== \n");
        #endif

        fwdPassDone = 0;
        if(backPassDone){
            //only implement serial backtracking line-search
            for (int alpha_index = 0; alpha_index < Op.alphaList.size(); alpha_index++){
                alpha = Op.alphaList[alpha_index];
                DoForwardPass();
                Op.dcost = accumulate(costList.begin(), costList.end(), 0.0) - accumulate(costListNew.begin(), costListNew.end(), 0.0);
                Op.expected = -alpha*(dV(0) + alpha*dV(1));

                double z;
                if (Op.expected > 0) {
                    z = Op.dcost / Op.expected;
                }
                else 
                {
                    z = static_cast<double>(-signbit(Op.dcost)); //[TODO:doublecheck]
                    PRINT("non-positive expected reduction: should not occur \n"); //warning
                }

                if (z > Op.zMin)
                { 
                    fwdPassDone = 1;
                    break;
                }
            }
            if (!fwdPassDone) alpha = sqrt(-1.0);
        }
        
        #ifdef DEBUG
        PRINT("======================================= STEP 4: accept step (or not), draw graphics, print status =====================================================")
        #endif

        if (Op.debug_level > 1 && Op.last_head == Op.print_head){
            Op.last_head = 0;
            PRINT("iteration,\t cost, \t reduction, \t expected, \t gradient, \t log10(lambda) \n");
        }
        
        if (fwdPassDone)
        {

            Op.dlambda = min(Op.dlambda / Op.lambdaFactor, 1.0 / Op.lambdaFactor);
            Op.lambda = Op.lambda * Op.dlambda * (Op.lambda > Op.lambdaMin);

            // accept changes
            xList = updatedxList;
            uList = updateduList;
            costList = costListNew;
            newDeriv = 1;

            // TODO: add constraint tolerance check
            if (Op.dcost < Op.tolFun) 
            {
                if (Op.debug_level >= 1)
                {
                    PRINT(("\nSUCCESS: cost change < tolFun\n"));
                }
                break;
            }
        }
        else 
        { 
            Op.dlambda= max(Op.dlambda * Op.lambdaFactor, Op.lambdaFactor);
            Op.lambda= max(Op.lambda * Op.dlambda, Op.lambdaMin);

            // terminate ?
            if (Op.lambda > Op.lambdaMax) 
            {
                if (Op.debug_level >= 1) 
                {
                    PRINT(("\nEXIT: lambda > lambdaMax\n"));
                }
                break;
            }
        }
    }

    Op.iterations = iter;

    if (!backPassDone) 
    {
        if (Op.debug_level >= 1) 
        {
            PRINT(("\nEXIT: no descent direction found.\n"));
        }
        
        return;    
    } else if (iter >= Op.max_iter) 
    {
        if (Op.debug_level >= 1)
        {
            PRINT(("\nEXIT: Maximum iterations reached.\n"));
        }
        return;
    }
}


bool ILQR::InitializeTrajectory()
{
    xList[0] = Op.xInit;
    commandVec_t zeroCommand;
    zeroCommand.setZero();

    // (low priority) TODO: implement control limit selection
    // (low priority) TODO: initialize trace data structure

    initFwdPassDone = 0;
    diverge = 1;
    
    for (int alpha_index = 0; alpha_index < Op.alphaList.size(); alpha_index++)
    {
        alpha = Op.alphaList[alpha_index]; // alpha = 1, 

        for (int i=0; i < N; i++)
        {
            uList[i] = initCommand[i];
        }

        DoForwardPass();

        //simplistic divergence test
        int diverge_element_flag = 0;
        for (int i = 0; i < xList.size(); i++)
        {
            for (int j = 0; j < xList[i].size(); j++){
                if (fabs(xList[i](j,0)) > 1e8) 
                {
                    diverge_element_flag = 1;
                }
            }
        }

        if (!diverge_element_flag)
        {
            cout << "Not Diverge?" << endl;
            diverge = 0;
            break;
        }
    }
    
    initFwdPassDone = 1;
    xList = updatedxList;

    //constants, timers, counters
    newDeriv = 1; 
    Op.lambda= Op.lambdaInit;

    // not used.
    Op.w_pen_l= Op.w_pen_init_l;
    Op.w_pen_f= Op.w_pen_init_f;
    Op.dcost = 0;
    Op.expected = 0;
    Op.print_head = 6;
    Op.last_head = Op.print_head;
    
    #ifdef
    PRINT("\n=========== begin iLQR ===========\n");
    #endif

    return true
}


void ILQR::DoBackwardPass()
{    
    if (Op.regType == 1)
        lambdaEye = Op.lambda*stateMat_t::Identity();
    else
        lambdaEye = Op.lambda*stateMat_t::Zero();

    diverge = 0;
    
    g_norm_sum = 0.0;
    Vx[N] = costFunction->getcx()[N];
    Vxx[N] = costFunction->getcxx()[N];
    dV.setZero();

    for (int i = N-1; i >= 0; i--)
    {
        Qx  = costFunction->getcx()[i]  + dynamicModel->getfxList()[i].transpose() * Vx[i+1];
        Qu  = costFunction->getcu()[i]  + dynamicModel->getfuList()[i].transpose() * Vx[i+1];
        Qxx = costFunction->getcxx()[i] + dynamicModel->getfxList()[i].transpose() * (Vxx[i+1])*dynamicModel->getfxList()[i];
        Quu = costFunction->getcuu()[i] + dynamicModel->getfuList()[i].transpose() * (Vxx[i+1])*dynamicModel->getfuList()[i];
        Qux = costFunction->getcux()[i] + dynamicModel->getfuList()[i].transpose() * (Vxx[i+1])*dynamicModel->getfxList()[i];

        if(Op.regType == 1)
            QuuF = Quu + Op.lambda*commandMat_t::Identity();
        else
            QuuF = Quu;
        
        QuuInv = QuuF.inverse();

        if(!isPositiveDefinite(Quu))
        {
            
            //To be Implemented : Regularization (is Quu definite positive ?)
            PRINT("============================= Quu is not positive definite ===============================\n");

            (Op.lambda==0.0) ? Op.lambda += 1e-4 : Op.lambda *= 10; 
            backPassDone = 0;
            break;
        }


        if (!enableQPBox)
        {
            //TRACE("Use Cholesky decomposition");
            Eigen::LLT<MatrixXd> lltOfQuuF(QuuF);
            Eigen::MatrixXd L = lltOfQuuF.matrixU(); 
            //assume QuuF is positive definite
            
            //A temporary solution: check the non-PD case
            if (lltOfQuuF.info() == Eigen::NumericalIssue)
            {
                diverge = i;
                PRINT("Possibly non semi-positive definitie matrix!");
                return;
            }

            Eigen::MatrixXd L_inverse = L.inverse();
            k = - L_inverse*L.transpose().inverse()*Qu;
            K = - L_inverse*L.transpose().inverse()*Qux;
        }

        // update cost-to-go approximation
        dV(0) += k.transpose()*Qu;
        scalar_t c_mat_to_scalar;
        c_mat_to_scalar = 0.5*k.transpose()*Quu*k;
        dV(1) += c_mat_to_scalar(0,0);
        Vx[i] = Qx + K.transpose()*Quu*k + K.transpose()*Qu + Qux.transpose()*k;
        Vxx[i] = Qxx + K.transpose()*Quu*K+ K.transpose()*Qux + Qux.transpose()*K;
        Vxx[i] = 0.5*(Vxx[i] + Vxx[i].transpose());

        kList[i] = k;
        KList[i] = K;

        g_norm_max= 0.0;
        for(unsigned int j=0; j<commandSize; j++) {
            g_norm_i = fabs(kList[i](j,0)) / (fabs(uList[i](j,0))+1.0);
            if(g_norm_i > g_norm_max) g_norm_max = g_norm_i;
        }
        g_norm_sum += g_norm_max;
    }
    Op.g_norm = g_norm_sum/(static_cast<double>(Op.n_hor));
}

void ILQR::DoForwardPass()
{
    updatedxList[0] = Op.xInit;
    int nargout = 2;

    stateVec_t x_unused;
    x_unused.setZero();
    commandVec_t u_NAN_loc;
    u_NAN_loc << sqrt(-1.0); // sqrt(-1)=NaN. After this line, u_nan has nan for [0] and garbage for rest
    isUNan = 0;

    //[TODO: to be optimized]
    if (!initFwdPassDone)
    {
        for (int i = 0; i < N; i++) {
            updateduList[i] = uList[i];
            updateduList[i] = uList[i] + alpha*kList[i] + KList[i]*(updatedxList[i]-xList[i]);
            // running cost, state, input
            ComputeDynamicsCost(nargout, dt, updatedxList[i], updateduList[i], isUNan, updatedxList[i+1], costFunction);
            costList[i] = costFunction->getc();
        }
        // getting final cost, state, input=NaN
        isUNan = 1;
        dynamicModel->kuka_arm_dyn_cst_min_output(nargout, dt, updatedxList[N], u_NAN_loc, isUNan, x_unused, costFunction);
        costList[N] = costFunction->getc();
    }
}


void ILQRSolver::ComputeDynamicsCost(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, 
                                CostFunctionKukaArm*& costFunction){
    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    int Nl = xList.size();
    
    costFunction->getc() = 0;
    AA.setZero();
    BB.setZero();

    scalar_t c_mat_to_scalar;


    if {
        const int nargout_update2 = 3;
        for (int k = 0; k < Nl; k++) {
            if(k == Nl-1) 
            {
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQf()*(xList[k] - xgoal);
                costFunction->getc() += c_mat_to_scalar(0,0);
            }
            else 
            {
                FList[k] = ForwardStep(nargout_update2, xList[k], uList[k]);    //assume three outputs, code needs to be optimized
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQ()*(xList[k] - xgoal);

                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunction->getR()*uList[k];
                costFunction->getc() += c_mat_to_scalar(0,0); // TODO: to be checked
    
                // AA = fx, BB = fu               
                A_temp[k] = AA;
                B_temp[k] = BB;
            }
        } 

    }

    #ifdef DEBUG
    PRINT("finish kuka_arm_dyn_cst\n");
    #endif
}


// Integrates the trajectoy for the forward pass using simple forward euler
State ILQRSolver::ForwardStep(const State& X, const Control){
    
    // output of kuka arm dynamics is xdot = f(x,u)
    Xdot1 = kuka_arm_dynamics(X, U);
    stateVec_t X_new;

    // simple Euler Integration (for debug)
    X_new = X + dt * Xdot1;

    return X_new;
}


// checks is the passed matrix is positive definite.
bool ILQR::IsPositiveDefinite(const commandMat_t & Quu_p)
{
    Eigen::VectorXcd singular_values = Quu_p.eigenvalues();
    for (int i = 0; i < Quu_p.cols(); ++i)
    {
        if (singular_values[i].real() < 0.)
        {
            PRINT("Matrix is not SDP...");
            return false;
        }
    }
    return true;
}

// trajectory
ILQRSolver::Trajectory ILQRSolver::getLastSolvedTrajectory()
{
    lastTraj.xList = xList;
    // for(unsigned int i=0;i<N+1;i++)lastTraj.xList[i] += xgoal;//retrieve original state with xgoal
    lastTraj.uList = uList;
    lastTraj.iter = iter;
    lastTraj.finalCost = accumulate(costList.begin(), costList.end(), 0.0);
    lastTraj.finalGrad = Op.g_norm;
    lastTraj.finalLambda = log10(Op.lambda);
    lastTraj.time_forward = Op.time_forward;
    lastTraj.time_backward = Op.time_backward;
    lastTraj.time_derivative = Op.time_derivative;
    return lastTraj;
}


}  // namespace kuka_iiwa_arm
