
/// @file
///
/// kuka_plan_runner is designed to wait for LCM messages contraining
/// a robot_plan_t message, and then execute the plan on an iiwa arm
/// (also communicating via LCM using the
/// lcmt_iiwa_command/lcmt_iiwa_status messages).
///
/// When a plan is received, it will immediately begin executing that
/// plan on the arm (replacing any plan in progress).
///
/// If a stop message is received, it will immediately discard the
/// current plan and wait until a new plan is received.

#include <iostream>
#include <memory>

#include "config.h"
#include "spline.h"
#include "ilqrsolver.h"
#include "kuka_arm.h"
#include "SoftContactModel.h"
#include "KukaModel.h"
#include "models.h"


/* DDP trajectory generation */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <string>
#include <list>
#include <chrono>


/* DDP trajectory generation */
static std::list< const char*> gs_fileName;
static std::list< std::string > gs_fileName_string;

/* ------------- Eigen print arguments ------------------- */
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
 /* ------------------------------------------------------- */


template <class DynamicsT, class PlantT, class costFunctionT, class OptimizerT, class OptimizerResultT>
class ModelPredictiveController
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Scalar                = double;            ///< Type of scalar used by the optimizer
    using Optimizer             = OptimizerT;
    using Dynamics 				= DynamicsT;
    using CostFunction   		= costFunctionT;
    using State                 = stateVec_t;            
    using Control               = commandVec_t;           
    using StateTrajectory       = stateVecTab_t;   
    using ControlTrajectory     = commandVecTab_t; 
    using Plant 				= PlantT;
    using Result                = OptimizerResultT;            ///< Type of result returned by the optimizer
    using TerminationCondition =
        std::function<bool(int, State)>;     ///< Type of termination condition function

    static const int MPC_BAD_CONTROL_TRAJECTORY = -1;
    Optimizer opt_;
    Dynamics& dynamics_;
    // Plant& plant_;
    CostFunction& cost_function_;
    bool verbose_;
    Scalar dt_;
    int H_;
    Logger* logger_;
    ControlTrajectory control_trajectory;
    StateTrajectory x_track_;

public:
    /**
     * @brief               Instantiate the receding horizon wrapper for the trajectory optimizer.
     * @param dt            Time step
     * @param time_steps    Number of time steps over which to optimize
     * @param iterations    The number of iterations to perform per time step
     * @param logger        util::Logger to use for informational, warning, and error messages
     * @param verbose       True if informational and warning messages should be passed to the logger; error messages are always passed
     * @param args          Arbitrary arguments to pass to the trajectory optimizer at initialization time
     */
    ModelPredictiveController(Scalar dt, int time_steps, int iterations, bool verbose, Logger *logger,
    		 Dynamics                       &dynamics,
	         CostFunction                   &cost_function,
	         Optimizer 						&opt,
	         const StateTrajectory  		&x_track)
    : dt_(dt), H_(time_steps), verbose_(verbose), dynamics_(dynamics), cost_function_(cost_function), opt_(opt), x_track_(x_track)  
    {
    	logger_ = logger;
    	control_trajectory.resize(H_);
    }

    /**
     * @brief                               Run the trajectory optimizer in MPC mode.
     * @param initial_state                 Initial state to pass to the optimizer
     * @param initial_control_trajectory    Initial control trajectory to pass to the optimizer
     * @param terminate                     Termination condition to check before each time step
     * @param dynamics                      Dynamics model to pass to the optimizer
     * @param plant                         Plant model
     * @param cost_function                 Running cost function L(x, u) to pass to the optimizer
     * @param terminal_cost_function        Terminal cost function V(xN) to pass to the optimizer
     * @param args                          Arbitrary arguments to pass to the trajectory optimizer at run time
     */
    template <typename TerminationCondition>
	void run(const Eigen::Ref<const State>  &initial_state,
	         ControlTrajectory              initial_control_trajectory,
	         Plant                          &plant_,
	         TerminationCondition           &terminate)
	         // TerminalCostFunction           &terminal_cost_function)
	{
	    if(initial_control_trajectory.size() != H_)
	    {
	        logger_->error("The size of the control trajectory does not match the number of time steps passed to the optimizer!");
	        std::exit(MPC_BAD_CONTROL_TRAJECTORY);
	    }

	    State x = initial_state, xold = initial_state;
	    Control u;
	    // Scalar true_cost = cost_function.c(xold, initial_control_trajectory[0]);
	    stateVec_t x_track;
	    x_track.setZero();
	    scalar_t true_cost = cost_function_.cost_func_expre(0, xold, initial_control_trajectory[0], x_track_[0]);
	    
	    // OptimizerResult<Dynamics> result;
	    Result result;
	    control_trajectory = initial_control_trajectory;
	    u = initial_control_trajectory[0];
	    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	    std::chrono::duration<float, std::milli> elapsed;

	    int64_t i = 0;
	    while(!terminate(i, x))
	    {
	        if(verbose_)
	        {
	            if(i > 0)
	            {
	                end = std::chrono::high_resolution_clock::now();
	                elapsed = end - start;
	                logger_->info("Completed MPC loop for time step %d in %d ms\n", i - 1, static_cast<int>(elapsed.count()));
	            }
	            logger_->info("Entered MPC loop for time step %d\n", i);
	            start = std::chrono::high_resolution_clock::now();
	        }

	        // Run the optimizer to obtain the next control
	        opt_.solve(xold, control_trajectory, x_track_);
	        result = opt_.getLastSolvedTrajectory();

	        u = result.uList[0];
	        if(verbose_)
	        {
	            logger_->info("Obtained control from optimizer: ");
	            for(int m = 0; m < u.rows(); ++m) { logger_->info("%f ", u(m)); }
	            logger_->info("\n");
	        }

	        // Apply the control to the plant and obtain the new state
	        x = plant_.f(xold, u);
	        if(verbose_)
	        {
	            logger_->info("Received new state from plant: ");
	            for(int n = 0; n < x.rows(); ++n) { logger_->info("%f ", x(n)); }
	            logger_->info("\n");
	        }

	        // Calculate the true cost for this time step
	        true_cost = cost_function_.cost_func_expre(0, xold, initial_control_trajectory[0], x_track_[0]);

	        if(verbose_) logger_->info("True cost for time step %d: %f\n", i, true_cost(0,0));

	        // Slide down the control trajectory
	        for (int i = 0; i < H_; i++)
	        {
	        	control_trajectory[i] = result.uList[i + 1];
	        }

	        // control_trajectory.leftCols(H_ - 1) = result.uList.rightCols(H_ - 1);
	        if(verbose_) logger_->info("Slide down the control trajectory\n");

	        xold = x;
	        ++i;
	    }
	}
};



