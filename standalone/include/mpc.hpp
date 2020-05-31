
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


/* DDP trajectory generation */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <string>
#include <list>

#include "mpc.h"

// using namespace std;
// using namespace Eigen;


/* DDP trajectory generation */
static std::list< const char*> gs_fileName;
static std::list< std::string > gs_fileName_string;

/* ------------- Eigen print arguments ------------------- */
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
 /* ------------------------------------------------------- */


template <class Dynamics, class OptimizerT>
class ModelPredictiveController
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Optimizer             = OptimizerT;
    // using RunningCost           = CostFunction;
    // using TerminalCost          = TerminalCostFunction;
    using Scalar                = double;            ///< Type of scalar used by the optimizer
    using State                 = stateVec_t;            
    using Control               = controlVec_t;           
    using StateTrajectory       = stateVecTab_t;   
    using ControlTrajectory     = controlVecTab_t; 
    // using Result                = OptimizerResult<Dynamics>;            ///< Type of result returned by the optimizer
    using TerminationCondition =
        std::function<bool(int,
                           OptimizerResult<Dynamics> &,
                           const Eigen::Ref<const State> &,
                           RunningCost &, TerminalCost &, Scalar)>;     ///< Type of termination condition function

    static const int MPC_BAD_CONTROL_TRAJECTORY = -1;

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
    ModelPredictiveController(Scalar dt, int time_steps, int iterations, bool verbose)
    : opt_(dt, time_steps, iterations,
      dt_(dt), H_(time_steps), logger_(logger), verbose_(verbose) {}

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
    template <typename TerminationCondition, typename Plant, typename CostFunction, typename TerminalCostFunction>
	void run(const Eigen::Ref<const State>  &initial_state,
	         Eigen::Ref<ControlTrajectory>  initial_control_trajectory,
	         TerminationCondition           &terminate,
	         Dynamics                       &dynamics,
	         Plant                          &plant,
	         CostFunction                   &cost_function)
	         // TerminalCostFunction           &terminal_cost_function)
	{
	    if(initial_control_trajectory.cols() != H_)
	    {
	        logger_->error("The size of the control trajectory does not match the number of time steps passed to the optimizer!");
	        std::exit(MPC_BAD_CONTROL_TRAJECTORY);
	    }

	    State x = initial_state, xold = initial_state;
	    Control u;
	    Scalar true_cost = cost_function.c(xold, initial_control_trajectory.col(0));
	    // OptimizerResult<Dynamics> result;
	    // result.control_trajectory = initial_control_trajectory;
	    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	    std::chrono::duration<float, std::milli> elapsed;

	    int64_t i = 0;
	    while(!terminate(i, result, x, cost_function, terminal_cost_function, true_cost))
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
	        result = opt_.run(xold, result.control_trajectory, dynamics, cost_function, terminal_cost_function, std::forward<ARGS>(args)...);
	        u = result.control_trajectory.col(0);
	        if(verbose_)
	        {
	            logger_->info("Obtained control from optimizer: ");
	            for(int m = 0; m < u.rows(); ++m) { logger_->info("%f ", u(m)); }
	            logger_->info("\n");
	        }

	        // Apply the control to the plant and obtain the new state
	        x = plant.f(xold, u);
	        if(verbose_)
	        {
	            logger_->info("Received new state from plant: ");
	            for(int n = 0; n < x.rows(); ++n) { logger_->info("%f ", x(n)); }
	            logger_->info("\n");
	        }

	        // Calculate the true cost for this time step
	        true_cost = cost_function.c(x, u);
	        if(verbose_) logger_->info("True cost for time step %d: %f\n", i, true_cost);

	        // Slide down the control trajectory
	        result.control_trajectory.leftCols(H_ - 1) = result.control_trajectory.rightCols(H_ - 1);
	        if(verbose_) logger_->info("Slide down the control trajectory\n");

	        xold = x;
	        ++i;
	    }
	}
}


