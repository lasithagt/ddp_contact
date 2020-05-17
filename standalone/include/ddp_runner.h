/// Repackaged from runddp.cc into a library to be called from another file

#include <iostream>
#include <memory>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <list>

#include "lcm/lcm-cpp.hpp"
#include "drake/lcmt_ddp_traj.hpp"

#include "drake/common/drake_assert.h"
#include "drake/common/find_resource.h"
#include "drake/common/trajectories/piecewise_polynomial.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"
#include "drake/multibody/joints/floating_base_types.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/manipulation/util/world_sim_tree_builder.h"
#include "drake/multibody/rigid_body_plant/rigid_body_plant.h"

#include "drake/lcmt_generic_string_msg.hpp"

#include "drake/DDP_traj_gen/config.h"
#include "drake/DDP_traj_gen/spline.h"
#include "drake/DDP_traj_gen/ilqrsolver.h"
// #include "drake/DDP_traj_gen/udpsolver.h"
#include "drake/DDP_traj_gen/kuka_arm.h"

using namespace std;
using namespace Eigen;

#define useILQRSolver 1
#define useUDPSolver 0

/* DDP trajectory generation */

static std::list< const char*> gs_fileName;
static std::list< std::string > gs_fileName_string;

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {


using trajectories::PiecewisePolynomial;
typedef PiecewisePolynomial<double> PPType;
typedef PPType::PolynomialType PPPoly;
typedef PPType::PolynomialMatrix PPMatrix;

using manipulation::kuka_iiwa::kIiwaArmNumJoints;
using manipulation::util::WorldSimTreeBuilder;
using manipulation::util::ModelInstanceInfo;
using systems::RigidBodyPlant;

class DDPRunner {
public:
void RunUDP(stateVec_t xinit, stateVec_t xgoal);
void saveVector(const Eigen::MatrixXd & _vec, const char * _name);
void saveValue(double _value, const char * _name);
void clean_file(const char * _file_name, std::string & _ret_file);

private:
lcm::LCM lcm_;
lcmt_ddp_traj ddp_traj_;

//UDP parameters
stateVecTab_t joint_state_traj;
commandVecTab_t torque_traj;
stateVecTab_t joint_state_traj_interp;
commandVecTab_t torque_traj_interp;
};

} // namespace kuka_iiwa_arm
} // examples
} // drake