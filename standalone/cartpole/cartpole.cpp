#include "cartpole.h"
#include "cartpolewindow.h"
#include <ddp/costs.hpp>
#include <ddp/ddp.hpp>
#include <ddp/mpc.hpp>
#include <ddp/util.hpp>
#include <QApplication>
#include <QtConcurrent/QtConcurrent>
#include <iostream>
#include <vector>

int main(int argc, char **argv)
{
    using Scalar = double;
    using Dynamics = CartPole<Scalar>;
    using Plant = CartPolePlant;
    using Cost = CartPoleCost<Scalar>;
    using TerminalCost = CartPoleTerminalCost<Scalar>;

    // Log output
    util::DefaultLogger logger;
    bool verbose = true;

    // Timing
    Scalar tf = 2.0;
    Scalar dt = 0.025;
    auto time_steps = util::time_steps(tf, dt);
    int max_iterations = 15;

    // Dynamics
    Dynamics cp_dynamics(9.80, 1.0, 0.1, 0.50, 10.0);

    // Noisy Plant
    Scalar ssigma = 0.01;
    Scalar csigma = 0.00;
    Plant cp_plant(cp_dynamics, dt, ssigma, csigma);

    // Initial state, desired state, initial control sequence
    Dynamics::State x0 = Dynamics::State::Zero();
    Dynamics::State xf; xf << 0.0, 0.0, M_PI, 0.0;
    Dynamics::ControlTrajectory u = Dynamics::ControlTrajectory::Zero(time_steps);

    // Costs
    Cost::StateHessian Q;
    Q.setZero();
    Q.diagonal() << 1.0, 1.0, 100.0, 10.0;

    Cost::ControlHessian R;
    R << 10.0;

    TerminalCost::Hessian Qf;
    Qf.setZero();
    Qf.diagonal() << 1.0, 1.0, 1000.0, 100.0;

    Cost cp_cost(xf, Q, R);
    TerminalCost cp_terminal_cost(xf, Qf);

    // Initialize receding horizon controller
    ModelPredictiveController<Dynamics, optimizer::DDP> mpc(dt, time_steps, max_iterations, &logger, verbose, /*DDP options*/ &logger, verbose);

    // Connect signals and configure plot window
    QApplication a(argc, argv);
    MainWindow win(xf(2), xf(0), dt);
    QObject::connect(&cp_plant, &Plant::update, &win, &MainWindow::update, Qt::QueuedConnection);
    win.show();

    // Define a termination condition
    using Result = ModelPredictiveController<Dynamics, optimizer::DDP>::Result;
    using StateRef = Eigen::Ref<const Dynamics::State>;
    auto termination =
        [&](int i, Result &r, const StateRef &x, Cost &rc, TerminalCost &tc, Scalar true_cost)
        {
            return (QThreadPool::globalInstance()->activeThreadCount() == 0);
        };

    // Run in receding horizon mode
    QtConcurrent::run([&](){ mpc.run(x0, u, termination, cp_dynamics, cp_plant, cp_cost, cp_terminal_cost); });

    return a.exec();
}
