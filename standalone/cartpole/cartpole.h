#ifndef DDP_CARTPOLE_H
#define DDP_CARTPOLE_H

#include <ddp/costs.hpp>
#include <ddp/dynamics.hpp>
#include <ddp/plant.hpp>
#include <ddp/eigenmvn.hpp>
#include <QtCore>
#include <QVector>

template <class T>
struct CartPole: public Dynamics<T, 4, 1>
{
    using State = typename Dynamics<T, 4, 1>::State;
    using Control = typename Dynamics<T, 4, 1>::Control;

    CartPole(T g, T mc, T mp, T l, T fm) : gravity(g), masscart(mc), masspole(mp), length(l), force_mag(fm) {}
    inline State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        using std::sin; using std::cos;

        // state transition = [xdot, xacc, thetadot, thetaacc]
        State dx;

        // state = [x, xdot, theta, thetadot]
        const T force = u.value() * force_mag;
        const T costheta = cos(x(2));
        const T sintheta = sin(x(2));

        const T temp = (force + polemass_length() * x(3) * x(3) * sintheta) / total_mass();
        const T thetaacc = (gravity * sintheta - costheta * temp) / (length * (4.0/3.0 - masspole * costheta * costheta / total_mass()));
        const T xacc = temp - polemass_length() * thetaacc * costheta / total_mass();

        dx(0) = x(1);
        dx(1) = xacc;
        dx(2) = x(3);
        dx(3) = thetaacc;

        return dx;
    }

    T gravity;
    T masscart;
    T masspole;
    T length;
    T force_mag;

private:
    T total_mass() const { return this->masspole + this->masscart; }
    T polemass_length() const { return this->masspole * this->length; }
};

template <class T>
struct CartPoleCost: public CostFunction<CartPole<T>>
{
    using Scalar = T;
    using Dynamics = CartPole<T>;
    using State = typename CostFunction<CartPole<T>>::State;
    using Control = typename CostFunction<CartPole<T>>::Control;
    using Gradient = typename CostFunction<CartPole<T>>::Gradient;
    using Hessian = typename CostFunction<CartPole<T>>::Hessian;
    using StateHessian = Eigen::Matrix<Scalar, CartPole<T>::StateSize, CartPole<T>::StateSize>;
    using ControlHessian = Eigen::Matrix<Scalar, CartPole<T>::ControlSize, CartPole<T>::ControlSize>;

    CartPoleCost(const Eigen::Ref<const State> &xf, const Eigen::Ref<const StateHessian> &Q, const Eigen::Ref<const ControlHessian> &R)
    : CostFunction<CartPole<T>>(xf), Q_(Q), R_(R)
    {
        QR_.setZero();
        QR_.topLeftCorner(Dynamics::StateSize, Dynamics::StateSize) = Q;
        QR_.bottomRightCorner(Dynamics::ControlSize, Dynamics::ControlSize) = R;
    }

    Scalar c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        State error = x - this->target();
        return (error.transpose() * Q_ * error).value() + (u.transpose() * R_ * u).value();
    }

    Gradient dc(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        Gradient g;
        g.head(Dynamics::StateSize) = Q_ * (x - this->target());
        g.tail(Dynamics::ControlSize) = R_ * u;
        return g;
    }

    Hessian d2c(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        return QR_;
    }

private:
    StateHessian Q_;
    ControlHessian R_;
    Hessian QR_;
};

template <class T>
struct CartPoleTerminalCost: public TerminalCostFunction<CartPole<T>>
{
    using Scalar = T;
    using Dynamics = CartPole<T>;
    using State = typename TerminalCostFunction<CartPole<T>>::State;
    using Gradient = typename TerminalCostFunction<CartPole<T>>::Gradient;
    using Hessian = typename TerminalCostFunction<CartPole<T>>::Hessian;

    CartPoleTerminalCost(const Eigen::Ref<const State> &xf, const Eigen::Ref<const Hessian> &Q)
    : TerminalCostFunction<CartPole<T>>(xf), Q_(Q) {}

    Scalar c(const Eigen::Ref<const State> &x)
    {
        return (x - TerminalCostFunction<CartPole<T>>::xf).transpose() *
                Q_ * (x - TerminalCostFunction<CartPole<T>>::xf);
    }

    Gradient dc(const Eigen::Ref<const State> &x)
    {
        return Q_ * (x - TerminalCostFunction<CartPole<T>>::xf);
    }

    Hessian d2c(const Eigen::Ref<const State> &x)
    {
        return Q_;
    }

private:
    Hessian Q_;
};

struct CartPolePlant: public QObject, public Plant<double, 4, 1>
{
    using Base                  = Plant<double, 4, 1>;
    using Scalar                = typename Base::Scalar;
    using State                 = typename Base::State;
    using Control               = typename Base::Control;
    using StateNoiseVariance    = Eigen::Matrix<Scalar, StateSize, StateSize>;
    using ControlNoiseVariance  = Eigen::Matrix<Scalar, ControlSize, ControlSize>;

    CartPolePlant(CartPole<Scalar> &cp, Scalar dt, Scalar state_var, Scalar control_var, QObject *parent = 0)
        : QObject(parent), cp_(cp), dt_(dt), sdist_(State::Zero(), state_var * StateNoiseVariance::Identity()),
          cdist_(Control::Zero(), control_var * ControlNoiseVariance::Identity()) {}

    inline State f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u)
    {
        static QVector<Scalar> sx(StateSize);
        static QVector<Scalar> su(ControlSize);

        Control u_noisy = u + cdist_.samples(1);

        State xnew = x + cp_.f(x, u_noisy) * dt_ + sdist_.samples(1);
        Eigen::Map<State>(sx.data(), StateSize) = xnew;
        Eigen::Map<Control>(su.data(), ControlSize) = u_noisy;
        update(sx, su);
        return xnew;
    }

signals:
    void update(const QVector<Scalar> &state, const QVector<Scalar> &control) const;

private:
    Q_OBJECT
    CartPole<Scalar> &cp_;
    Scalar dt_;
    Eigen::EigenMultivariateNormal<Scalar, StateSize> sdist_;
    Eigen::EigenMultivariateNormal<Scalar, ControlSize> cdist_;
};

#endif //DDP_CARTPOLE_H
