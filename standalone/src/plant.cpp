
#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "KukaModel.h"
#include "models.h"


stateVec_t KukaPlant::f(const Eigen::Ref<const State> &x, const Eigen::Ref<const Control> &u) 
{
    controlVec_t u_noisy = u + cdist_.samples(1);

    stateVec_t xnew = x + cp_.f(x, u_noisy) * dt_ + sdist_.samples(1);
    Eigen::Map<State>(sx.data(), StateSize) = xnew;
    Eigen::Map<Control>(su.data(), ControlSize) = u_noisy;
    return xnew;

}