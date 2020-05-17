#ifndef KUKA_ANALYTICAL_SOLUTIONS_H
#define KUKA_ANALYTICAL_SOLUTIONS_H

class KUKAAnalyticalSolutions
{
public:

	KUKAAnalyticalSolutions() = default;

	~KUKAAnalyticalSolutions() = default;

	void FK( double* FK, const double* q);

	void Jacobian( double* jac, const double* q);

	void MassMatrix( double* M, const double* parms, const double* q);

	void Coriolis( double* C, const double* parms, const double* q, const double* dq);

	void Gravity(double* G, const double* parms, const double* q);

	void InverseDynamics( double* INV, const double* parms, const double* q, const double* dq, const double* ddq );

};

#endif //KUKA_H