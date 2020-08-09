#include "curvature.hpp"
#include <cmath>
#include <Eigen/Dense>

int main() {

	Curvature curve = Curvature();

	int N = 100;
	double dt = 10.0 / (N - 1);
	Eigen::MatrixXd X(N, 3);
	double t = 0.0;

	for(int i = 0; i < N; i++) {
		X(i, 0) = std::sin(3*t);
		X(i, 1) = std::cos(t);
		X(i, 2) = 0.8;
		t += dt;
	}

	Eigen::VectorXd L(N);
	Eigen::VectorXd R(N);
	Eigen::MatrixXd k(N ,3);
	curve.curvature(X, L, R, k);

	return 0;
}