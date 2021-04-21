#pragma once
#include "CommonMethods.h"
#include <vector>
void adaptiveNestedRK(vector<double>& nodes, vector<double>& values,
	const methodRK& innerMethod, const methodRK& outerMethod) {
	double* k = new double[outerMethod.s];
	nodes = vector<double>();
	values = vector<double>();

	nodes.push_back(X0);
	double xCurr = X0;
	values.push_back(Y0);
	double yCurr = Y0;
	double totalTol = EPS / abs(X1 - X0);
	double step = (X1 - X0) / 2;
	double localError = 0;
	bool lastStep = false;
	while (abs(xCurr - X1) > DBL_EPSILON) {
		kFromStep(xCurr, yCurr, step, k, outerMethod, f);
		localError = 0;
		for (int i = 0; i < outerMethod.s; i++) {
			localError += (innerMethod.b[i] - outerMethod.b[i]) * k[i];
		}
		localError *= step;

		if (lastStep || (totalTol * 0.5 < abs(localError) / step && abs(localError) / step < totalTol)) {
			xCurr += step;
			yCurr = nextValue(xCurr, yCurr, step, k, outerMethod, f, true);
			nodes.push_back(xCurr);
			values.push_back(yCurr);
		}
		else {

			if (abs(localError) / step <= totalTol * 0.5) {
				if (localError == 0.0) {
					step *= deltaMin;
				}
				else {
					step *= max(pow(totalTol * step / abs(localError) / 2, 1. / degree), deltaMin);
				}
			}
			else {
				step *= min(pow(totalTol * step / abs(localError), 1. / degree), deltaMax);
			}
		}
		if (abs(step) > abs(X1 - xCurr)) {
			step = X1 - xCurr;
			lastStep = true;
		}
	}
	delete[]k;
}

void constRK(vector<double>& nodes, vector<double>& values, int N, const methodRK& method) {
	double* k = new double[method.s];
	nodes = vector<double>();
	values = vector<double>();

	double step = (X1 - X0) / N;
	double xCurr = X0;
	double yCurr = Y0;
	nodes.push_back(xCurr);
	values.push_back(yCurr);
	for (int i = 1; i <= N; i++) {
		yCurr = nextValue(xCurr, yCurr, step, k, method, f);
		xCurr += step;
		nodes.push_back(xCurr);
		values.push_back(yCurr);
	}

	delete[]k;
}
