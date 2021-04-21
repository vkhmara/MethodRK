#pragma once
#include "CommonMethods.h"
#include <vector>
void fff() {}
void adaptiveSimpleRK(vector<double>& nodes, vector<double>& values,
	const methodRK& method) {
	double* k = new double[method.s];
	nodes = vector<double>();
	values = vector<double>();

	nodes.push_back(X0);
	double xCurr = X0;
	values.push_back(Y0);
	double totalTol = EPS / abs(X1 - X0);
	double step = (X1 - X0) / 2;
	double localError = 0;
	double yCurr = Y0, yNext1, yNext2;
	bool lastStep = false;
	while (abs(xCurr - X1) > DBL_EPSILON) {
		if (step < 1e-5) {
			fff();
		}
		yNext1 = nextValue(xCurr, yCurr, step, k, method, f);
		yNext2 = nextValue(xCurr, yCurr, step / 2, k, method, f);
		yNext2 = nextValue(xCurr + step / 2, yNext2, step / 2, k, method, f);
		
		localError = abs(yNext1 - yNext2) / 15;

		if (lastStep || (totalTol * 0.5 * step < abs(localError) && abs(localError) < totalTol * step)) {
			xCurr += step;
			yCurr = yNext1;
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
			/*if (lastStep || (abs(localError) < totalTol * step)) {
			xCurr += step;
			yCurr = yNext1;
			nodes.push_back(xCurr);
			values.push_back(yCurr);
			step = X1 - xCurr;
		}
		else {
			step *= min(pow(totalTol * step / abs(localError), 1. / degree), deltaMax);
		}*/
		if (abs(step) > abs(X1 - xCurr)) {
			step = X1 - xCurr;
			lastStep = true;
		}
	}
	delete[]k;
}