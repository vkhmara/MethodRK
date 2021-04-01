#pragma once
#include "Consts.h"
struct methodRK {
	double** A;
	double* b;
	double* c;
	int s;
	methodRK():A(nullptr), b(nullptr), c(nullptr), s(0) {}
};

void kFromStep(double xValue, double yValue, double step, double *k, const methodRK& method,
	double func(double, double)) {
	double sum;
	for (int i = 0; i < method.s; i++) {
		sum = 0;
		for (int j = 0; j < i; j++) {
			sum += method.A[i][j] * k[j];
		}
		k[i] = func(xValue + step * method.c[i], yValue + step * sum);
	}
}

double nextValue(double xValue, double yValue, double step, double* k,
	const methodRK& method, double func(double, double), bool calculatedK=false) {
	if (!calculatedK) {
		kFromStep(xValue, yValue, step, k, method, func);
	}
	double sum = 0;
	for (int i = 0; i < method.s; i++) {
		sum += k[i] * method.b[i];
	}
	return yValue + step * sum;
}