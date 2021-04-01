#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "CommonMethods.h"

double** A;
double* b1;
double* b2;
double* c;

int s;
methodRK method1, method2;

double getNumberFromString(string str) {
	int sleshIndex = str.find('/');
	if (sleshIndex == string::npos) {
		return atof(str.c_str());
	}
	return atof(str.substr(0, sleshIndex).c_str()) /
		atof(str.substr(sleshIndex + 1).c_str());
}

void downloadArray(string filename, double* arr) {
	ifstream file(filename);
	char ch[200];
	int prevPos = -1;
	file.getline(ch, 200);
	string content = ch;
	content += " ";
	for (int i = 0; i < s; i++) {
		arr[i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
		prevPos = content.find(' ', prevPos + 1);
	}
	file.close();
}

void downloadMatrix(string filename, double** matrix) {
	string content;
	int prevPos = -1;

	ifstream file(filename);
	char ch[200];
	for (int i = 0; i < s; i++) {
		matrix[0][i] = 0;
	}
	for (int iter = 1; iter < s; iter++) {
		prevPos = -1;
		file.getline(ch, 200);
		content = ch;
		content += " ";
		for (int i = 0; i < iter; i++) {
			matrix[iter][i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
			prevPos = content.find(' ', prevPos + 1);
		}
		for (int i = iter; i < s; i++) {
			matrix[iter][i] = 0;
		}
	}
	file.close();
}

void printAll() {
	cout << "A\n";
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < s; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << "b1\n";
	for (int i = 0; i < s; i++) {
		cout << b1[i] << " ";
	}
	cout << "\nb2\n";
	for (int i = 0; i < s; i++) {
		cout << b2[i] << " ";
	}
	cout << "\nc\n";
	for (int i = 0; i < s; i++) {
		cout << c[i] << " ";
	}
}

void initAll() {
	ifstream file(dimenFile);
	file >> s;
	file.close();

	A = new double* [s];
	for (int i = 0; i < s; i++) {
		A[i] = new double[s];
	}
	b1 = new double[s];
	b2 = new double[s];
	c = new double[s];

	method1.A = A;
	method1.b = b1;
	method1.c = c;
	method1.s = s;

	method2.A = A;
	method2.b = b2;
	method2.c = c;
	method2.s = s;
}

void downloadAll() {
	downloadMatrix(AFile, A);
	downloadArray(b1File, b1);
	downloadArray(b2File, b2);
	downloadArray(cFile, c);
}

void clearAll() {
	for (int i = 0; i < s; i++) {
		delete[]A[i];
	}
	delete[]A;

	delete[]b1;
	delete[]b2;
	delete[]c;
}

void adaptiveRK(vector<double>& nodes, vector<double>& values) {
	double* k = new double[s];
	nodes = vector<double>();
	values = vector<double>();

	nodes.push_back(X0);
	double xCurr = X0;
	values.push_back(Y0);
	double solut = Y0;
	double totalTol = EPS / abs(X1 - X0);
	double step = (X1- X0) / 2;
	double localError = 0;
	bool lastStep = false;
	while (abs(xCurr - X1) > DBL_EPSILON) {
		kFromStep(xCurr, solut, step, k, method2, f);
		localError = 0;
		for (int i = 0; i < s; i++) {
			localError += (b1[i] - b2[i]) * k[i];
		}
		localError *= step;
		
		if (lastStep || (totalTol * 0.5 < abs(localError) / step && abs(localError) / step < totalTol)) {
			xCurr += step;
			solut = nextValue(xCurr, solut, step, k, method2, f, true);
			nodes.push_back(xCurr);
			values.push_back(solut);
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

void constRK(vector<double>& nodes, vector<double>& values, int N) {
	double* k = new double[s];
	nodes = vector<double>();
	values = vector<double>();

	double step = (X1 - X0) / N;
	double xCurr = X0;
	double yCurr = Y0;
	nodes.push_back(xCurr);
	values.push_back(yCurr);
	for (int i = 1; i <= N; i++) {
		yCurr = nextValue(xCurr, yCurr, step, k, method2, f);
		xCurr += step;
		nodes.push_back(xCurr);
		values.push_back(yCurr);
	}

	delete[]k;
}

void outputResults(vector<double>& x, vector<double>& y, string filename) {
	ofstream out(filename);

	int N = x.size();
	for (int i = 0; i < N; i++) {
		out << setprecision(16) <<  x[i] << " " << setprecision(16) << y[i] << endl;
	}

	out.close();
}

int main() {
	initAll();
	downloadAll();

	vector<double> nodes;
	vector<double> values;

	adaptiveRK(nodes, values);
	outputResults(nodes, values, whereOutput + "AdaptiveRK.txt");
//	constRK(nodes, values, 314);
//	outputResults(nodes, values, whereOutput + "ConstRK.txt");
	
	clearAll();
	system("pause");
	return 0;
}