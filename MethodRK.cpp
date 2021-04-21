#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "CommonMethods.h"
#include "NestedMethodRK.h"
#include "SimpleRK.h"

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

void downloadArray(string filename, double* arr, int dim) {
	ifstream file(filename);
	char ch[200];
	int prevPos = -1;
	file.getline(ch, 200);
	string content = ch;
	content += " ";
	for (int i = 0; i < dim; i++) {
		arr[i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
		prevPos = content.find(' ', prevPos + 1);
	}
	file.close();
}

void downloadMatrix(string filename, double** matrix, int dim) {
	string content;
	int prevPos = -1;

	ifstream file(filename);
	char ch[200];
	for (int i = 0; i < dim; i++) {
		matrix[0][i] = 0;
	}
	for (int iter = 1; iter < dim; iter++) {
		prevPos = -1;
		file.getline(ch, 200);
		content = ch;
		content += " ";
		for (int i = 0; i < iter; i++) {
			matrix[iter][i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
			prevPos = content.find(' ', prevPos + 1);
		}
		for (int i = iter; i < dim; i++) {
			matrix[iter][i] = 0;
		}
	}
	file.close();
}

void printAll(const methodRK& method) {
	cout << "A\n";
	for (int i = 0; i < method.s; i++) {
		for (int j = 0; j < method.s; j++) {
			cout << method.A[i][j] << " ";
		}
		cout << endl;
	}
	cout << "b1\n";
	for (int i = 0; i < method.s; i++) {
		cout << method.b[i] << " ";
	}
	cout << "\nc\n";
	for (int i = 0; i < method.s; i++) {
		cout << method.c[i] << " ";
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
	downloadMatrix(AFile, A, s);
	downloadArray(b1File, b1, s);
	downloadArray(b2File, b2, s);
	downloadArray(cFile, c, s);
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

	adaptiveNestedRK(nodes, values, method1, method2);
	outputResults(nodes, values, whereOutput + "AdaptiveRK.txt");
//	constRK(nodes, values, 314, method2);
//	outputResults(nodes, values, whereOutput + "ConstRK.txt");
	adaptiveSimpleRK(nodes, values, method2);
	outputResults(nodes, values, whereOutput + "AdaptiveSimpleRK.txt");
	
	clearAll();
	system("pause");
	return 0;
}