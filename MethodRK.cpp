#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

const string AFile = "A.txt";
const string b1File = "b1.txt"; // more accurate
const string b2File = "b2.txt";
const string cFile = "c.txt";
const string dimenFile = "Dimen.txt";
const string whereOutput = "D:/”чЄба/ урсач/6 сем/Graphics/";
const double E = 2.718281828459045;
const double PI = 3.14159265358979;
const double deltaMax = 0.95;

const double X0 = 0.0;
const double Y0 = 3.0;
const double X1 = 4;
//const double expected = pow(E, 3.0) / 2 - 3.0;
const double EPS = 1e-5;
const int degree = 4;

double** A;
double* b1;
double* b2;
double* c;

int s;

double f(double x, double y) {
	return -y + x * x;
//	return y;
}

double getNumberFromString(string str) {
	int sleshIndex = str.find('/');
	if (sleshIndex == string::npos) {
		return atof(str.c_str());
	}
	return atof(str.substr(0, sleshIndex).c_str()) /
		atof(str.substr(sleshIndex + 1).c_str());
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
}

void downloadAll() {
	string content;
	int prevPos = -1;

	ifstream file(AFile);
	char ch[200];
	for (int i = 0; i < s; i++) {
		A[0][i] = 0;
	}
	for (int iter = 1; iter < s; iter++) {
		prevPos = -1;
		file.getline(ch, 200);
		content = ch;
		content += " ";
		for (int i = 0; i < iter; i++) {
			A[iter][i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
			prevPos = content.find(' ', prevPos + 1);
		}
		for (int i = iter; i < s; i++) {
			A[iter][i] = 0;
		}
	}
	file.close();

	prevPos = -1;
	file.open(b1File);
	file.getline(ch, 200);
	content = ch;
	content += " ";
	for (int i = 0; i < s; i++) {
		b1[i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
		prevPos = content.find(' ', prevPos + 1);
	}
	file.close();

	prevPos = -1;
	file.open(b2File);
	file.getline(ch, 200);
	content = ch;
	content += " ";
	for (int i = 0; i < s; i++) {
		b2[i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
		prevPos = content.find(' ', prevPos + 1);
	}
	file.close();

	prevPos = -1;
	file.open(cFile);
	file.getline(ch, 200);
	content = ch;
	content += " ";
	for (int i = 0; i < s; i++) {
		c[i] = getNumberFromString(content.substr(prevPos + 1, content.find(' ', prevPos + 1)));
		prevPos = content.find(' ', prevPos + 1);
	}
	file.close();
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

void kFromStep(double xValue, double yValue, double step, double* k) {
	double sum = 0;
	for (int i = 0; i < s; i++) {
		sum = 0;
		for (int t = 0; t < i; t++) {
			sum += A[i][t] * k[t];
		}
		k[i] = f(xValue + c[i] * step, yValue + step * sum);
	}
}

double nextValue(double xValue, double yValue, double step, double* b, double* k, bool calculatedK=false) {
	if (!calculatedK) {
		kFromStep(xValue, yValue, step, k);
	}
	double sum = 0;
	for (int i = 0; i < s; i++) {
		sum += b[i] * k[i];
	}
	return yValue + step * sum;
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
	while (abs(xCurr - X1) > DBL_EPSILON) {
		kFromStep(xCurr, solut, step, k);
		localError = 0;
		for (int i = 0; i < s; i++) {
			localError += (b1[i] - b2[i]) * k[i];
		}
		localError *= step;
		if (step < 1e-5) {
			cout << "Stop" << endl;
		}
		if (abs(localError) / step < totalTol) {
			xCurr += step;
			solut = nextValue(xCurr, solut, step, b2, k, true);
			nodes.push_back(xCurr);
			values.push_back(solut);
		}
		else {
			step *= min(pow(totalTol * step /abs(localError), 1. / degree), deltaMax);
		}
		if (abs(step) > abs(X1 - xCurr)) {
			step = X1 - xCurr;
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
		yCurr = nextValue(xCurr, yCurr, step, b2, k);
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
	constRK(nodes, values, 400);
	outputResults(nodes, values, whereOutput + "ConstRK.txt");
	
	clearAll();
	system("pause");
	return 0;
}