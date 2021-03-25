#include <iostream>
#include <fstream>
#include <string>
using namespace std;

const string AFile = "A.txt";
const string b1File = "b1.txt"; // more accurate
const string b2File = "b2.txt";
const string cFile = "c.txt";
const string dimenFile = "Dimen.txt";
const double E = 2.718281828459045;
const double PI = 3.14159265358979;
class myEx : exception{
public:
	myEx() {}
};

const double X0 = 0.0;
const double Y0 = 0.0;
const double X1 = PI / 4;
const double expected = 1;//pow(E, 0.5);//2;
const double EPS = 1e-10;
const int degree = 4;

double hMin = 1e-3;

double** A;
double* b1;
double* b2;
double* c;
double* kk;

int s;

double f(double x, double y) {
	return y * y + 1;
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
	kk = new double[s];
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
	delete[]kk;
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

void solution(double& solut, double& accurSolut) {
	double xCurr = X0;
	solut = Y0;
	const double deltaMax = 0.9;
	double totalTol = EPS / abs(X1 - X0);
	double step = (X1- X0) / 2;
	double localError = 0;
	while (abs(xCurr - X1) > DBL_EPSILON) {
		kFromStep(xCurr, solut, step, kk);
		localError = 0;
		for (int i = 0; i < s; i++) {
			localError += (b1[i] - b2[i]) * kk[i];
		}
// Here must be localError=step*sum(...), but as totalTol is also multiplyed by step
// then I can reduce it
//		localError *= step;
		if (abs(localError) < totalTol) {
			xCurr += step;
			solut = nextValue(xCurr, solut, step, b2, kk, true);
		}
		else {
			step *= min(pow(totalTol * step /abs(localError), 1. / degree), deltaMax);
		}
		if (abs(step) > abs(X1 - xCurr)) {
			step = X1 - xCurr;
		}
	}
	accurSolut = solut + localError;
}


int main() {
	initAll();
	downloadAll();
	double solut, accurSolut;
	solution(solut, accurSolut);
	cout << endl << abs(expected - solut) << endl << abs(expected - accurSolut) << endl;
	clearAll();
	system("pause");
	return 0;
}