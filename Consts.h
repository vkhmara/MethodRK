#pragma once
#include <string>
using namespace std;
const string AFile = "A.txt";
const string b1File = "b1.txt"; // more accurate
const string b2File = "b2.txt";
const string cFile = "c.txt";
const string dimenFile = "Dimen.txt";
const string whereOutput = "D:/”чЄба/ урсач/6 сем/Graphics/";
const double E = 2.718281828459045;
const double PI = 3.14159265358979;
const double deltaMax = 0.8;
const double deltaMin = 1.2;

const double X0 = 0.0;
const double Y0 = 0.0;
const double X1 = 50;
//const double expected = pow(E, 3.0) / 2 - 3.0;
const double EPS = 1e-13;
const int degree = 4;
const double a = 2;

double f(double x, double y) {
	return -y + cos(x / a) + sin(x / a);
}