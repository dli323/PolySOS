#include <iostream>
#include <vector>
#include "matrixcpp.h"
#include "utility.h"
#include "mex.h"

using namespace std;

class ConInfo {
public:
	ConInfo(vector<double> cfv, double csv, Matrix pwv) : cf(cfv), cs(csv), pw(pwv) {} /* Initialization */
	vector<double> cf;
	double cs;
	Matrix pw;
};

int add_cstr(Node *head1, Node *head2, Node *head3, int *locpow, int size);
void func_pwlist(Matrix pw, Node **pwlist);
double max_degree(Matrix pw);
void con_moment_generator(int nvar, int degree, ConInfo cstr, vector<int> &row, vector<int> &col, vector<double> &val);

