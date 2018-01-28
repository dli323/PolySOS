#include <iostream>
#include <vector>
#include "matrixcpp.h"
#include "utility.h"
#include "mex.h"

using namespace std;

bool isMon(Node *head1, Node *head2);
void func_pwlist(Matrix pw, Node **pwlist);
void vector_generator(int degree, double *coeff, Matrix pw, vector<double> &vec);

