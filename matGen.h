#include <iostream>
#include <vector>
#include "utility.h"
#include "mex.h"

using namespace std;

void moment_generator(int nvar, int degree, vector<int> &row, vector<int> &col);
void sparse_moment_generator(int nvar, int degree, vector<int> &varset, vector<int> &row, vector<int> &col);

