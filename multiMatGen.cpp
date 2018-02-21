#include <iostream>
#include "multiMatGen.h"
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n, degree, nbasis;
	double *basis, *idx;

	vector<int> monb;
	Matrix monidx;

	n = (int)mxGetScalar(prhs[0]);
	degree = (int)mxGetScalar(prhs[1]);
	basis = mxGetPr(prhs[2]);
	nbasis = (int)mxGetN(prhs[2]);

	for (int i = 0; i < nbasis; i++)
	{
		monb.push_back((int)basis[i]);
	}

	multiplication_matrix_generator(n, degree, monb, monidx);

	plhs[0] = mxCreateDoubleMatrix(n, nbasis, mxREAL);
	idx = mxGetPr(plhs[0]);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < nbasis; j++)
			idx[i + n*j] = monidx(i + 1, j + 1);
	}
}


void multiplication_matrix_generator(int nvar, int degree, vector<int> monb, Matrix &monidx)
{
	int Nmat = nchoosek(nvar + degree, degree);
	Node **pwlist = new Node*[Nmat];
	pwlist[0] = NULL;

	int idx = 1;
	for (int i = 1; i <= degree; i++)
	{
		hmg_pwlist(nvar, i, &pwlist[idx]);
		idx += nchoosek(nvar + i - 1, i);
	}

	int Nb = monb.size();
	Node **blist = new Node*[Nb];

	for (int i = 0; i < Nb; i++)
		blist[i] = pwlist[monb.at(i) - 1];  /* -1 is for c index */

	int *locpow = new int[nvar];
	int sumlocp = 0;
	monidx = Matrix(nvar, Nb);
	for (int i = 1; i <= nvar; i++)
	{
		for (int j = 1; j <= Nb; j++)
		{
			sumlocp = add_dense(pwlist[i], blist[j - 1], locpow, nvar);
			monidx(i, j) = 1;
			if (sumlocp != 0)
				monidx(i, j) = nchoosek(nvar + sumlocp - 1, sumlocp - 1) + hmg_pow_index(locpow, nvar, sumlocp);
		}
	}
}