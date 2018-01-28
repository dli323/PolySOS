#include <iostream>
#include <algorithm>
#include "matGen.h"
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n, degree, size;
	double *row, *col;
    
    vector<int> Arow;
	vector<int> Acol;

	n = (int)mxGetScalar(prhs[0]);
	degree = (int)mxGetScalar(prhs[1]);

	if(nrhs == 2)
    {
        moment_generator(n, degree, Arow, Acol);
    } 
    else 
    {
        double *vars;
        int nvars;
        vector<int> varset;
        vars = mxGetPr(prhs[2]);
        nvars = (int)mxGetN(prhs[2]);
        
        for(int i = 0;i < nvars; i++)
	    {
            varset.push_back((int)vars[i]);
	    }
        sparse_moment_generator(n, degree, varset, Arow, Acol);
    }

	size = Arow.size();

	plhs[0] = mxCreateDoubleMatrix(size,1,mxREAL);
	row = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(size,1,mxREAL);
	col = mxGetPr(plhs[1]);

	for(int i = 0;i<size;i++)
	{
		row[i] = Arow[i];
		col[i] = Acol[i];
	}
}

void moment_generator(int nvar, int degree, vector<int> &row, vector<int> &col)
{
	int Nmat = nchoosek(nvar + degree, degree);

	Node ** pwlist = new Node*[Nmat];
	pwlist[0] = NULL;

	int idx = 1;
	for (int i = 1; i <= degree; i++)
	{
		hmg_pwlist(nvar, i, &pwlist[idx]);
		idx += nchoosek(nvar + i - 1, i);
	}

	int Nvec = nchoosek(nvar + 2 * degree, 2 * degree);

	int *locpow = new int[nvar];
	int sumlocp = 0;
	for (int i = 1; i <= Nmat; i++)
	{
		for (int j = i; j <= Nmat; j++)
		{
			sumlocp = add_dense(pwlist[i - 1], pwlist[j - 1], locpow, nvar);
			int mon_2d_idx = 1;
			if (sumlocp != 0)
				mon_2d_idx = nchoosek(nvar + sumlocp - 1, sumlocp - 1) + hmg_pow_index(locpow, nvar, sumlocp);

			/* assign the (row, col, 1) to sparse matrix*/
			row.push_back(mon_2d_idx);
			col.push_back(i + (j - 1)*Nmat);

			if (i < j)
			{
				row.push_back(mon_2d_idx);
				col.push_back(j + (i - 1)*Nmat);
			}
		}
	}

	delete[]locpow;
	for (int i = 0; i<Nmat; i++) release_memory(pwlist[i]);
	delete[]pwlist;
}

void sparse_moment_generator(int nvar, int degree, vector<int> &varset, vector<int> &row, vector<int> &col)
{
	/* Dense monomial */
	int Nvec = nchoosek(nvar + 2 * degree, 2 * degree);
	int Nmat = nchoosek(nvar + degree, degree);
	Node **pwlist = new Node*[Nmat];
	pwlist[0] = NULL;

	int idx = 1;
	for (int i = 1; i <= degree; i++)
	{
		hmg_pwlist(nvar, i, &pwlist[idx]);
		idx += nchoosek(nvar + i - 1, i);
	}

	/* monomial location in sparse matrix */
	int Nvarsp = varset.size();
	int Nmatsp = nchoosek(Nvarsp + degree, degree);

	vector<int> locsp; /* location in sparse vector */
	locsp.assign(Nmat, -1);
	int loc = 1;
	locsp.at(0) = loc; /* constant term */

	for (int i = 1; i < Nmat; i++)
	{
		Node *onelist = pwlist[i];
		bool flag = 1;
		while (onelist != NULL)
		{
			if (find(varset.begin(), varset.end(), onelist->index) == varset.end())
			{
				flag = 0;
				break;
			}
			onelist = onelist->next;
		}
		if (flag)
		{
			loc++;
			locsp.at(i) = loc;
		}
	}

	/* construct the 'sparse' selection matrix */
	int *locpow = new int[nvar];
	int sumlocp = 0;
	for (int i = 1; i <= Nmat; i++)
	{
		for (int j = i; j <= Nmat; j++)
		{
			sumlocp = add_dense(pwlist[i - 1], pwlist[j - 1], locpow, nvar);
			int mon_2d_idx = 1;
			if (sumlocp != 0)
				mon_2d_idx = nchoosek(nvar + sumlocp - 1, sumlocp - 1) + hmg_pow_index(locpow, nvar, sumlocp);

			/* assign the (row, col, 1) to sparse matrix*/
			if ((locsp.at(i - 1) > 0) && (locsp.at(j - 1) > 0))
			{
				row.push_back(mon_2d_idx);
				col.push_back(locsp.at(i - 1) + (locsp.at(j - 1) - 1)*Nmatsp);

				if (i < j)
				{
					row.push_back(mon_2d_idx);
					col.push_back(locsp.at(j - 1) + (locsp.at(i - 1) - 1)*Nmatsp);
				}
			}
			
		}
	}

	delete[]locpow;
	for (int i = 0; i<Nmat; i++) release_memory(pwlist[i]);
	delete[]pwlist;
}