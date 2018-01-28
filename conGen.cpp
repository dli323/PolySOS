#include <iostream>
#include <math.h>
#include "conGen.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nvar;
	int degree;
  int mpw, npw;
  double *pw;       
  double *cf;
	vector<double> coeff;
	double cs;
    
	vector<int> Brow;
	vector<int> Bcol;
	vector<double> Bval;
  double *row, *col, *val;

	nvar = (int)mxGetScalar(prhs[0]);
	degree = (int)mxGetScalar(prhs[1]);
  pw =  mxGetPr(prhs[2]);
  mpw = (int)mxGetM(prhs[2]);
  npw = (int)mxGetN(prhs[2]);
  cf = mxGetPr(prhs[3]);
	cs = mxGetScalar(prhs[4]);

	for (int i = 0; i < mpw; i++)
		coeff.push_back(cf[i]);
    
    Matrix pw_in_c = Matrix(mpw, npw);
    
    for(int i = 0; i < mpw; i++)
    {
        for(int j = 0; j < npw; j++)
        {
            pw_in_c(i+1, j+1) = pw[i + mpw*j];
        }
    }
    
	ConInfo cstr(coeff, cs, pw_in_c);
	con_moment_generator(nvar, degree, cstr, Brow, Bcol, Bval);
       
	int size = Brow.size();

	plhs[0] = mxCreateDoubleMatrix(size, 1, mxREAL);
	row = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(size, 1, mxREAL);
	col = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(size, 1, mxREAL);
	val = mxGetPr(plhs[2]);

	for (int i = 0; i < size; i++)
	{
		row[i] = Brow[i];
		col[i] = Bcol[i];
		val[i] = Bval[i];
	}
      
}

int add_cstr(Node *head1, Node *head2, Node *head3, int *locpow, int size)
{
	int degree = 0;

	for (int i = 0; i < size; i++) locpow[i] = 0;

	Node*p = head1;
	while (p != NULL)
	{
		locpow[p->index] += p->degree;
		degree += p->degree;
		p = p->next;
	}

	p = head2;
	while (p != NULL)
	{
		locpow[p->index] += p->degree;
		degree += p->degree;
		p = p->next;
	}

	p = head3;
	while (p != NULL)
	{
		locpow[p->index] += p->degree;
		degree += p->degree;
		p = p->next;
	}

	return degree;
}

void func_pwlist(Matrix pw, Node **pwlist)
{
	int nmon = pw.GetRows();
	int nvar = pw.GetCols();

	int nonzero = 0;
	Node *node;
	for (int i = 0; i < nmon; i++)
	{
		nonzero = 0;
		for (int j = 0; j < nvar; j++)
		{
			if (pw(i + 1, j + 1) > 0)
			{
				nonzero++;
				if (nonzero == 1) { pwlist[i] = new Node(j, pw(i + 1, j + 1)); }
				else
				{
					node = find_create_idx(pwlist[i], j);
					node->degree = pw(i + 1, j + 1);
				}
			}
		}
	}
}

double max_degree(Matrix pw)
{
	int nmon = pw.GetRows();
	int nvar = pw.GetCols();
	double maxdeg, mondeg;

	maxdeg = 0;
	for (int i = 1; i <= nmon; i++)
	{
		mondeg = 0;
		for (int j = 1; j <= nvar; j++)
		{
			mondeg = mondeg + pw(i, j);
		}
		if (mondeg >= maxdeg)
		{
			maxdeg = mondeg;
		}
	}

	return maxdeg;
}

void con_moment_generator(int nvar, int degree, ConInfo cstr, vector<int> &row, vector<int> &col, vector<double> &val)
{
	/* Dense monomial */
	int Nvec = nchoosek(nvar + 2 * degree, 2 * degree);

	double maxdeg;
	maxdeg = max_degree(cstr.pw);
	int redeg = degree - ceil(maxdeg / 2);

	/* Remaining monomial */
	int Nmat = nchoosek(nvar + redeg, redeg);
	Node ** pwlist = new Node*[Nmat];
	pwlist[0] = NULL;
	int idx = 1;
	for (int i = 1; i <= redeg; i++)
	{
		hmg_pwlist(nvar, i, &pwlist[idx]);
		idx += nchoosek(nvar + i - 1, i);
	}

	/* Constraint's monomial */
	int Nmon = cstr.pw.GetRows();
	Node **cstrpwlist = new Node*[Nmon];
	func_pwlist(cstr.pw, cstrpwlist);

	int *locpow = new int[nvar];
	int sumlocp = 0;
	for (int i = 1; i <= Nmat; i++)
	{
		for (int j = i; j <= Nmat; j++)
		{
			for (int k = 0; k <= Nmon; k++)
			{
				if (k == 0 && cstr.cs != 0)
				{
					sumlocp = add_dense(pwlist[i - 1], pwlist[j - 1], locpow, nvar);
					int mon_2d_idx = 1;
					if (sumlocp != 0)
						mon_2d_idx = nchoosek(nvar + sumlocp - 1, sumlocp - 1) + hmg_pow_index(locpow, nvar, sumlocp);

					/* assign the (row, col, cstr.cs) to sparse matrix*/
					row.push_back(mon_2d_idx);
					col.push_back(i + (j - 1)*Nmat);
					val.push_back(cstr.cs);

					if (i < j)
					{
						row.push_back(mon_2d_idx);
						col.push_back(j + (i - 1)*Nmat);
						val.push_back(cstr.cs);
					}
				}
				else if (k > 0)
				{
					sumlocp = add_cstr(pwlist[i - 1], pwlist[j - 1], cstrpwlist[k - 1], locpow, nvar);
					int mon_2d_idx = 1;
					if (sumlocp != 0)
						mon_2d_idx = nchoosek(nvar + sumlocp - 1, sumlocp - 1) + hmg_pow_index(locpow, nvar, sumlocp);

					/* assign the (row, col, cstr.cf) to sparse matrix*/
					row.push_back(mon_2d_idx);
					col.push_back(i + (j - 1)*Nmat);
					val.push_back(cstr.cf.at(k - 1));

					if (i < j)
					{
						row.push_back(mon_2d_idx);
						col.push_back(j + (i - 1)*Nmat);
						val.push_back(cstr.cf.at(k - 1));
					}
				}
			}
		}
	}

	free(locpow);
	for (int i = 0; i < Nmat; i++)
	{
		release_memory(pwlist[i]);
	}
	delete[]pwlist;
	for (int i = 0; i < Nmon; i++)
	{
		release_memory(cstrpwlist[i]);
	}
	delete[]cstrpwlist;

}