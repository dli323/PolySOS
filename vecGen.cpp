#include <iostream>
#include "vecGen.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int degree;
    int mpw, npw;
    double *pw;       
    double *coeff;
    int len;
    
	double *bvec; 
    vector<double> vec;

	degree = (int)mxGetScalar(prhs[0]);
    pw =  mxGetPr(prhs[1]);
    mpw = (int)mxGetM(prhs[1]);
    npw = (int)mxGetN(prhs[1]);
    coeff = mxGetPr(prhs[2]);
    
    Matrix pw_in_c = Matrix(mpw, npw);
    
    for(int i = 0; i < mpw; i++)
    {
        for(int j = 0; j < npw; j++)
        {
            pw_in_c(i+1, j+1) = pw[i + mpw*j];
        }
    }
    
    vector_generator(degree, coeff, pw_in_c, vec);
    
    len = vec.size();
    
	plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL);
	bvec = mxGetPr(plhs[0]);
    
    for(int i = 0; i<len; i++)
	{
		bvec[i] = vec[i];
	}
    
}

bool isMon(Node *head1, Node *head2)
{
	Node *p1 = head1;
	Node *p2 = head2;
	
	while (p1 != NULL && p2 != NULL && p1->index == p2->index && p1->degree == p2->degree)
	{
		p1 = p1->next;
		p2 = p2->next;
	}
	return p1 == p2 ? 1:0;
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
			if (pw(i+1,j+1) > 0) 
			{ 
				nonzero++;
				if (nonzero == 1) { pwlist[i] = new Node(j, pw(i+1,j+1)); }
				else
				{
					node = find_create_idx(pwlist[i], j);
					node->degree = pw(i+1,j+1);
				}
			}
		}
	}
}


void vector_generator(int degree, double *coeff, Matrix pw, vector<double> &vec)
{
	int nmon = pw.GetRows();
	int nvar = pw.GetCols();

	Node **pwlist = new Node*[nmon];
	func_pwlist(pw, pwlist);
	
	int Nvec = nchoosek(nvar +  degree,  degree);
	Node **densepwlist = new Node*[Nvec];
	densepwlist[0] = NULL;

	int idx = 1;
	for (int i = 1; i <= degree; i++)
	{
		hmg_pwlist(nvar, i, &densepwlist[idx]);
		idx += nchoosek(nvar + i - 1, i);
	}

	for (int i = 1; i<Nvec; i++)
	{
		Node *onelist = densepwlist[i];
	}

	int flag;
	int i, j;

	for (i = 1; i<Nvec; i++)
	{
		flag = 0;
		for (j = 0; j < nmon; j++)
		{
			if (isMon(pwlist[j], densepwlist[i]))
			{
				flag = 1;
				break;
			}
		}
		vec.push_back(coeff[j]*flag);
	}

	for (int i = 0; i<Nvec; i++)
	{
		release_memory(densepwlist[i]);
	}

	delete[]densepwlist;

	for (int i = 0; i<nmon; i++)
	{
		release_memory(pwlist[i]);
	}

	delete[]pwlist;

}