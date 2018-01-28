#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

// Declarations
class Node;
int nchoosek(int n, int k);
void hmg_pwlist(int nvar, int degree, Node **pwlist);
void hmg_pwlist_recursive(Node *init_vec, int size, int pos, Node **&pwlist, int &pwlist_size);
Node* find_idx(Node *head, int pos);
Node* find_create_idx(Node *head, int pos);
void insert_after(Node *current, Node *after);
Node* copy_node_list(Node *current);
void print_link_struct(Node *head);
void release_memory(Node *current);
int add_dense(Node *head1, Node *head2, int *locpow, int size);
int hmg_pow_index(int *monpower, int n, int deg);

// Functions

class Node {
public:
	Node(int idx, int deg) : index(idx), degree(deg), next(NULL) {} /* Initialization */
	int index;
	int degree;
	Node *next;
};


int nchoosek(int n, int k)
{
	if (k == 0 || k == n)
		return 1;

	if (n < k)
	{
		cout << "error: n=" << n << "k=" << k << "\n" << endl;
		exit(0);
	}

	double numerator = 1;
	for (int i = n - k + 1; i <= n; i++) numerator *= i;

	double denominator = 1;
	for (int i = 1; i <= k; i++) denominator *= i;

	if (numerator < 0)
	{
		cout << "overflow\n" << endl;
		exit(1);
	}

	return (int)(numerator / denominator);
}

void hmg_pwlist(int nvar, int degree, Node **pwlist)
{
	Node *init_vec = new Node(0, degree);
	pwlist[0] = init_vec;

	int pwlist_size = 1;
	hmg_pwlist_recursive(init_vec, nvar, 0, pwlist, pwlist_size);
}

void hmg_pwlist_recursive(Node *init_vec, int size, int pos, Node **&pwlist, int &pwlist_size)
{
	if (pos < size - 1)
	{
		Node *vec = find_idx(init_vec, pos);
		int ele = 0;
		if (vec->index == pos) ele = vec->degree;

		for (int i = ele - 1; i >= 0; i--)
		{
			Node *temp_vec = copy_node_list(init_vec);

			Node *node_pos = find_idx(temp_vec, pos);
			node_pos->degree = i;

			if (node_pos->degree == 0)
			{
				Node * pre_node = find_idx(temp_vec, pos - 1);
				if (pre_node == NULL)
					temp_vec = node_pos->next;
				else
					pre_node->next = node_pos->next;
				delete node_pos;
			}

			if (temp_vec != NULL)
			{
				Node *node_pos_plus_1 = find_create_idx(temp_vec, pos + 1);
				node_pos_plus_1->degree += (ele - i);
			}
			else
			{
				Node *node_pos_plus_1 = new Node(pos + 1, 0);
				node_pos_plus_1->degree += (ele - i);
				temp_vec = node_pos_plus_1;
			}

			pwlist[pwlist_size] = temp_vec;
			pwlist_size++;

			hmg_pwlist_recursive(temp_vec, size, pos + 1, pwlist, pwlist_size);
		}
	}
}


Node* find_idx(Node *head, int pos)
{
	Node *list = head;
	Node *pre_node = NULL;
	while (list != NULL && list->index <= pos)
	{
		pre_node = list;
		list = list->next;
	}

	return pre_node;
}

Node* find_create_idx(Node *head, int pos)
{
	Node *list = head;
	Node *node = find_idx(list, pos);
	if (node->index == pos)
		return node;
	else
	{
		Node *new_node = new Node(pos, 0);
		insert_after(node, new_node);
		return new_node;
	}
}

void insert_after(Node *current, Node *after)
{
	after->next = current->next;
	current->next = after;
}

Node* copy_node_list(Node *current)
{
	Node *list = current;
	Node *node = new Node(list->index, list->degree);
	Node *head = node;

	Node *tail = node;
	while (list->next != NULL)
	{
		list = list->next;
		node = new Node(list->index, list->degree);
		tail->next = node;
		tail = tail->next;
	}

	return head;
}

void release_memory(Node *current)
{
	while (current != NULL)
	{
		Node *p = current;
		current = current->next;
		delete p;
	}
}

void print_link_struct(Node *head)
{
	Node *p = head;
	while (p != NULL)
	{
		printf("%d:%d ", p->index, p->degree);
		p = p->next;
	}
	printf("\n");
}

int add_dense(Node *head1, Node *head2, int *locpow, int size)
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
	return degree;
}


int hmg_pow_index(int *monpower, int n, int deg)
{
	int *temp_monpower = new int[n + 1]; // to satisfy matlab code from 1 to n
	temp_monpower[0] = 0;
	for (int i = 1; i <= n; i++) temp_monpower[i] = monpower[i - 1];

	int indx = 0;
	vector<int> I;
	I.push_back(-1); // to satisfy matlab code from 1 to n

	for (int i = 1; i <= n; i++)
	{
		if (temp_monpower[i] > 0)
			I.push_back(i);
	}

	if (I.size() == 1) // since add a new entry to fill up value in the index 0 
	{
		indx = 1;
		return indx;
	}

	int kdg = deg;
	indx = 1;
	int Isize = I.size();
	for (int k = 1; k <= Isize - 1; k++)
	{
		if (k == 1 && I[1] > 1)
		{
			for (int i = 1; i <= kdg; i++) indx = indx + nchoosek(I[1] - 1 + i - 1, i) * nchoosek(n - I[1] + 1 + kdg - i - 1, kdg - i);
		}
		else if (k > 1 && I[k] > I[k - 1] + 1)
		{
			for (int i = 1; i <= kdg; i++) indx = indx + nchoosek(I[k] - I[k - 1] - 1 + i - 1, i) * nchoosek(n - I[k] + 1 + kdg - i - 1, kdg - i);
		}

		for (int i = temp_monpower[I[k]] + 1; i <= kdg; i++)
			indx = indx + nchoosek(n - I[k] + kdg - i - 1, kdg - i);

		kdg = kdg - temp_monpower[I[k]];
	}
	delete[]temp_monpower;

	return indx;
}