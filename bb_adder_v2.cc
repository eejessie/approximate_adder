#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <sys/timeb.h>
#include <time.h>
#include "header/helper.h"

using namespace std;


int compute_G_v2(int size, unsigned long long add1, unsigned long long add2)
{
	int G = add1 & add2;
	if (G == 0)
		return 0;
	int P = add1 ^ add2;
	vector<int> G_vec, P_vec;
	d2b(G, G_vec, size);
	d2b(P, P_vec, size);
/*
	cout << "G_vec: ";
	for (int i = 0; i < G_vec.size(); i++)
		cout << G_vec[i];
	cout << endl << "P_vec: ";
	for (int i = 0; i < P_vec.size(); i++)
		cout << P_vec[i];
	cout << endl;
*/	
	int group_G = 0;
	for (int i = size-1; i > 0; i--)
	{
		int prod_gp = 1;
		int prod_p = 1;
		for (int j = i-1; j >= 0; j--)
			prod_p *= P_vec[j];
		prod_gp = G_vec[i] * prod_p;
		group_G += prod_gp;
	} 
	group_G += G_vec[0];
	if (group_G)
		return 1;
	else
		return 0;
}


int compute_block_pgk_v2(int m, int k, int l, int index, unsigned long long add1_low, unsigned long long add1_high, 
						 unsigned long long add2_low, unsigned long long add2_high)
{
//	cout << endl << "in compute_block_pgk, index = " << index << endl;

	unsigned long long add1_tmp, add2_tmp;
	
	int size = l;
	unsigned int bit_num;
		
	if (index-1 <= m/2-1)
	{	
	//	cout << "case1. " << endl;
		int num_shift = index * k - l;
		if (num_shift < 0)
		{
			num_shift = 0;
			size = index * k;
		}
		add1_tmp = add1_low >> num_shift;
		add2_tmp = add2_low >> num_shift;
		
		bit_num = pow(2.0, size) - 1;
		add1_tmp &= bit_num;
		add2_tmp &= bit_num;		
	}
	else if ( (index-1 > m/2-1) && (l - (index-1-(m/2-1))*k) > 0)
	{
	//	cout << "case2. " << endl;
		int num_bit_high = (index-m/2)*k;
		int num_bit_low = l - num_bit_high;
		
		unsigned long long add1_low_tmp, add1_high_tmp, add2_low_tmp, add2_high_tmp;
		
		unsigned int bit_num_high = pow(2.0, num_bit_high) - 1;
		unsigned int bit_num_low = pow(2.0, num_bit_low) - 1;
		
	//	cout << "num_bit_high = " << num_bit_high << ", num_bit_low = " << num_bit_low << endl;
		
		add1_low_tmp = add1_low >> (m*k/2 - num_bit_low);	
		add1_low_tmp &= bit_num_low;	
		add1_high_tmp = add1_high & bit_num_high;
		add1_tmp = add1_high_tmp * pow(2.0, num_bit_low) + add1_low_tmp;
		
		add2_low_tmp = add2_low >> (m*k/2 - num_bit_low);
		add2_low_tmp &= bit_num_low;
		add2_high_tmp = add2_high & bit_num_high;
		add2_tmp = add2_high_tmp * pow(2.0, num_bit_low) + add2_low_tmp;
		
	//	cout << "add1_tmp = " << add1_tmp << ", add2_tmp = " << add2_tmp << endl;
		
		bit_num = pow(2.0, size) - 1;
		
	}
	else
	{
	//	cout << "case3. " << endl;
		int num_shift = (index-m/2)*k - l;
		add1_tmp = add1_high >> num_shift;
		add2_tmp = add2_high >> num_shift;
		
		bit_num = pow(2.0, size) - 1;
		add1_tmp &= bit_num;
		add2_tmp &= bit_num;
	}
	
	unsigned long long p_block = add1_tmp ^ add2_tmp;
	if (p_block == bit_num)                             //group_P = 1
		return 1;                
	else
	{
		int res = compute_G_v2(size, add1_tmp, add2_tmp);
		if (res)
			return 2;                                  //group_G = 1
		else
			return 3;                                  //group_K = 1
	}
	
}


void bb_adder_v2(int n, int k, int l, unsigned long long add1_low, unsigned long long add1_high, unsigned long long add2_low, unsigned long long add2_high, unsigned long long &sum_est_low, unsigned long long &sum_est_high)
{
//	cout << endl << "in bb_adder: " << endl;
	int m = n/k;
	int t = l/k;
	unsigned int *sum_block = new unsigned int[m];
	int m1 = n/(2*k);
	int m2 = m - m1;
	int ci = 0, co;

	unsigned int bit_num = pow(2.0, k)-1;
	unsigned long long add1_block, add2_block;
	
	//step1. compute the sum for each block in the lower part
	for (int i = 0; i < m/2; i++)
	{	
		int shift_bit = i*k;
		add1_block = add1_low >> shift_bit;
		add2_block = add2_low >> shift_bit;
		add1_block &= bit_num;
		add2_block &= bit_num;
		
		//compute sum 
		sum_block[i] = add1_block + add2_block + ci;
		if (sum_block[i] > bit_num) sum_block[i] -= pow(2.0, k);
		
	//	cout << "i = " << i << endl;
	//	cout << "add1_block = " << add1_block << ", add2_block = " << add2_block << ", ci = " << ci << ", sum_block = " << sum_block[i] << endl;

		//compute co
		int res = compute_block_pgk_v2(m, k, l, i+1, add1_low, add1_high, add2_low, add2_high);
		if (res == 1)      co = 0;
		else if (res == 2) co = 1;
		else               co = 0;
		
		ci = co;
	//	cout << "co = " << co << endl;
	}
	
	//step2. compute the sum for each block in the higher part
	for (int i = m/2; i < m; i++)
	{	
		int shift_bit = (i-m/2)*k;
		add1_block = add1_high >> shift_bit;
		add2_block = add2_high >> shift_bit;
		add1_block &= bit_num;
		add2_block &= bit_num;
		
		//compute sum 
		sum_block[i] = add1_block + add2_block + ci;
		if (i < m-1 && sum_block[i] > bit_num) sum_block[i] -= pow(2.0, k);
		
	//	cout << "i = " << i << endl;
	//	cout << "add1_block = " << add1_block << ", add2_block = " << add2_block << ", ci = " << ci << ", sum_block = " << sum_block[i] << endl;

		//compute co
		int res = compute_block_pgk_v2(m, k, l, i+1, add1_low, add1_high, add2_low, add2_high);
		if (res == 1)      co = 0;
		else if (res == 2) co = 1;
		else               co = 0;
		
		ci = co;
	//	cout << "co = " << co << endl;
	}
		
	cout << "sum_block: " << endl;
	for (int i = m-1; i >= 0; i--)
		cout << sum_block[i] << " ";
	cout << endl;	
		
	unsigned int max_low = 1;
	for (int i = 0; i < n/2-1; i++)
		max_low = max_low * 2 + 1;
		
	for (int i = 0; i < m/2; i++)
		sum_est_low += sum_block[i] * pow(2.0, i*k);
	if (sum_est_low > max_low)
		sum_est_low -= max_low;

	for (int i = m/2; i < m; i++)
		sum_est_high += sum_block[i] * pow(2.0, (i-m/2)*k);

	delete []sum_block;

}
