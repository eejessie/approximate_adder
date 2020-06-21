#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <sys/timeb.h>
#include <time.h>

using namespace std;


void d2b(int d, vector<int> &bin, int numBit)
{
	int mod = 0;
	vector<int> tmpstr;
	
	if(d == 0)
	{
		for(int i = 0; i < numBit; i++)
			bin.push_back(0);
		return;	
	}

	while(d > 0)
	{
		mod = d%2;
		d /= 2;
		tmpstr.push_back(mod);
	}
	unsigned int len = tmpstr.size();
	int minus = numBit - len;
	if(minus > 0)
		for(int i = 0; i < minus; i++)
			bin.push_back(0);	
	for(int i = len - 1; i >= 0; i--)
		bin.push_back(tmpstr[i]);
}


void d2b_ull(unsigned long long d, vector<int> &bin, int numBit)
{
	int mod = 0;
	vector<int> tmpstr;
	
	if(d == 0)
	{
		for(int i = 0; i < numBit; i++)
			bin.push_back(0);
		return;	
	}

	while(d > 0)
	{
		mod = d%2;
		d /= 2;
		tmpstr.push_back(mod);
	}
	unsigned int len = tmpstr.size();
	int minus = numBit - len;
	if(minus > 0)
		for(int i = 0; i < minus; i++)
			bin.push_back(0);	
	for(int i = len - 1; i >= 0; i--)
		bin.push_back(tmpstr[i]);
}



void obtain_block_dist(int n, int k, vector<double> &input_dist, vector<double> &P_prob, vector<double> &G_prob, vector<double> &K_prob)
{
	map<int, double>::iterator itrm_id;
	struct timeb st, et;
	ftime(&st);

	//step2. get the distribution for each block from input_dist
	int m = n/k;
        vector<map<int, double> > block_dist;
        for (int i = 0; i < m; i++)
	{
		cout << endl << "#block " << i << endl;
		map<int, double> this_block_dist;
		for (int j = 0; j < pow(2.0, k); j++)
		{
			cout << "j = " << j << endl;
			double prob = 0;
			if (i == 0)
			{
				for (unsigned int r = 0; r < pow(2.0, n-k); r++)
				{
					unsigned rr = r << k;
					unsigned int value = rr + j; 
			//		cout << "value = " << value << ", ";
					double this_prob = input_dist[value];
		       //		cout << "this_prob = " << this_prob << endl;
					prob += this_prob; 
				}
		//		cout << "j = " << j << ", prob = " << prob << endl;
			}
			else
			{
				double this_prob;
				for (int s = 0; s < pow(2.0, i*k); s++)
				{
					int value;
					if (m-1-i > 0)
					{
						unsigned int mid_value = j << (k*i);	
						for (unsigned int q = 0; q < pow(2.0, k*(m-1-i)); q++)
						{
							int shift_bits = k * (i+1);
							unsigned int qq = q << shift_bits;	
						 	value = s + mid_value + qq; 		
							this_prob = input_dist[value];
							prob += this_prob; 
						}
					}
					else if (m-1-i == 0)
					{
						unsigned int left_value = j << (k*i);
						value = s + left_value;
						this_prob = input_dist[value];
						prob += this_prob; 
					}
				}
			}
			this_block_dist.insert(pair<int, double>(j, prob));
		}
       	//	cout << "this_block_dist: " << endl;
		double sum_block = 0;
		for (itrm_id = this_block_dist.begin(); itrm_id != this_block_dist.end(); itrm_id++)
		{
	//		cout << itrm_id->first << ": " <<  itrm_id->second  << endl;
			sum_block += itrm_id->second;
		}
		cout << "sum_block = " << sum_block << endl;
		block_dist.push_back(this_block_dist);
	}		
	
	//step3. obtain the distribution for the sum of each block by computing convolution
	vector<map<int, double> > sum_block_dist; 
	int full_value = pow(2.0, k) - 1;
	for (int i = 0; i < m; i++)
	{
		map<int, double> &cur_block_dist = block_dist[i];
		cout << "block = " << i << endl;
		double prob_P = 0;
		for (int r = 0; r <= full_value; r++)
		{
			itrm_id = cur_block_dist.find(r);
			double prob1 = itrm_id->second;
			itrm_id = cur_block_dist.find(full_value-r);
			double prob2 = itrm_id->second;
			prob_P += prob1 * prob2;
		}		
		P_prob.push_back(prob_P);

		double prob_G = 0;
		double prob, prob1, prob2;
		for (int j = pow(2.0, k); j < pow(2.0, k+1)-1; j++)
		{
			prob = 0;
			for (int r = 0; r <= j; r++)
			{
				itrm_id = cur_block_dist.find(r);
				if (itrm_id == cur_block_dist.end())
				{
					prob1 = 0;
					continue;
				}
				else
					prob1 = itrm_id->second;

				itrm_id = cur_block_dist.find(j-r);
				if (itrm_id == cur_block_dist.end())
				{
					prob2 = 0;
					continue;
				}
				else
					prob2 = itrm_id->second;
				prob += prob1 * prob2;
			}		
			prob_G += prob;
		}
		G_prob.push_back(prob_G);
		K_prob.push_back(1-prob_P-prob_G);
        } 

	ftime(&et);
    	double runtime = ((et.time - st.time)*1000 + (et.millitm - st.millitm))/1000.0;
	cout << "runtime for obtaining block distribution: " << runtime << endl;
}



void obtain_block_dist_GL(int n, int k, int kk, vector<double> &input_dist, vector<double> &P_prob, vector<double> &G_prob, vector<double> &K_prob, vector<double> &GL_prob)
{
	int m = n/k;
	obtain_block_dist(n, k, input_dist, P_prob, G_prob, K_prob);

        vector<map<int, double> > block_dist_left;
	map<int, double>::iterator itrm_id;
        for (int i = 0; i < m; i++)
	{
		cout << endl << "#block " << i << endl;
		map<int, double> this_block_dist_left;
		//obtain this_block_dist_left
		int kr = k - kk;
		for (int j = 0; j < pow(2.0, kk); j++)
		{
	//		cout << "j = " << j << endl;
			double prob = 0, this_prob;
			for (int s = 0; s < pow(2.0, k*i+kr); s++)
			{
				int value;
				if (m-1-i > 0)
				{
					unsigned int mid_value = j << (k*i+kr);	
					for (unsigned int q = 0; q < pow(2.0, k*(m-1-i)); q++)
					{
						int shift_bits = k * (i+1);
						unsigned int qq = q << shift_bits;	
					 	value = s + mid_value + qq; 		
						this_prob = input_dist[value];
				//		cout << "value = " << value << ", this_prob = " << this_prob << endl;
						prob += this_prob; 
					}
				}
				else if (m-1-i == 0)
				{
					unsigned int left_value = j << (k*i+kr);
					value = s + left_value;
					this_prob = input_dist[value];
			//		cout << "value = " << value << ", this_prob = " << this_prob << endl;
					prob += this_prob; 
				}
			}
			this_block_dist_left.insert(pair<int, double>(j, prob));
		}
		cout << "this_block_dist_left: " << endl;
		double sum_block_left = 0;
		for (itrm_id = this_block_dist_left.begin(); itrm_id != this_block_dist_left.end(); itrm_id++)
		{
			cout << itrm_id->first << ": " <<  itrm_id->second  << endl;
			sum_block_left += itrm_id->second;
		}
//		cout << "sum_block_left = " << sum_block_left << endl;
		block_dist_left.push_back(this_block_dist_left);
	}		
	
	//step2. obtain the distribution for the sum of each block by computing convolution
	for (int i = 0; i < m; i++)
	{
		//obtain GL_prob
		map<int, double> cur_block_dist_left = block_dist_left[i];
		double prob_GL = 0;
		for (int j = 0; j < pow(2.0, kk+1)-1; j++)
		{
			int value = j;
			double prob = 0;
			for (int r = 0; r <= j; r++)
			{
				int x = r;
				int y = j - r;
				itrm_id = cur_block_dist_left.find(x);
				double prob1 = itrm_id->second;
				itrm_id = cur_block_dist_left.find(y);
				double prob2 = itrm_id->second;
				prob += prob1 * prob2;
			}		
			if (value > pow(2.0, kk) - 1)
				prob_GL += prob;
		}
		GL_prob.push_back(prob_GL);
        } 
}


/*void obtain_block_sep_dist(int n, int k, int kk, vector<double> &input_dist, vector<double> &PL_prob, vector<double> &PR_prob, vector<double> &GL_prob, vector<double> &GR_prob, vector<double> &KL_prob, vector<double> &KR_prob)
{
	map<int, double>::iterator itrm_id;

	//step2. get the distribution for each block from input_dist
	int m = n/k;
        vector<map<int, double> > block_dist_left, block_dist_right;
        for (int i = 0; i < m; i++)
	{
		cout << endl << "#block " << i << endl;
		map<int, double> this_block_dist_left, this_block_dist_right;
		int kr = k - kk;
		//obtain this_block_dist_left
		for (int j = 0; j < pow(2.0, kk); j++)
		{
			cout << "j = " << j << endl;
			double prob = 0;
			double this_prob;
			for (int s = 0; s < pow(2.0, i*k+kr); s++)
			{
				int value;
				if (m-1-i > 0)
				{
					unsigned int mid_value = j << (k*i+kr);	
					for (unsigned int q = 0; q < pow(2.0, k*(m-1-i)); q++)
					{
						int shift_bits = k * (i+1);
						unsigned int qq = q << shift_bits;	
					 	value = s + mid_value + qq; 		
						this_prob = input_dist[value];
						cout << "value = " << value << ", this_prob = " << this_prob << endl;
						prob += this_prob; 
					}
				}
				else if (m-1-i == 0)
				{
					unsigned int left_value = j << (k*i+kr);
					value = s + left_value;
					this_prob = input_dist[value];
					cout << "value = " << value << ", this_prob = " << this_prob << endl;
					prob += this_prob; 
				}
			}
			this_block_dist_left.insert(pair<int, double>(j, prob));
		}
		cout << "this_block_dist_left: " << endl;
		double sum_block_left = 0;
		for (itrm_id = this_block_dist_left.begin(); itrm_id != this_block_dist_left.end(); itrm_id++)
		{
			cout << itrm_id->first << ": " <<  itrm_id->second  << endl;
			sum_block_left += itrm_id->second;
		}
		cout << "sum_block_left = " << sum_block_left << endl;
		block_dist_left.push_back(this_block_dist_left);

		
		//obtain this_block_dist_right
		for (int j = 0; j < pow(2.0, kr); j++)
		{
			cout << "j = " << j << endl;
			double prob = 0;
			double this_prob;
			for (int s = 0; s < pow(2.0, i*k); s++)
			{
				int value;
				unsigned int mid_value = j << (k*i);	
				for (unsigned int q = 0; q < pow(2.0, k*(m-1-i)+kk); q++)
				{
					int shift_bits = k * i + kr;
					unsigned int qq = q << shift_bits;	
				 	value = s + mid_value + qq; 		
					this_prob = input_dist[value];
					cout << "value = " << value << ", this_prob = " << this_prob << endl;
					prob += this_prob; 
				}
			}
			this_block_dist_right.insert(pair<int, double>(j, prob));
		}
		cout << "this_block_dist_right: " << endl;
		double sum_block_right = 0;
		for (itrm_id = this_block_dist_right.begin(); itrm_id != this_block_dist_right.end(); itrm_id++)
		{
			cout << itrm_id->first << ": " <<  itrm_id->second  << endl;
			sum_block_right += itrm_id->second;
		}
		cout << "sum_block_right = " << sum_block_right << endl;
		block_dist_right.push_back(this_block_dist_right);
	}		
	
	//step3. obtain the distribution for the sum of each block by computing convolution
	vector<map<int, double> > sum_block_dist; 
	for (int i = 0; i < m; i++)
	{
		//obtain PL_prob, GL_prob, KL_prob
		map<int, double> cur_block_dist_left = block_dist_left[i];
		map<int, double> cur_sum_block_dist_left;
		double prob_GL = 0, prob_KL = 0;
		for (int j = 0; j < pow(2.0, kk+1)-1; j++)
		{
			int value = j;
			double prob = 0;
			for (int r = 0; r <= j; r++)
			{
				int x = r;
				int y = j - r;
				itrm_id = cur_block_dist_left.find(x);
				double prob1 = itrm_id->second;
				itrm_id = cur_block_dist_left.find(y);
				double prob2 = itrm_id->second;
				prob += prob1 * prob2;
			}		
		        cur_sum_block_dist_left.insert(pair<int, double>(value, prob));	
			if (value == pow(2.0, kk) - 1)
				PL_prob.push_back(prob);
			if (value > pow(2.0, kk) - 1)
				prob_GL += prob;
			if (value < pow(2.0, kk) - 1)
				prob_KL += prob;
		}
		GL_prob.push_back(prob_GL);
		KL_prob.push_back(prob_KL);
		sum_block_dist_left.push_back(cur_sum_block_dist_left);

		//obtain PR_prob, GR_prob, KR_prob
		map<int, double> cur_block_dist_right = block_dist_right[i];
		map<int, double> cur_sum_block_dist_right;
		double prob_GR = 0, prob_KR = 0;
		for (int j = 0; j < pow(2.0, kr+1)-1; j++)
		{
			int value = j;
			double prob = 0;
			for (int r = 0; r <= j; r++)
			{
				int x = r;
				int y = j - r;
				itrm_id = cur_block_dist_right.find(x);
				double prob1 = itrm_id->second;
				itrm_id = cur_block_dist_right.find(y);
				double prob2 = itrm_id->second;
				prob += prob1 * prob2;
			}		
		        cur_sum_block_dist_right.insert(pair<int, double>(value, prob));	
			if (value == pow(2.0, kr) - 1)
				PR_prob.push_back(prob);
			if (value > pow(2.0, kr) - 1)
				prob_GR += prob;
			if (value < pow(2.0, kr) - 1)
				prob_KR += prob;
		}
		GR_prob.push_back(prob_GR);
		KR_prob.push_back(prob_KR);
		sum_block_dist_right.push_back(cur_sum_block_dist_right);
        } 
}
*/
