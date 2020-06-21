#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <sys/timeb.h>
#include <time.h>
#include "header/helper.h"
#include "header/proposed.h"
#include "header/enumerate.h"
#include "header/bb_adder.h"
#include "header/global_var.h"

using namespace std;


void MC(int n, int k, int l, double &ER, double &MSE, map<unsigned long long, int> &error_dist)
{
	cout << "MC starts: " << endl;
	map<unsigned long long, int>::iterator itrm;
	srand((unsigned)time(NULL));
	
	int m = n/k;
	unsigned long long add1_high = 0, add1_low = 0, add2_high = 0, add2_low = 0;   
	unsigned long long sum_true_low = 0, sum_true_high = 0, sum_est_low = 0, sum_est_high = 0, dev_low = 0, dev_high = 0;
	int ci = 0;
	
	multimap<unsigned long long, struct inout_pair> results;
	multimap<unsigned long long, struct inout_pair>::iterator itr;
	
	for (int i = 0; i < sample_size; i++)
	{
		int seed = rand();
		add1_low = 0, add1_high = 0, add2_low = 0, add2_high = 0;
		sum_true_low = 0, sum_true_high = 0, sum_est_low = 0, sum_est_high = 0;

		vector<int> add1_low_vec, add1_high_vec, add2_low_vec, add2_high_vec;
		for (int j = 0; j < n/2; j++)
		{
			int c1 = rand() % 2;
			add1_low = add1_low * 2 + c1;
			if (c1) add1_low_vec.push_back(1);
			else add1_low_vec.push_back(0);
			
			int c2 = rand() % 2;
			add2_low = add2_low * 2 + c2;
			if (c2) add2_low_vec.push_back(1);
			else add2_low_vec.push_back(0);
		}
		
		for (int j = 0; j < n/2; j++)
		{
			int c1 = rand() % 2;
			add1_high = add1_high * 2 + c1;
			if (c1) add1_high_vec.push_back(1);
			else add1_high_vec.push_back(0);
			
			int c2 = rand() % 2;
			add2_high = add2_high * 2 + c2;
			if (c2) add2_high_vec.push_back(1);
			else add2_high_vec.push_back(0);
		}
		cout << endl;
		cout << "add1:" << endl;
		for (int i = 0; i < n/2; i++)
			cout << add1_high_vec[i];
		cout << " ";
		for (int i = 0; i < n/2; i++)
			cout << add1_low_vec[i];
		cout << endl;
		cout << "add2:" << endl;
		for (int i = 0; i < n/2; i++)
			cout << add2_high_vec[i];
		cout << " ";
		for (int i = 0; i < n/2; i++)
			cout << add2_low_vec[i];
		cout << endl;
		cout << "add1_high = " << add1_high << ", add1_low = " << add1_low << endl;
		cout << "add2_high = " << add2_high << ", add2_low = " << add2_low << endl;
	
		//compute sum_true
		sum_true_low = add1_low + add2_low;
		unsigned int max_low = 1;
		for (int i = 0; i < n/2-1; i++)
			max_low = max_low * 2 + 1;
		cout << "max_low = " << max_low << endl;
		ci = 0;
		if (sum_true_low > max_low) 
		{
			ci = 1;
			sum_true_low = sum_true_low - pow(2.0, n/2);
		}
		sum_true_high = add1_high + add2_high + ci;
		cout << "sum_true_high = " << sum_true_high << ", sum_true_low = " << sum_true_low << endl;
				
		//compute sum_est
		bb_adder(n, k, l, add1_low, add1_high, add2_low, add2_high, sum_est_low, sum_est_high);
		cout << "sum_est_high = " << sum_est_high << ", sum_est_low = " << sum_est_low << endl;
		
		unsigned long long dev = (sum_true_low - sum_est_low) + ((sum_true_high - sum_est_high) * pow(2.0, n/2));
		cout << "dev = " << dev << endl;
		if (dev == 0) continue;
		if (dev < 0)
		{
		 	cout << "dev < 0, error !!!" << endl;
		 	exit(1);
		}		
		
		struct inout_pair ip;
		ip.add1_high = add1_high;
		ip.add1_low = add1_low;
		ip.add2_high = add2_high;
		ip.add2_low = add2_low;
		ip.sum_high = sum_est_high;
		ip.sum_low = sum_est_low;
		results.insert(make_pair(dev, ip));
				
		itrm = error_dist.find(dev);
		if(itrm == error_dist.end())	
			error_dist.insert(make_pair(dev, 1));
		else
			itrm->second++;

	/*	int dev_int = (int)log2(dev);
		if (log2(dev) - dev_int >= 0.5)
			dev_int += 1;
		itrm = error_dist.find(dev_int);
		if(itrm == error_dist.end())	
			error_dist.insert(make_pair(dev_int, 1));
		else
			itrm->second++;
	*/
	}
	cout << endl << endl;
	
	cout << "results: " << endl;
	for (itr = results.begin(); itr != results.end(); itr++)
	{
		struct inout_pair ip = itr->second;
		cout << itr->first << ": " << ip.add1_high << "," << ip.add1_low << " + " << ip.add2_high << "," << ip.add2_low << " = " <<  ip.sum_high << "," << ip.sum_low << endl; 
	}
	
	unsigned long long num_errors = 0;
	MSE = 0;
	cout << "error distribution: " << endl;
	for(itrm = error_dist.begin(); itrm != error_dist.end(); itrm++)
	{
		cout << itrm->first << ", " << itrm->second << endl;
		if (itrm->first != 0)
			num_errors += itrm->second;
		unsigned long long ed = itrm->first;
		int num_ed = itrm->second;
		double prob = num_ed/(double)sample_size;
	//	cout << "prob = " << prob << endl;
		MSE += ed * ed * prob;
	//	cout << "0. MSE = " << MSE << endl; 
	}
		
	cout << "num_errors = " << num_errors << endl;
	ER = num_errors/(double)sample_size;
	cout << "MC: ER = " << ER << ", MSE = " << MSE << endl;
}


void MC_simulation(int n, int k, int l)
{
	map<unsigned long long, int> error_dist;
	map<unsigned long long, int>::iterator itrm;
	map<int, double> ed_prob;
	map<int, double>::iterator itrm_id;
	vector<double> max_prob;
	vector<double> min_prob;
	for(int i = 0; i < sample_times; i++)
	{
		double ER = 0, MSE = 0;
		error_dist.clear();
		MC(n, k, l, ER, MSE, error_dist);
		cout << "i = " << i << ", ER = " << ER << ", MSE = " << MSE << endl;

		if (max_prob.empty())
		{
			for (int j = 0; j < error_dist.size(); j++)
			{
				max_prob.push_back(0);
				min_prob.push_back(1);
			}
		}
	//	cout << "probs: " << endl;
		int j = 0;
		for (itrm = error_dist.begin(); itrm != error_dist.end(); itrm++, j++)
		{
			double prob = itrm->second/(double)sample_size;
			if (prob > max_prob[j])
				max_prob[j] = prob;
			if (prob < min_prob[j])
				min_prob[j] = prob;
		}
	//	cout << endl;
	}
	
/*	cout << "max_prob: " << endl;
	for (int i = 0; i < max_prob.size(); i++)
		cout << max_prob[i] << " ";
	cout << endl;
	cout << "min_prob: " << endl;
	for (int i = 0; i < max_prob.size(); i++)
		cout << min_prob[i] << " ";
	cout << endl;
*/	
/*	cout << "errors: " << ed_prob.size() << endl;
	for (itrm_id = ed_prob.begin(); itrm_id != ed_prob.end(); itrm_id++)
		cout << itrm_id->first << " ";
	cout << endl;
	cout << "probs: " << ed_prob.size() << endl;
	for (itrm_id = ed_prob.begin(); itrm_id != ed_prob.end(); itrm_id++)
	{
		int num = itrm_id->second;
		double prob = num/((double)sample_size * sample_times);
		cout << prob << " ";
	}
	cout << endl;
*/	
}
