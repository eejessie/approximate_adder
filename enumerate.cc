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
#include "header/bb_adder.h"
#include "header/global_var.h"

using namespace std;

void enumerate(int n, int k, int l)
{
	cout << "enumerate starts: " << endl;
	map<unsigned long long, int> error_dist;
	map<unsigned long long, int>::iterator itrm;
	
	multimap<unsigned long long, struct inout_pair> results;
	multimap<unsigned long long, struct inout_pair>::iterator itr;

	int num = pow(2.0, n);
	int bit_num = pow(2.0, n/2) - 1;
	unsigned long long add1_low, add1_high, add2_low, add2_high;
	unsigned long long sum_accurate, sum_est_low, sum_est_high, dev;
	for (unsigned long long i = 0; i < num; i++) {
		for (unsigned long long j = 0;  j < num; j++)
		{
			sum_accurate = i + j;
			
			add1_low = i & bit_num;
			add1_high = i >> n/2;
			add1_high &= bit_num;
			
			add2_low = j & bit_num;
			add2_high = j >> n/2;
			add2_high &= bit_num;
			
		//	cout << endl << "$add1 = " << i << ", add2 = " << j << endl;
		/*	cout << "add1_low = " << add1_low << endl;
			cout << "add1_high = " << add1_high << endl;
			cout << "add2_low = " << add2_low << endl;
			cout << "add2_high = " << add2_high << endl;
		*/
			
			sum_est_low = 0; sum_est_high = 0;
			bb_adder(n, k, l, add1_low, add1_high, add2_low, add2_high, sum_est_low, sum_est_high);
			dev = sum_accurate - (sum_est_high * pow(2.0, n/2) + sum_est_low); 
		//	cout << "sum_acc = " << sum_accurate << ", sum_est_high = " << sum_est_high << 
		//		", sum_est_low = " << sum_est_low << endl;
		//	cout << "dev = " << dev << endl;
		//	cout << "add1 = " << i << ", add2 = " << j << ", sum_acc = " << sum_accurate << 
	   //			", sum_est = " << sum_estimate << ", dev = " << dev << endl;
			if (dev > 0)
			{
			//	cout << "add1 = " << i << ", add2 = " << j << ", sum_acc = " << sum_accurate << 
		//			", sum_est = " << (sum_est_high << n/2 + sum_est_low) << ", dev = " << dev << endl;
				itrm = error_dist.find(dev);
				if(itrm == error_dist.end())	
					error_dist.insert(make_pair(dev, 1));
				else
					itrm->second++;
			
			/*	struct inout_pair ip;
				ip.add1_high = add1_high;
				ip.add1_low = add1_low;
				ip.add2_high = add2_high;
				ip.add2_low = add2_low;
				ip.sum_high = sum_est_high;
				ip.sum_low = sum_est_low;
				results.insert(make_pair(dev, ip));
			*/
			}
		}
	}		
		
/*	cout << "results: " << endl;
	for (itr = results.begin(); itr != results.end(); itr++)
	{
		struct inout_pair ip = itr->second;
		cout << itr->first << ": " << ip.add1_high << "," << ip.add1_low << " + " << ip.add2_high << 
			"," << ip.add2_low << " = " <<  ip.sum_high << "," << ip.sum_low << endl; 
	}
*/
				
	double ER, MSE = 0;
	unsigned long long num_errors = 0;
//	cout << "error distribution: " << endl;
	for(itrm = error_dist.begin(); itrm != error_dist.end(); itrm++)
	{
	//	cout << itrm->first << ", " << itrm->second << endl;
		unsigned long long ed = itrm->first;
		int num_ed = itrm->second;
		if (itrm->first != 0)
			num_errors += num_ed;
		double prob = num_ed/pow(2.0, 2*n);
		MSE += ed * ed * prob;
	}
//	ER = num_errors/pow(2.0, 2*n) * 2;
	ER = num_errors/pow(2.0, 2*n);
	cout << "ER = " << ER << ", MSE = " << MSE << endl;
}

