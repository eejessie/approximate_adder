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
#include "header/global_var.h"

using namespace std;

void ED(int i, int j, long long ePar, double pPar, int m, int k, int t, double *d, double **e, 
		vector<long long> &ed_vec, vector<double> &prob_vec, int &count)
{
//	cout << endl << "Entering into ED!" << endl;
//	cout << "i = " << i << ", j = " << j << endl; 
	count++;
	if (i >= m)
	{
		pPar *= d[m-j-1];
	//	setprecision(10);
//		cout << "ePar = " << setprecision(12) << ePar << ", pPar = " << pPar << endl;
	//	cout << "log2(ePar) = " << setprecision(7) << log2(ePar) << ", log2(pPar) = " << log2(pPar) << endl;
		if (ePar < 1e-20 && ePar >= 0);
		else 
		{
		//	ed_vec.push_back(log2(ePar));
		//	prob_vec.push_back(log2(pPar));
			ed_vec.push_back(ePar);
			prob_vec.push_back(pPar);
		}
		return;
	}
	long long new_ePar = ePar + pow(2.0, i*k);
	double new_pPar = pPar * e[i][j];
	ED(i+t+1, i, new_ePar, new_pPar, m, k, t, d, e, ed_vec, prob_vec, count);
	ED(i+1, j, ePar, pPar, m, k, t, d, e, ed_vec, prob_vec, count);
	return;
}

//proposed method
void proposed(int n, int k, int l)
{	 
	int m = n/k;
	int t = l/k;
	double *d = new double[m];
	if (d == NULL)
	{
		cout << "allocating memory for d fails!" << endl;
		exit(2);
	}
	double **e = new double* [m];
	if (e == NULL)
	{
		cout << "allocating memory for e fails!" << endl;
		exit(2);
	}
	for (int i = 0; i < m; i++)
	{
		e[i] = new double[m];
		if (e[i] == NULL)
		{
			cout << "allocating memory for e[i] fails!" << endl;
			exit(3);
		}
	}

	double ER = 0;
	vector<double> input_dist;
	if (dist_class == 2 || dist_class == 3)  //Gaussian distribution or Geometric distribution  
	{
		//generate input distributions and then obtain P_prob, G_prob, K_prob 
		if (dist_class == 2)
		{
			double coefficient = 1.0/(sqrt(2*M_PI)*sigma);
			for (int i = 0; i < pow(2.0, n); i++)
			{
				double prob = coefficient * exp(-pow(i-miu, 2)/(2*sigma*sigma));
				input_dist.push_back(prob);
			}
		}
		else if (dist_class == 3)
		{
			for (int i = 0; i < pow(2.0, n); i++)
			{
				double prob = pow((1 - p_geo), i) * p_geo;	
				input_dist.push_back(prob);
			}
		}
		cout << "input_dist: " << endl;
		double sum = 0;
		for (int i = 0; i < input_dist.size(); i++)
			sum += input_dist[i];
		cout << "sum = " << sum << endl;

	    vector<double> P_prob, G_prob, K_prob, GL_prob;
		if (l - t*k < 1e-6 && l - t*k >= 0) // l is a multiple of k
		{
	                obtain_block_dist(n, k, input_dist, P_prob, G_prob, K_prob);
			cout << "P_prob: ";
			for (int i = 0; i < P_prob.size(); i++) cout << P_prob[i] << " ";	
			cout << endl;
			cout << "G_prob: ";
			for (int i = 0; i < G_prob.size(); i++) cout << G_prob[i] << " ";
			cout << endl;
			cout << "K_prob: ";
			for (int i = 0; i < K_prob.size(); i++) cout << K_prob[i] << " ";
			cout<< endl;
			//calculate di
			for (int i = 0; i <= t; i++) d[i] = 1;
			for (int i = t+1; i < m; i++)
			{
				double part1 = 1, part2 = 0, part3 = 0;
				for (int j = 0; j < i; j++)
					part1 = part1 * P_prob[j];
				for (int j = 1; j <= t; j++)
				{
					double part2_each = G_prob[i-j] * d[i-j];
					for (int s = i-j+1; s <= i-1; s++)
						part2_each = part2_each * P_prob[s];	
					part2 = part2 + part2_each;
				}
				for (int j = 1; j <= i; j++)
				{
					double part3_each = K_prob[i-j] * d[i-j];
					for (int s = i-j+1; s <= i-1; s++)
						part3_each = part3_each * P_prob[s];
					part3 = part3 + part3_each;
				}
				d[i] = part1 + part2 + part3;
			}
			cout << "values of d: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << d[i] << " ";
				if (i == m-1)	ER = 1 - d[i];
			}

		        //calculate e[i][j]
			for (int i = 0; i < m; i++)
				for (int j = 0; j < m; j++)
				{
					if (i-j <= t)	e[i][j] = 0;
					else
					{
						e[i][j] = G_prob[i-t-1] * d[i-t-j-1];
						for (int s = i-t; s <= i-1; s++)
							e[i][j] = e[i][j] * P_prob[s];
					}
				}
				
			cout << "values of e: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << "i = " << i << endl;	
				for (int j = 0; j < m; j++)
					cout << e[i][j] << " "; 
				cout << endl;
			}
		}
		else   // l is not a multiple of k
		{
			int kk = l - t * k;
	                obtain_block_dist_GL(n, k, kk, input_dist, P_prob, G_prob, K_prob, GL_prob);
			cout << "P_prob: ";
			for (int i = 0; i < P_prob.size(); i++) cout << P_prob[i] << " ";	
			cout << endl;
			cout << "G_prob: ";
			for (int i = 0; i < G_prob.size(); i++) cout << G_prob[i] << " ";
			cout << endl;
			cout << "K_prob: ";
			for (int i = 0; i < K_prob.size(); i++) cout << K_prob[i] << " ";
			cout<< endl;
			cout << "GL_prob: ";
			for (int i = 0; i < GL_prob.size(); i++) cout << GL_prob[i] << " ";
			cout << endl;
			//calculate di
			for (int i = 0; i <= t; i++)    d[i] = 1;
			for (int i = t+1; i < m; i++)
			{
				double part1 = 1, part2 = 1, part3 = 0, part4 = 0;
				for (int j = 0; j <= i-1; j++)
					part1 *= P_prob[j];
				part2 = GL_prob[i-t-1] * d[i-t-1]; 				
				for (int j = i-t; j <= i-1; j++)
					part2 *= P_prob[j];
				for (int j = 1; j <= i; j++)
				{
					double part3_each = K_prob[i-j] * d[i-j];
					for (int s = i-j+1; s <= i-1; s++)
						part3_each *= P_prob[s];
					part3 += part3_each;
				}
				for (int j = 1; j <= t; j++)
				{
					double part4_each = G_prob[i-j] * d[i-j];
					for (int s = i-j+1; s <= i-1; s++)
						part4_each *= P_prob[s];
					part4 += part4_each;
				}
				d[i] = part1 + part2 + part3 + part4;
			}
			cout << "values of d: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << d[i] << " ";
				if (i == m-1)	{
					ER = 1 - d[i];
					cout <<  "\n ER = " << ER << endl;
				}					
			}
		}
	} else if (dist_class == 1) {
		//calculate p(Pi), p(Gi), p(Ki)
		double pp = 1/pow(2.0, k);
		double pg = 0.5 - 1/pow(2.0, k+1);
		double pk = pg;
		cout << "pp = " << pp << ", pg = " << pg << ", pk = " << pk << endl;
			
		double ER;
		if (l - t*k < 1e-6 && l - t*k >= 0)  // == 0: l is a multiple of k
		{		
			cout << "case 1: l is a multiple of k!" << endl;
			//calculate di
			for (int i = 0; i <= t; i++)
				d[i] = 1;
			for (int i = t+1; i < m; i++)
			{
				double v_part1 = pow(pp, i);
				
				double v_part2 = 0;
				for (int j = 1; j <= t; j++)
					v_part2 += pow(pp, j-1)*pg*d[i-j];
					
				double v_part3 = 0;
				for (int j = 1; j <= i; j++)
					v_part3 += pow(pp, j-1)*pk*d[i-j];
					
				d[i] = v_part1 + v_part2 + v_part3;
			}
			cout << "values of d: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << d[i] << " ";
				if (i == m-1) {
					ER = 1 - d[i];
					cout << "\n ER = " << ER << endl;
				}							
			}
			cout << endl;
		
			//calculate e_{i,j}
			for (int i = 0; i < m; i++)
			{
				if (i < t+1)
					e[i][0] = 0;
				else
					e[i][0] = pow(pp, t) * pg * d[i-t-1];
			}
				
			for (int i = 0; i < m; i++)
				for (int j = 1; j < m; j++)
				{
					if (i-j <= t)
						e[i][j] = 0;
					else
						e[i][j] = pow(pp, t) *pg * d[i-t-j-1];
				}
				
			cout << "values of e: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << "i = " << i << endl;	
				for (int j = 0; j < m; j++)
					cout << e[i][j] << " "; 
				cout << endl;
			}
		} else { // != 0: l is not a multiple of k 
			cout << "case 1: l is not a multiple of k!" << endl;
			int kk = l - t*k;
		//	cout << "kk = " << kk << endl;
			//calculate p(PLi), p(GLi), p(KLi), p(PRi), p(GRi), p(KRi)
			double ppl = 1/pow(2.0, kk);
			double pgl = 0.5 - 1/pow(2.0, kk+1);
			double pkl = pgl;		
			double ppr = 1/pow(2.0, k - kk);
			double pgr = 0.5 - 1/pow(2.0, k-kk+1);
			double pkr = pgr;
			cout << "ppl = " << ppl << ", pgl = " << pgl << ", pkl = " << pkl << endl;
			cout << "ppr = " << ppr << ", pgr = " << pgr << ", pkr = " << pkr << endl;
		
			//calculate di
			for (int i = 0; i <= t; i++)
				d[i] = 1;
			for (int i = t+1; i < m; i++)
			{
				double v_part1 = pow(pp, i);
				
				double v_part2 = 0;
				for (int j = 1; j <= t; j++)
					v_part2 += pow(pp, j-1)*pg*d[i-j];
					
				double v_part3 = 0;
				for (int j = 1; j <= i; j++)
					v_part3 += pow(pp, j-1)*pk*d[i-j];
					
				double v_part4 = pow(pp, t)*pgl*d[i-t-1];
					
				d[i] = v_part1 + v_part2 + v_part3 + v_part4;
			}
			
			cout << "values of d: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << d[i] << " ";
				if (i == m-1) {
					ER = 1 - d[i];
					cout << "\n ER = " << ER << endl;
				}					
			}
			cout << endl;
		
			//calculate e_{i,j}
			for (int i = 0; i < m; i++)
			{
				if (i < t+1)
					e[i][0] = 0;
				else if (i == t+1)
					e[i][0] = pow(pp, t) * ppl * pgr;
				else
					e[i][0] = pow(pp, t) * ppl * pgr * d[i-t-1] + pow(pp, t+1) * pgl * d[i-t-2];
			}
			
			for (int i = 0; i < m; i++) {
				for (int j = 1; j < m; j++)
				{
					if (i-j < t+1)
						e[i][j] = 0;
					else if (i-j == t+1)
						e[i][j] = pow(pp, t) * ppl * pgr;
					else
						e[i][j] = pow(pp, t) *ppl * pgr * d[i-t-j-1] + pow(pp, t+1) * pgl* d[i-t-j-2];
				}
			}
			cout << "values of e: " << endl;
			for (int i = 0; i < m; i++)
			{
				cout << "i = " << i << endl;	
				for (int j = 0; j < m; j++)
					cout << e[i][j] << " "; 
				cout << endl;
			}
		}
	}
	
	//calculate error distritution
	vector<long long> ed_vec;
	vector<double> prob_vec;
	int count = 0;
	ED(t+1, 0, 0, 1, m, k, t, d, e, ed_vec, prob_vec, count);
	cout << "count = " << count << endl;
	cout << "ed_vec: " << ed_vec.size() << endl;
	map<long long, double> ed_prob;
	map<long long, double>::iterator itrm_id;
	for(int i = 0; i < ed_vec.size(); i++)
	{
	//	cout << ed_vec[i] << " " << prob_vec[i] << endl;;
		
		int integ_part = (int)log2(ed_vec[i]);
		if (log2(ed_vec[i]) - integ_part >= 0.5)
			integ_part++;
		itrm_id = ed_prob.find(integ_part);
		if (itrm_id == ed_prob.end())
			ed_prob.insert(make_pair(integ_part, prob_vec[i]));
		else
			itrm_id->second += prob_vec[i];	
	
	/*	itrm_id = ed_prob.find(ed_vec[i]);
		if (itrm_id == ed_prob.end())
			ed_prob.insert(make_pair(ed_vec[i], prob_vec[i]));
		else
			itrm_id->second += prob_vec[i];	
	*/
	}
	cout << endl;

	cout << endl << "ed_prob_0: " << endl;
	for (itrm_id = ed_prob.begin(); itrm_id != ed_prob.end(); itrm_id++)
		cout << "(" << itrm_id->first << ", " << itrm_id->second << ")" << endl;
	cout << endl;

	double MSE = 0;
	cout << endl << "prob_vec: " << prob_vec.size() << endl;
	for(int i = 0; i < prob_vec.size(); i++)
	{
    //	cout << prob_vec[i] << " ";
		MSE += (double)ed_vec[i] * (double)ed_vec[i] * prob_vec[i];
	}
	cout << endl;

//	cout << "proposed method: ER = " << ER << ", MSE = " << MSE << endl;
	cout << "proposed method: MSE = " << MSE << endl;

	
//	cout << "prob in map: " << endl;
	for (itrm_id = ed_prob.begin(); itrm_id != ed_prob.end(); itrm_id++)
	{
	//	cout << itrm_id->second << " ";
		double prob = itrm_id->second;
		int ir = itrm_id->first/k;
		itrm_id->second = prob/d[m-1-ir];
	}
	cout << endl;
	
	cout << "ed_prob_1: " << endl;
	for (itrm_id = ed_prob.begin(); itrm_id != ed_prob.end(); itrm_id++)
		cout << "(" << itrm_id->first << ", " << itrm_id->second << ")" << endl;
	cout << endl;

	
	//free dynamic-allocated memory
	delete []d;
	for (int i = 0; i < m; i++)
		delete []e[i];
	delete []e;
}
