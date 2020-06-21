#ifndef HELPER_H
#define HELPER_H 
 
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
#include <sys/time.h>

using namespace std;

void d2b(int d, vector<int> &bin, int numBit);

void d2b_ull(unsigned long long d, vector<int> &bin, int numBit);

void obtain_block_dist(int n, int k, vector<double> &input_dist, vector<double> &P_prob, vector<double> &G_prob, vector<double> &K_prob);

void obtain_block_dist_GL(int n, int k, int kk, vector<double> &input_dist, vector<double> &P_prob, vector<double> &G_prob, vector<double> &K_prob, vector<double> &GL_prob);


void obtain_block_sep_dist(int n, int k, int kk, vector<double> &input_dist, vector<double> &PL_prob, vector<double> &PR_prob, vector<double> &GL_prob, vector<double> &GR_prob, vector<double> &KL_prob, vector<double> &KR_prob);

#endif 

