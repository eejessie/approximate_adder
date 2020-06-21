#ifndef GLOBAL_VAR_H
#define GLOBAL_VAR_H

#define sample_size 10000
#define sample_times 1
#define dist_class 1 
#define miu 5e+5 
#define sigma 7e+4 
#define p_geo 0.00001
//int sample_size = 10000;
//int sample_times = 1;
//int dist_class = 1; 
//double miu = 5e+5;
//int sigma = 7e+4;
//double p_geo = 0.00001;

struct inout_pair
{
	unsigned long long add1_high, add1_low;
	unsigned long long add2_high, add2_low;
	unsigned long long sum_high, sum_low;
	unsigned dev;
};

#endif 
