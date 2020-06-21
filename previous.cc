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

using namespace std;

void previous(int n, int k, int l)
{
	int kk = k+l;
	int B = k;
	double p = (pow(2.0, B) - 1)/pow(2.0, kk+1);
	int power = ceil((n-kk)/B);
	double ER = 1 - pow((1-p), power);


//	double ER = 1 - pow((1-(pow(2.0, B)-1)/(2.0, kk+1)), ceil((n-kk)/B));	
	double MSE = pow(2.0, 2*n-kk-1)/(pow(2.0, B) + 1);
	
	cout << "previous method: ER = " << ER << ", MSE = " << MSE << endl;
}
