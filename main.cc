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
#include "header/previous.h"
#include "header/MC.h"
#include "header/global_var.h"

using namespace std;

int main(int argc, char *argv[]) 
{

    if(argc != 5)
    {
    	//n: total number of bits of the adder
    	//k: block size
    	//l: carry generator length
        cout << "Correct usage: ./main n k l sel > record_file " << endl;
        exit(1);
    }
    
	struct timeb startTime, endTime;                         //Record the computing time.
	ftime(&startTime);
	
	struct timeval st, et;
	gettimeofday(&st, NULL);
 	
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	int l = atoi(argv[3]);
	int sel = atoi(argv[4]);
	int m = n/k;
	int t = l/k;
	cout << "n = " << n << ", m = " << m << ", k = " << k << ", l = " << l << ", t = " << t << endl;
	
	if (sel == 1)
		proposed(n, k, l); //proposed method
	else if (sel == 2)
		enumerate(n, k, l); //enumerate
	else if (sel == 3) //Monte-Carlo simulation
	{
		MC_simulation(n, k, l);
	}
	else if (sel == 4) //previous method
		previous(n, k, l);

    ftime(&endTime);
        
    double runtime = ((endTime.time - startTime.time)*1000 + (endTime.millitm - startTime.millitm))/1000.0;
    cout << "total runtime: " << runtime << endl;
    gettimeofday(&et, NULL);
//  cout << "total time: " << ((et.tv_sec - st.tv_sec) * 1000000 + et.tv_usec - st.tv_usec) << endl;
	
    return 0;
}

