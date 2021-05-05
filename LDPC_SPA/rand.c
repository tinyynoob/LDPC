
/************************************************************************
*                E X T E R N A L   R E F E R E N C E S
*************************************************************************
*/
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/timeb.h>

//#include "main.h"
#include "rand.h"


double 
drand(void)
{
    long a;
    double ran,x;
    a=0x7fff;
#ifdef WORK_STATION 
    x=2147483647.;
    ran=(double)((double)(random())/x);
#endif

#ifndef WORK_STATION
    x=(double) a;
    ran=(double)((double)(rand())/x);
#endif
    return(ran);
}


void rand_init()
{


struct timeb timeptr; 

unsigned a;

#ifdef WORK_STATION
    long now;
    srandom(time(&now));
#endif
#ifndef WORK_STATION    
	ftime(&timeptr);
	a = (unsigned) (timeptr.millitm  + timeptr.time);
	//printf("rand seed = %d\n", a); 
    srand( a );
#endif
}


double gaussian(void)
{
double V1,V2,U1,U2,S,X;
static int cnt=0;
static double Y;
    if (cnt==1) {
        cnt=0; 
        return (Y);
    }
    else {
        do{
            U1=drand();
            U2=drand();
            V1=2*U1-1;
            V2=2*U2-1;
            S=pow(V1,2.0)+pow(V2,2.0);
        }while (S>=1);
         X=V1*pow(-2.0*log(S)/S,0.5);
        Y=V2*pow(-2.0*log(S)/S,0.5);
        cnt=1;
        return(X);
    }
}




void rand_info(char *infoBits, long infoLen)
{
	long i;
	

	for (i = 0; i < infoLen; ++i) {
	    infoBits[i] = rand() % 2;
	}

}



void randPos(long *bias_pos, long rang, long n)
{
	long i;
	long j;
	long r;
	long rep;

	for (i = 0; i < n; ++i) {
	    
		do{
		r = rand()%rang;
		rep = 0;
		for (j = 0; j < i; ++j) {
		    if(bias_pos[j] == r) {
                rep = 1;
			    break;
			}
		}
		if(rep == 0) {
		    bias_pos[i] = r;
		}
		}while(rep == 1);
	}
}

