#include "sqrt6.h"
double sqrt6(double y)
{
	double x, z, tempf;
	unsigned long *tfptr = ((unsigned long *)&tempf) + 1;
	tempf = y;
	*tfptr = (0xbfcdd90a - *tfptr) >> 1;
	x = tempf;
	z = y*0.5;
	x = (1.5*x) - (x*x)*(x*z);    //The more you make replicates of this statement 
	//the higher the accuracy, here only 2 replicates are used  
	x = (1.5*x) - (x*x)*(x*z);
	return x*y;
}