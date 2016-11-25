#include <iostream>
using namespace std;

//#define  tolerance  0.000001
#define  tolerance  0.000001
#define  maxsteps  4000000


float F(float x)
{	
	//f(x) = 2x -1;
	return 2*x*x - 1; 
}
float dF(float x)
{
    //f'(x)
    return 4*x;
}

float newton(float xs)
{
	float x, xprim;
	float t1, t2;
	x = xs;
	xprim = xs + 2*tolerance;
	int iterCount = 0;
	while(((x - xprim  >= tolerance) || (x  - xprim  <= -tolerance))  && (iterCount < maxsteps))
        {
		xprim = x;
		t1 = F(x);
		t2 = dF(x);
		x = x  - t1 / t2;
		iterCount++;
		cout << "x= " << x << endl;
	}
	if(!((x  - xprim  <=  tolerance) && (x  - xprim  >= -tolerance))) 
        {
		//x = INFTY;
		x = 999999;
	}
	
	cout << "Newton's method result " << x << ", IterCount = " << iterCount << endl;
	return x;
}

int main()
{	float xs = 0.001;
	float x = newton(xs);
	return 0;
}
