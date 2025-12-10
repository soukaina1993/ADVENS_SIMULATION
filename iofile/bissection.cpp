#include <math.h>
//Ajoute par Danilo:
#include <iostream>

#define JMAX 40

float rtbis(float (*func)(float), float Xmin, float Xmax, float DeltaX)
{

	int j;
	float dx,f,fmid,xmid,rtb;

	f=(*func)(Xmin);
	fmid=(*func)(Xmax);
	if (f*fmid >= 0.0)
		std::cout<<"Root must be bracketed for bisection in rtbis"<<std::endl;
	rtb = f < 0.0 ? (dx=Xmax-Xmin,Xmin) : (dx=Xmin-Xmax,Xmax);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < DeltaX || fmid == 0.0) return rtb;
	}
	std::cout<<"Too many bisections in rtbis"<<std::endl;
	return 0.0;
}
#undef JMAX
