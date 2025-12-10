#ifndef _MATH_
#define _MATH_

#define _USE_MATH_DEFINES

#include <iostream>
using namespace std;  //introduces namespace std

#include <math.h>
//#include "fp.h"
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <string.h>
#if !defined(__STRICT_ANSI__) || defined(_POSIX_C_SOURCE) || defined(_POSIX_SOURCE) || defined(_XOPEN_SOURCE) || defined(_GNU_SOURCE) || defined(_BSD_SOURCE) || defined(_USE_MATH_DEFINES)

#define M_PI		3.14159265358979323846

#endif

#define NDEF	HUGE_VAL
#define NIL		(0)

//double sqr(const double& x);

//double D2A(const double& D);

//double Rad2Deg (const double& Rad);

//double Deg2Rad (const double& Deg);

//void checkmalloc(const char*	message,
//				void*	ptrname);


double sqr (double const &x);
void Solve2ndPolynome (double const &a,double const &b,	double const &c,double s1,double s2);

void SolveBinome (double const &a,double const &b,double const &c,double const &DeltaPrec,double s1,double s2);

double XPwrY(double const &x,double const &y);
inline double XpwrY(const double  &x, const double  &y)
{	return pow(x,y);
}

double XPwrI(double const &x,long int &I);

inline double XpwrI(const double  &x, int i)
{	return pow(x,i);
}

double Atan2Pi (double const &x,
                double const & y);

double AtanPi (double const &x,
               double const & y);

double arcsin (double const &a);

double ModReal (double const &a,
                double const &d);

double DivReal (double const &a,
                double const &d);

double RoundReal (double const &a,
                  double const &d);

double Rad2Deg (double const &Rad);

double Deg2Rad (double const &Deg);

double StringTodouble (const char word);

const char doubleToString (double const &Num,
                           int dec);
inline double sign(const double &x)
{	return x>=0 ? 1:-1;
}
inline double min (const double &x1,const double &x2)
{	return x1<x2 ? x1:x2;
}
inline double pow3 (const double  &x)
{	double temp=x;
    temp*=x;
    temp*=x;
    return temp;
}

inline int round2int(const double &x)
{	return round(x);
}

double exp10(const double &d);

//----------------------------------------------
inline double round(const double &x,
                    const double &i)
{	double aux = exp10(i);
    return floor(x*aux+0.5)/aux;
}




//----------------------------------------------
// cut the fractionnal part of x
inline double roundinf(const double &x,
                       const double &i)
{	double aux = sign(x)*exp10(i);
    return floor(x*aux)/aux;
}
//----------------------------------------------
inline double roundsup(const double &x,
                       const double &i)
{	double aux = sign(x)*exp10(i);
    return ceil(x*aux)/aux;
}
//----------------------------------------------
inline int trunc2int(const double &x)
{	return floor(x);
}
//----------------------------------------------
inline double trunc(const double &x,
                    const double &i)
{	double aux = exp10(i);
    return floor(x*aux)/aux;
}
//------------------------------------
inline double frac(const double &x, const double &NbDec)
{ 	double DecPart = x-floor(x);
    return  round(DecPart, NbDec);
}

//----------------------------------------------
#pragma mark       === UNIT CONVERSION ===
//----------------------------------------------
// conversion rad -> deg

//----------------------------------------------
#pragma mark       === TRIGO ===
//----------------------------------------------

//----------------------------------------------
/*inline double AtanPi (const double &x,
				const double & y)
//-�<=Atan < �
{	arc double;
	arc = ArcTan(y / x);
	if (x >= 0)	return arc;
	else 		if (y >= 0)	return arc + pi;
				else		return arc - pi;
}*/
//----------------------------------------------

//----------------------------------------------
inline double arccos (const double &a)
// result is in [0,pi]
{	return acos(a);
}
//----------------------------------------------
#pragma mark       === STRING ===
//----------------------------------------------
inline char *strdup(char *s)
{	char *temp =(char*)malloc((strlen(s)+1)*sizeof(char));
    strcpy(temp, s);
    return temp;
}
//----------------------------------------------
inline double StringTodouble (const char* word)
{	return atof(word);
}
//----------------------------------------------
inline void StringTodouble (const char* word, double num)
{	num = atof(word);
}
//------------------------------------
const char dig2Char[] = {'0','1','2','3','4','5','6','7','8','9'};
//------------------------------------
inline char* IntegerToString (const int Num, char* st)
{	char* scan = st;
    if (Num!=0)
    {	double cNum=Num;
        if (cNum<0)
        {	cNum=-cNum;
            *scan++='-';
        }
        int NbDig = floor(log10(cNum))+1;
        cNum/=pow(10.0,NbDig);
        double dig;
        while (NbDig--)
        {	cNum*=10;
            cNum = modf(round(cNum,NbDig), &dig);
            *scan++=dig2Char[(int)dig];
        }
    }
    else *scan++='0';
    *scan='\0';
    return scan;
}
//------------------------------------
inline void DoubleToString (const double Num, int NbDec, char *st)
{	int intP = floor(Num);
    char* scan = IntegerToString(intP,st);
    if (NbDec!=0)
    {	*scan++='.';
        int fracP=frac(fabs(Num),NbDec)*pow(10.0,NbDec);
        IntegerToString(fracP,scan);
    }
}
/*----------------------------------------------
inline void Fcc2Str(FourCharCode tag, char* outStr)
{	char* src = (char*)&tag;
	char* dst = outStr;
	for (unsigned i = 4 ; i ; i-- )
		*dst++ = *src++;
	*dst = 0;
//	*(FourCharCode*)outStr = tag;
//	outStr[4] = 0;
}
//----------------------------------------------
inline void Str2Fcc(const char* word, FourCharCode& outTag)
{	const char* src = word;
	char* dst = (char*)&outTag;
	for (unsigned i = 4 ; i ; i-- )
		*dst++ = *src++;
//	outTag = *(FourCharCode*)word;
}

//----------------------------------------------
#pragma mark       === RECT ===
//----------------------------------------------
inline void InsetRect (Rect &r, short dleft, short dtop, short dright,short dbottom)
{	r.left -= dleft;
	r.top  -= dtop;
	r.right += dright;
	r.bottom += dbottom;
}
*/
//----------------------------------------------
#pragma mark       === MATH SOLVER ===
//----------------------------------------------
inline void Solve2ndPolynome (  const double &a,
                                const double &b,
                                const double &c,
                                double &s1,
                                double &s2)
//a x^2+b x+c=0, return INF if there are no real solutions
{	double Delta, Delta2;
    Delta2 = sqr(b) - 4.0 * a * c;
    if (Delta2 < 0.0)
    {	s1 = HUGE_VAL;
        s2 = HUGE_VAL;
    }
    else
    {	Delta = sqrt(Delta2);
        s1 = (-b + Delta) / 2.0 / a;
        s2 = (-b - Delta) / 2.0 / a;
    }
}
//----------------------------------------------

//--------------------------------------------------------
inline void
QuadraticRoots(	const double& c,
                   const double& b,
                   const double& a,
                   double& r1,
                   double& r2)
{	// precise quadratic solution formulation
    // ax^2+bx+c=0
    double q = -0.5*(b+sign(b)*sqrt(sqr(b)-4.0*a*c));
    r1 = q/a;
    r2 = c/q;
}


double A(const double &x, const double &y);

//----------------------------------------------
inline void CubicRoots (	const double& a0,
                            const double& a1,
                            const double& a2,
                            const double& a3,
                            double& Sr1,
                            double& Sr2,
                            double& Si2,
                            double& Sr3,
                            double& Si3)
{	double a03 = a0 / a3;
    double a13 = a1 / a3;
    double a23 = a2 / a3;
    double q = a13 / 3.0 - sqr(a23) / 9.0;
    double r = (a13 * a23 - 3.0 * a03) / 6.0 - sqr(a23) * a23 / 27.0;
    double q3s2 = sqr(q) * q + sqr(r);
    if (q3s2 <= 0.0)
    {	// 3 real solutions
        double q2 = sqrt(-q);
        double temp = r/q/q2;
        double phi = acos(temp);
        Sr1 = -2.0 * q2 * cos(phi / 3.0) - a23 / 3.0;
        Sr3 = -2.0 * q2 * cos((phi + 2.0 * M_PI) / 3.0) - a23 / 3.0;
        Sr2 = -2.0 * q2 * cos((phi + 4.0 * M_PI) / 3.0) - a23 / 3.0;
        Si2 = 0.0;
        Si3 = 0.0;
    }
    else
    {	// 1 real solution & 2 complexes solutions
        double s1 = A(r + sqrt(q3s2), 1.0 / 3.0);
        double s2 = A(r - sqrt(q3s2), 1.0 / 3.0);
        Sr1 =  s1 + s2 - a23 / 3.0;
        Sr2 = -(s1 + s2) / 2.0 - a23 / 3.0;
        Si2 =  (s1 - s2) / 2.0 * sqrt(3.0);
        Sr3 =  Sr2;
        Si3 = -Si2;
    }
}




inline double XpY(double &x, double&y);


//----------------------------------------------
inline void CubicRoots2 (	const double& a0,
                             const double& a1,
                             const double& a2,
                             const double& a3,
                             double& Sr1,
                             double& Sr2,
                             double& Si2,
                             double& Sr3,
                             double& Si3)
{	double a13 = a1 / a3;
    double a23 = a2 / a3;

    double q = a13 / 3.0 - sqr(a23) / 9.0;
    double r = (a13 * a23 - 3.0 * a0) / 6.0 - sqr(a23) * a23 / 27.0;

    if (q <= 0.0)
    { 	// 3 real solutions
        double q2 = sqrt(-q);
        double temp = r/q/q2;

        //double phi = pi/2.0-atan(temp/sqrt(1.0-sqr(temp)));

        double phi = acos(temp);//(-r / sqrt(-sqr(q) * q));
        Sr1 = -2.0 * q2 * cos(phi / 3.0) 			  - a23 / 3.0;
        Sr3 = -2.0 * q2 * cos((phi + 2.0 * M_PI) / 3.0) - a23 / 3.0;
        Sr2 = -2.0 * q2 * cos((phi + 4.0 * M_PI) / 3.0) - a23 / 3.0;
        Si2 = 0.0;
        Si3 = 0.0;
    }
    else
    {	// 1 real solution & 2 complexes solutions
        double q3s2 = sqrt(sqr(q) * q + sqr(r));
        double s1 = pow(r + q3s2, 1.0 / 3.0);
        double s2 = A(r - q3s2, 1.0 / 3.0);
        Sr1 = s1 + s2 - a23 / 3.0;
        Sr2 = -(s1 + s2) / 2.0 - a23 / 3.0;
        Si2 =  (s1 - s2) / 2.0 * sqrt(3.0);
        Sr3 = Sr2;
        Si3 = -Si2;
    }
}


//----------------------------------------------
// return 10^i
inline double exp10(const int &i)
{	return pow(10.0,i);
}
//----------------------------------------------
// return natural log of x
inline double ln (const double  &x)
{	return log(x);
}
//----------------------------------------------
inline double D2A(const double &D)
// diameter -> disc area
{	return sqr(D)/4*M_PI;
}

#pragma mark       === VECTOR ===

double sum(vector<double> var);
vector<double> sum(vector<double> var1,vector<double> var2);
vector<double> vminus(vector<double> var1,vector<double> var2);
int sum(vector<int> var);
double mean(vector<double> var);
double mean(vector<int> var);
vector<double> vN(vector<double> var);
vector<double> vN(vector<double> var,double tot);
double max(vector<double> var);
double min(vector<double> var);
int min(vector<int> var);
int min_pos(vector<int> var);
vector<double> multi(vector<double> &vec, int var);
vector<double> multi(vector<double> &vec, double var);
vector<double>Syear_month(vector<double> var); // sum for one month
vector<double>Myear_month(vector<double> var); // mean for one month
vector<double> Mhour_year(vector<double> var); // from one day to year in hour
vector<double> Mmonth_year(vector<double> var);
#endif
