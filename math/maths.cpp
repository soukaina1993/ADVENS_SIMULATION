
#include "maths.h"

using namespace std;
//----------------------------------------------
double sqr (double const &x)
{	return pow(x,2);
}
//----------------------------------------------
void Solve2ndPolynome (double const &a,
                       double const &b,
                       double const &c,
                       double s1,
                       double s2)
//a x^2+b x+c=0, return INF if there are no real solutions
{	double Delta, Delta2;
    Delta2 = sqr(b) - 4 * a * c;
    if (Delta2 < 0.0)
    {	s1 = NAN;
        s2 = NAN;
    }
    else
    {	Delta = sqrt(Delta2);
        s1 = (-b + Delta) / 2 / a;
        s2 = (-b - Delta) / 2 / a;
    }
}
//----------------------------------------------
void SolveBinome (double const &a,
                  double const &b,
                  double const &c,
                  double const &DeltaPrec,
                  double s1,
                  double s2)
//a x^2+b x+c=0, return INF if there are no real solutions,if �=b^2-4ac<DeltaPrec -> �=0
{	double Delta, Delta2;
    Delta2 = sqr(b) - 4 * a * c;
    if (Delta2 < -DeltaPrec)
    {	s1 = NAN;
        s2 = NAN;
    }
    else
    {	if (Delta2 < DeltaPrec)	Delta = 0.0;
        else					Delta = sqrt(Delta2);
        s1 = (-b + Delta) / 2 / a;
        s2 = (-b - Delta) / 2 / a;
    }
}
//----------------------------------------------
double XPwrY(double const &x,
             double const &y)
{	return pow(x,y);
}
//----------------------------------------------
double XPwrI(double const &x,
             long int &I)
{	double cur,ans ;
    if (I==0)  return 1.0;
    if (I<0)  { I=(-I); cur=1/x; }
    else cur=x;
    ans=1;
    while (I>1)
    {	if (!I == 1)  ans = ans * cur;
        cur = cur*cur;
        I = I / 2;
    }
    return ans*cur;
}
//----------------------------------------------
double Atan2Pi (double const &x,
                double const & y)
{	double arc;
    arc = atan(y / x);
    if (x >= 0)	if (y >= 0)	return arc;
        else		return arc + 2 * M_PI;
    else		return arc + M_PI;
}
//----------------------------------------------
/*double AtanPi (double const &x,
				double const & y)
//-�<=Atan < �
{	arc double;
	arc = ArcTan(y / x);
	if (x >= 0)	return arc;
	else 		if (y >= 0)	return arc + pi;
				else		return arc - pi;
}*/
//----------------------------------------------
double AtanPi (double const &x,
               double const & y)
//-�<=Atan < �
{	return atan2(y,x);
}
//----------------------------------------------
double arcsin (double const &a)
{	return atan(a / sqrt(1 - sqr(a)));
}
//----------------------------------------------
double ModReal (double const &a,
                double const &d)
// reste= a modulo d, avec d fractionnaire
{	return (a / d - trunc(a / d)) * d;
}
//----------------------------------------------
double DivReal (double const &a,
                double const &d)
// division enti�re= a div d, avec d fractionnaire
{	return trunc(a / d);
}
//----------------------------------------------
double RoundReal (double const &a,
                  double const &d)
// arrondi a en base d,  avec d fractionnaire
{	return d * trunc(a / d);
}

//----------------------------------------------
double Rad2Deg (double const &Rad)
{	return Rad / M_PI * 180;
}
//----------------------------------------------
double Deg2Rad (double const &Deg)
{	return Deg * M_PI / 180;
}
//----------------------------------------------
double StringTodouble (const std::string word)
{
    //short Index,ValidPrefix;
    //decimal dec ;
    //Index = 0;
    //void str2dec ( const char *s, short *ix, decimal *d, short *vp );
    //str2dec(&word, &Index, &dec, &ValidPrefix);
    //if (ValidPrefix == 0) ;//do what is necessary
    //return dec2num(&dec);
    return std::stod(word);
}
//----------------------------------------------
/*const char doubleToString (double const &Num,
						int dec)
{	char word;
	decform format;
	decimal d;
	format.style = FIXEDDECIMAL; //FLOATDECIMAL or FIXEDDECIMAL
	format.digits = dec;
	dec2str(&format,&d, &word);
	if (word=='?')	word='NAN';
	return word;
}*/
//----------------------------------------------
inline double XpY(double & x, double& y){
    return sign(x) * exp(y * log(fabs(x)));
}
double A(const double &x, const double &y) {
    return sign(x) * exp(y * log(fabs(x)));
}
double sum(vector<double> var)
{return std::accumulate(var.begin(),var.end(),0.0);}

vector<double> sum(vector<double> a,vector<double> b){
    std::transform (a.begin(), a.end(), b.begin(), a.begin(), std::plus<double>());
    return a;
}

vector<double> vminus(vector<double> a,vector<double> b){
    std::transform (a.begin(), a.end(), b.begin(), a.begin(), std::minus<double>());
    return a;
}

int sum(vector<int> var){
    return accumulate(var.begin(),var.end(),0);
}
double mean(vector<double> var){
    return sum(var)/var.size();
}

double mean(vector<int> var){
    return sum(var)/var.size();
}
vector<double> vN(vector<double> var){

    double suma=sum(var);

    for(int i(0); i<var.size(); ++i)
        if (suma!=0)
            var[i]=var[i]/suma;

    return var;
}

vector<double> vN(vector<double> var,double tot){

    double suma=sum(var);

    for(int i(0); i<var.size(); ++i)
        if (suma!=0)
            var[i]=var[i]/suma*tot;

    return var;
}

vector<double> multi(vector<double> &vec, int var){
    vector<double> mvec(vec.size());
    for(int i(0); i<vec.size(); ++i)
        mvec[i]=vec[i]*var;

    return mvec;
}

vector<double> multi(vector<double> &vec, double var){
    vector<double> mvec(vec.size());
    for(int i(0); i<vec.size(); ++i)
        mvec[i]=vec[i]*var;

    return mvec;
}

double exp10(const double &i) {
    return pow(10.0,i);
}

double max(vector<double> var){
    return *max_element(var.begin(), var.end());
}

double min(vector<double> var){
    return *min_element(var.begin(), var.end());
}

int min(vector<int> var){
    return *min_element(var.begin(), var.end());
}

int min_pos(vector<int> var){
    return distance(var.begin(),min_element(var.begin(),var.end()));
}


vector<double> Syear_month(vector<double> var){
    vector<int> dat= {0, 744, 1416, 2160,2880, 3624,4344,5088,5832,6552,7296,8016,8760};
    vector<double> month(12);
    for(int j(0);j<12;j++) {
        for (int i=dat[j];i < dat[j+1]; i++){
            month[j] +=var[i];
        }
    }
    return month;
}

vector<double> Myear_month(vector<double> var){
    vector<int> dat= { 0, 744, 1416, 2160,2880, 3624,4344,5088,5832,6552,7296,8016,8760};
    vector<double> month(12);
    for(int j(0);j<12;j++) {
        for (int i=dat[j]; i < dat[j + 1]; i++){
            month[j]+=var[i];}
        month[j]=month[j]/(dat[j+1]-dat[j]);
    }
    return month;
}

vector<double> Mhour_year(vector<double> var){
    vector<double> year;
    for(int j(0);j<365;j++) {
        for (int i(0); i < 24; i++){
            year.push_back(var[i]);}
    }
    return year;
}

vector<double> Mmonth_year(vector<double> var){
    vector<double> year;
    for (int i(0); i < 31*24; i++){
        year.push_back(var[1]/31*24);}
    for (int i(0); i < 28*24; i++){
        year.push_back(var[2]/28*24);}
    for (int i(0); i < 31*24; i++){
        year.push_back(var[3]/31*24);}
    for (int i(0); i < 30*24; i++){
        year.push_back(var[4]/30*24);}
    for (int i(0); i < 31*24; i++){
        year.push_back(var[5]/31*24);}
    for (int i(0); i < 30*24; i++){
        year.push_back(var[6]/30*24);}
    for (int i(0); i < 31*24; i++){
        year.push_back(var[7]/31*24);}
    for (int i(0); i < 31*24; i++){
        year.push_back(var[8]/31*24);}
    for (int i(0); i < 30*24; i++){
        year.push_back(var[9]/30*24);}
    for (int i(0); i < 31*24; i++){
        year.push_back(var[10]/31*24);}
    for (int i(0); i < 30*24; i++){
        year.push_back(var[11]/30*24);}
    for (int i(0); i < 31*24; i++){
        year.push_back(var[12]/31*24);}
    return year;
}

