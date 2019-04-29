#include <Python.h>
#include <math.h>
#include "nrutil.c"

// Convert phred score to probability
static double Q2P(int e) {
	if (e == 0) return 1-1e-16;
	else return pow(10.0,-e/10.);
} 


// static int A2Q[94] = {'!', '"', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '{', '|', '}', '~'};

static int aA2Q(char x) {
	if (x=='!') return 0;
	else if (x=='"') return 1;
	else if (x=='#') return 2;
	else if (x=='$') return 3;
	else if (x=='%') return 4;
	else if (x=='&') return 5;
	else if (x=='\'') return 6;
	else if (x=='(') return 7;
	else if (x==')') return 8;
	else if (x=='*') return 9;
	else if (x=='+') return 10;
	else if (x==',') return 11;
	else if (x=='-') return 12;
	else if (x=='.') return 13;
	else if (x=='/') return 14;
	else if (x=='0') return 15;
	else if (x=='1') return 16;
	else if (x=='2') return 17;
	else if (x=='3') return 18;
	else if (x=='4') return 19;
	else if (x=='5') return 20; 
	else if (x=='6') return 21;
	else if (x=='7') return 22;
	else if (x=='8') return 23;
	else if (x=='9') return 24;
	else if (x==':') return 25;
	else if (x== ';') return 26; 
	else if (x== '<') return 27; 
	else if (x== '=') return 28; 
	else if (x== '>') return 29; 
	else if (x== '?') return 30; 
	else if (x== '@') return 31; 
	else if (x== 'A') return 32; 
	else if (x== 'B') return 33; 
	else if (x== 'C') return 34; 
	else if (x== 'D') return 35; 
	else if (x== 'E') return 36; 
	else if (x== 'F') return 37; 
	else if (x== 'G') return 38; 
	else if (x== 'H') return 39; 
	else if (x== 'I') return 40; 
	else if (x== 'J') return 41; 
	else if (x== 'K') return 42; 
	else if (x== 'L') return 43; 
	else if (x== 'M') return 44; 
	else if (x== 'N') return 45; 
	else if (x== 'O') return 46; 
	else if (x== 'P') return 47; 
	else if (x== 'Q') return 48; 
	else if (x== 'R') return 49; 
	else if (x== 'S') return 50; 
	else if (x== 'T') return 51; 
	else if (x== 'U') return 52; 
	else if (x== 'V') return 53; 
	else if (x== 'W') return 54; 
	else if (x== 'X') return 55; 
	else if (x== 'Y') return 56; 
	else if (x== 'Z') return 57; 
	else if (x== '[') return 58; 
	else if (x== '\\') return 59; 
	else if (x== ']') return 60; 
	else if (x== '^') return 61; 
	else if (x== '_') return 62; 
	else if (x== '`') return 63; 
	else if (x== 'a') return 64; 
	else if (x== 'b') return 65; 
	else if (x== 'c') return 66; 
	else if (x== 'd') return 67; 
	else if (x== 'e') return 68; 
	else if (x== 'f') return 69; 
	else if (x== 'g') return 70; 
	else if (x== 'h') return 71; 
	else if (x== 'i') return 72; 
	else if (x== 'j') return 73; 
	else if (x== 'k') return 74; 
	else if (x== 'l') return 75; 
	else if (x== 'm') return 76; 
	else if (x== 'n') return 77; 
	else if (x== 'o') return 78; 
	else if (x== 'p') return 79; 
	else if (x== 'q') return 80; 
	else if (x== 'r') return 81; 
	else if (x== 's') return 82; 
	else if (x== 't') return 83; 
	else if (x== 'u') return 84; 
	else if (x== 'v') return 85; 
	else if (x== 'w') return 86; 
	else if (x== 'x') return 87; 
	else if (x== 'y') return 88; 
	else if (x== 'z') return 89; 
	else if (x== '{') return 90; 
	else if (x== '|') return 91; 
	else if (x== '}') return 92; 
	else if (x== '~') return 93; 
}

// This is fast enough
static int A2Q(char x) {
	switch (x) {
		case '!': return 0;
		case '"': return 1;
		case '#': return 2;
		case '$': return 3;
		case '%': return 4;
		case '&': return 5;
		case '\'': return 6;
		case '(': return 7;
		case ')': return 8;
		case '*': return 9;
		case '+': return 10;
		case ',': return 11;
		case '-': return 12;
		case '.': return 13;
		case '/': return 14;
		case '0': return 15;
		case '1': return 16;
		case '2': return 17;
		case '3': return 18;
		case '4': return 19;
		case '5': return 20; 
		case '6': return 21;
		case '7': return 22;
		case '8': return 23;
		case '9': return 24;
		case ':': return 25;
		case ';': return 26; 
		case '<': return 27; 
		case '=': return 28; 
		case '>': return 29; 
		case '?': return 30; 
		case '@': return 31; 
		case 'A': return 32; 
		case 'B': return 33; 
		case 'C': return 34; 
		case 'D': return 35; 
		case 'E': return 36; 
		case 'F': return 37; 
		case 'G': return 38; 
		case 'H': return 39; 
		case 'I': return 40; 
		case 'J': return 41; 
		case 'K': return 42; 
		case 'L': return 43; 
		case 'M': return 44; 
		case 'N': return 45; 
		case 'O': return 46; 
		case 'P': return 47; 
		case 'Q': return 48; 
		case 'R': return 49; 
		case 'S': return 50; 
		case 'T': return 51; 
		case 'U': return 52; 
		case 'V': return 53; 
		case 'W': return 54; 
		case 'X': return 55; 
		case 'Y': return 56; 
		case 'Z': return 57; 
		case '[': return 58; 
		case '\\': return 59; 
		case ']': return 60; 
		case '^': return 61; 
		case '_': return 62; 
		case '`': return 63; 
		case 'a': return 64; 
		case 'b': return 65; 
		case 'c': return 66; 
		case 'd': return 67; 
		case 'e': return 68; 
		case 'f': return 69; 
		case 'g': return 70; 
		case 'h': return 71; 
		case 'i': return 72; 
		case 'j': return 73; 
		case 'k': return 74; 
		case 'l': return 75; 
		case 'm': return 76; 
		case 'n': return 77; 
		case 'o': return 78; 
		case 'p': return 79; 
		case 'q': return 80; 
		case 'r': return 81; 
		case 's': return 82; 
		case 't': return 83; 
		case 'u': return 84; 
		case 'v': return 85; 
		case 'w': return 86; 
		case 'x': return 87; 
		case 'y': return 88; 
		case 'z': return 89; 
		case '{': return 90; 
		case '|': return 91; 
		case '}': return 92; 
		case '~': return 93;
	}
}

// Returns the value ln[gamma(xx)] for xx > 0. 
// Internal arithmetic will be done in double precision, a nicety 
// that you can omit if ﬁve-ﬁgure accuracy is good enough.
static double 
gammln(double xx) { 
	double x,y,tmp,ser; 
	static double cof[6] = {76.18009172947146,-86.50532032941677, 
	                        24.01409824083091,-1.231739572450155, 
	                        0.1208650973866179e-2,-0.5395239384953e-5}; 
	int j; 
	y= x = xx; 
	tmp = x + 5.5; 
	tmp -= (x+0.5)*log(tmp); 
	ser = 1.000000000190015; 
	for (j=0;j<=5;j++) ser += cof[j]/++y; 
	return -tmp+log(2.5066282746310005*ser/x); 
} 

// Returns the incomplete gamma function P(a,x).
double
gammp(double a, double x) { 
	void gcf(double *gammcf, double a, double x, double *gln); 
	void gser(double *gamser, double a, double x, double *gln); 
	void nrerror(char error_text[]); 
	double gamser,gammcf,gln; 
	
	if (x <0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x <(a+1.0)) { // Use the series representation.
		gser(&gamser,a,x,&gln); 
		return gamser; 
	} else { // Use the continued fraction representation.
		gcf(&gammcf,a,x,&gln); 
		return 1.0-gammcf; // and take its complement. 
	} 
} 
// helper functions for gammp
#define ITMAX 1000 // Maximum allowed number of iterations. 
#define EPS 1.0e-100 // Relative accuracy. 
#define FPMIN 1.0e-100 // Number near the smallest representable FP number

// Returns the incomplete gamma function P(a,x) evaluated by its 
// series representation as gamser. 
// Also returns ln gamma(a) as gln. 
void gser(double *gamser, double a, double x, double *gln) { 
	double gammln(double xx);
	void nrerror(char error_text[]); 
	int n; 
	double sum,del,ap; 
	
	*gln=gammln(a); 
	if (x <=0.0) { 
		if (x< 0.0) nrerror("x less than 0 in routine gser"); 
		*gamser=0.0; 
		return; 
	} else { 
		ap=a; 
		del=sum=1.0/a; 
		for (n=1;n<=ITMAX;n++) { 
			++ap; 
			del *=x/ap; 
			sum +=del; 
			if(fabs(del) < fabs(sum)*EPS) { 
				*gamser=sum*exp(-x+a*log(x)-(*gln)); 
				return; 
			} 
		} 
		nrerror("a too large, ITMAX too small in routine gser"); 
		return; 
	} 
}
// Returns the incomplete gamma function Q(a,x) evaluated by its continued 
// fraction representation as gammcf. Also returns ln gamma(a) as gln. 
void gcf(double *gammcf, double a, double x, double *gln) { 
	double gammln(double xx); 
	void nrerror(char error_text[]); 
	int i; 
	float an,b,c,d,del,h; 
	*gln=gammln(a); 
	b=x+1.0-a; // Setup for evaluating continued fraction by modiﬁed Lentz’s method(§5.2) with b0=0. 
	c=1.0/FPMIN; 
	d=1.0/b; 
	h=d; 
	for (i=1;i<=ITMAX;i++) { // Iterate to convergence. 
		an =-i*(i-a); 
		b+= 2.0; 
		d=an*d+b; 
		if (fabs(d) < FPMIN) d=FPMIN; 
		c=b+an/c; 
		if (fabs(c) < FPMIN) c=FPMIN; 
		d=1.0/d; 
		del=d*c; 
		h*= del; 
		if (fabs(del-1.0) <EPS) break; 
	} 
	if (i >ITMAX) nrerror("a too large, ITMAX too small in gcf"); 
	*gammcf=exp(-x+a*log(x)-(*gln))*h; // Put factors in front. 
}


//Returns the error functionerf(x).  
double
erfff(double x) {
	double gammp(double a, double x); 
	return x< 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x); 
} 

double 
gaussiancdf(double x, double m, double s) {
	double erfff(double x);
	// double y = 
	return .5*(1.0 + erfff( (x-m)/(s*sqrt(2.0)) ));
}

static double 
incompleteBeta(double x, double a, double b) {
	double beta = 0.0;
	double dt = 0.000003;
	double t = dt;
	while (t < x) {
		beta += pow(t, a-1) * pow(1-t, b-1) * dt;
		t += dt;
	}
	return beta;
}


static double 
incompleteBetaFast(double x, double a, double b, float dt) {
	double beta = 0.0;
	double t = dt;
	while (t < x) {
		beta += pow(t, a-1) * pow(1-t, b-1) * dt;
		t += dt;
	}
	return beta;
}


static float betaDFS(float a, float b) {
	return incompleteBeta(1.0, a, b);
}

static double 
beta(double z, double w) { 
	// Returns the value of the beta function B(z,w). 
	double gammln(double xx); 
	return exp(gammln(z)+gammln(w)-gammln(z+w)); 
} 

static double
regularizedIncompleteBeta(double x, double a, double b) {
	return incompleteBeta(x,a,b)/beta(a,b);
}

static double
regularizedIncompleteBetaFast(double x, double a, double b, float dt) {
	return incompleteBetaFast(x,a,b,dt)/beta(a,b);
}


static double
regularizedIncompleteGamma(double a, double x) {
	double gammp(double aa, double xx);
	return gammp(a,x);
}

static double
gaussianpdf(double x, double m, double s) {
	return 1/(2*s)*exp(-(pow(x-m,2)/(2*pow(s,2))));
}

static double
tdistcdf(double t, double v) {
	return 1 - regularizedIncompleteBeta(v/(v+(t*t)), 0.5*v, 0.5);
}

static double
factrl(int n) {
	double gammln(double xx);
	void nchar(char error_text[]);
	static int ntop = 4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;
	
	if (n<0) nrerror("Negative factorial in routine factrl");
	if (n>32) return exp(gammln(n+1.0));
	while (ntop<n) {
		j=ntop++;
		a[ntop] = a[j]*ntop;
	}
	return a[n];
}

// Iterative implementation of factorial function.
static double factit(int x) {
	if (x <= 1) return 1;
	double a[x+1];
	a[0]=1; a[1]=1;
	int i;
	for (i=2;i<=x;i++) a[i] = i*a[i-1];
	return a[x];
}

static double prefact[201];

static void initPrefact(int nana) {
	int j;
	double factit(int j);
	// double prefact[maxn];
	for (j=0;j<201;j++) {
		prefact[j] = factit(j);
	}
}

// static double Power(double a, int b) { 
//     if (b==0) return 1; 
//     if (a==0) return 0; 
//     if (b%2==0) { 
//         return Power(a*a, b/2); 
//     } else if (b%2==1) { 
//         return a*Power(a*a,b/2); 
//     } 
//     return 0;
// }

static double Pbin(int n, int k, double p) {
	// double factit(int n);
	// double nchoosek = factit(n)/(factit(k)*factit(n-k));
	// double Power(double a, int b);
	double nchoosek = prefact[n]/(prefact[k]*prefact[n-k]);
	return nchoosek * pow(p, k) * pow(1-p, n-k);
	// return nchoosek * Power(p,k) * Power(1-p, n-k);
}


///////////////////////////////////////////////////////////

static PyObject *
maths_GaussianPDF(PyObject *self, PyObject *args) {
	const float x, m, s;
	if (!PyArg_ParseTuple(args, "fff", &x, &m, &s))
		return NULL;
	return Py_BuildValue("f", gaussianpdf(x, m, s));
}

static PyObject *
maths_GaussianCDF(PyObject *self, PyObject *args) {
	const float x, m, s;
	if (!PyArg_ParseTuple(args, "fff", &x, &m, &s))
		return NULL;
	return Py_BuildValue("f", gaussiancdf(x, m, s));
}


static PyObject *
maths_tdistPDF_old(PyObject *self, PyObject *args) {
	const float t; const float s; const float m;
	float zdistpdf(float tt, float mm, float ss);
	if (!PyArg_ParseTuple(args, "fff", &t, &m, &s))
		return NULL;
	return Py_BuildValue("f", gaussianpdf(t,m,s));
}

static PyObject *
maths_tdistCDF(PyObject *self, PyObject *args) {
	const float t; const float df;
	double tdistcdf(double tt, double dd);
	if (!PyArg_ParseTuple(args, "ff", &t, &df))
		return NULL;
	return Py_BuildValue("f", tdistcdf(t,df));
}


static PyObject *
maths_gamma(PyObject *self, PyObject *args) {
	const float a;
	double gammln(double xx);
	if (!PyArg_ParseTuple(args, "f", &a))
		return NULL;
	return Py_BuildValue("f", exp(gammln(a)));
}

static PyObject *
maths_beta(PyObject *self, PyObject *args) {
	// Integral approximation to complete Beta distribution
	const float a, b;
	if (!PyArg_ParseTuple(args, "ff", &a, &b))
		return NULL;
	return Py_BuildValue("f", beta(a,b));
}

static PyObject *
maths_incompleteBeta(PyObject *self, PyObject *args) {
	// Integral approximation to incomplete Beta distribution
	const float x, a, b;
	if (!PyArg_ParseTuple(args, "fff", &x, &a, &b))
		return NULL;
	float beta = incompleteBeta(x,a,b);
	return Py_BuildValue("f", beta);
}


static PyObject *
maths_gammaRegularized(PyObject *self, PyObject *args) {
	const float a, z;
	if (!PyArg_ParseTuple(args, "ff", &a, &z))
		return NULL;
	float p = regularizedIncompleteGamma(a, z);
	return Py_BuildValue("f", p);
}


static PyObject *
maths_betaCDF(PyObject *self, PyObject *args) {
	const float x, a, b;
	if (!PyArg_ParseTuple(args, "fff", &x, &a, &b))
		return NULL;
	float p = regularizedIncompleteBeta(x,a,b);
	return Py_BuildValue("d", p);
}

static PyObject *
maths_betaCDFfast(PyObject *self, PyObject *args) {
	const float x, a, b, dt;
	if (!PyArg_ParseTuple(args, "ffff", &x, &a, &b, &dt))
		return NULL;
	float p = regularizedIncompleteBetaFast(x,a,b, dt);
	float Z = regularizedIncompleteBetaFast(1,a,b, dt);
	return Py_BuildValue("d", p/Z);
}



static PyObject *
maths_FdistCDF(PyObject *self, PyObject *args) {
	const float f, df1, df2, dt;
	if (!PyArg_ParseTuple(args, "ffff", &f, &df1, &df2, &dt))
		return NULL;
	float x = df1*f/(df1*f+df2);
	float p = regularizedIncompleteBetaFast(x,df1/2,df2/2, dt);
	return Py_BuildValue("d", p);
}

static PyObject *
maths_ChiSquareCDF(PyObject *self, PyObject *args) {
	const float x, k;
	if (!PyArg_ParseTuple(args, "ff", &x, &k))
		return NULL;
	float p = regularizedIncompleteGamma(k/2.0, x/2.0);
	return Py_BuildValue("f", p);
}


// this has issues for n > 10
static PyObject *
maths_factorial(PyObject *self, PyObject *args) {
	const int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	int fn = factit(n);
	return Py_BuildValue("i", fn);
}


static PyObject *
maths_prishelper(PyObject *self, PyObject *args) {
	char* S;
	char* R;
	char* Q;
	const int offset;
	const int length;
	int j = 0;
	int A2Q(char x);
	
	if (!PyArg_ParseTuple(args, "sssii", &S, &R, &Q, &offset, &length))
		return NULL;
	
	double p = 1.0;
	int d = 0;
	for (j=0;j<length;j++) {
		if (j != offset) { 
			if (S[j] == R[j]) {
				p *= 1.0-Q2P(A2Q(Q[j]));
				d += 0;
			} else {
				p *= Q2P(A2Q(Q[j]));
				d += 1;
			}
		}
	}
	
	return Py_BuildValue("di", p,d);
}


static PyObject *
maths_ntPriS(PyObject *self, PyObject *args) {
	char* S;
	char* R;
	char* Q;
	const int offset;
	const int length;
	double lambda;
	double err;
	
	if (!PyArg_ParseTuple(args, "sssiidd", &S, &R, &Q, &offset, &length, &lambda, &err))
		return NULL;
	
	int A2Q(char x);
	double Pbin(int n, int k, double p);
	double p = 1.0;
	int d = 0;
	int j = 0;
	for (j=0;j<length;j++) {
		if (j != offset) { 
			if (S[j] == R[j]) {
				p *= 1.0-Q2P(A2Q(Q[j]));
				d += 0;
			} else {
				p *= Q2P(A2Q(Q[j]));
				d += 1;
			}
		}
	}
	
	double qaqo = Q2P(A2Q(Q[offset]));
	double olp = (1.0-lambda)*p;
	double primero = lambda*Pbin(length, d+1, err) + olp*     qaqo;
	double segundo = lambda*Pbin(length, d,   err) + olp*(1.0-qaqo);
	double a,t,c,g;
	
	// this doesn't work properly
	// switch(S[offset]) {
	// 	case 'A': a = segundo; t = c = g = primero;
	// 	case 'T': t = segundo; a = c = g = primero;
	// 	case 'C': c = segundo; a = t = g = primero;
	// 	case 'G': g = segundo; a = t = c = primero;
	// }
	
	if (S[offset]!='A') a = primero;
	else                a = segundo;
	if (S[offset]!='T') t = primero;
	else                t = segundo;
	if (S[offset]!='C') c = primero;
	else                c = segundo;
	if (S[offset]!='G') g = primero;
	else                g = segundo;
	
	return Py_BuildValue("{s:d,s:d,s:d,s:d}", "A",a, "T",t, "C",c, "G",g);
}

static PyObject *
maths_statePriS(PyObject *self, PyObject *args) {
	char* S;
	char* R;
	char* Q;
	const int length;
	double lambda;
	double err;
	
	if (!PyArg_ParseTuple(args, "sssidd", &S, &R, &Q, &length, &lambda, &err))
		return NULL;
	
	int A2Q(char x);
	double Pbin(int n, int k, double p);
	int j;
	
	double yp = 1.0;
	int yd = 0;
	for (j=0;j<length;j++) {
		if (S[j] == R[j]) {
			yp *= 1.0-Q2P(A2Q(Q[j]));
			yd += 0;
		} else {
			yp *= Q2P(A2Q(Q[j]));
			yd += 1;
		}
	}
	// // no case
	// double np = 1.0;
	// int nd = 0;
	// for (j=0;j<length;j++) {
	// 	if (S[j] == R[j]) {
	// 		np *= Q2P(A2Q(Q[j]));
	// 		nd += 0;
	// 	} else {
	// 		np *= 1.0-Q2P(A2Q(Q[j]));
	// 		nd += 1;
	// 	}
	// }
	
	double y = lambda*Pbin(length, yd, err) + (1.0-lambda)*yp;
	// double n = lambda*Pbin(length, nd, err) + (1.0-lambda)*np;
	
	// return Py_BuildValue("{s:d,s:d}", "Y",y, "N",n);
	return Py_BuildValue("{s:d}", "Y",y);
}


static PyObject *
maths_initprefact(PyObject *self, PyObject *args) {
	void initPrefact(int x);
	initPrefact(0);
	return Py_BuildValue("i", 1);
}

////////////////////////////
static PyMethodDef MathsMethods[] = {
	{"ntPriS", maths_ntPriS, METH_VARARGS,
	"ntPriS for Sniper."},
	{"statePriS", maths_statePriS, METH_VARARGS,
	"statePriS for Glimmr."},
	{"prishelper", maths_prishelper, METH_VARARGS,
	"ntPriS for Sniper."},
	{"initprefact", maths_initprefact, METH_VARARGS,
	"TESTING"},
	{"factit", maths_factorial, METH_VARARGS,
	"Iterative factorial function."},
	{"GaussianPDF", maths_GaussianPDF, METH_VARARGS,
	"Compute p-value from a standard normal distribution."},
	{"GaussianCDF", maths_GaussianCDF, METH_VARARGS,
	"Compute cumulative probability of a standard normal distribution."},
	{"gamma", maths_gamma, METH_VARARGS,
	"Compute gamma function."},
	{"beta", maths_beta, METH_VARARGS,
	"Compute complete Beta(a,b)."},
	{"incompleteBeta", maths_incompleteBeta, METH_VARARGS,
	"Compute incomplete Beta(x,a,b)."},
	{"gammaRegularized", maths_gammaRegularized, METH_VARARGS,
	"Compute regularized 2p gamma function P(a,z)."},
	{"FdistCDF", maths_FdistCDF, METH_VARARGS,
	"Compute cumulative probability for F distribution."},
	{"betaCDF", maths_betaCDF, METH_VARARGS,
	"Compute cumulative probability for Beta distribution."},
	{"betaCDFfast", maths_betaCDFfast, METH_VARARGS,
	"Compute cumulative probability for Beta distribution."},
	{"ChiSquareCDF", maths_ChiSquareCDF, METH_VARARGS,
	"Compute cumulative probability for Chi Square distribution."},
	{"tdistCDF", maths_tdistCDF, METH_VARARGS,
	"Compute the p-value of a t distribution."},
	// {"factorial", maths_factorial, METH_VARARGS,
	// "Compute factorial function of long int n."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initmaths(void) {
    (void) Py_InitModule("maths", MathsMethods);
}