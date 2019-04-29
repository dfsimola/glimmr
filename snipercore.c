// Version: 1.0
// 
// Daniel F. Simola (simola@mail.med.upenn.edu)
// Laboratory of Junhyong Kim
// University of Pennsylvania
// May 2011
// 
// Copyright (c) 2011, Daniel F. Simola and Junhyong Kim, University of
// Pennsylvania.  All Rights Reserved.
// 
// You may not use this file except in compliance with the terms of our
// License. You may obtain a copy of the License at http://kim.bio.upenn.edu/software/
// 
// Unless required by applicable law or agreed to in writing, this
// software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied.  See the License
// for the specific language governing permissions and limitations
// under the License.


#include <Python.h>
#include <math.h>

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

static double Pbin(int n, int k, double p) {
	double nchoosek = prefact[n]/(prefact[k]*prefact[n-k]);
	return nchoosek * pow(p, k) * pow(1-p, n-k);
}


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


// attempt at a speed up
// static PyObject *
// maths_statePriS(PyObject *self, PyObject *args) {
// 	char* S;
// 	char* R;
// 	char* Q;
// 	const int length;
// 	double lambda;
// 	double err;
// 	
// 	double y;
// 	
// 	if (!PyArg_ParseTuple(args, "sssidd", &S, &R, &Q, &length, &lambda, &err))
// 		return NULL;
// 	
// 	if (lambda < 0) {
// 		
// 		y = 0.99;
// 		
// 	} else {
// 	
// 		int A2Q(char x);
// 		double Pbin(int n, int k, double p);
// 		int j;
// 	
// 		double yp = 1.0;
// 		int yd = 0;
// 		for (j=0;j<length;j++) {
// 			if (S[j] == R[j]) {
// 				yp *= 1.0-Q2P(A2Q(Q[j]));
// 				yd += 0;
// 			} else {
// 				yp *= Q2P(A2Q(Q[j]));
// 				yd += 1;
// 			}
// 		}
// 		// // no case
// 		// double np = 1.0;
// 		// int nd = 0;
// 		// for (j=0;j<length;j++) {
// 		// 	if (S[j] == R[j]) {
// 		// 		np *= Q2P(A2Q(Q[j]));
// 		// 		nd += 0;
// 		// 	} else {
// 		// 		np *= 1.0-Q2P(A2Q(Q[j]));
// 		// 		nd += 1;
// 		// 	}
// 		// }
// 		
// 		y = lambda*Pbin(length, yd, err) + (1.0-lambda)*yp;
// 		// double n = lambda*Pbin(length, nd, err) + (1.0-lambda)*np;
// 	}
// 	
// 	// return Py_BuildValue("{s:d,s:d}", "Y",y, "N",n);
// 	return Py_BuildValue("{s:d}", "Y",y);
// }



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
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initsnipercore(void) {
    (void) Py_InitModule("snipercore", MathsMethods);
}