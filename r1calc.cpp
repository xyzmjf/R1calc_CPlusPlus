/*  	
	Program r1calc 
	This code uses the Gauss Newton non linear regression algorithm 
	to provide the least squares fit to data from the NMR 
	inversion recovery Fourier Transform (IRFT) experiment for a single isolated spin 
 
	There are 3 optimisation parameters for which initial esimates must be provided.
	These are (effectively) peak area at zero time, peak area at long time and relaxation rate R1
	Given these estimates the program iterates until the final and optimum least squares solution is found.
	Note that if the initial estimates are very poor then the method may not converge

	Note that this code is easily modified to fit other functional forms.
	Change the function called func to specify the mathematical form of the function to use for fitting. 
	Change the function called deriv to specify the gradient of func with respect to the optimisation parameters. 
	The optimisation parameters are held in the floating point vector called a

	Compile code with command line:
	g++ r1calc.cpp -lm -o r1calc

	Execute code with command line 
	./r1calc 

	Enter number of x and y data points when prompted
	Enter x (time) and y (peak area) data points when prompted 
	Enter initial estimates when prompted 
	Code should then iterate reducing least squares residuals 
	until convergence criterion < 0.001 
	The x y(experimental) and y(calculated) values will then be printed
	after the optimisation parameter final values. 

	Example input data - 6 data points
	0 	-10
	1 	-5
	2 	0
	3 	3
	5 	4
	10 	-10

	Example for initial parameter estimates
	x or peak area at zero time = -10
	y or peak area at long time = 10
	relaxation rate r1 = 0.5
	

*/



using namespace std;
#include <iostream>
#include <cmath>

/* 
-----------------------------------------------
*/

/* generate space for array of floating point variables 
   using malloc, access the array elements using pointers */
float *fvector(int nl,int nh)
{
float *v;
v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
if (!v) {
	cout << "allocation error in function fvector" << endl; 
}
return (v-nl);
}

/* allocate space for 2D array of floats
   with rows nrl..nrh and columns ncl..nch */
float **fmatrix(int nrl,int nrh,int ncl,int nch)
{
int i;
float **m;
/* allocate array of pointers to rows */
m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
if (!m) { 
	cout << "allocation failure 1 in fmatrix" << endl;
}
m-=nrl;
/* allocate rows and set pointers to them */
for (i=nrl;i<=nrh;i++) {
        m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
        if (!m[i]) {
		cout << "allocation failure 2 in fmatrix" << endl;
	}
        m[i]-=ncl; }; /* end for i */
/* return pointer to array of pointers to rows */
return(m);
} /* end fmatrix */

 

float func(float *a,float tau)
{
float f;
f=a[1]+a[2]*exp(-a[3]*tau);
return(f);
}



float deriv(float *a,float x,int i)
/* derivative function in numeric form */
/* derivative of function wrt parameter i */
{
float dela,f,f2,f1;
dela=a[i]*0.001;
if (dela==0.0) dela=1.0e-8;
a[i]=a[i]+dela;
f2=func(a,x);
a[i]=a[i]-dela;
f1=func(a,x);
f=(f2-f1)/dela;
return(f);
}

void getxy(int npoints,float *x,float *y)
{
int i;
cout << "enter time value and peak area value separated by space " << endl;
for (i=1;i<=npoints;i++) {
	cin >> x[i] >> y[i];
	}
}



void printxy(int npoints,float *x,float *y)
{
int i;
char c;
cout << "tau value, peak area data entered" << endl;
for (i=1;i<=npoints;i++) {
	cout << x[i] << "\t" << y[i] << endl;
	}
}



void holdit()
{
char c;
c=getchar();c=getchar();
}


void zcalc(int nparam,int npoints,float *z,float *x,float *y,float *a)
{
int j,k,pt;
for (j=1;j<=nparam;j++) {
        z[j]=0;
        for (pt=1;pt<=npoints;pt++)
        z[j]=z[j]+(y[pt]-func(a,x[pt]))*deriv(a,x[pt],j);
        } /* end for j */
} /* end zcalc */



void pcalc(int nparam,int npoints,float *x,float *a,float **p)
{
int j,k,pt;
float dj,dk,sum;
for (j=1;j<=nparam;j++) {
        for (k=1;k<=nparam;k++) {
                sum=0;
                for (pt=1;pt<=npoints;pt++) {
                        dj=deriv(a,x[pt],j);
                        dk=deriv(a,x[pt],k);
                        sum=sum+dj*dk;
                        }; /* end for pt */
                p[j][k]=sum;
        }; /* end for k */
}; /* end for j */
} /* end pcalc */


void copymatp(int nparam,float **p,float **pnv)
/* make pnv a copy of matrix p
   prior to inversion in place */
{
int i,j;
for (i=1;i<=nparam;i++) {
        for (j=1;j<=nparam;j++) pnv[i][j]=p[i][j];
        }; /* end for i */
} /* end copymatp */




float sumsq(int npoints,int nparam,float *x,float *y,float *a)
{
float sum,delta;
int i;
sum=0;
for (i=1;i<=npoints;i++) {
        delta=(y[i]-func(a,x[i]));
        sum=sum+delta*delta;
        }; /* end for i */
return(sum);
} /* end sumsq */


float abval(float x)
/* return absolute value of x */
{
float y;
y=x;
if (x<0) y=(-1.0)*x;
return(y);
}

void matinv(int nparam,float **p,float **pnv)
/* calculate pnv as inverse of matrix p */
{
int i,j,k,n;
float small;
copymatp(nparam,p,pnv);
n=nparam;
small=1e-8;
for (i=1;i<=n;i++) {
        if (abval(pnv[i][i])<small) {
		cout << "singularity in row " << i << endl;
	}
        pnv[i][i]=1.0/pnv[i][i];
        for (j=1;j<=n;j++) {
                if (i!=j) pnv[i][j]=pnv[i][j]*pnv[i][i];
                }; /* end for j */
        for (k=1;k<=n;k++) {
                if (i!=k) {
                        for (j=1;j<=n;j++) {
                        if (i!=j) pnv[k][j]=pnv[k][j]-pnv[i][j]*pnv[k][i];
                        }; /* end for j */
                }; /* end if */
        }; /* end for k */
        for (k=1;k<=n;k++) {
                if (i!=k) pnv[k][i]=(-1.0)*pnv[k][i]*pnv[i][i];
        }; /* end for k */
}; /* end for i */
} /* end matinv */


/* function to optimise parameters stored in array a
   so as to minimise sum of errors squared between data
   held in arrays x,y and function func.
   Gauss-Newton linearisation algorithm used. */
/* dynamically allocated arrays x,y hold data points
   a is vector of optimisation parameters */
void iterate(int nparam,int npoints,float *x,float *y,float *a)
{
int k,j;
float conv,*z,*acor;
float **p,**pnv;
/* dynamically allocate space for array elements
   z and acor 1..nparam, arrays acor and z local to function iterate()
   matrices p and pnv also local to iterate()  */
z=fvector(1,nparam);
acor=fvector(1,nparam);
p=fmatrix(1,nparam,1,nparam);
pnv=fmatrix(1,nparam,1,nparam);
cout << "iterating" << endl;
cout << "exit when convergence criterion<0.001" << endl;
cout << "conv\t\tsumsq" << endl;
do {
pcalc(nparam,npoints,x,a,p);
zcalc(nparam,npoints,z,x,y,a);
matinv(nparam,p,pnv);
/* calculate increments */
for (k=1;k<=nparam;k++) {
        acor[k]=0;
        for (j=1;j<=nparam;j++) acor[k]=acor[k]+pnv[k][j]*z[j];
        }; /* end for k */
/* add corrections and calculate conv */
conv=0;
for (k=1;k<=nparam;k++) {
        a[k]=a[k]+acor[k];
        conv=conv+abval(acor[k]/a[k]);
        }; /* end for k */
cout << conv << "\t" << sumsq(npoints,nparam,x,y,a) << endl;
} while (conv>0.001);
} /*end iterate */





/* -------------------------------------------------------- */ 

int main()
{
int i,npoints,nparam;
int one;
float mzero,tau;
/* the function fvector dynamically allocates memory
   one dimensional arrays x,y,a,acor required */
float *x,*y,*a,*acor;
nparam=3;

cout << "Program r1calc" << endl;
cout << "Non linear fit to single exponential recovery curve " << endl;
cout << "Enter number of data points ";
cin >> npoints;
cout << "Number of data points entered = " << npoints << endl;
/* dynamically allocate space for arrays x and y */
one=1;
x=fvector(one,npoints);
y=fvector(one,npoints);
a=fvector(1,nparam);
getxy(npoints,x,y);
printxy(npoints,x,y);

/* get initial estimates of optimisation parameters */
cout << "enter peak area at long time ";
cin >> a[1];
cout << "enter peak area at zero time " ;
cin >> mzero;
a[2]=mzero-a[1];
cout << "enter estimate for longitudinal relaxation rate R1 " ;
cin >> a[3];

/* optimise parameters */
iterate(nparam,npoints,x,y,a);
cout << "iteration finished " << endl;
cout << "Peak area at long time=\t" << a[1] << endl;
cout << "Peak area at zero time=\t" << a[1]+a[2] << endl;
cout << "R1 in inverse seconds\t" << a[3] << endl;
cout << "x(obs)\ty(obs)\ty(calc)" << endl;
for (i=1;i<=npoints;i++) { 
	cout << x[i] << "\t" << y[i] << "\t" << func(a,x[i]) << endl;
}


return 0;
}



