/********* Misc functions **********/
double cos_pi(double nu)
{
	// Detect if nu is an integer. If |nu|>1e14, the significand is saturated
	// with numbers before the decimal point, and we cannot make the difference between an
	// integer and a real number.
	double nup5 = nu + 0.5;
	if (floor(nup5) == nup5 && fabs(nu) < 1e14)
		return 0.0;

	return cos(M_PI*nu);
}

/// We compute the value of sin(pi*nu), paying particular
/// attention to the case where nu is an integer.
double sin_pi(double nu)
{
	// Detect if nu is an integer. Same comment as above if |nu|>1e14.
	if (floor(nu) == nu && fabs(nu) < 1e14)
		return 0.0;

	return sin(M_PI*nu);
}

// Print complex number
void printcomplex(double complex z)
{
	double r = creal(z), i = cimag(z);
	printf("%e %c i%e",r,(i>=0.0f)? '+':'-',fabs(i));
}

// Print result
void printr(double complex z)
{
	printf("result is: ");
	printcomplex(z);
	printf("\n");
}

// Print solution
void prints(double r, double i)
{
	printf("should be: ");
	printcomplex(r + i * I);
	printf("\n");
}

/********* Bessel functions *************/

/* FORTRAN functions in library */
/* Bessel function of the first kind. */
extern void zbesj_wrap(double,double,double,int,int,double*,double*,int*,int*);
/* Bessel function of the second kind. */
extern void zbesy_wrap(double,double,double,int,int,double*,double*,int*, double*,double*,int*);

/* Computes the Bessel functions of the first kind
For negative order, we use the reflection formula
J_{-\nu}(z) = (-1)^\nu J_\nu(z) */
double complex bessel(double order, int kind, double complex z)
{
	// Input values for Fortran subroutines.
	double zr = creal(z);
	double zi = cimag(z);
	double nu = fabs(order);
	int kode = 1;
	int N = 1;

	// Output values.
	double cyr, cyi, cwrkr, cwrki;
	int nz, ierr;

	// External function call
	if(kind == 1)
	{
		// First kind
		zbesj_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&ierr);
	}
	else if(kind == 2)
	{
		// Second kind
		zbesy_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&cwrkr,&cwrki,&ierr);
	}
	else
	{
		// Not supported
		return 0.0;
	}

	// If the input is real, then the output should be real as well.
	if (zi == 0.0 && zr >= 0.0)
	{
		cyi = 0.0;
	}
	double complex answer = cyr + cyi *I;  // Placeholder for output.

	// If the return code is not normal, we print the error code.
	if (ierr!=0)
		printf("bessel(%f, %i, %e + i%e): Error code %i\n",
			order, kind, creal(z), cimag(z), ierr);

	// If order is negative, then we must apply the reflection formula.
	if (order < 0.0)
	{
		// We prepare the rotation coefficients.
		double c = cos_pi(nu);
		double s = sin_pi(nu);

		// External function call
		if(kind == 1)
		{
			// We compute the Bessel Y function.
			double cyrY, cyiY, cwrkr, cwrki;
			int nzY, ierrY, kodeY = 1, NY = 1;
			zbesy_wrap(zr,zi,nu,kodeY,NY,&cyrY,&cyiY,&nzY,&cwrkr,&cwrki,&ierrY);
			double complex answerY = cyrY + cyiY * I;
			answer = c*answer - s*answerY;
		}
		else if(kind == 2)
		{
			// We compute the Bessel J function.
			double cyrJ,cyiJ;
			int nzJ, ierrJ, kodeJ = 1, NJ = 1;
			zbesj_wrap(zr,zi,nu,kodeJ,NJ,&cyrJ,&cyiJ,&nzJ,&ierrJ);
			double complex answerJ = cyrJ + cyiJ * I;
			answer = s*answerJ + c*answer;
		}
	}

	return answer;
}

// Implement short handles for besselJ
double complex besselJ(double order, double complex z)
{
	return bessel(order, 1, z);
}

// Implement short handles for besselY
double complex besselY(double order, double complex z)
{
	return bessel(order, 2, z);
}

// Computes n-th derivative of the Bessel function of first or second kind
double complex besselDiff(double order, int kind, int n, double complex z)
{
	double complex s = bessel(order-n, kind, z), q;

	// Compute series
	double p = 1.0;
	for(int i=1; i<=n; i++)
	{
		p = -1.0 * (p*(n-i+1)) / i;
		q = bessel(order-n+2*i, kind, z);
		s += p*q;
	}

	return s/pow(2.0,n);
}

// Implement short handles for besselJdiff
double complex besselJdiff(double order, int n, double complex z)
{
	return besselDiff(order, 1, n, z);
}

// Implement short handles for besselYdiff
double complex besselYdiff(double order, int n, double complex z)
{
	return besselDiff(order, 2, n, z);
}


/*********** Hankel functions ************/
/* FORTRAN functions in library complex_bessel */
/*! Hankel function of both kinds. Kind determined by integer argument. */
extern void zbesh_wrap(double,double,double,int,int,int,double*,double*,int*,int*);

/* Computes the Hankel function of requested kind.
   For negative order, we use the reflection formula
   H^{(1)}_{-\nu}(z) = H^{(1)}_\nu(z)\exp(\pi \nu I) */
double complex hankel(double order, int kind, double complex z)
{
	// Input values for Fortran subroutines.
	double zr = creal(z);
	double zi = cimag(z);
	double nu = fabs(order);
	int kode = 1;
	int N = 1;

	// Output values
	double cyr,cyi;
	int nz,ierr;

	// External function call.
	zbesh_wrap(zr,zi,nu,kode,kind,N,&cyr,&cyi,&nz,&ierr);
	double complex answer = cyr + cyi *I;

	// Reflection formula if order is negative.
	if (order < 0.0)
		answer *= cexp(M_PI*nu*I);

	return answer;
}

// Implement short handles for hankelH1
double complex hankelH1(double order, double complex z)
{
	return hankel(order, 1, z);
}

// Implement short handles for hankelH2
double complex hankelH2(double order, double complex z)
{
	return hankel(order, 2, z);
}

/* Computes n-th derivative of the Hankel function of first or second kind
   The derivative of H_n^k(z) is given by (H_{n-1}^k(z) - H_{n+1}^k(z))/2 */
double complex hankelDiff(double order, int kind, int n, double complex z)
{
	double complex s = hankel(order-n, kind, z);

	// Compute series
	double p = 1.0;
	for(int i=1; i<=n; i++)
	{
		p = -1.0 * (p*(n-i+1)) / i;
		s += p*hankel(order-n+2*i, kind, z);
	}

	return s/pow(2.0,n);
}

// Implement short handles for hankelH1diff
double complex hankelH1diff(double order, int n, double complex z)
{
	return hankelDiff(order, 1, n, z);
}

// Implement short handles for hankelH2diff
double complex hankelH2diff(double order, int n, double complex z)
{
	return hankelDiff(order, 2, n, z);
}
