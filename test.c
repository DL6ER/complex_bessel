#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "cbessel.h"

int main(int argc, char const *argv[])
{
	// Chose some parameters
	double a = 1.0;
	double b = 1.0;
	double complex in = a + b * I, r;

	printf("Test program for Bessel / Hankel function computation\nusing complex_bessel library\n\n");

	printf("************ Testing zeroth order ************\n\n");
	printf("Test 1: Bessel function of first kind, zeroth order\n");
	r = besselJ(0,in);
	printr(r);
	// BesselJ(0, (1+i1))
	prints(9.376085e-01, -4.965299e-01);

	printf("Test 2: Bessel function of second kind, zeroth order\n");
	r = besselY(0,in);
	printr(r);
	// BesselY(0, (1+i1))
	prints(4.454745e-01, 7.101586e-01);

	printf("Test 3: Hankel function of first kind, zeroth order\n");
	r = hankelH1(0,in);
	printr(r);
	// hankelH1(0, (1+i1))
	prints(2.274499e-01, -5.105546e-02);

	printf("Test 4: Hankel function of second kind, zeroth order\n");
	r = hankelH2(0,in);
	printr(r);
	// hankelH2(0, (1+i1))
	prints(1.647767e+00, -9.420044e-01);

	printf("\n************ Testing negative order ************\n\n");

	printf("Test 5: Bessel function of first kind, -1 order\n");
	r = besselJ(-1,in);
	printr(r);
	// BesselJ(-1, (1+i1))
	prints(-6.141603e-01, -3.650280e-01);

	printf("Test 6: Bessel function of second kind, -1 order\n");
	r = besselY(-1,in);
	printr(r);
	// BesselY(-1, (1+i1))
	prints(6.576945e-01, -6.298010e-01);

	printf("Test 7: Hankel function of first kind, -1 order\n");
	r = hankelH1(-1,in);
	printr(r);
	// hankelH1(-1, (1+i1))
	prints(1.564067e-02, 2.926665e-01);

	printf("Test 8: Hankel function of second kind, -1 order\n");
	r = hankelH2(-1,in);
	printr(r);
	// hankelH1(-2, (1+i1))
	prints(-1.243961e+0, -1.022723e+00);

	printf("\n************ Testing derivatives (1st order) ************\n\n");

	printf("Test 9: Bessel function of first kind, zeroth order\n");
	r = besselJdiff(0,1,in);
	printr(r);
	// 1/2 (BesselJ(-1, (1+i1)) - BesselJ(1, (1+i1)))
	prints(-6.141603e-01, -3.650280e-01);

	printf("Test 10: Bessel function of second kind, zeroth order\n");
	r = besselYdiff(0,1,in);
	printr(r);
	// 1/2 (BesselJ(-1, (1+i1)) - BesselJ(1, (1+i1)))
	prints(6.576945e-01, -6.298010e-01);

	printf("Test 11: Hankel function of first kind, zeroth order\n");
	r = hankelH1diff(0,1,in);
	printr(r);
	// 1/2 (HankelH1(-1, (1+i1)) - HankelH1(1, (1+i1)))
	prints(1.564067e-02, 2.926665e-01);

	printf("Test 12: Hankel function of second kind, zeroth order\n");
	r = hankelH2diff(0,1,in);
	printr(r);
	// 1/2 (HankelH2(-1, (1+i1)) - HankelH2(1, (1+i1)))
	prints(-1.243961e+0, -1.022723e+00);

	printf("\n************ Testing derivatives (2nd order) ************\n\n");

	printf("Test 13: Bessel function of first kind, zeroth order\n");
	r = besselJdiff(0,2,in);
	printr(r);
	// (-(1+i1) BesselJ(-1, (1+i1)) + (- (1+i1)^2) BesselJ(0, (1+i1)))/(1+i1)^2
	prints(-4.480143e-01, 3.719638e-01);

	printf("Test 14: Bessel function of second kind, zeroth order\n");
	// (-(1+i1) BesselY(-1, (1+i1)) + (- (1+i1)^2) BesselY(0, (1+i1)))/(1+i1)^2
	r = besselYdiff(0,2,in);
	printr(r);
	prints(-4.594213e-01, -6.641081e-02);

	printf("Test 15: Hankel function of first kind, zeroth order\n");
	// (-(1+i1) HankelH1(-1, (1+i1)) + (- (1+i1)^2) HankelH1(0, (1+i1)))/(1+i1)^2
	r = hankelH1diff(0,2,in);
	printr(r);
	prints(-3.816035e-01, -8.745746e-02);

	printf("Test 16: Hankel function of second kind, zeroth order\n");
	// (-(1+i1) HankelH2(-1, (1+i1)) + (- (1+i1)^2) HankelH2(0, (1+i1)))/(1+i1)^2
	r = hankelH2diff(0,2,in);
	printr(r);
	prints(-5.144251e-01, 8.313850e-01);

	return 0;
}
