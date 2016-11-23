#include <stdio.h>
#include "KernelUtils.h"

using namespace kernels;

// This small standalone program tests the ExaHyPE::kernel::idx and similar
// utilities

void idxtest() {
	int m1, m2, m3, m4, m5, m6;
	m1 = 1;
	m2 = 2;
	m3 = 2;
	m4 = 1;
	m5 = 2;
	m6 = 1;
	
	printf("Running idxtest\n");
	
	idx6 b(m1, m2, m3, m4, m5, m6);
	idx a(m1, m2, m3, m4, m5, m6);
	for(int i=0; i < m1; i++) {
	for(int j=0; j < m2; j++) {
	for(int k=0; k < m3; k++) {
	for(int l=0; l < m4; l++) {
	for(int m=0; m < m5; m++) {
	for(int n=0; n < m6; n++) {
		int s = b(i,j,k,l,m,n);
		int ri,rj,rk,rl,rm,rn;
		a.rev(s, &ri,&rj,&rk,&rl,&rm,&rn);
		printf("b(%d,%d,%d,%d,%d,%d)=rev(%d,%d,%d,%d,%d,%d)=(%s,%s,%s,%s,%s,%s) ist pos %d=%d\n",
		       i, j, k, l, m, n,
		       ri, rj, rk, rl, rm, rn,
		       (ri==i ? "Y" : " "),
		       (rj==j ? "Y" : " "),
		       (rk==k ? "Y" : " "),
		       (rl==l ? "Y" : " "),
		       (rm==m ? "Y" : " "),
		       (rn==n ? "Y" : " "),
		       s, a(i,j,k,l,m,n) );
	}}}}}}
	
	/* The above output can easily be checked for correctness on the command line:
	 * Just type
	 *    ./kernelIndexTest | grep -v 'Y'
	 * to see if there are wrong lines. In case of correctness, nothing is printed.
	 **/
	
	for(int i=0; i < a.size; i++) {
		int i0, i1, i2, i3, i4, i5;
		a.rev(i, &i0,&i1,&i2,&i3,&i4,&i5);
		printf("a(%d,%d,%d,%d,%d,%d)=%d=%d=%d\n", i0,i1,i2,i3,i4,i5, i, a(i0,i1,i2,i3,i4,i5), b(i0,i1,i2,i3,i4,i5));
	}
}

void printat(const char* name, int i, int j, int k, double val) {
	printf("%s(%d,%d,%d) = %f\n", name, i, j, k, val);
}

void arraytest() {
	printf("Running arraytest\n");
	
	darray L(2,3);
	// use as setter
	L(1,2) = 7;
	// or as getter:
	double val = L(1,2);
	printat("owned", 1, 2, 0, L(1,2));
	
	// or with not owned data:
	double storage[5*7*9];
	dshadow convenient(storage, 5, 7, 9);
	convenient(2,3,8) = 10;
	printat("notowned", 2, 3, 8, convenient(2,3,8));
}

int main() {
	idxtest();
	arraytest();
}