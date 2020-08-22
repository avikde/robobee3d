
#include "wlcontroller.h"
#include <stdio.h>

int main() {
	float A[] = {1,2,3,4,5,6,7,8,9};
	float b[] = {2,43,324};
	float C[3];

	matMult(C, A, b, 3, 1, 3, 1, 1, 0);
	for (int i = 0; i < 3; ++i)
		printf("%f,", C[i]);
	printf("\n");
	
	// Simple test
	float u0[4] = {140.0, 0., 0., 0.};
	float h0[6] = {0, 0, 0.981, 0, 0, 0};
	float pdotdes[6] = {0, 0, 10, 0, 0, 0};
	float u[4];

	wlControllerUpdate(u, u0, h0, pdotdes);

	for (int i = 0; i < 4; ++i)
		printf("%f,", u[i]);
	printf("\n");

	return 0;
}
