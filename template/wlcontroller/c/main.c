
// #include "funapprox.h"
#include "eigenc.h"
#include "wlcontroller.h"
#include <stdio.h>

void wrenchMap(float *w, const float *u) {

}
void wrenchJacMap(float *dw_du, const float *u) {
	
}

int main() {
	// FunApprox_t fa;
	// float popts[] = {-0.41498716334913877, 0.004161530470390393, -0.04921896208897451, -0.08114084916912893, 1.0, -2.1498519749986776e-05, 0.0007452193674740687, -0.000592201745189127, 1.0, 2.489837075232938, -0.8220389384978928, 1.0, -0.6424779320800695, 1.0, 1.0, 0.04853925386307669, -0.0005410195533586718, 0.03325717298367986, 0.11402108349065512, 1.0, 2.9591070779816484e-06, -0.0003656669558602931, -0.0006235124178038529, 1.0, -0.15022060357201647, 0.0003189736297552722, 1.0, 0.09522768136462195, 1.0, 1.0, -1.0443154298902466, 0.008730954429386921, 3.774256296882532, -6.565041851772019, 1.0, 7.255944919680254e-05, -0.024440375440139876, 0.03221142824195457, 1.0, -23.321197610424026, 3.865481128472532, 1.0, 1.7626691411728337, 1.0, 1.0, -39.162566526689496, 0.4347272529556003, 35.20801123124533, 8.960105550238962, 1.0, -0.002179040590683337, -0.20825719312770322, -0.2710167102094739, 1.0, 98.4856348629078, 34.454898506387536, 1.0, -68.47304863790173, 1.0, 1.0, 3.86697945239634, -0.054406534814573505, 4.136310537263415, -2.40430210592086, 1.0, 0.00034641731566390497, -0.05418716099028047, 0.011121246757241961, 1.0, 5.142757475687521, -10.274310705200136, 1.0, 11.418953851491452, 1.0, 1.0, 1.5868629373519454, -0.013584604586031352, 8.897152221477027, 6.856427471376807, 1.0, 1.2282136287999407e-05, -0.05776369789309159, -0.0360868798609794, 1.0, 8.468979086282319, 0.7365484310116243, 1.0, 0.44877875839694553, 1.0, 1.0};
	// funApproxInit(&fa, popts);

	float A[] = {1,2,3,4,5,6,7,8,9};
	float b[] = {2,43,324};
	float C[3];

	matMult(C, A, b, 3, 1, 3, 1, true, false);
	for (int i = 0; i < 3; ++i)
		printf("%f,", C[i]);
	
	WLQP_t wlqp;
	wlqpInit(&wlqp);
	float u0[4] = {140.0, 0., 0., 0.};
	const float mb = 100;
	const float g = 9.81e-3f;
	float h0[6] = {0, 0, mb * g, 0, 0, 0};
	float pdotdes[6] = {0, 0, 10, 0, 0, 0};
	float u[4];
	wlqpUpdate(&wlqp, u, u0, h0, pdotdes);

	return 0;
}
