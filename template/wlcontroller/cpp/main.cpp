
#include <iostream>
#include "wlcontroller.hpp"

int main(int argc, char **argv) {
	WLController wlc;
	
	u_t u0;
	w_t h0, pdotdes;
	pdotdes << 0, 0, 10, 0, 0, 0;
	const float mb = 100;
	const float g = 9.81e-3f;
	h0 << 0, 0, mb * g, 0, 0, 0;
	u0 << 140.0f, 0., 0., 0.;

	auto u = wlc.update(u0, h0, pdotdes);

	std::cout << u.transpose() <<  " hello\n";
	return 0;
}
