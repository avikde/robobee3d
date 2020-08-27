
#include <iostream>
#include "wlcontroller.hpp"

int main(int argc, char **argv) {
	WLController wlc;
	
	u_t u0;
	w_t p0, h0, pdes, Qdiag;
	Qdiag.head<3>().fill(1.0f);
	Qdiag.tail<3>().fill(0.1f);
	pdes << 0, 0, 20, 0, 0, 0;
	p0.setZero();
	const float mb = 100;
	const float g = 9.81e-3f;
	h0 << 0, 0, mb * g, 0, 0, 0;
	u0 << 140.0, 0., 0., 0.;

	wlc.setLimits(u_t(90., -0.5, -0.2, -0.1), u_t(160., 0.5, 0.2, 0.1), u_t(5., 0.01, 0.01, 0.01));
	wlc.setWeight(Qdiag);

	auto u = wlc.update(u0, p0, h0, pdes);

	std::cout << u.transpose() <<  " hello\n";
	return 0;
}
