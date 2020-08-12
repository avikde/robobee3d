
#include <iostream>
#include "wlcontroller.hpp"

int main(int argc, char **argv) {
  float popts[] = {0};
  WLController wlc(popts);
  
  u_t u0;
  w_t p0, h0, pdes, Qdiag;

  std::cout << wlc.update(u0, p0, h0, pdes, Qdiag) <<  "hello\n";
  return 0;
}
