
#include <iostream>
#include "wlqp.hpp"

int main(int argc, char **argv) {
  std::cout << wlqp(u_t::Zero(), dw_du_t::Zero(), w_t::Zero(), 2).transpose() <<  "hello\n";
  return 0;
}
