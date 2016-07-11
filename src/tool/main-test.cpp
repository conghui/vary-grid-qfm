#include <cmath>
#include <vector>
#include <iostream>
#include <cmath>
#include "logger.h"
extern "C" {
#include <rsf.h>
}

int main(int argc, char* argv[]) {
  bool verb;

  sf_init(argc, argv);

  if (!sf_getbool("verbose", &verb)) verb = 0;
  INFO() << "program exit normally";
  return 0;
}

