#include <cassert>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <utility>

#include "EMViterbiPackage/BasicHelper.h"
#include "EMViterbiPackage/Notation.h"

using namespace std;

int main(int argc, char *argv[]) {
  Notation n1("P", {"A"}, Notation::GIVEN_DELIM, {"B"});
  map<Notation, double> data;

  cout << -DBL_MAX << endl;
  cout << -DBL_MAX + -DBL_MAX << endl;

  return 0;
}
