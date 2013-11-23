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
  data[n1] = 1.5;
  Notation n2("P", {"A"}, Notation::GIVEN_DELIM, {"B"});
  data[n2] = .75;

  for (auto it = data.begin(); it != data.end(); ++it) {
    cout << it->first << " " << it->second << endl;
  }

  return 0;
}
