// Takes learned_probabilities_for_best_run.txt output from Main and reformats
// in a variety of ways.  E.g., log probs -> normal probs
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <vector>

bool REVERT_PROBS = true;  // Convert log probs to regular probs.
// Command line parameters.
bool FEWER_DECIMAL_PLACES = false;
int NUM_DEC_PLACES = 3;

using namespace std;

void OrganizeIntoList(const vector<pair<double, string> > &pr_list,
    map<string, vector<pair<double, string> > > *organized_list) {
  // ASSUMPTION: | and ) will never be part of a tag/word.
  for (auto aPair : pr_list) {
    string probStr = aPair.second;
    int start = probStr.find("|");
    int end = probStr.find(")");
    string tag = aPair.second.substr(start + 1, -1 + end - start);

    auto it = organized_list->find(tag);
    if (it == organized_list->end()) {
      vector<pair<double, string> > inner;
      inner.push_back(aPair);
      (*organized_list)[tag] = inner;
    } else {
      (*organized_list)[tag].push_back(aPair);
    }
  }
}

void PrintOrganizedData(map<string, vector<pair<double, string> > >
                        organized_list) {
  if (FEWER_DECIMAL_PLACES)
    cout << fixed << setprecision(NUM_DEC_PLACES);
  for (auto x : organized_list) {
    cout << "Tag: " << x.first << endl;
    for (auto y : x.second) {
      if (FEWER_DECIMAL_PLACES) {
        if (REVERT_PROBS){
          if (exp(y.first) >= pow(10, -NUM_DEC_PLACES))
            cout << "\t" << y.second << " " << exp(y.first) << endl;
        } else {
          if (y.first >= pow(10, -NUM_DEC_PLACES))
            cout << "\t" << y.second << " " << y.first << endl;
        }
      } else {
        if (REVERT_PROBS)
          cout << "\t" << y.second << " " << exp(y.first) << endl;
        else
          cout << "\t" << y.second << " " << y.first << endl;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << "Usage: ./<exec> learned_probabilities_for_best_run.txt [num decimal places]" << endl;
    return 0;
  } else if (argc == 3) {
    stringstream ss;
    ss << argv[2];
    ss >> NUM_DEC_PLACES;
    FEWER_DECIMAL_PLACES = true;
  }
  string filename = argv[1];
  ifstream fin;
  fin.open(filename.c_str());
  if (fin.fail()) {
    cerr << "Could not open file " << filename << endl;
    return 1;
  } else {
    vector<pair<double, string> > ch_pr_list;
    vector<pair<double, string> > lm_pr_list;
    bool channel_probs = true;
    while (true) {
      string line;
      getline(fin, line);
      if (fin.eof())
        break;
      // Ignore title strings like "Channel probabilities"
      if (line.find("probabilities") == string::npos) {
        if (channel_probs) {
          double prob;
          string notation;
          stringstream ss(line);
          ss >> notation >> prob;
          ch_pr_list.push_back(make_pair(prob, notation));
        } else {
          double prob;
          string notation;
          stringstream ss(line);
          ss >> notation >> prob;
          lm_pr_list.push_back(make_pair(prob, notation));
        }
      } else if (line.find("Language model probabilities") != string::npos) {
        channel_probs = false;
      }
    }
    fin.close();
    // Reorganize and print out.
    // Key: tag, Value: list of pairs of notations and associated probs.
    map<string, vector<pair<double, string> > > organized_list;
    OrganizeIntoList(ch_pr_list, &organized_list);
    cout << "Channel probabilities:" << endl;
    PrintOrganizedData(organized_list);
    organized_list.clear();
    OrganizeIntoList(lm_pr_list, &organized_list);
    cout << "Language model probabilities:" << endl;
    PrintOrganizedData(organized_list);
  }
  return 0;
}
