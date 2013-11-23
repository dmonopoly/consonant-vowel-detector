#include <cassert>
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
#include "EMViterbiPackage/Node.h"
#include "EMViterbiPackage/Edge.h"
#include "EMViterbiPackage/NotationConstants.h"
#include "EMViterbiPackage/NotationHelper.h"
#include "EMViterbiPackage/TrellisAid.h"

#include "CypherReader.h"
#include "TagGrammarFinder.h"

using namespace std;

#define NUMBER_ITERATIONS 50
#define EXTRA_PRINTING false
#define SHOW_PROBS_BEFORE_EM false

void DisambiguateDuplicates(const set<string> &obs_symbols,
                            vector<string> *tag_list,
                            map<Notation, double> *data) {
  // Alters all tags that have overlapping meanings with observed symbols.
  // Disambiguates by adding an '. Updates data accordingly.
  for (auto i = obs_symbols.begin(); i != obs_symbols.end(); ++i) {
    for (auto j = tag_list->begin(); j != tag_list->end(); ++j) {
      if (*i == *j) {
        stringstream ss;
        ss << *j << "'";
        string old_val = *j;
        *j = ss.str();
        string new_val = ss.str();
        vector<pair<Notation, Notation> > values_to_replace; // old key, new key
        Notation key;
        for (auto data_pair = data->begin(); data_pair != data->end();
             ++data_pair) {
          key = data_pair->first;
          size_t pos = key.repr().find(old_val);
          if (pos != string::npos) {
            Notation old_key = key;
            NotationHelper::ReplaceSymbol(old_val, new_val, &key);
            values_to_replace.push_back(make_pair(old_key, key));
          }
        }
        for (auto i = values_to_replace.begin(); i != values_to_replace.end();
             ++i) {
          Notation old_key = i->first;
          Notation new_key = i->second;
          if (old_key.repr() != new_key.repr()) {
            (*data)[new_key] = (*data)[old_key];
            data->erase(old_key);
          }
        }
      }
    }
  }
}

void PrepareObsTagProbs(const vector<string> &observed_data,
                        const vector<string> &tag_list,
                        const set<string> &obs_symbols,
                        map<Notation, double> *data) {
  // Sets initial observed char / tag probabilities. Can seed to uniform
  // probability (1/# of unique observed symbols) or randomize.
  for (auto obs = observed_data.begin(); obs != observed_data.end(); ++obs) {
    for (auto tag = tag_list.begin(); tag != tag_list.end(); ++tag) {
      Notation nObsTagProb("P", {*obs}, Notation::GIVEN_DELIM, {*tag});
      (*data)[nObsTagProb] = (double) 1/obs_symbols.size();
    }
  }
  // Deal with spaces in the substitute table: set
  // P(any obs except space |space tag) to 0 and P(space obs|space tag) to 1.
  // "_" means "space" in the ciphertext letter sequence.
  // "_" means "pause" in the spoken spanish plaintext. Altered to "_'" for
  // uniqueness in DisambiguateDuplicates.
  for (auto obs = observed_data.begin(); obs != observed_data.end(); ++obs) {
    Notation nAnythingGivenSpaceTag("P", {*obs}, Notation::GIVEN_DELIM, {"_'"});
    (*data)[nAnythingGivenSpaceTag] = 0;
  }
  Notation nSpaceTagGivenSpaceObs("P", {"_"}, Notation::GIVEN_DELIM, {"_'"});
  (*data)[nSpaceTagGivenSpaceObs] = 1;
}

void SeedNotationConstants(map<Notation, double> *data) {
  // Also seed NotationConstants.
  (*data)[NotationConstants::p1] = 1;
}

void ChangeAbsoluteProbsToLogProbs(map<Notation, double> *data) {
  for (auto it = data->begin(); it != data->end(); ++it) {
    if (it->first.is_probability()) {
      it->second = log(it->second);
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    // TODO: allow trigam counts file + probabilities file too.
    cerr << "Usage: ./<exec> <bigram-counts-file> <cyphertext>" << endl;
    return 0;
  }
  string filename_for_bigrams = argv[1];
  map<Notation, double> data;  // Storage for log probabilities and counts.
  vector<string> tag_list;
  bool found = TagGrammarFinder::GetTagGrammarFromOrganizedRows(
      filename_for_bigrams, &data, &tag_list);
  if (!found)
    return 0;
  else if (EXTRA_PRINTING)
    cout << "Found tag grammar.\n";

  string filename_for_cypher = argv[2];
  vector<string> observed_data;
  set<string> obs_symbols;
  vector<Node *> nodes;
  vector<Edge *> edges_to_update;
  vector<Edge *> all_edges; // for deletion later
  bool got_obs_data = CypherReader::GetObservedData(filename_for_cypher,
                                                    &observed_data,
                                                    &obs_symbols);
  if (!got_obs_data)
    return 0;
  else if (EXTRA_PRINTING)
    cout << "Found cyphertext.\n";

  DisambiguateDuplicates(obs_symbols, &tag_list, &data);
  PrepareObsTagProbs(observed_data, tag_list, obs_symbols, &data);
  SeedNotationConstants(&data);
  ChangeAbsoluteProbsToLogProbs(&data);
  clock_t t = clock();
  TrellisAid::BuildTrellis(&nodes, &edges_to_update, &all_edges, observed_data,
                           tag_list);
  if (EXTRA_PRINTING) {
    cout << "Built trellis.\n";
  }

  if (SHOW_PROBS_BEFORE_EM) {
    cout << "Printing probs..." << endl;
    for (auto it = data.begin(); it != data.end(); ++it) {
      cout << it->first << " " << it->second << endl;
    }
  }

  cout << NUMBER_ITERATIONS << " iterations:" << endl;

  Notation nObsSeq("P", observed_data, Notation::SEQ_DELIM);
  TrellisAid::ForwardBackwardAndViterbi(NUMBER_ITERATIONS, nodes,
                                        edges_to_update, all_edges, &data,
                                        observed_data);
  TrellisAid::DestroyTrellis(&nodes, &all_edges);
  t = clock() - t;
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);
  return 0;
}
