#ifndef TRELLIS_H_
#define TRELLIS_H_

#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <map>
#include <set>
#include <unordered_map>
#include <sstream>
#include <vector>

#include "BasicHelper.h"
#include "Notation.h"
#include "NotationConstants.h"
#include "Node.h"
#include "Edge.h"

using namespace std;

struct Node;
struct Edge;

namespace TrellisAid {
  // WARNING: Creates data on heap. Call DestroyTrellis when done. Post:
  // 'nodes' points to a vector where front() is the start node, back() is the
  // end, and the vector lists the nodes in topological order. Each node has a
  // unique name. 'edges' points to a vector of corresponding edges with
  // representations like P(A|X). **Nodes and edges are in topological order!**.
  void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                    vector<Edge *> *all_edges, const vector<string>
                    &observed_data, const vector<string> &tag_list);
  void DestroyTrellis(vector<Node *> *nodes, vector<Edge *> *all_edges);

  // Traverses the trellis with the Viterbi algorithm to find the best matching
  // tag sequence and prints the results.
  // Returns the best matching string.
  string Viterbi(const map<Notation, double> &data,
                        const vector<Node *> &nodes,
                        const vector<string> observed_data);
  // Runs ForwardBackward 'num_iterations' times to determine best probabilities
  // along edges select_edges and then calls Viterbi(). Updates data's count
  // keys (e.g., C(X,A)) in the process.
  void ForwardBackwardAndViterbi(const int num_iterations,
                                 const set<string> &obs_symbols, // word_list
                                 const vector<string> &tag_list,
                                 const vector<Node *> &nodes,
                                 const vector<Edge *> &select_edges,
                                 const vector<Edge *> &all_edges,
                                 map<Notation, double> *data,
                                 vector<double> *increasing_probs,
                                 string *best_match,
                                 const vector<string> observed_data);
}

#endif  // TRELLIS_H_
