#include <fstream>

#include "TrellisAid.h"

#define EXTRA_PRINTING false
#define PRINT_VITERBI_RESULTS_OFTEN false

namespace TrellisAid {
  void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                    vector<Edge *> *all_edges,
                    const vector<string> &observed_data,
                    const vector<string> &tag_list) {
    if (EXTRA_PRINTING)
      cout << "Building trellis." << endl;
    // Store the last column of nodes to add links to in prev_nodes. Accumulate
    // the next set of prev_nodes in future_prev_nodes.
    vector<Node *> prev_nodes, future_prev_nodes;

    Node *start_node = new Node("START_NODE", 0);
    nodes->push_back(start_node);

    int topol_index = 1;
    for (int i = 0; i < observed_data.size(); ++i)  {
      future_prev_nodes.clear();
      if (EXTRA_PRINTING)
        cout << "obs " << observed_data.at(i) << "--\n";
      for (int j = 0; j < tag_list.size(); ++j) {
        if (EXTRA_PRINTING)
          cout << Basic::Tab(1) << "tag " << tag_list.at(j) << "--\n";
        // Encode unique names for each node. Note that 'first' and 'second' are
        // useful for retrieving the 'sister' node if you only have one. Also
        // note the hack that enables us to extract the tag and the word.
        string the_tag = tag_list.at(j);
        string the_word = observed_data.at(i);
        stringstream ss;
        ss << the_tag << "#TAG#" << the_word << "#WORD#" << i << j << "first";
        Node *n1 = new Node(ss.str(), topol_index, the_tag, the_word);
        ss.clear(); ss.str("");
        ss << the_tag << "#TAG#" << the_word << "#WORD#" << i << j << "second";
        Node *n2 = new Node(ss.str(), topol_index + 1, the_tag, the_word);
        if (EXTRA_PRINTING) {
          cout << Basic::Tab(2) << "node name: " << n1->repr() << endl;
          cout << Basic::Tab(2) << "node name: " << n2->repr() << endl;
        }
        nodes->push_back(n1);
        nodes->push_back(n2);
        if (topol_index == 1) {
          Notation notation_obj("P", {the_tag});
          Edge *e = new Edge(notation_obj, start_node, n1);
          e->set_type(Edge::EdgeType::SINGLE_TAG);
          all_edges->push_back(e);
        } else {
          for (Node *p : prev_nodes) {
            // E.g., P(t2|t1)
            Notation notation_obj("P", {the_tag}, Notation::GIVEN_DELIM,
                                  {p->tag});
            if (EXTRA_PRINTING)
              cout << Basic::Tab(2) << "new edge: " << notation_obj << endl;
            Edge *e = new Edge(notation_obj, p, n1);
            e->set_type(Edge::EdgeType::LANGUAGE_MODEL);
            select_edges->push_back(e);  // Make LM probs also updatable.
            all_edges->push_back(e);
          }
        }
        future_prev_nodes.push_back(n2);
        Notation notation_obj("P", {the_word}, Notation::GIVEN_DELIM,
                              {the_tag});
        Edge *e = new Edge(notation_obj, n1, n2);
        e->set_type(Edge::EdgeType::CHANNEL);
        select_edges->push_back(e);
        all_edges->push_back(e);
      }
      topol_index += 2;
      prev_nodes = future_prev_nodes;
    }

    // To the end point.
    Node *end_node = new Node("END_NODE", topol_index);
    nodes->push_back(end_node);
    for (Node *p : prev_nodes) {
      Edge *e = new Edge(NotationConstants::p1, p, end_node);
      e->set_type(Edge::EdgeType::CONSTANT);
      all_edges->push_back(e);
    }

    if (EXTRA_PRINTING)
      cout << "Done building trellis." << endl;
  }

  void DestroyTrellis(vector<Node *> *nodes, vector<Edge *> *all_edges) {
    // Deletes nodes and edges.
    for (Node *n : *nodes) {
      delete n;
    }
    for (Edge *e : *all_edges) {
      delete e;
    }
  }

  void Viterbi(const map<Notation, double> &data, const vector<Node *> &nodes,
               const vector<string> observed_data) {
    // Key: string representation of node; Value: best value of P(t, w) so far
    // to that node. Best P(t, w) is stored in opt.at(last node).
    map<string, double> opt;
    // Key: string representation of node; Value: previous node string repr.
    map<string, string> best_path;

    // Set all nodes except start to default low value.
    if (EXTRA_PRINTING) {
      cout << "Initializing optimal values." << endl;
    }
    for (Node *n : nodes) {
      opt[n->repr()] = -DBL_MAX;
    }
    opt[nodes.front()->repr()] = log(1); // Max value probability.

    // Run through trellis. Topological order assumed.
    if (EXTRA_PRINTING) {
      cout << "About to run through trellis." << endl;
    }
    for (int i = 0; i < nodes.size(); ++i) {
      Node *current_node = nodes.at(i);
      for (Edge *e : current_node->child_edges) {
        Node *next = e->dest;
        try {
          double new_val = opt.at(current_node->repr()) + data.at(e->repr());
          // Underflow error originally happened here.
          if (new_val > opt.at(next->repr())) {
            opt[next->repr()] = new_val;
            best_path[next->repr()] = current_node->repr();
          }
        } catch (out_of_range &e) {
          cerr << "Out of range error in Viterbi while going through trellis: "
            << e.what() << endl;
          exit(0);
        }
      }
    }
    if (EXTRA_PRINTING) {
      cout << "Done running through trellis." << endl;
    }

    // Output best result following backpointer map.
    if (EXTRA_PRINTING) {
      cout << "Printing result from backpointer map." << endl;
    }
    vector<string> best_tag_seq;
    vector<string> assoc_word_seq;
    string next_node_repr = nodes.back()->repr();

    // Follow the backpointers of best_path len(observed_data) times.
    for (int i = 0; i < observed_data.size(); ++i) {
      string name;
      try {
        name = best_path.at(next_node_repr);
      } catch (out_of_range &e) {
        // This can happen when best_path doesn't set the node, which can happen
        // when, in the main Viterbi loop above, new_val > opt.at(next.repr())
        // isn't satisfied when it should be. In particular, if new_val isn't
        // actually improved upon - the + data.at(e->repr()) doesn't add
        // anything - then you have this issue!
        cerr << "Out of range error in Viterbi while getting name: " <<
          next_node_repr << " | " << e.what() << endl;
      }
      string tag = name.substr(0, name.find("#TAG#"));
      string word = name.substr(name.find("#TAG#") + 5,
                                name.find("#WORD#") - (name.find("#TAG#") + 5));
      best_tag_seq.push_back(tag);
      assoc_word_seq.push_back(word);
      try {
        next_node_repr = best_path.at(name); // Skip sister node.
      } catch (out_of_range &e) {
        cerr << "Out of range error in Viterbi while getting next " <<
          "node from best path: " << name << " | " << e.what() << endl;
      }
    }
    // Reverse these because they were retrieved backwards.
    reverse(best_tag_seq.begin(), best_tag_seq.end());
    reverse(assoc_word_seq.begin(), assoc_word_seq.end());

    double best_prob_pTAndW;
    best_prob_pTAndW = opt.at(nodes.back()->repr());
    Notation n_best_match_pTAndW("P", assoc_word_seq, Notation::AND_DELIM,
        best_tag_seq);
    Notation n_best_match_pTGivenW("P", best_tag_seq, Notation::GIVEN_DELIM,
        assoc_word_seq);
    cout << "\n--Viterbi results--\n";
    stringstream ss;
    for (int i = 0; i < best_tag_seq.size(); ++i) {
      ss << best_tag_seq.at(i);
    }
    string best_match_pTAndW_str = ss.str();
    cout << "The highest probability found belongs to " << n_best_match_pTAndW
        << ": " << best_prob_pTAndW << endl;
    cout << "Best matching tag sequence: " << best_match_pTAndW_str << endl;
  } // End Viterbi

  void ForwardBackwardAndViterbi(const int num_iterations,
                                 const set<string> &word_list, 
                                 const vector<string> &tag_list,
                                 const vector<Node *> &nodes,
                                 const vector<Edge *> &select_edges,
                                 const vector<Edge *> &all_edges,
                                 map<Notation, double> *data,
                                 const vector<string> observed_data) {
    Notation nObsSeq("P", observed_data, Notation::SEQ_DELIM);
    ofstream fout;
    if (PRINT_VITERBI_RESULTS_OFTEN) {
      fout.open("observed_data_probabilities.txt");
    }
    if (EXTRA_PRINTING)
      cout << "Beginning Forward-Backward." << endl;

    // PRECONDITION: The order of nodes/edges is already in topological order.
    map<string, double> alpha;  // Sum of all paths from start to this node.
    map<string, double> beta;  // Sum of all paths from this node to final.
    alpha[nodes.front()->repr()] = log(1);
    beta[nodes.back()->repr()] = log(1);
    for (int iter_count = 0; iter_count < num_iterations; ++iter_count) {
      if (EXTRA_PRINTING)
        cout << "Forward pass... ";
      // Forward pass. Assumes start node is at i = 0.
      try {
        for (int i = 1; i < nodes.size(); ++i) {
          double sum = -DBL_MAX;
          for (Edge *e : nodes[i]->parent_edges) {
            if (EXTRA_PRINTING) {
              cout << Basic::Tab(1) << "Accumulating new alpha value by adding " 
                << sum << " (previous sum) to \n[" << alpha[e->src->repr()] << 
                " (alpha of src) and " << data->at(e->repr()) << 
                " (data at edge)]: " << "[" << 
                  alpha[e->src->repr()] + data->at(e->repr()) << "]" << endl;
            }
            sum = Basic::AddLogs(sum,
                alpha[e->src->repr()] + data->at(e->repr()));
            if (EXTRA_PRINTING) {
              cout << Basic::Tab(1) << "Sum is now " << sum << endl;
            }
          }
          if (EXTRA_PRINTING) {
            cout << Basic::Tab(2) << "Resulting alpha value for " << nodes[i]->repr() <<
                ": " << sum << endl;
          }
          alpha[nodes[i]->repr()] = sum;
        }
      } catch (out_of_range &e) {
        cerr << "Out of range error in forward pass: " << e.what() << endl;
        exit(1);
      } catch (exception &e) {
        cerr << "Issue in forward pass: " << e.what() << endl;
      }

      if (EXTRA_PRINTING) {
        cout << "Backward pass... ";
      }
      // Backward pass. Assumes end node is at i = size - 1.
      for (int i = nodes.size() - 2; i >= 0; --i) {
        double sum = -DBL_MAX;
        for (Edge *e : nodes[i]->child_edges) {
          sum = Basic::AddLogs(sum,
              beta[e->dest->repr()] + data->at(e->repr()));
        }
        if (EXTRA_PRINTING) {
          cout << Basic::Tab(1) << "Beta value for " << nodes[i]->repr() << ": "
            << sum << endl;
        }
        beta[nodes[i]->repr()] = sum;
      }

      if (EXTRA_PRINTING) {
        cout << "Counting pass... " << endl;
      }

      // Counting pass. e->notation.first[0] gets the word if select edge
      // represents a channel probability; if LM probability, then tag. This is
      // because in BuildTrellis edges are constructed with a Notation object
      // with this form: P(w|t) or P(t|t).

      // First reset the counts.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->notation.second[0]}, Notation::AND_DELIM,
                             {e->notation.first[0]});
        (*data)[n_count_key] = -DBL_MAX; // log(0)
      }

      // Key: tag. Value: total fractional count associated with that tag.
      // Example: When normalizing probabilities - i.e., computing P(b|x) from
      // C(x, b), we do C(x, b)/sum of all C(x, w) over all w. The denominator
      // is the total fractional count for the tag x. This is for _tw. For _tt,
      // it is the same idea, except sum C(x, t) over all tags t.
      // Tag is accessed via e->notation.second[0].
      unordered_map<string, double> total_fract_counts_tw; // Sum C(t, w_i)
      unordered_map<string, double> total_fract_counts_tt; // Sum C(t, t_i)

      // Iterate over select edges to update count_keys, used for updating
      // probabilities later.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation count_key("C", {e->notation.second[0]}, Notation::AND_DELIM,
                           {e->notation.first[0]});
        if (EXTRA_PRINTING) {
          cout << Basic::Tab(1) << "Getting count key from edge " << e->repr()
              << ": " << count_key << endl;
          cout << Basic::Tab(2) << alpha[e->src->repr()] <<
              "<-alpha, edge prob->" << data->at(e->repr()) << endl;
          cout << Basic::Tab(2) << beta[e->dest->repr()] << "<-beta, end->" <<
              alpha[nodes.back()->repr()] << endl;
        }
        if (EXTRA_PRINTING) {
          cout << Basic::Tab(2) << "Data at count key before: " << 
            (*data)[count_key] << endl;
          double val = Basic::AddLogs((*data)[count_key],
                         alpha[e->src->repr()] +
                         data->at(e->repr()) +
                         beta[e->dest->repr()] -
                         alpha[nodes.back()->repr()]); // positive... BAD? or okay? TODO: Issue #1
          cout << Basic::Tab(2) << "Value about to add -> " << val << endl;
          if (val > 0) {
            cout << "Bad AddLog: result was positive." << endl;
          }
          cout << Basic::Tab(2) << "Final data at count key: " << 
            (*data)[count_key] << endl;
          (*data)[count_key] = val;
        } else {
          (*data)[count_key] = Basic::AddLogs((*data)[count_key],
                                              alpha[e->src->repr()] +
                                              data->at(e->repr()) +
                                              beta[e->dest->repr()] -
                                              alpha[nodes.back()->repr()]);
        }
      }
      // Compute total_fract_counts by summing C(tag, tag_i).
      for (string tag : tag_list) {
        double tag_total_frac_count = -DBL_MAX;
        for (string tag_i : tag_list) {
          Notation count_of_tag_and_tag_i("C", {tag}, Notation::AND_DELIM,
                                          {tag_i});
          tag_total_frac_count = Basic::AddLogs(tag_total_frac_count,
              (*data)[count_of_tag_and_tag_i]);
          // TODO: Debugging for Issue #1. Remove.
//           cout << "  just added " << (*data)[count_of_tag_and_tag_i] << endl;
        }
        total_fract_counts_tt[tag] = tag_total_frac_count;
        // TODO: Debugging for Issue #1. Remove.
//         cout << "tt tag_total_frac_count for " << tag << ": " <<
//           tag_total_frac_count << endl;
      }
      // Compute total_fract_counts by summing C(tag, word_i).
      for (string tag : tag_list) {
        double tag_total_frac_count = -DBL_MAX;
        for (string word_i : word_list) {
          Notation count_of_tag_and_word_i("C", {tag}, Notation::AND_DELIM,
                                          {word_i});
          tag_total_frac_count = Basic::AddLogs(tag_total_frac_count,
              (*data)[count_of_tag_and_word_i]);
          // TODO: Debugging for Issue #1. Remove.
//           cout << "  just added " << (*data)[count_of_tag_and_word_i] << endl;
        }
        total_fract_counts_tw[tag] = tag_total_frac_count;
        // TODO: Debugging for Issue #1. Remove.
//         cout << "tw tag_total_frac_count for " << tag << ": " <<
//           tag_total_frac_count << endl;
      }

      // Normalization step: Update the unknown probabilities that we want to
      // find. Use them in the next iteration.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->notation.second[0]}, Notation::AND_DELIM,
                             {e->notation.first[0]});
        if (e->type == Edge::EdgeType::CHANNEL) {
          (*data)[e->repr()] = (*data)[n_count_key] -
            total_fract_counts_tw.at(e->notation.second[0]);
        } else if (e->type == Edge::EdgeType::LANGUAGE_MODEL) {
          (*data)[e->repr()] = (*data)[n_count_key] -
            total_fract_counts_tt.at(e->notation.second[0]);
        } else {
          cerr << "Unknown edge type in normalization step of Forward-Backward:"
            << " " << e->type << endl;
          exit(0);
        }
      }

      // Update probability of observed data sequence. This should increase
      // with each iteration.
      (*data)[nObsSeq] = alpha[nodes.back()->repr()];

      if (PRINT_VITERBI_RESULTS_OFTEN) {
        // Print viterbi.
        cout << iter_count + 1 << ": \n";
        TrellisAid::Viterbi(*data, nodes, observed_data);
        // Print P(obs).
        cout << "P(observed sequence): " << alpha[nodes.back()->repr()] <<
                                            endl << endl;
        fout << alpha[nodes.back()->repr()] << endl;
      }
    }
    if (EXTRA_PRINTING) {
      cout << "Done with Forward-Backward. Proceeding to Viterbi." << endl;
    }

    // By this point, P(A|X), P(A|Y), etc. have been maximized thanks to alpha
    // and beta values in the forward-backward passes. The counting pass
    // collected fractional counts of e.g. C(X|A), which were then used to
    // update the "Given" probabilities (P(A|X), P(A|Y), etc.). Now we use
    // Viterbi to find the highest-probability path based on the collected
    // probabilities and print the result.
    if (PRINT_VITERBI_RESULTS_OFTEN) {
      cout << "Final results----" << endl;
    }
    TrellisAid::Viterbi(*data, nodes, observed_data);
  }
} // end namespace TrellisAid

