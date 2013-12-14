#ifndef EDGE_H_
#define EDGE_H_

#include <iostream>
#include <vector>

#include "Node.h"
#include "Notation.h"
#include "GraphAid.h"

using namespace std;

struct Node;

struct Edge {
  enum EdgeType {
    CHANNEL,  // P(w|t)
    LANGUAGE_MODEL,  // P(t_2|t_1)
    SINGLE_TAG,  // P(t)
    CONSTANT,  // P(1)
  };
  EdgeType type;
  Notation notation;
  Node *src, *dest;
  Edge(Notation n, Node *src, Node *dest);
  void set_type(EdgeType type) {
    this->type = type;
  }
  Notation repr() {
    return notation;
  }
};

#endif  // EDGE_H_
