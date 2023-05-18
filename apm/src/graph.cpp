#include "graph.h"

#include <iostream>

#include "general.h"
#include "rng.h"

namespace apm {

Graph::Graph(int n, double alpha, Rng &rng) : rng(rng) {
  //   this->rng = rng;
  this->n = n;
  this->alpha = alpha;
  // probability for an edge
  this->p = this->alpha / (this->n - 1);
}

Graph &Graph::generate() {
  vvInt neighbors(this->n);
  // Rng rng;
  // vvInt edgeList;

  for (int i = 0; i < n; i++) {
    for (int j = (i + 1); j < n; j++) {
      // draw edge with probability p
      if (this->rng.real() < p) {
        neighbors.at(i).push_back(j);
        neighbors.at(j).push_back(i);
        // edgeList.push_back({i, j});
      }
    }
  }

  this->neighbors = move(neighbors);
  // this->edgeList = move(edgeList);

  return *this;
}

Graph &Graph::print() {
  std::cout << "[Verbose]\t[graph]\tfull graph (p=" << this->p << "):\n";
  for (size_t i = 0; i < this->neighbors.size(); i++) {
    std::cout << i << ": ";
    for (int n : this->neighbors.at(i)) {
      std::cout << n << ' ';
    }
    std::cout << std::endl;
  }

  // if (this->edgeList.size() > 0) {
  //   std::cout << "[Verbose]\t[graph]\tedge list: ";
  //   for (size_t i = 0; i < this->edgeList.size(); i++) {
  //     std::cout << "(" << this->edgeList.data()[i][0] << ", "
  //               << this->edgeList.data()[i][1] << "), ";
  //   }
  //   std::cout << '\n';
  // }

  return *this;
}

Graph &Graph::printGiantComponent() {
  std::cout << "[Verbose]\t[graph]\tfull graph (p=" << this->p << "):\n";
  for (size_t i = 0; i < this->giantComponent.size(); i++) {
    std::cout << i << ": ";
    for (int n : this->giantComponent.at(i)) {
      std::cout << n << ' ';
    }
    std::cout << std::endl;
  }

  // if (this->edgeList.size() > 0) {
  //   std::cout << "[Verbose]\t[graph]\tedge list: ";
  //   for (size_t i = 0; i < this->edgeList.size(); i++) {
  //     std::cout << "(" << this->edgeList.data()[i][0] << ", "
  //               << this->edgeList.data()[i][1] << "), ";
  //   }
  //   std::cout << '\n';
  // }

  return *this;
}

Graph &Graph::findGiantComponent(bool logging) {
  std::vector<sInt> components;
  sInt visited;
  sInt c;
  size_t giantSize = 0;
  int giantIndex = 0;

  // loop over all nodes
  for (size_t i = 0; i < this->neighbors.size(); i++) {
    c.clear();
    // check if we haven't seen the node before
    if (!visited.contains(i)) {
      // collect all nodes that are in the component of i
      Graph::addToComponent(this->neighbors, i, c, visited);
      components.push_back(c);
    }
  }

  // find the largest component
  for (size_t i = 0; i < components.size(); i++) {
    if (components.at(i).size() > giantSize) {
      giantSize = components.at(i).size();
      giantIndex = i;
    }
  }

  this->giantComponentSet = components.at(giantIndex);
  this->giantComponent = Graph::extractComponent(this->neighbors,
                                                 this->giantComponentSet);

  if (logging) {
    std::cout << "[Info]\t[graph]\tFound giant component with "
              << this->giantComponent.size()
              << " nodes.\n";
  }

  return *this;
}

Graph &Graph::computeEdgeList() {
  vvInt edgeList;

  for (size_t i = 0; i < this->giantComponent.size(); i++) {
    for (size_t j = 0; j < this->giantComponent.at(i).size(); j++) {
      edgeList.push_back({static_cast<int>(i), giantComponent.at(i).at(j)});
    }
  }

  this->edgeList = move(edgeList);

  return *this;
}

Graph &Graph::printGiantComponentSet() {
  std::cout << "[Debug]\t[graph]\tgiant component: { ";
  for (int e : this->giantComponentSet) {
    std::cout << e << ", ";
  }
  std::cout << "}\n";

  return *this;
}

void Graph::addToComponent(vvInt &graph,
                           int node,
                           sInt &component,
                           sInt &visited) {
  // check if we haven't visited the node before
  if (!visited.contains(node)) {
    // now we visited the node
    visited.insert(node);
    component.insert(node);
    // loop through the neighbors
    for (int i : graph.at(node)) {
      addToComponent(graph, i, component, visited);
    }
  }
}

void Graph::addToChainList(vvInt &graph,
                           int node,
                           vInt &component,
                           sInt &visited,
                           int &gateway) {
  // check if we haven't visited the node before
  if (!visited.contains(node)) {
    // now we visited the node
    visited.insert(node);
    component.push_back(node);
    // loop through the neighbors
    for (int i : graph.at(node)) {
      // check if the neighbors were visited before
      if (!visited.contains(i)) {
        gateway = i;
        // only continue if the node is in a chain, i.e. degree is less than 3
        if (graph.at(i).size() < 3) {
          addToChainList(graph, i, component, visited, gateway);
        }
      }
    }
  }
}

vvInt Graph::extractComponent(vvInt &graph, sInt component) {
  vvInt subgraph;
  std::vector<int> missing;
  std::vector<int> *neighbors;
  size_t numMissing = 0;
  int indexOld;

  // find the nodes to exclude
  for (size_t i = 0; i < graph.size(); i++) {
    if (!component.contains(i)) {
      missing.push_back(i);
    }
  }
  numMissing = missing.size();

  // copy the neighbor lists for all the nodes in the giant component
  for (int node : component) {
    subgraph.push_back(graph.at(node));
  }

  // shift down the indices in the neighbor lists
  for (size_t i = 0; i < subgraph.size(); i++) {
    neighbors = &(subgraph.at(i));
    // loop through the neighbors
    for (size_t j = 0; j < neighbors->size(); j++) {
      indexOld = neighbors->at(j);
      // loop over all the missing values
      for (size_t k = 0; k < numMissing; k++) {
        // reduce the index of the neighbor if it is larger than one
        // of the missing nodes because they all got shifted down
        if (indexOld > missing.at(k)) {
          neighbors->at(j)--;
        }
      }
    }
  }

  return subgraph;
}

matrixSparse Graph::adjacencyMatrix() {
  matrixSparse adj;
  for (auto neighbors : this->giantComponent) {
    VectorSparse row(neighbors, 1.);
    adj.push_back(row);
  }
  return adj;
}

Graph &Graph::computeDegreeSequence() {
  vInt k(this->giantComponent.size());
  for (size_t i = 0; i < k.size(); i++) {
    k.at(i) = giantComponent.at(i).size();
  }
  this->degreeSequence = move(k);
  return *this;
}

Graph &Graph::printDegreeSeqence() {
  std::cout << "[Debug]\t[graph]\tdegree sequence:\n";
  for (size_t i = 0; i < (this->degreeSequence).size(); i++) {
    std::cout << i << ": " << degreeSequence.at(i) << '\n';
  }
  return *this;
}

vInt Graph::findDanglingChains(bool logging) {
  vInt endPoints;  // nodes of degree 1
  vInt danglingChains;
  int neigh;
  int nNeigh;

  // find all endpoints that are nodes of degree 1
  for (size_t i = 0; i < this->degreeSequence.size(); i++) {
    if (this->degreeSequence.at(i) == 1) {
      endPoints.push_back(i);
    }
  }

  // std::cout << "Found " << endPoints.size() << " endpoints: ";
  // for (int i : endPoints) {
  //   std::cout << i << ", ";
  // }
  // std::cout << std::endl;

  // check the number of neighbors of the endpoints
  for (size_t i = 0; i < endPoints.size(); i++) {
    neigh = this->giantComponent.at(endPoints.at(i)).at(0);
    nNeigh = this->giantComponent.at(neigh).size();
    // std::cout << "Node " << endPoints.at(i)
    //           << " has " << nNeigh << " neighbors\n";
    if (nNeigh == 2) {
      danglingChains.push_back(endPoints.at(i));
    }
  }

  if (logging) {
    std::cout << "[Info]\t[graph]\tFound " << danglingChains.size()
              << " dangling chains at nodes: ";
    for (int i : danglingChains) {
      std::cout << i << ", ";
    }
    std::cout << std::endl;
  }

  return danglingChains;
}

void Graph::getDanglingChainList(vInt &danglingChains,
                                 vInt &dcList,
                                 vInt &gwList) {
  // vInt dcList;
  // vInt gwList;
  sInt visited;
  int gateway;
  // size_t last = 0;

  for (int endPoint : danglingChains) {
    addToChainList(this->giantComponent, endPoint, dcList, visited, gateway);
    gwList.push_back(gateway);

    // std::cout << "Chain: ";
    // for (size_t idx = last; idx < dcList.size(); idx++) {
    //   std::cout << dcList.at(idx) << ", ";
    // }
    // std::cout << "gateway: " << gateway << "\n";

    // remember previous ending of the list before
    // last = dcList.size();
  }

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[graph]\tFound " << dcList.size()
            << " nodes in the all dangling chains and "
            << gwList.size() << " gateways\n";
#endif
  // return dcList;
}

void Graph::getDanglingChainSet(vInt &danglingChains,
                                sInt &dcSet,
                                sInt &gwSet) {
  vInt dcList;
  vInt gwList;
  this->getDanglingChainList(danglingChains, dcList, gwList);
  sInt dcSetLocal(dcList.begin(), dcList.end());
  sInt gwSetLocal(gwList.begin(), gwList.end());

  dcSet = std::move(dcSetLocal);
  gwSet = std::move(gwSetLocal);
}

matrixSparse Graph::transitionMatrix() {
  matrixSparse pi;
  for (size_t i = 0; i < (this->giantComponent).size(); i++) {
    VectorSparse row((this->giantComponent).at(i), 1.);
    row.divide((this->degreeSequence).at(i));
    pi.push_back(row);
  }
  return pi;
}

Graph &Graph::generateGiantComponent(bool logging) {
  generate();
  findGiantComponent(logging);
  computeDegreeSequence();
  computeEdgeList();

  return *this;
}

}  // namespace apm
