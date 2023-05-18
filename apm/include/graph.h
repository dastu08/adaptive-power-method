#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "general.h"
#include "rng.h"

namespace apm {

/* Construct and manage an ER graph

*/
class Graph {
 private:
  Rng &rng;                // reference to the rng
  int n;                   // number of nodes
  double alpha;            // mean degree of the node
  double p;                // probability for an edge
  vvInt neighbors;         // neighbor list
  sInt giantComponentSet;  // index set of the giant component
  vvInt giantComponent;    // neighbor list of the giant component
  vInt degreeSequence;     // list of the node degrees
  vvInt edgeList;          // list of the edges of the giant component;

  /* Extract the subgraph given by the nodes in the component.

  **Parameter**
      - graph: neighbor list of the whole graph
      - component: set of nodes that are to extracted

  **Return**
      Neighbor list of the subgraph.

  **Description**
      Extract the neighbor lists of the nodes that make up the component.
      Relabel the indices to make them be consistent with the new subgraph.
  */
  static vvInt extractComponent(vvInt &graph, sInt component);

  /* Recursively add the neighbors of the node to the component.

  **Parameter**
      - graph: neighbor list of the graph
      - node: current node to consider
      - component: set of nodes that already belong the component
      - visited: set of nodes that were visited by previous iterations

  **Description**
      Check if the node was already visited before. If not add the node to
      the component and call the function for the neighbors of the node again.
      The set `component` will be updated over the recursion. After the first
      call of the function finishes all the nodes connected to the starting
      node are in `component`.
  */
  static void addToComponent(vvInt &graph,
                             int node,
                             sInt &component,
                             sInt &visited);

  /* Recursively add the neighbors of the node to the chain list.

  **Parameter**
      - graph: neighbor list of the graph
      - node: current node to consider
      - component: list of nodes that already belong the chain list
      - visited: set of nodes that were visited by previous iterations
      - gateway: keep the current guess for the gateway

  **Description**
    Check if the node was already visited before. If not add the node to the
    component. Check if the neighbors of the node have degree less than 3,
    meaning that they are in a chain. Then call the function for those neighbors
    of the node again if not visited before. The set `component` will be updated
    over the recursion. After the first call of the function finishes all the
    nodes in the chain connected to the initial node are in `component`.
    The last value of gateway is the valid gateway.
  */
  static void addToChainList(vvInt &graph,
                             int node,
                             vInt &component,
                             sInt &visited,
                             int &gateway);

 public:
  /* Specify the graph size and mean degree.

  **Parameter**
      - n: number of nodes
      - alpha: mean degree of the nodes
      - rng: random number generator to use

  **Description**
      Initialize the class properties. Compute the probability p=alpha/n.
  */
  Graph(int n, double alpha, Rng &rng);

  vvInt getEdgeList() {
    return this->edgeList;
  }

  vvInt getNeighbors() {
    return this->neighbors;
  }

  vvInt getGiantComponent() {
    return this->giantComponent;
  }

  sInt getGiantComponentSet() {
    return this->giantComponentSet;
  }

  vInt getDegreeSequence() {
    return this->degreeSequence;
  }

  /* Generate the random graph.

  **Return**
      Reference to `this`.

  **Description**
      Connect two nodes with this probability `p` and save the neighbor list.
  */
  Graph &generate();

  /* Construct the edge list from the neighbor list of the giant component

  **Return**
      Reference to `this`.

  **Prerequisites**
      Need to have called `findGiantComponent` to have found the giant
      component.

  **Description**
      Loop through the neighbor list and add all pairs to the edge list.
      Note that each edge will appear twice but with the nodes reversed.

  */
  Graph &computeEdgeList();

  /* Print the neighbor list of the graph.

  **Return**
      Reference to `this`.

  **Description**
      Print in each line the node index and its neighbors.
  */
  Graph &print();

  /* Print the neighbor list of the giant component of the graph.

  **Prerequisites**
    Need to have called `findGiantComponent` to have found the giant
    component.

  **Return**
    Reference to `this`.

  **Description**
    Print in each line the node index and its neighbors.
  */
  Graph &printGiantComponent();

  /* Find the giant connected component of the graph.

  **Parameter**
    - logging: flag to disable info logging

  **Return**
      Reference to `this`.

  **Description**
      Find all the connected components of the graph. Find of these components
      the one with the most elements. Save the set of nodes contained in the
      giant component. Save the subgraph of the giant component.
  */
  Graph &findGiantComponent(bool logging = true);

  /* Print the set of nodes that are the giant component.

  **Prerequisites**
    Need to have called `findGiantComponent` to have found the giant
    component.

  **Return**
    Reference to `this`.

  **Description**
    Print in one line the nodes that are in the giant component. The
    numbers refer to the original graph.
  */
  Graph &printGiantComponentSet();

  /* Compute the adjacency matrix of the giant component.

  **Prerequisites**
      Need to have called `findGiantComponent` to have found the giant
      component.

  **Return**
      Sparse matrix were each row is a sparse vector with elements where the
      node is connected to its neighbors.

  **Description**
      Convert the neighbor lists into a vector of sparse vectors which is a
      sparse matrix.
  */
  matrixSparse adjacencyMatrix();

  /* Compute the degree sequence of the giant component.

  **Prerequisites**
      Need to have called `findGiantComponent` to have found the giant
      component.

  **Description**
      Count the number of neighbors of each node. This is the degree of the
      node.
  */
  Graph &computeDegreeSequence();

  /* Print the degree sequence of the giant component.

 **Prerequisites**
      Need to have called `computeDegreeeSequence`.

  **Description**
      Print in each row the node index and its degree.
  */
  Graph &printDegreeSeqence();

  /* Find the number of dangling chains in the giant component.

  **Prerequisites**
    Need to have called `computeDegreeeSequence`.

  **Parameter**
    - logging: flag to disable info logging

  **Return**
    List of node labels that are the endpoints of dangling chains.

  **Description**
    Go through the degree sequence and find all endpoints that are nodes of
    degree 1. Of the endpoints select the one neighbor and get its number of
    neibhbors. If this is 2 then we have found a dangling chain.
  */
  vInt findDanglingChains(bool logging = false);

  /* Find the nodes that are in the dangling chains and their gateways.

  **Parameter**
    - danglingChains: list of dangling chain endpoints, e.g. from
      `findDanglingChains()`
    - dcList: list that will contain all the nodes in the chains
    - gwList: list of all the gateways

  **Description**
    Go through the endpoints and recursively add the neighbor which have a
    degree of less than 3 to the dcList. Add the last node with more than 2
    neighbors to the gwList.
  */
  void getDanglingChainList(vInt &danglingChains,
                            vInt &dcList,
                            vInt &gwList);

  /* Find the nodes that are in the dangling chains and their gateways.
  
  Same as `getDanglingChainList` just with type `set` instead of `vector`.
  */
  void getDanglingChainSet(vInt &danglingChains,
                           sInt &dcSet,
                           sInt &gwSet);

  /* Compute the transition matrix of the giant component.

  **Prerequisites**
      Need to have called `findGiantComponent` and `computeDegreeSequence`
      to have found the giant component and to have the node degrees.

  **Return**
      Sparse matrix were each row is a sparse vector with transiton
      probabilities to its connected neighbors.

  **Description**
      Convert the neighbor lists into a vector of sparse vectors which is a
      sparse matrix. Then divide by the node degree to get a stochastic
      matrix.
  */
  matrixSparse transitionMatrix();

  /* Generate the random graph and find its giant component

  **Parameter**
    - logging: flag that is passed to the subroutines

  **Return**
      Reference to `this`.

  **Description**
      Call `generate`, `findConnectedComponent`, `computeDegreeSequence` and
      `computeEdgeList`.
  */
  Graph &generateGiantComponent(bool logging = true);
};

}  // namespace apm

#endif  // _GRAPH_H_