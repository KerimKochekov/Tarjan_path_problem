# Fast algorithm for solving path problems
Implementation of [Tarjan's 1981 paper](https://dl.acm.org/doi/pdf/10.1145/322261.322273) for "Computing Path Expressions for Reducible Flow Graphs"

Abstract: Let G=(V, E) be a directed graph with a distinguished source vertex s. The single-source path expression problem is to find, for each vertex v, a regular expression P(s, v) which represents the set of all paths in G from s to v.

Definitions: 
- A flow graph G = (V, E, r) is a directed graph with a distinguished start vertex r such that every vertex in G Is reachable from r.
- A reducible flow graph G = (V, E, r) is a flow graph that can be reduced to the graph consisting of the single vertex r and no edges by means of the following transformations: 
  - T1 (Remove a loop): If e is an edge such that h(e) = t(e), delete edge e. 
  - T2 (Remove a vertex): If w # r is a vertex such that all edges e with t(e) = w have h(e) = v for some vertex v, contract w into v by deleting w and all edges entering w, and converting any edge e with h(e) = w into an edge e' with h(e') = v and t(e') = t(e).
![Example graph](https://github.com/KerimKochekov/Tarjan_path_problems/blob/main/bin/example_graph.png)

Notes from paper:
- For dense graphs the time bound is O(n^3 + m) and the space bound is O(n^2). (Note that m, the number of edges, is bounded by n 2 unless the graph contains multiple edges.)
- Flow graph is reducible if every cycle has a single entry from the start vertex.
