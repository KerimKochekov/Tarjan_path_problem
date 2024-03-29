# Fast algorithm for solving path problems
Implementation of [Tarjan's 1981 paper](https://dl.acm.org/doi/pdf/10.1145/322261.322273) for *"Computing Path Expressions for Reducible Flow Graphs"*.

## **Abstract**

Let $G=(V, E)$ be a directed graph with a distinguished source vertex $s$. The single-source path expression problem is to find, for each vertex $v$, a regular expression $P(s, v)$ which represents the set of all paths in $G$ from $s$ to $v$.

## **How to run?**

```bash
$ g++ algorithm.cpp -std=c++17 -o algorithm
$ ./algorithm
```

### **Input**

The `graph.txt` file contains a directed flow graph $G$ with $n$ vertices and $m$ edges with the source vertex $s$. The first line of the file given three integers of $n, m, s$. Following $m$ lines given directed edges in form $u \ v$, where $0 \le u \neq v < n$.

### **Definitions**
- A flow graph $G = (V, E, r)$ is a directed graph with a distinguished start vertex $r$ such that every vertex in $G$ is reachable from $r$.
- A reducible flow graph $G = (V, E, r)$ is a flow graph that can be reduced to the graph consisting of the single vertex $r$ and no edges by means of the following transformations: 
  - $T_1$ *(Remove a loop)*: If $e$ is an edge such that $h(e) = t(e)$, delete edge $e$. 
  - $T_2$ *(Remove a vertex)*: If $w \ne r$ is a vertex such that all edges $e$ with $t(e) = w$ have $h(e) = v$ for some vertex $v$, contract $w$ into $v$ by deleting $w$ and all edges entering $w$, and converting any edge $e$ with $h(e) = w$ into an edge $e'$ with $h(e') = v$ and $t(e') = t(e)$.

<div align=center>
    <img src="https://github.com/KerimKochekov/Tarjan_path_problems/blob/main/bin/example_graph.png" width="50%" height="auto" alt="Example graph">
</div>

### **Notes from paper**
- Flow graph is reducible if every cycle has a single entry from the start vertex or the flow graph's all dominator strong components are single vertices.

Here is the derived graph of $G$ with the help of dominator decomposition.

<div align=center>
    <img src="https://github.com/KerimKochekov/Tarjan_path_problems/blob/main/bin/derived_graph.png" width="50%" height="auto" alt="Derived graph">
</div>
