# On Directed Densest Subgraph Detection

This project aims to detect directed dense subgraphs in directed graphs.

dds is the executable, and is compiled on Ubuntu 18.04.5, with -O3 optimization.

TwitterList.txt is an example directed graph from [KONECT](http://konect.cc/networks/)

## Running Format

./dds [1]input_graph [2]H/A/M/G

**Running example for harmonic meand-based denese subgraph detection**

./dds ./TwitterList.txt H

Similarly, you can run instances for arithmetic, minimum, and geometric mean-based dense subgraph detection by replacing 'H' with 'A', 'M', or 'G' respectively.

## Graph Format

number of vertices \t number of edges \n

v0 \t v1 \n

v0 \t v2 \n

...
