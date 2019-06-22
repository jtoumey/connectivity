# Connectivity

Toy code to generate an edge list based on a discretized bounding box. 

## Method

1. Compute vertices and cell centers
3. Determine the four vertices for each quad cell
4. Determine the four (maximum) cells sharing each vertex
5. Build the edge list
  * For each cell, build the four edges and determine if they are unique
