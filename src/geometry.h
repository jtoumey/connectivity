#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

struct Parameters
{
    double x, y;
    int nx, ny, np;
    double dx, dy;
    int num_verts;

    Parameters(double x_, double y_, int nx_, int ny_) : 
        x(x_), y(y_), nx(nx_), ny(ny_) {}
};

struct Node 
{
    double x, y;
    int neighbors[4] = {-1, -1, -1, -1};
};

struct Edge
{
    int start_vertex;
    int end_vertex;
};

struct UniqueEdgeStatus
{
    int unique_neighbor;
    int neighbor_count;
};

class Geometry
{
private:
    Parameters geometry_params;

    std::vector<Node> vertices;
    std::vector<Node> cell_centers;
    std::vector<Edge> edge_list;

public:
    Geometry(Parameters inputs_) : geometry_params(inputs_) {};

    void calculateVertices();
    void calculateCellCenters();
    void calculateConnectivity();
    void generateEdgeList();
    UniqueEdgeStatus findUniqueFaceNeighbors(int, int, int);
    void determineCellEdgeAssociation();

    void writeEdgeList();
};

#endif // GEOMETRY_H
