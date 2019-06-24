#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>

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

struct UniqueEdge
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
    int findFaceNeighbors(int, int, int);
    int getFaceNeighbors(int, int, int);
    UniqueEdge findUniqueFaceNeighbors(int, int, int);

    void writeEdgeList();
};

int main (void)
{
    Parameters grid_inputs(0.3, 0.4, 3, 4);

    grid_inputs.np = grid_inputs.nx*grid_inputs.ny;
    grid_inputs.num_verts = (grid_inputs.nx + 1)*(grid_inputs.ny + 1);
    grid_inputs.dx = grid_inputs.x/(float)(grid_inputs.nx);
    grid_inputs.dy = grid_inputs.y/(float)(grid_inputs.ny);

    Geometry geom(grid_inputs);

    geom.calculateVertices();
    geom.calculateCellCenters();
    geom.calculateConnectivity();

    geom.generateEdgeList();
    geom.writeEdgeList();

    return(0);
}

void Geometry::calculateVertices()
{
    double xc = 0.0; 
    double yc = 0.0;
    
    for (int ii = 0; ii < (geometry_params.nx + 1); ++ii)
    {
        for (int jj = 0; jj < (geometry_params.ny + 1); ++jj)
        {
            Node vc;
            vc.x = xc; 
            vc.y = yc;

            vertices.push_back(vc);

            yc += geometry_params.dy;
        }
        xc += geometry_params.dx;
        yc = 0.0;
    }
}

void Geometry::calculateCellCenters()
{
    double xc = 0.0 + geometry_params.dx/2.0;
    double yc = 0.0 + geometry_params.dy/2.0;

    for (int ii = 0; ii < geometry_params.nx; ++ii)
    {
        for (int jj = 0; jj < geometry_params.ny; ++jj)
        {
            Node cel;

            cel.x = xc;
            cel.y = yc;

            cell_centers.push_back(cel);

            yc += geometry_params.dy;
        }
        xc += geometry_params.dx;
        yc = 0.0 + geometry_params.dy/2.0;
    }
    for (int ii = 0; ii < geometry_params.np; ++ii)
    {
        std::cout << ii << ": " << cell_centers[ii].x << ", " << cell_centers[ii].y << std::endl;
    }
    std::cout << std::endl;
}

void Geometry::calculateConnectivity()
{
    double eps = 1e-8;

    for (int ii = 0; ii < geometry_params.np; ++ii)
    {
        Node cc = cell_centers[ii];

        for (int jj = 0; jj < geometry_params.num_verts; ++jj)
        {
            Node vc = vertices[jj];

            if (std::fabs((cc.x - geometry_params.dx/2.0) - vc.x) < eps)
            {
                if (std::fabs((cc.y - geometry_params.dy/2.0) - vc.y) < eps)
                {
                    std::cout << "Cell " << ii << " match with SW vert " << jj << " in x" << std::endl;
                    cell_centers[ii].neighbors[0] = jj;
                    vertices[jj].neighbors[1] = ii;
                }
                else if (std::fabs((cc.y + geometry_params.dy/2.0) - vc.y) < eps)
                {
                    std::cout << "Cell " << ii << " match with NW vert " << jj << " in x" << std::endl;
                    cell_centers[ii].neighbors[3] = jj;
                    vertices[jj].neighbors[0] = ii;
                }
            }
            else if (std::fabs((cc.x + geometry_params.dx/2.0) - vc.x) < eps)
            {
                if (std::fabs((cc.y - geometry_params.dy/2.0) - vc.y) < eps)
                {
                    std::cout << "Cell " << ii << " match with SE vert " << jj << " in x" << std::endl;
                    cell_centers[ii].neighbors[1] = jj;
                    vertices[jj].neighbors[2] = ii;
                }
                else if (std::fabs((cc.y + geometry_params.dy/2.0) - vc.y) < eps)
                {
                    std::cout << "Cell " << ii << " match with NE vert " << jj << " in x" << std::endl;
                    cell_centers[ii].neighbors[2] = jj;
                    vertices[jj].neighbors[3] = ii;
                }
            }
        }
        std::cout << std::endl;
    }
}

UniqueEdge Geometry::findUniqueFaceNeighbors(int current_cell_index, int edge_start_vert, int edge_end_vert)
{
    UniqueEdge current_edge;

    int neighbor_count = 0;
    int unique_neighbor = 0;
    // Loop over the four (at maximum) cells connected to the edge start vertex
    for (int jj = 0; jj < 4; ++jj)
    {
        int neighbor_cell_index = vertices[edge_start_vert].neighbors[jj];

        // Some vertices don't have neighbors in particular quadrants
        if (neighbor_cell_index == -1)
        {
            //continue;
        }
        else
        {
            // Loop over vertices of this shared cell 
            for (int nn = 0; nn < 4; ++nn)
            {
                // If this shared cell also possesses the end end vertex and is NOT the cell of interest itself, then the face
                // has a neighbor
                if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && neighbor_cell_index != current_cell_index)
                {
                    neighbor_count++;

                    if (current_cell_index < neighbor_cell_index)
                    {
                        std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is unique." << std::endl; 
                        unique_neighbor = 1;
                    }
                }
                else if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && current_cell_index < neighbor_cell_index)
                {
                    //continue;
                    unique_neighbor = 1;
                }
            }
        }
    }

    current_edge.neighbor_count = neighbor_count;
    current_edge.unique_neighbor = unique_neighbor;

    return(current_edge);
}

int Geometry::findFaceNeighbors(int current_cell_index, int edge_start_vert, int edge_end_vert)
{
    int neighbor_count = 0;
    // Loop over the four (at maximum) cells connected to the edge start vertex
    for (int jj = 0; jj < 4; ++jj)
    {
        int neighbor_cell_index = vertices[edge_start_vert].neighbors[jj];

        // Some vertices don't have neighbors in particular quadrants
        if (neighbor_cell_index == -1)
        {
            //continue;
        }
        else
        {
            // Loop over vertices of this shared cell 
            for (int nn = 0; nn < 4; ++nn)
            {
                // If this shared cell also possesses the end end vertex and is NOT the cell of interest itself, then the face
                // has a neighbor
                if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && neighbor_cell_index != current_cell_index)
                {
                    neighbor_count++;

                    if (current_cell_index < neighbor_cell_index)
                    {
                        std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is unique." << std::endl; 
                    }
                }
                else
                {
                    //continue;
                }
            }
        }
    }
    return(neighbor_count);
}

int Geometry::getFaceNeighbors(int current_cell_index, int edge_start_vert, int edge_end_vert)
{
    int unique_neighbor = 0;

    // Loop over the four (at maximum) cells connected to the edge start vertex
    for (int jj = 0; jj < 4; ++jj)
    {
        int neighbor_cell_index = vertices[edge_start_vert].neighbors[jj];

        // Some vertices don't have neighbors in particular quadrants
        if (neighbor_cell_index == -1)
        {
            // continue;
        }
        else
        {
            // Loop over vertices of this shared cell 
            for (int nn = 0; nn < 4; ++nn)
            {
                // If this shared cell also possesses the end end vertex and is NOT the cell of interest itself, then the face
                // has a neighbor
                if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && current_cell_index < neighbor_cell_index)
                {
                    std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is unique." << std::endl; 

                    unique_neighbor = 1;
                }
                else
                {
                    //continue;
                }
            }
        }
    }
    return(unique_neighbor);
}

void Geometry::generateEdgeList()
{
    // Loop over all cells
    for (int ii = 0; ii < geometry_params.np; ++ii)
    {
        int neighbor_count = 0;
        // For a South face of the cell, collect the SW and SE vertices as the start and end
        int sw_vertex_index = cell_centers[ii].neighbors[0];
        int se_vertex_index = cell_centers[ii].neighbors[1];

        //create the edge here
        Edge edge_c;
        edge_c.start_vertex = sw_vertex_index;
        edge_c.end_vertex = se_vertex_index;

        // return a boolean for if edge is unique
        UniqueEdge current_edge = findUniqueFaceNeighbors(ii, sw_vertex_index, se_vertex_index);

        if (current_edge.neighbor_count == 0 || current_edge.unique_neighbor == 1)
        {
            edge_list.push_back(edge_c);
        }
        
        // For an East face of the cell, collect the SE and NE vertices as the start and end
        int ne_vertex_index = cell_centers[ii].neighbors[2];

        //create the edge here
        edge_c.start_vertex = se_vertex_index;
        edge_c.end_vertex = ne_vertex_index;

        // return a boolean for if edge is unique
        current_edge = findUniqueFaceNeighbors(ii, se_vertex_index, ne_vertex_index);

        if (current_edge.neighbor_count == 0 || current_edge.unique_neighbor == 1)
        {
            edge_list.push_back(edge_c);
        }

        // For a North face of the cell, collect the NE and NW vertices as the start and end
        int nw_vertex_index = cell_centers[ii].neighbors[3];

        //create the edge here
        edge_c.start_vertex = ne_vertex_index;
        edge_c.end_vertex = nw_vertex_index;

        // return a boolean for if edge is unique
        neighbor_count = findFaceNeighbors(ii, ne_vertex_index, nw_vertex_index);

        if (neighbor_count == 0)
        {
            edge_list.push_back(edge_c);
        }
        else 
        {
            int unique_neighbor_flag = getFaceNeighbors(ii, ne_vertex_index, nw_vertex_index);

            if (unique_neighbor_flag == 1)
            {
                edge_list.push_back(edge_c);
            }
        }

        // For a West face of the cell, collect the NW and SW vertices as the start and end

        //create the edge here
        edge_c.start_vertex = nw_vertex_index;
        edge_c.end_vertex = sw_vertex_index;

        // return a boolean for if edge is unique
        neighbor_count = findFaceNeighbors(ii, nw_vertex_index, sw_vertex_index);

        if (neighbor_count == 0)
        {
            edge_list.push_back(edge_c);
        }
        else 
        {
            int unique_neighbor_flag = getFaceNeighbors(ii, nw_vertex_index, sw_vertex_index);

            if (unique_neighbor_flag == 1)
            {
                edge_list.push_back(edge_c);
            }
        }
    }
}

void Geometry::writeEdgeList()
{
    int num_edges = edge_list.size();
    // Loop over all cells
    for (int ii = 0; ii < num_edges; ++ii)
    {
        std::cout << "Edge #: " << ii << "; Vertices: (" << edge_list[ii].start_vertex << ", " << edge_list[ii].end_vertex << ")" << std::endl;
    }

    std::ofstream edge_list_output;
    edge_list_output.open("edge_list.dat");

    for (size_t ii = 0; ii < num_edges; ++ii)
    {
        edge_list_output << "Edge #: " << ii << "; Vertices: (" << edge_list[ii].start_vertex << ", " << edge_list[ii].end_vertex << ")" << std::endl;
    }
    edge_list_output.close();
} 
