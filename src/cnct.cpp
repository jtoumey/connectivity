#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

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
    int e_start;
    int e_end;
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

int Geometry::findFaceNeighbors(int current_cell_index, int edge_start_vert, int edge_end_vert)
{
    int edge_validity;
    int neighbor_count = 0;
    // Loop over the four (at maximum) cells connected to the edge start vertex
    for (int jj = 0; jj < 4; ++jj)
    {
        int neighbor_cell_index = vertices[edge_start_vert].neighbors[jj];

        // Some vertices don't have neighbors in particular quadrants
        if (neighbor_cell_index == -1)
        {
            continue;
        }
        else
        {
            if (neighbor_cell_index != current_cell_index)
            {
                neighbor_count++;
                
            }
            // // Loop over vertices of this shared cell 
            // for (int nn = 0; nn < 4; ++nn)
            // {

            //     if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && neighbor_cell_index > current_cell_index)
            //     {
            //         std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is not valid." << std::endl;

            //         edge_validity = 0;
            //     }
            //     else
            //     {
            //         edge_validity = 1;
            //     }

                //     neighbor_count++;

                //     if (current_cell_index < neighbor_cell_index)
                //     {
                //         std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is unique." << std::endl;
                //     }
                //     else if (current_cell_index > neighbor_cell_index)
                //     {
                //         std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") was already stored." << std::endl;
                //     }
                //     else if (current_cell_index == neighbor_cell_index)
                //     {
                //         std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is unique and on a boundary." << std::endl;
                //     }
                // }
                // else if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && neighbor_cell_index == current_cell_index )
                // {
                //     if (current_cell_index )
                //     {
                //         std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") is unique and on a boundary." << std::endl;
                //     }
                // }             
            // }
        }
    }
    return(edge_validity);
    // std::cout << "Edge with vertices (" << edge_start_vert << ", " << edge_end_vert << ") has " << neighbor_count << " neighbors." << std::endl;
}

void Geometry::generateEdgeList()
{
    // Loop over all cells
    for (int ii = 0; ii < geometry_params.np; ++ii)
    {
        // For a South face of the cell, collect the SW and SE vertices as the start and end
        int sw_vertex_index = cell_centers[ii].neighbors[0];
        int se_vertex_index = cell_centers[ii].neighbors[1];

        //create the edge here
        // return a boolean for if edge is unique
        int s_edge_validity = findFaceNeighbors(ii, sw_vertex_index, se_vertex_index);

        if (s_edge_validity == 1)
        {
            std::cout << "Edge with vertices (" << sw_vertex_index << ", " << se_vertex_index << ") is unique." << std::endl;
            // if yes, push onto edge list
        }


        int ne_vertex_index = cell_centers[ii].neighbors[2];

        findFaceNeighbors(ii, se_vertex_index, ne_vertex_index);

        int nw_vertex_index = cell_centers[ii].neighbors[3];

        findFaceNeighbors(ii, ne_vertex_index, nw_vertex_index);

        findFaceNeighbors(ii, nw_vertex_index, sw_vertex_index);

    }
}

    //     int ne_vertex_index = cell_centers[ii].neighbors[2];

    //     // Loop over cells that use the sw_vertex_index 
    //     for (int jj = 0; jj < 4; ++jj)
    //     {
    //         int neighbor_cell_index = vertices[se_vertex_index].neighbors[jj];

    //         // some vertices don't have neighbors in particular quadrants
    //         if (neighbor_cell_index == -1)
    //         {
    //             continue;
    //         }
    //         else
    //         {
    //             // Loop over associated vertices of cells that contain sw_vertex_index
    //             for (int nn = 0; nn < 4; ++nn)
    //             {
    //                 if (cell_centers[neighbor_cell_index].neighbors[nn] == ne_vertex_index && neighbor_cell_index != ii)
    //                 {
    //                     std::cout << "Found a n-s edge with a neighbor: verts ( " << se_vertex_index << ", " << ne_vertex_index << ")" << std::endl;
    //                     std::cout << "  Neighbor:" << ii << "(parent cell), with " << neighbor_cell_index << std::endl;

    //                     if (ii < neighbor_cell_index)
    //                     {
    //                         std::cout << "Found a unique edge:" << std::endl;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     continue;
    //                 }
    //             }
    //         }
    //     }

    //     int nw_vertex_index = cell_centers[ii].neighbors[3];
    //     // Loop over cells that use the sw_vertex_index 
    //     for (int jj = 0; jj < 4; ++jj)
    //     {
    //         int neighbor_cell_index = vertices[ne_vertex_index].neighbors[jj];

    //         // some vertices don't have neighbors in particular quadrants
    //         if (neighbor_cell_index == -1)
    //         {
    //             continue;
    //         }
    //         else
    //         {
    //             // Loop over associated vertices of cells that contain sw_vertex_index
    //             for (int nn = 0; nn < 4; ++nn)
    //             {
    //                 if (cell_centers[neighbor_cell_index].neighbors[nn] == nw_vertex_index && neighbor_cell_index != ii)
    //                 {
    //                     std::cout << "Found a n-s edge with a neighbor: verts ( " << ne_vertex_index << ", " << nw_vertex_index << ")" << std::endl;
    //                     std::cout << "  Neighbor:" << ii << "(parent cell), with " << neighbor_cell_index << std::endl;

    //                     if (ii < neighbor_cell_index)
    //                     {
    //                         std::cout << "Found a unique edge:" << std::endl;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     continue;
    //                 }
    //             }
    //         }
    //     }


    //     // Loop over cells that use the sw_vertex_index 
    //     for (int jj = 0; jj < 4; ++jj)
    //     {
    //         int neighbor_cell_index = vertices[nw_vertex_index].neighbors[jj];

    //         // some vertices don't have neighbors in particular quadrants
    //         if (neighbor_cell_index == -1)
    //         {
    //             continue;
    //         }
    //         else
    //         {
    //             // Loop over associated vertices of cells that contain sw_vertex_index
    //             for (int nn = 0; nn < 4; ++nn)
    //             {
    //                 if (cell_centers[neighbor_cell_index].neighbors[nn] == sw_vertex_index && neighbor_cell_index != ii)
    //                 {
    //                     std::cout << "Found a n-s edge with a neighbor: verts ( " << nw_vertex_index << ", " << sw_vertex_index << ")" << std::endl;
    //                     std::cout << "  Neighbor:" << ii << "(parent cell), with " << neighbor_cell_index << std::endl;

    //                     if (ii < neighbor_cell_index)
    //                     {
    //                         std::cout << "Found a unique edge:" << std::endl;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     continue;
    //                 }
    //             }
    //         }
    //     }

    // }
