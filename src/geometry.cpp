#include "geometry.h"

#include <cmath>
#include <fstream>

#include <iostream> // TODO: print to screen, so remove when done

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
                    cell_centers[ii].neighbors[0] = jj;
                    vertices[jj].neighbors[1] = ii;
                }
                else if (std::fabs((cc.y + geometry_params.dy/2.0) - vc.y) < eps)
                {
                    cell_centers[ii].neighbors[3] = jj;
                    vertices[jj].neighbors[0] = ii;
                }
            }
            else if (std::fabs((cc.x + geometry_params.dx/2.0) - vc.x) < eps)
            {
                if (std::fabs((cc.y - geometry_params.dy/2.0) - vc.y) < eps)
                {
                    cell_centers[ii].neighbors[1] = jj;
                    vertices[jj].neighbors[2] = ii;
                }
                else if (std::fabs((cc.y + geometry_params.dy/2.0) - vc.y) < eps)
                {
                    cell_centers[ii].neighbors[2] = jj;
                    vertices[jj].neighbors[3] = ii;
                }
            }
        }
    }
}

UniqueEdgeStatus Geometry::findUniqueFaceNeighbors(int current_cell_index, int edge_start_vert, int edge_end_vert)
{
    UniqueEdgeStatus current_edge_status;

    int neighbor_count = 0;
    int unique_neighbor = 0;

    // Loop over the four (at maximum) cells connected to the edge start vertex
    for (int jj = 0; jj < 4; ++jj)
    {
        int neighbor_cell_index = vertices[edge_start_vert].neighbors[jj];

        
        if (neighbor_cell_index == -1)
        {
            // Some vertices don't have neighbors in particular quadrants, so we skip these
        }
        else
        {
            // Loop over vertices of this shared cell 
            for (int nn = 0; nn < 4; ++nn)
            {
                // If this shared cell also possesses the edge end vertex and is NOT the cell of interest itself, then the current face
                // has a neighbor
                if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && neighbor_cell_index != current_cell_index)
                {
                    // Counts the neighbors that are not the cell of interest itself
                    neighbor_count++;

                    if (current_cell_index < neighbor_cell_index)
                    {
                        // Further, if our current cell index is less than this neighbor index, we know we found a unique edge
                        // and thus set the flag s.t. it will get pushed onto the list
                        unique_neighbor = 1;
                    }
                }
                else if (cell_centers[neighbor_cell_index].neighbors[nn] == edge_end_vert && current_cell_index < neighbor_cell_index)
                {
                    // If the edge has no neighbors beyond its parent cell, it is on a boundary and hence has not been counted yet
                    unique_neighbor = 1;
                }
            }
        }
    }

    current_edge_status.neighbor_count = neighbor_count;
    current_edge_status.unique_neighbor = unique_neighbor;

    return(current_edge_status);
}

void Geometry::generateEdgeList()
{
    // Loop over all cells
    for (int ii = 0; ii < geometry_params.np; ++ii)
    {
        for (int jj = 0; jj < 4; ++jj)
        {
            int e_start_vertex = cell_centers[ii].neighbors[jj];
            int e_end_vertex = cell_centers[ii].neighbors[(jj+1)%4];

            Edge edge_c;
            edge_c.start_vertex = e_start_vertex;
            edge_c.end_vertex = e_end_vertex;

            UniqueEdgeStatus edge_status = findUniqueFaceNeighbors(ii, e_start_vertex, e_end_vertex);

            if (edge_status.neighbor_count == 0 || edge_status.unique_neighbor == 1)
            {
                edge_list.push_back(edge_c);
            }
        }
    }
}

void Geometry::writeEdgeList()
{
    int num_edges = edge_list.size();

    std::ofstream edge_list_output;
    edge_list_output.open("edge_list.dat");

    for (size_t ii = 0; ii < num_edges; ++ii)
    {
        edge_list_output << "Edge #: " << ii << "; Vertices: (" << edge_list[ii].start_vertex << ", " << edge_list[ii].end_vertex << ")" << std::endl;
    }
    edge_list_output.close();
} 

void Geometry::determineCellEdgeAssociation()
{
    for (size_t cell_id = 0; cell_id < geometry_params.np; ++cell_id)
    {

        for (size_t face_dir = 0; face_dir < 4; ++face_dir)
        {
            Edge current_edge;
            current_edge.start_vertex = cell_centers[cell_id].neighbors[face_dir];
            current_edge.end_vertex = cell_centers[cell_id].neighbors[(face_dir+1)%4];
            size_t num_edges = edge_list.size();

            for (size_t edge_id = 0; edge_id < num_edges; ++edge_id)
            {
                // Find an edge with the same start and end vertices (test both orderings)
                if ((current_edge.start_vertex == edge_list[edge_id].start_vertex && current_edge.end_vertex == edge_list[edge_id].end_vertex) || (current_edge.start_vertex == edge_list[edge_id].end_vertex && current_edge.end_vertex == edge_list[edge_id].start_vertex))
                {   
                    cell_centers[cell_id].edge_neighbors[face_dir] = edge_id;
                }
            }

        }
    }

    size_t num_edges = edge_list.size();
    for (size_t edge_id = 0; edge_id < num_edges; ++edge_id)
    {
        for (size_t cell_id = 0; cell_id < geometry_params.np; ++cell_id)
        {
            for (size_t face_dir = 0; face_dir < 4; ++face_dir)
            {
                if (cell_centers[cell_id].edge_neighbors[face_dir] == edge_id)
                {
                    if (edge_list[edge_id].cell_neighbors[0] == -1)
                    {
                        edge_list[edge_id].cell_neighbors[0] = cell_id;
                    }
                    else
                    {
                        edge_list[edge_id].cell_neighbors[1] = cell_id;
                    }
                }
            }
        }
    }

}
