#include "geometry.h"

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
