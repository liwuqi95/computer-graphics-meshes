#include "per_vertex_normals.h"
#include "triangle_area_normal.h"
#include "vertex_triangle_adjacency.h"

#include <vector>

void per_vertex_normals(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &N) {

    ////////////////////////////////////////////////////////////////////////////
    N.resize(V.rows(), V.cols());

    std::vector<std::vector<int>> surfaces;
    vertex_triangle_adjacency(F, V.rows(), surfaces);

    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < surfaces[i].size(); j++) {
            int face = surfaces[i][j];
            N.row(i) = N.row(i) + triangle_area_normal(V.row(F(face, 0)), V.row(F(face, 1)), V.row(F(face, 2)));
        }

        N.row(i) = N.row(i).normalized();
    }
    ////////////////////////////////////////////////////////////////////////////
}
