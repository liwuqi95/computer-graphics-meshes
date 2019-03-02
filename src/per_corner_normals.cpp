#include "per_corner_normals.h"
#include "triangle_area_normal.h"
// Hint:
#include "vertex_triangle_adjacency.h"
#include <iostream>

void per_corner_normals(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const double corner_threshold,
        Eigen::MatrixXd &N) {
    ////////////////////////////////////////////////////////////////////////////
    N = Eigen::MatrixXd::Zero(F.rows() * 3, 3);

    std::vector<std::vector<int>> VF;
    vertex_triangle_adjacency(F, V.rows(), VF);


    for (int i = 0; i < F.rows(); i++) {

        Eigen::RowVector3d face_normal = triangle_area_normal(V.row(F(i, 0)), V.row(F(i, 1)),
                                                              V.row(F(i, 2))).normalized();
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < VF[F(i, j)].size(); k++) {
                int face = VF[F(i, j)][k];
                Eigen::RowVector3d n = triangle_area_normal(V.row(F(face, 0)), V.row(F(face, 1)), V.row(F(face, 2)));

                if (face_normal.dot(n.normalized()) > cos(corner_threshold * M_PI / 180.0)) {
                    N.row(i * 3 + j) += n;
                }
            }
            N.row(3 * i + j) = N.row(3 * i + j).normalized();
        }
    }
    ////////////////////////////////////////////////////////////////////////////
}
