#include "cube.h"

void cube(
        Eigen::MatrixXd &V,
        Eigen::MatrixXi &F,
        Eigen::MatrixXd &UV,
        Eigen::MatrixXi &UF,
        Eigen::MatrixXd &NV,
        Eigen::MatrixXi &NF) {

    V.resize(8, 3);
    F.resize(6, 4);
    UV.resize(14, 2);
    UF.resize(6, 4);
    NV.resize(6, 3);
    NF.resize(6, 4);

    V << 0, 0, 1,
            0, 0, 0,
            1, 0, 1,
            1, 0, 0,
            0, 1, 1,
            0, 1, 0,
            1, 1, 1,
            1, 1, 0;

    F << 0, 1, 3, 2,
            1, 3, 7, 5,
            2, 3, 7, 6,
            0, 2, 6, 4,
            4, 5, 7, 6,
            0, 1, 5, 4;

    UV << 0, 0.25,
            0, 0.5,
            0.25, 0,
            0.25, 0.25,
            0.25, 0.5,
            0.25, 0.75,
            0.5, 0,
            0.5, 0.25,
            0.5, 0.5,
            0.5, 0.75,
            0.75, 0.25,
            0.75, 0.5,
            1, 0.25,
            1, 0.5;

    UF << 0, 1, 4, 3,
            2, 3, 7, 6,
            3, 4, 8, 7,
            4, 5, 9, 8,
            7, 8, 11, 10,
            10, 11, 13, 12;

    NV << 0, -1, 0,
            0, 0, -1,
            1, 0, 0,
            0, 0, 1,
            0, 1, 0,
            -1, 0, 0;

    NF << 0, 0, 0, 0,
            1, 1, 1, 1,
            2, 2, 2, 2,
            3, 3, 3, 3,
            4, 4, 4, 4,
            5, 5, 5, 5;


}
