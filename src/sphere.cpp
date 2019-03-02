#include "sphere.h"
#include <iostream>

void sphere(
        const int num_faces_u,
        const int num_faces_v,
        Eigen::MatrixXd &V,
        Eigen::MatrixXi &F,
        Eigen::MatrixXd &UV,
        Eigen::MatrixXi &UF,
        Eigen::MatrixXd &NV,
        Eigen::MatrixXi &NF) {
    ////////////////////////////////////////////////////////////////////////////


    int num_face = num_faces_u * num_faces_v;
    int num_vertice = (num_faces_u + 1) * (num_faces_v + 1);

    V.resize(num_vertice, 3);
    F.resize(num_face, 4);
    UV.resize(num_vertice, 2);
    UF.resize(num_face, 4);
    NV.resize(num_vertice, 3);
    NF.resize(num_face, 4);


    double increament_theta = 2.0 * M_PI / num_faces_u;
    double increament_phi = M_PI / num_faces_v;

    double theta, phi;
    int index, next, face_row;


    for (int i = 0; i <= num_faces_v; i++) {
        for (int j = 0; j <= num_faces_u; j++) {

            index = i * (num_faces_u + 1) + j;
            phi = increament_phi * i;
            theta = increament_theta * j;

            V.row(index) = Eigen::RowVector3d(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
            UV.row(index) = Eigen::RowVector2d((double) j / (double) num_faces_u, (double) i / (double) num_faces_v);
            NV.row(index) = V.row(index);

        }
    }


    for (int i = 0; i < num_faces_v; i++) {
        for (int j = 0; j < num_faces_u; j++) {

            face_row = i * num_faces_u + j;
            index = i * (num_faces_u + 1) + j;
            next = index + num_faces_u + 1;

            F.row(face_row) = Eigen::RowVector4i(index, index + 1, next + 1, next);
            UF.row(face_row) = Eigen::RowVector4i(index, index + 1, next + 1, next);
            NF.row(face_row) = Eigen::RowVector4i(index, index + 1, next + 1, next);
        }
    }

    ////////////////////////////////////////////////////////////////////////////
}
