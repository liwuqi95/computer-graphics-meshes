#include "catmull_clark.h"
#include <unordered_map>
#include <vector>
#include <utility>
#include <functional>
#include "vertex_triangle_adjacency.h"
#include "iostream"

void catmull_clark(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const int num_iters,
        Eigen::MatrixXd &SV,
        Eigen::MatrixXi &SF) {
    ////////////////////////////////////////////////////////////////////////////


    if (num_iters < 1)
        return;

    SV = V;
    SF.resize(4 * F.rows(), 4);

    // Step 1 -> Face Points

    Eigen::MatrixXd FP;
    FP.resize(F.rows(), 3);

    for (int j = 0; j < F.rows(); j++) {
        FP.row(j) = (V.row(F(j, 0)) + V.row(F(j, 1)) + V.row(F(j, 2)) + V.row(F(j, 3))) / 4.0;

        SV.conservativeResize(SV.rows() + 1, SV.cols());
        SV.row(SV.rows() - 1) = FP.row(j);
    }



    // Step 2 -> Edge Points

    // Get all edges with adjance face

    std::unordered_map<std::string, std::vector<int>> Edge;
    for (int j = 0; j < F.rows(); j++) {
        for (int k = 0; k < F.cols(); k++) {

            int v1 = F(j, k);
            int v2 = F(j, (k + 1) % F.cols());
            std::string key = std::to_string(fmin(v1, v2)) + "-" + std::to_string(fmax(v1, v2));

            if (Edge.find(key) == Edge.end())
                Edge.insert(std::make_pair(key, std::vector<int>()));

            Edge.find(key)->second.push_back(j);
        }
    }

    // Insert all edge points
    std::unordered_map<std::string, int> EP;
    for (int j = 0; j < F.rows(); j++) {
        for (int k = 0; k < F.cols(); k++) {
            int v1 = F(j, k);
            int v2 = F(j, (k + 1) % F.cols());
            std::string key = std::to_string(fmin(v1, v2)) + "-" + std::to_string(fmax(v1, v2));

            if (EP.find(key) == EP.end()) {

                int index = SV.rows();
                SV.conservativeResize(SV.rows() + 1, SV.cols());
                SV.row(index) = V.row(v1) + V.row(v2);

                for (int z = 0; z < Edge.find(key)->second.size(); z++)
                    SV.row(index) += FP.row(Edge.find(key)->second[z]);

                SV.row(index) /= (2 + Edge.find(key)->second.size());

                EP.insert(std::make_pair(key, index));
            }
        }
    }


    // Step 3 -> Find new vertix

    std::vector<std::vector<int>> surfaces;
    vertex_triangle_adjacency(F, V.rows(), surfaces);


    std::vector<std::vector<int>> vertexes;
    vertexes.resize(V.rows());

    for (int j = 0; j < F.rows(); j++) {
        for (int k = 0; k < F.cols(); k++) {
            int v1 = F(j, k);
            int v2 = F(j, (k + 1) % F.cols());

            if (std::find(vertexes[v1].begin(), vertexes[v1].end(), v2) == vertexes[v1].end())
                vertexes[v1].push_back(v2);

            if (std::find(vertexes[v2].begin(), vertexes[v2].end(), v1) == vertexes[v2].end())
                vertexes[v2].push_back(v1);
        }
    }

    for (int j = 0; j < V.rows(); j++) {

        int n = surfaces[j].size();

        // calculating F
        Eigen::RowVector3d cF(0, 0, 0);
        for (int k = 0; k < surfaces[j].size(); k++)
            cF += FP.row(surfaces[j][k]);
        cF /= surfaces[j].size();

        // calculating R
        Eigen::RowVector3d R(0, 0, 0);
        for (int k = 0; k < vertexes[j].size(); k++)
            R += (V.row(j) + V.row(vertexes[j][k]));
        R /= vertexes[j].size();

        Eigen::RowVector3d result = (cF + R + (n - 3) * V.row(j)) / n;

        SV.row(j) = result;
    }

    // Create New Face
    for (int j = 0; j < F.rows(); j++) {
        for (int k = 0; k < F.cols(); k++) {
            int v = F(j, k);

            int v1 = F(j, (k + 1) % F.cols());
            int v2 = F(j, (k + 3) % F.cols());


            std::string key1 = std::to_string(fmin(v, v1)) + "-" + std::to_string(fmax(v, v1));
            std::string key2 = std::to_string(fmin(v, v2)) + "-" + std::to_string(fmax(v, v2));

            int p1 = EP.find(key1)->second;
            int p2 = EP.find(key2)->second;

            SF.row(4 * j + k) = Eigen::RowVector4i(v, p1, V.rows() + j, p2);
        }
    }


    catmull_clark(Eigen::MatrixXd(SV), Eigen::MatrixXi(SF), num_iters - 1, SV, SF);


    ////////////////////////////////////////////////////////////////////////////
}
