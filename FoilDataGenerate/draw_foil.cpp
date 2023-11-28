#define _USE_MATH_DEFINES
#include <cmath>
#include "matplotlibcpp.h"
#include <vector>
#include <fstream>
#include <cmath>
//default output path $(SolutionDir)$(Platform)\$(Configuration)
namespace plt = matplotlibcpp;

extern "C" void draw() {

    std::vector<double> x_u, y_u, x_l, y_l;
    std::ifstream input_upper("foil_u.txt");
    std::ifstream input_lower("foil_l.txt");
    if (!input_lower.is_open() || !input_upper.is_open()) {
        std::cerr << "Unable to open the file." << std::endl;
        return;
    }
    double tmp_x, tmp_y;
    while (input_upper >> tmp_x >> tmp_y) {
        x_u.push_back(tmp_x);
        y_u.push_back(tmp_y);
    }
    while (input_lower >> tmp_x >> tmp_y) {
        x_l.push_back(tmp_x);
        y_l.push_back(tmp_y);
    }
    input_upper.close();
    input_lower.close();
    plt::figure_size(1000, 250);
    plt::xlim(0, 1);
    plt::ylim(-0.125, 0.125);
    plt::plot(x_u, y_u, { {"color", "red"}, {"linestyle", "--"}, {"linewidth", "0.5"} });
    plt::plot(x_l, y_l, { {"color", "blue"}, {"linestyle", "--"}, {"linewidth", "0.5"} });
    //plt::legend();
    plt::show();
    //return 0;
}