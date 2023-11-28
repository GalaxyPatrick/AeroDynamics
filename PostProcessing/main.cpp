#define _USE_MATH_DEFINES
#include <cmath>
#include "matplotlibcpp.h"
#include <vector>
#include <fstream>
#include <cmath>

namespace plt = matplotlibcpp;

int main() {

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> x_r, y_r;

    /* Open result.txt file */
    std::ifstream inputFile("result.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open the file." << std::endl;
        return 1;
    }
    double xValue , yValue;
    while (inputFile >> xValue >> yValue) {
        x.push_back(xValue);
        y.push_back(yValue);
    }
    inputFile.close();
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "Read: x = " << x[i] << ", y = " << y[i] << std::endl;
    }
    for (size_t i = 0; i < 101; i++) {
        x_r.push_back(i * 2 * M_PI / 100);
        y_r.push_back(1 - 4 * pow(sin(i * 2 * M_PI / 100), 2));
    }
    //plt::xkcd();
    plt::named_plot("solution_result", x, y, "b+-");
    plt::named_plot("theory_result", x_r, y_r, "r--");
    plt::title("AN ORDINARY SIN WAVE");
    plt::legend();
    plt::show();
    return 0;
}