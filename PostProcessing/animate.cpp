#include "matplotlibcpp.h"
#include <cmath>

namespace plt = matplotlibcpp;

void update_plot(int frame) {
    // Clear the current plot
    plt::clf();

    // Generate x values
    std::vector<double> x;
    for (double i = -5.0; i <= 5.0; i += 0.1) {
        x.push_back(i);
    }

    // Calculate y values using a sine function with a changing phase
    std::vector<double> y;
    for (double val : x) {
        y.push_back(std::sin(val + 0.1 * frame));
    }

    // Plot the function
    plt::plot(x, y);

    // Show the plot
    plt::pause(0.01);
}

int main() {
    // Set the interactive mode to true
    plt::ion();

    // Create an animation with 100 frames
    for (int i = 0; i < 100; ++i) {
        update_plot(i);
    }

    // Keep the plot window open
    plt::show();

    return 0;
}
