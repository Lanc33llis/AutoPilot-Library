#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    Path path = {Waypoint(0, 0, 30), Waypoint(1, 1, 45), Waypoint(2, 2, 45)};
    auto curve = curveGenerator(path);
    TankConfig drive(curve, 2, 2);
    drive.testTrajectory();
    return 0;
} 