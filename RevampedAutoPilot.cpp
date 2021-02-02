#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    Path path = {Waypoint(0, 0, 0), Waypoint(1, 1, 0), Waypoint(2, 2, 0)};
    TankConfig drive(curveGenerator(path), 2, 2);
    drive.testTrajectory();
    return 0;
} 