#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    cout.precision(6);
    Path path = { Waypoint(0, 0, 180), Waypoint(-1, 1, 230), Waypoint(-3, -2, 250)};
    auto curve = curveGenerator(path);
    TankConfig drive(curve, 2, 2);
    drive.testTrajectory();
    createDesmosGraph(drive, "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\graph.html");
    return 0;
}

//probably partion curve into many different segments to maintain width between wheels
