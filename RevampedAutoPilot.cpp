#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    cout.precision(6);
    Path path = {Waypoint(0, 0, 0), Waypoint(1, 1, 0) };
    auto curve = curveGenerator(path);
    // TankConfig drive(curve, 2, 2);
    // drive.testTrajectory();
    createDesmosGraph(curve, "", "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\graph.html");
    cout << "ArcLengthTest: " << ArcLengthDistance(curve[0]) << "\n";
    return 0;

    //
}

//probably partion curve into many different segments to maintain width between wheels
//don't know what coordinate funtions r
