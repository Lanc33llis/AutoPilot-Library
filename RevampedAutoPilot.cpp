#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    cout.precision(6);
    Path path = {Waypoint(1, 1, 0), Waypoint(2, 3, 0) };
    auto curve = curveGenerator(path);
    // TankConfig drive(curve, 2, 2);
    // drive.testTrajectory();
    createDesmosGraph(curve, "", "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\graph.html");
    // cout << "ArcLengthTest: " << ArcLengthDistance(curve[0]) << "\n";
    Spline XFunction = HermiteFinder(Waypoint( path[0].X, path[0].X, path[0].Angle ), Waypoint( path[1].X, path[1].X, path[1].Angle ));
    Spline YFunction = HermiteFinder(Waypoint( path[0].Y, path[0].Y, path[0].Angle ), Waypoint( path[1].Y, path[1].Y, path[1].Angle ));
    createDesmosGraph(Curve{XFunction}, "", "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\graph1.html");
    createDesmosGraph(Curve{YFunction}, "", "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\graph2.html");
    return 0;
}

//probably partion curve into many different segments to maintain width between wheels
//don't know what coordinate funtions r
