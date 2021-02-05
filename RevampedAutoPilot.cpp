#include "AutoPilot.hpp"

using namespace std;

string desktop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";
string laptop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";

int main() {
    cout << "Tesing AP\n";
    cout.precision(10);
    Path path = {Waypoint(0, 0, 0), Waypoint(1, 1, 0)};
    auto curve = curveGenerator(path);
    TankConfig drive(curve, 2, 2);
    // drive.testTrajectory();
    createDesmosGraph(drive, "", desktop + "graph.html");
    // Spline XSpline = HermiteFinder(Waypoint( 0, path[0].X, 45), Waypoint( path[1].X - path[0].X, path[1].X, 45));
    // Spline YSpline = HermiteFinder(Waypoint( 0, path[0].Y, path[0].Angle ), Waypoint( path[1].X - path[0].X, path[1].Y, path[1].Angle ));
    // createDesmosGraph(Curve{XSpline}, "", desktop + "graph1.html");
    // createDesmosGraph(Curve{YSpline}, "", desktop + "graph2.html");
    // cout << "ArcLengthTest: " << ArcLengthDistance(curve[0], 1, -1) << "\n";
    // cout << "ArcLengthToXValue Test: " << ArcLengthToXValue(curve[0], 1, -.25) << "\n";
    // cout << "Distance test: " << distance(0, 0, 2, 3    ) << "\n";
    return 0;

    //
}

//probably partion curve into many different segments to maintain width between wheels
//don't know what coordinate funtions r
