#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    cout.precision(10);
    Path path = {Waypoint(0, 0, 0), Waypoint(1, 1, 0)};
    auto curve = curveGenerator(path);
    // TankConfig drive(curve, 2, 2);
    // drive.testTrajectory();
    createDesmosGraph(curve, "", "C:\\Users\\scdel\\Desktop\\Coding\\AutoPilot-Library\\graph.html");
    Spline XSpline = HermiteFinder(Waypoint( 0, path[0].X, 45), Waypoint( path[1].X - path[0].X, path[1].X, 45));
    Spline YSpline = HermiteFinder(Waypoint( 0, path[0].Y, path[0].Angle ), Waypoint( path[1].X - path[0].X, path[1].Y, path[1].Angle ));
    createDesmosGraph(Curve{XSpline}, "", "C:\\Users\\scdel\\Desktop\\Coding\\AutoPilot-Library\\graph1.html");
    createDesmosGraph(Curve{YSpline}, "", "C:\\Users\\scdel\\Desktop\\Coding\\AutoPilot-Library\\graph2.html");
    cout << "ArcLengthTest: " << ArcLengthDistance(curve[0], .5, 1) << "\n";
    // cout << "Distance test: " << distance(0, 0, 2, 3    ) << "\n";
    return 0;

    //
}

//probably partion curve into many different segments to maintain width between wheels
//don't know what coordinate funtions r
