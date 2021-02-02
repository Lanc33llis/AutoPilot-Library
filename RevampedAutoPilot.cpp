#include "AutoPilot.hpp"

using namespace std;

int main() {
    cout << "Tesing AP\n";
    auto s = HermiteFinder(Waypoint(0, 0, 0), Waypoint(1, 1, 0));
    // cout << "Values: " << s.function.A << " " << s.function.B << " " << s.function.C << " " << s.function.D;
    auto ss = Segment(s, 2);
    cout << "Values: " << ss.;
}