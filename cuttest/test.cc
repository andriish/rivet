#include <Rivet/Cuts.hh>
#include <iostream>

using namespace std;


int main() {
    cout << "start" << endl;
    Rivet::FourMomentum test(1000,1000,1000,1000);
    Rivet::CutPtr ptr = Rivet::PtGtr(6.);
    Rivet::MakeCuttable<Rivet::FourMomentum> mom(test);

    cout << ptr->cut(mom) << endl;
    return 0;
}
