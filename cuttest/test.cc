#include <Rivet/Cuts.hh>
#include <iostream>

using namespace std;


int main() {
    cout << "FourMomentum:" << endl;

    /// Particle and Jet both inherit FM

    Rivet::FourMomentum fm(10,10,10,10); /// pT ~ 14.14
    Rivet::MakeCuttable<Rivet::FourMomentum> mom(fm);

    for(int i = 0; i < 10; i++) {
        Rivet::CutPtr ptr = Rivet::PtIn(i*2.,i*6.);
        cout << i << " ";
        cout << ptr->cut(mom) << endl;
    }

    cout << endl;

    cout << "PseudoJet:" << endl;

    fastjet::PseudoJet pjet(10,10,10,10); /// pT ~ 14.14
    Rivet::MakeCuttable<fastjet::PseudoJet> jet(pjet);

    for(int i = 0; i < 10; i++) {
        Rivet::CutPtr ptr = Rivet::PtIn(i*2.,i*6.);
        cout << i << " ";
        cout << ptr->cut(jet) << endl;
    }

    return 0;
}
