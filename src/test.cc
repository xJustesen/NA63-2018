// EXT. LIBRARYS
#include <iostream>
#include <fstream>

// ROOT LIBRARYS
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVersionCheck.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TObject.h"
#include "TStorage.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraph.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TNtuple.h"
#include <armadillo>

using namespace std;

int main (int argc, char const *argv[]) {

        /* Open root file and get tree */
        gROOT->ProcessLine("#include <vector>");
        TFile *f = TFile::Open("test.root", "READ");
        TTree* T1 = (TTree *)f->Get("T");
        TLeaf *NHits = (TLeaf *)T1->GetLeaf("NHits");
        TLeaf *x =  (TLeaf *)T1->GetLeaf("x");
        TLeaf *y =  (TLeaf *)T1->GetLeaf("y");
        TLeaf *plane =  (TLeaf *)T1->GetLeaf("plane");

        /* Extract data */
        int N_events = T1->GetEntries();
        cerr << N_events << endl;

        for (int i  = 0; i < N_events; i++) {

                T1->GetEntry(i);
                int N = NHits->GetValue();
                cout << x->GetValue() << endl;

        }

        // string filename = "test.out";
        // ofstream output (filename);
        //
        // for (size_t i = 0; i < Events.size(); i++) { // events
        //
        //         for (size_t j = 0; j < Events[i].size(); j++) { // hits
        //
        //                 output << Events[i][j][0] << ' ' << Events[i][j][1] << ' ' << Events[i][j][2] << "\n";
        //
        //         }
        //
        // }
        //
        return 0;
}
