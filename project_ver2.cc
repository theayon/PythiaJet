#include <iostream>
#include <vector>
#include <fastjet/ClusterSequence.hh>
#include <Pythia8/Pythia.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TApplication.h>

using namespace std;
using namespace Pythia8;
using namespace fastjet;

int main() {
    // Initialize Pythia
    Pythia pythia;
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 5000.");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 50.");
    pythia.init();

    // Define FastJet parameters
    double R = 0.4;
    JetDefinition jetDef(antikt_algorithm, R, E_scheme, Best);

    // Vectors to store pT of jets and partons
    std::vector<double> jetPts;
    std::vector<double> partonPts;

    // Generate events
    for (int iEvent = 0; iEvent < 5000; ++iEvent) {
        if (!pythia.next()) continue;

        // Vector to store final state particles
        std::vector<PseudoJet> particles;
        double partonPt = 0.0;

        // Loop over particles in the event
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal() && pythia.event[i].isHadron()) {
                particles.push_back(PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
            }
            if (pythia.event[i].statusAbs() == 23) { // Hard parton
                partonPt = pythia.event[i].pT();
            }
        }

        // Cluster particles into jets
        ClusterSequence cs(particles, jetDef);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

        // Store the pT of the leading jet and the corresponding parton
        if (!jets.empty()) {
            jetPts.push_back(jets[0].pt());
            partonPts.push_back(partonPt);
        }
    }

    // Create the histograms
    TH2F *hist = new TH2F("hist", "Jet pT vs Parton pT; Jet p_{T} (GeV/c); Parton p_{T} (GeV/c)", 100, 0, 100, 50, 50, 100);

    // Fill the histogram
    for (size_t i = 0; i < jetPts.size(); ++i) {
        hist->Fill(jetPts[i], partonPts[i]);
    }

    // Create the canvas
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetRightMargin(0.13);

    // Set the style
    gStyle->SetOptStat(0);

    // Draw the histogram
    canvas->cd();
    canvas->SetLogz();
    hist->Draw("colz");

    // Save the canvas
    canvas->Print("jet_vs_parton_pt.pdf");
}