#include <iostream>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include "TPave.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

double pTmin_jet = 0;
double pTmin_hadron = 0;

// Style format. Colours used by various drawn markers.
int colHS = kBlack, colPos = kRed, colNeg = kBlue;
int colNeut = kGreen + 3, colPU = kGray + 1;

using namespace std;
using namespace Pythia8;
using namespace fastjet;

// Draws text to the canvas using non-direct coordinates:
// (x, y)=(0, 0) -> lower-left, (1, 1) -> upper-right.

void drawText(double x, double y, TString txt, int col= kBlack, double tsize = 0.032, int align = 11) {
    static auto tex = new TLatex();
    tex->SetTextColor(col); // Text color
    tex->SetTextSize(tsize); // Text size
    tex->SetTextFont(42); // Text font
    tex->SetNDC(); // Use Normalized Device Coordinates
    tex->SetTextAlign(align); // Text alignment, 11 = center
    tex->DrawLatex(x, y, txt); // Draw the text
}

//==========================================================================

// Draws right-justified text, as above.
void drawTextR(double x, double y, TString txt, int col = kBlack, double tsize = 0.032) {
  drawText(x, y, txt, col, tsize, 31);
}

//==========================================================================

// Text to draw a marker at the (y, phi) coordinates of a particle.
// Absolute coordinates.

void drawParticleMarker(const Particle &p, int style, int col, double size= 1.0) {
  static auto m = new TMarker();
  m->SetMarkerStyle(style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->DrawMarker(p.eta(), p.phi());
}

void drawJetMarker(const PseudoJet &j, int style, int col, double size= 1.0) {
  static auto m = new TMarker();
  m->SetMarkerStyle(style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->DrawMarker(j.eta(), j.phi());
}

//==========================================================================

// Draw a line.

void drawLine(double x1, double y1, double x2, double y2, int col, int style) {
  static auto line = new TLine();
  line->SetLineColor(col);
  line->SetLineStyle(style);
  line->DrawLine(x1, y1, x2, y2);
}

//==========================================================================

//==========================================================================

// Draws a box for text to appear.

void drawLegendBox(double x1, double y1, double x2, double y2) {
  static auto *box = new TPave(x1, y1, x2, y2, 1, "ndc");
  box->SetFillColor(kWhite);
  box->Draw();
}

//==========================================================================

// Draw a marker for legend.

void drawMarker(double x, double y, int style, int col, double size= 1.0) {
  auto m = new TMarker(x, y, style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->SetNDC(true);
  m->Draw();
}

int main() {
    // Adjust ROOTs default style.
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    // Tick marks on top and RHS.
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02, "x");
    gStyle->SetTickLength(0.015, "y");
    // Good with SetMax higher. 57, 91 and 104 also OK.
    gStyle->SetPalette(57);

    // Define the canvas.
    auto can = new TCanvas();
    double x = 0.06, y = 0.96;
    // Left-right-bottom-top
    can->SetMargin(x, 0.02, 0.08, 0.06);

    // Define the energy-flow histogram.
    int NetaBins = 500/2, NphiBins = 314/2;
    double etaMax = 5, phiMax = TMath::Pi();
    auto axis = new TH2F("",";Pseudo-rapidity #it{#eta};Azimuth #it{#phi};Jet #it{p}_{T} [GeV]", NetaBins, -etaMax, etaMax, NphiBins, -phiMax, phiMax);
    auto pTflow = new TH2F("",";Pseudo-rapidity #it{#eta};Azimuth #it{#phi};Jet #it{p}_{T} [GeV]", NetaBins, -etaMax, etaMax, NphiBins, -phiMax, phiMax);
    axis->GetYaxis()->SetTitleOffset(0.8);
    axis->GetZaxis()->SetTitleOffset(1.1);

    // Name of output pdf file + open canvas for printing pages to it.
    TString pdf = "figures.pdf";
    can->Print(pdf + "[");
    // Turn off loads of TCanvas::Print Info messages.
    // Current canvas added to pdf file.
    gErrorIgnoreLevel = 4000;

    // Generator. Process selection. LHC initialization.
    Pythia pythia;
    // Description of the process (using ROOT's TLatex notation).
    TString desc = "#it{p+p} #rightarrow final state hadrons" " #sqrt{#it{s}} = 5 TeV";

    // PYTHIA setup. Process selection. LHC initialization.
    pythia.readString("Beams:idA = 2212"); //sets the first beam to be a proton
    pythia.readString("Beams:idB = 2212"); //sets the second beam to be a proton
    pythia.readString("Beams:eCM = 5000."); // Set center-of-mass energy to 5 TeV
    pythia.readString("HardQCD:all = on");   // Enable QCD processes
    pythia.readString("PhaseSpace:pTHatMin = 50."); // Set minimum pT of hard process to 20 GeV

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Setup fasjet. Create map with (key, value) = (descriptive text, jetDef).
    std::map<TString, fastjet::JetDefinition> jetDefs;
    jetDefs["Anti-#it{k_{t}} jets, #it{R} = 0.4"] = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
    jetDefs["#it{k_{t}} jets, #it{R} = 0.4"] = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
    jetDefs["Cambridge-Aachen jets, #it{R} = 0.4"] = fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);

    // Event loop.
    auto &event = pythia.event;
    const int nEvents = 25;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        // Generate event. (Skip to next if pythia.next() returns false = error.)
        if (!pythia.next()) continue;

        // Clear and draw an empty canvas to paint on.
        axis->Reset();
        axis->Draw();
        drawLine( -etaMax, TMath::Pi(), etaMax, TMath::Pi(), kGray, 7);
        drawLine( -etaMax, -TMath::Pi(), etaMax, -TMath::Pi(), kGray, 7);

        // Pseudojets for the fastjet clustering.
        std::vector<Particle> ptcls_hs;
        std::vector<fastjet::PseudoJet> stbl_ptcls;

        // Draw the final state hadrons
        for (int i = 0; i < event.size(); ++i) {
            auto &p = event[i];
            if (not p.isFinal() && not p.isHadron()) continue;
            if (std::abs(p.eta())>etaMax) continue;
            stbl_ptcls.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
            ptcls_hs.push_back(p);
        }

        // Add in ghost particles on the grid defined by the axis histogram.
        fastjet::PseudoJet ghost;
        double pTghost = 1e-100;
        for (int ieta = 1;ieta <= NetaBins; ++ieta) {
            for (int iphi = 1;iphi <= NphiBins; ++iphi) {
                double eta = axis->GetXaxis()->GetBinCenter(ieta);
                double phi = axis->GetYaxis()->GetBinCenter(iphi);
                ghost.reset_momentum_PtYPhiM(pTghost, eta, phi, 0);
                stbl_ptcls.push_back(ghost);
            }
        }

        can->SetLogz();
        can->SetRightMargin(0.13);
        bool first = true;
        for (auto jetDef:jetDefs) {
            fastjet::ClusterSequence clustSeq(stbl_ptcls, jetDef.second);
            auto jets = sorted_by_pt( clustSeq.inclusive_jets(pTmin_jet) );
            // Fill the pT flow.
            pTflow->Reset();
            axis->Reset();
            // For each jet:
            for (auto jet:jets) {
                // For each particle:
                for (auto c:jet.constituents()) {
                    if (c.pt() > 1e-50) continue;
                    axis->Fill(c.eta(), c.phi_std(), jet.pt());
                }
                pTflow->Fill(jet.eta(), jet.phi_std(), jet.pt());
            }
            axis->GetZaxis()->SetRangeUser(pTmin_jet/2, axis->GetBinContent(axis->GetMaximumBin())*2);
            // pTflow->GetZaxis()->SetRangeUser(8, 1100);
            // pTflow->GetZaxis()->SetMoreLogLabels();
            axis->Draw("colz");
            pTflow->SetMarkerSize(0.4);
            pTflow->SetMarkerStyle(21);
            pTflow->SetMarkerColor(kBlack);
            pTflow->Draw("scat Same0");

            // Draw the stable particles.
            for (auto &p:ptcls_hs) {
                if ( std::abs(p.eta()) < etaMax && p.pT() > pTmin_hadron) {
                    if (p.charge()>0) {
                        drawParticleMarker(p, 5, colPos, 0.8);
                    } else if (p.charge()<0) {
                        drawParticleMarker(p, 5, colNeg, 0.8);
                    } else {
                        drawParticleMarker(p, 21, colNeut, 0.4);
                        drawParticleMarker(p, 5, colNeut, 0.8);
                    }
                }
            }

            // Draw legend box
            drawLegendBox(0.7, 0.7, 0.9, 0.9);
            //drawText(0.71, 0.88, "Legend", kBlack, 0.03);

            // Draw jet marker in legend
            drawMarker(0.72, 0.88, 21, kBlack, 1.0);
            drawText(0.74, 0.88, "Jet-axis", kBlack, 0.03);

            // Draw positive hadron marker in legend
            drawMarker(0.72, 0.82, 5, colPos, 0.8);
            drawText(0.74, 0.82, "Positive hadron", colPos, 0.03);

            // Draw negative hadron marker in legend
            drawMarker(0.72, 0.76, 5, colNeg, 0.8);
            drawText(0.74, 0.76, "Negative hadron", colNeg, 0.03);

            // Draw neutral hadron marker in legend
            drawMarker(0.72, 0.72, 21, colNeut, 0.4);
            drawMarker(0.72, 0.72, 5, colNeut, 0.8);
            drawText(0.74, 0.72, "Neutral hadron", colNeut, 0.03);

            drawText( x, y, desc);
            drawTextR(0.87, y, jetDef.first + Form(", #it{p}_{T} > %.0f GeV", pTmin_jet), 31);
            can->Print(pdf);
        }
        break;
    }
    // Close the pdf
    can->Print(pdf + "]");
    printf( "\nProduced %s\n\n", pdf.Data());

    // Done.
    return 0;
}