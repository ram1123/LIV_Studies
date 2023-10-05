/*
 * =====================================================================================
 *
 *       Filename:  LHEanalyzer.cpp
 *
 *    Description:  This macro reades the lhe file and convert it into root file.
 *
 *        Version:  1.0
 *        Created:  Friday 02 January 2015 04:45:09  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Nhan Viet Tran
 *	Edited By:  Ramkrishna Sharma
 *   Organization:  CERN
 *
 * =====================================================================================
 */

#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <cmath>
#include "TTree.h"
#include "TClonesArray.h"
#include "TApplication.h"
// #include "LHEanalyzer.h"
// #include "hFactory.h"
// #include "h2Factory.h"
// #include "hFunctions.h"

/*

 TO COMPILE:

 export ROOTSYS=~/Desktop/root
 export PATH=$ROOTSYS/bin:$PATH
 export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
 export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH

 c++ -o analysis01 `root-config --glibs --cflags` -lm hFactory.cc hChain.cc h2Factory.cc h2Chain.cc analysis01.cpp

 TO RUN:

 ./analysis01 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe

 numerazione dei plot:

 2 -- segnale selezionati
 4 -- segnale persi
 0 -- segnale totale
 3 -- bkg irr selezionati
 5 -- bkg irr persi
 1 -- bkg irr totale
 7 -- bkg qcd selezionati
 8 -- bkg qcd persi
 6 -- bkg qcd totale

    // Name of variables:
    // costheta1 :
    // costheta2 :
    // costhetastar :
    // phi :
    // phi1 :
    // costhetaV1 :
    // costhetaV2 :
    //

 */

void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double &costheta1, double &costheta2, double &Phi, double &costhetastar, double &Phi1);

double deltaPhi(double phi1, double phi2)
{
    double deltaphi = fabs(phi1 - phi2);
    if (deltaphi > 6.283185308)
        deltaphi -= 6.283185308;
    if (deltaphi > 3.141592654)
        deltaphi = 6.283185308 - deltaphi;
    return deltaphi;
}

//! ========================================================================================

int main(int argc, char **argv)
{

    // gSystem->Load("libPhysics");

    // TApplication a("a", 0, 0); // just to make sure that the autoloading of ROOT libraries works

    TFile file(argv[2], "RECREATE");
    file.SetCompressionLevel(2);
    //   TFile file("phantom.root","RECREATE");
    TTree *tree = new TTree("tree", "Particles Info");

    // int FS_Leptons = 0;
    std::vector<int> initialQuarks_;
    std::vector<int> finalParticles_;
    std::vector<int> intermediateParticles_;

    //    tree->Branch( "",	&	,	"/F" );
    // tree->Branch("FS_Leptons", &FS_Leptons, "/I");
    tree->Branch("initialQuarks", &initialQuarks_);
    tree->Branch("finalParticles", &finalParticles_);
    tree->Branch("intermediateParticles", &intermediateParticles_);

    // PG loop over bkg
    // PG -------------
    int BKGnumber = 0;
    int BKGnumberWithObj = 0;
    int VBFBKGnumber = 0;
    int count = 0;
    int NSignal = 0;
    int NTotal = 0;
    double SM_Weight = 0.0;

    std::ifstream ifsbkg(argv[1]);
    // Create the Reader object
    LHEF::Reader bkgReader(ifsbkg);

    // PG loop over BKG
    while (bkgReader.readEvent())
    {
        // std::cout << "========================" << std::endl;
        ++BKGnumber;
        if (BKGnumber % 1000 == 0)
            std::cout << "BKG event " << BKGnumber << "\n";
        // if (BKGnumber > 5) break;
        // std::cout<< "Number of particles = "<<bkgReader.hepeup.NUP<<std::endl;
        // std::cout<< "Event weight = "<<bkgReader.hepeup.XWGTUP<<std::endl;
        // std::cout<<"rwgt size = "<<bkgReader.hepeup.namedweights.size()<<std::endl;
        // std::cout<<"rwgt size = "<<bkgReader.hepeup.weights.size()<<std::endl;
        // for (int iPart = 0 ; iPart < bkgReader.hepeup.namedweights.size(); ++iPart){
        // if (bkgReader.hepeup.namedweights[iPart].name == "fs0_0p0")
        // {
        // //std::cout<<iPart+1<<"\tWeight info : "<<std::scientific<<bkgReader.hepeup.namedweights[iPart].weights[0]<<"\t"<< bkgReader.hepeup.namedweights[iPart].name <<std::endl;
        // //std::cout<<iPart+1<<"\tWeight info : "<<std::setprecision(10)<<bkgReader.hepeup.namedweights[iPart].weights[0]<<"\t"<< bkgReader.hepeup.namedweights[iPart].name <<std::endl;
        // SM_Weight = bkgReader.hepeup.namedweights[iPart].weights[0];
        // }
        // }

        std::vector<int> leptons;
        std::vector<int> finalQuarks;
        std::vector<int> intermediates;
        std::vector<int> tops;
        TLorentzVector Is_Iqrk1, Is_Iqrk0;
        // PG loop over particles in the event
        int incomingPart = 0;
        for (int iPart = 0; iPart < bkgReader.hepeup.IDUP.size(); ++iPart)
        {

            int mother1 = bkgReader.hepeup.MOTHUP.at(iPart).first;

            // PG incoming particle
            if (bkgReader.hepeup.ISTUP.at(iPart) == -1)
            {
                initialQuarks_.push_back(bkgReader.hepeup.IDUP.at(iPart));
                incomingPart++;
                if (incomingPart == 1)
                {
                    Is_Iqrk0.SetPxPyPzE(
                        bkgReader.hepeup.PUP.at(incomingPart).at(0), // PG px
                        bkgReader.hepeup.PUP.at(incomingPart).at(1), // PG py
                        bkgReader.hepeup.PUP.at(incomingPart).at(2), // PG pz
                        bkgReader.hepeup.PUP.at(incomingPart).at(3)  // PG E
                    );
                }
                if (incomingPart == 2)
                {
                    Is_Iqrk1.SetPxPyPzE(
                        bkgReader.hepeup.PUP.at(incomingPart).at(0), // PG px
                        bkgReader.hepeup.PUP.at(incomingPart).at(1), // PG py
                        bkgReader.hepeup.PUP.at(incomingPart).at(2), // PG pz
                        bkgReader.hepeup.PUP.at(incomingPart).at(3)  // PG E
                    );
                }
                // std::cout<<"incoming particle = "<<incomingPart<<std::endl;
                count++;
            }

            // PG outgoing particles
            if (bkgReader.hepeup.ISTUP.at(iPart) == 1)
            {
                finalParticles_.push_back(bkgReader.hepeup.IDUP.at(iPart));

                // PG leptons
                if (abs(bkgReader.hepeup.IDUP.at(iPart)) == 11 || // PG electron
                    abs(bkgReader.hepeup.IDUP.at(iPart)) == 13 || // PG muon
                    abs(bkgReader.hepeup.IDUP.at(iPart)) == 15 || // PG tau
                    abs(bkgReader.hepeup.IDUP.at(iPart)) == 12 || // PG neutrino
                    abs(bkgReader.hepeup.IDUP.at(iPart)) == 14 || // PG neutrino
                    abs(bkgReader.hepeup.IDUP.at(iPart)) == 16)   // PG neutrino
                {
                    leptons.push_back(iPart);
                } // PG leptons
                else
                {
                    finalQuarks.push_back(iPart);
                }
            }

            // PG intermediates
            if (bkgReader.hepeup.ISTUP.at(iPart) == 2)
            {
                intermediateParticles_.push_back(bkgReader.hepeup.IDUP.at(iPart));
                intermediates.push_back(iPart);
            }

            // PG tops
            if (abs(bkgReader.hepeup.IDUP.at(iPart)) == 6)
            {
                tops.push_back(iPart);
            }
        } // PG loop over particles in the event

        // --------------- Indices for final state particles  -------------------
        int i_olep_part = -1;
        int i_olep_anti = -1;
        int i_wqrk_1 = -1;
        int i_wqrk_2 = -1;
        int i_iqrk_1 = -1;
        int i_iqrk_2 = -1;

        int tmpfill1 = 0;
        int signalFlag = 0;
        int signalWCtr = 0;




        tree->Fill();
        initialQuarks_.clear();
        finalParticles_.clear();
        intermediateParticles_.clear();

        //	if (BKGnumber > 24000)  break;
    }

    std::cout << "BKGnumber = " << BKGnumber << ", and NSignal = " << NSignal << std::endl;

    file.cd();
    tree->Write();
    file.Close();

    std::cout << "Count = " << count << std::endl;

    // Now we are done.
    return 0;
}

//////////////////////////////////
//// P A P E R   4 - V E C T O R   D E F I N I T I O N   O F   P H I   A N D   P H I 1
//////////////////////////////////
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double &costheta1, double &costheta2, double &Phi, double &costhetastar, double &Phi1)
{

    ///////////////////////////////////////////////
    // check for z1/z2 convention, redefine all 4 vectors with convention
    ///////////////////////////////////////////////
    TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
    p4H = thep4H;

    p4Z1 = thep4Z1;
    p4M11 = thep4M11;
    p4M12 = thep4M12;
    p4Z2 = thep4Z2;
    p4M21 = thep4M21;
    p4M22 = thep4M22;
    //// costhetastar
    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector thep4Z1inXFrame(p4Z1);
    TLorentzVector thep4Z2inXFrame(p4Z2);
    thep4Z1inXFrame.Boost(boostX);
    thep4Z2inXFrame.Boost(boostX);
    TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
    TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());
    costhetastar = theZ1X_p3.CosTheta();

    //// --------------------------- costheta1
    TVector3 boostV1 = -(thep4Z1.BoostVector());
    TLorentzVector p4M11_BV1(p4M11);
    TLorentzVector p4M12_BV1(p4M12);
    TLorentzVector p4M21_BV1(p4M21);
    TLorentzVector p4M22_BV1(p4M22);
    p4M11_BV1.Boost(boostV1);
    p4M12_BV1.Boost(boostV1);
    p4M21_BV1.Boost(boostV1);
    p4M22_BV1.Boost(boostV1);

    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Dot(p4M11_BV1.Vect()) / p4V2_BV1.Vect().Mag() / p4M11_BV1.Vect().Mag();

    //// --------------------------- costheta2
    TVector3 boostV2 = -(thep4Z2.BoostVector());
    TLorentzVector p4M11_BV2(p4M11);
    TLorentzVector p4M12_BV2(p4M12);
    TLorentzVector p4M21_BV2(p4M21);
    TLorentzVector p4M22_BV2(p4M22);
    p4M11_BV2.Boost(boostV2);
    p4M12_BV2.Boost(boostV2);
    p4M21_BV2.Boost(boostV2);
    p4M22_BV2.Boost(boostV2);

    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Dot(p4M21_BV2.Vect()) / p4V1_BV2.Vect().Mag() / p4M21_BV2.Vect().Mag();

    //// --------------------------- Phi and Phi1
    //    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector p4M11_BX(p4M11);
    TLorentzVector p4M12_BX(p4M12);
    TLorentzVector p4M21_BX(p4M21);
    TLorentzVector p4M22_BX(p4M22);

    p4M11_BX.Boost(boostX);
    p4M12_BX.Boost(boostX);
    p4M21_BX.Boost(boostX);
    p4M22_BX.Boost(boostX);

    TVector3 tmp1 = p4M11_BX.Vect().Cross(p4M12_BX.Vect());
    TVector3 tmp2 = p4M21_BX.Vect().Cross(p4M22_BX.Vect());

    TVector3 normal1_BX(tmp1.X() / tmp1.Mag(), tmp1.Y() / tmp1.Mag(), tmp1.Z() / tmp1.Mag());
    TVector3 normal2_BX(tmp2.X() / tmp2.Mag(), tmp2.Y() / tmp2.Mag(), tmp2.Z() / tmp2.Mag());

    //// Phi
    TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;
    double tmpSgnPhi = p4Z1_BX.Vect().Dot(normal1_BX.Cross(normal2_BX));
    double sgnPhi = tmpSgnPhi / fabs(tmpSgnPhi);
    Phi = sgnPhi * acos(-1. * normal1_BX.Dot(normal2_BX));

    //////////////

    TVector3 beamAxis(0, 0, 1);
    TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();

    TVector3 p3V1_BX(tmp3.X() / tmp3.Mag(), tmp3.Y() / tmp3.Mag(), tmp3.Z() / tmp3.Mag());
    TVector3 tmp4 = beamAxis.Cross(p3V1_BX);
    TVector3 normalSC_BX(tmp4.X() / tmp4.Mag(), tmp4.Y() / tmp4.Mag(), tmp4.Z() / tmp4.Mag());

    //// Phi1
    double tmpSgnPhi1 = p4Z1_BX.Vect().Dot(normal1_BX.Cross(normalSC_BX));
    double sgnPhi1 = tmpSgnPhi1 / fabs(tmpSgnPhi1);
    Phi1 = sgnPhi1 * acos(normal1_BX.Dot(normalSC_BX));

    //    std::cout << "extractAngles: " << std::endl;
    //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
    //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;
}
