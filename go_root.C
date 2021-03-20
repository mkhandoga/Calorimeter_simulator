//
//  main.cpp
//  calo_sim
//
//  Created by Mykola Khandoga on 11/03/2021.
//  Copyright © 2021 Mykola Khandoga. All rights reserved.
//
//#include "TH1D.h"
//#include "TDatime.h"
//#include "TMatrixD.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "CaloClusterGenerator.cpp"
#include "inv_functions.cpp"

#include <tgmath.h>
//#include "VarsConstants.cpp"

//#include <stdio.h>
using namespace std;

int go_root() {
    vector<Double_t> Signal,                    // Cell signal
    dSignal,                   // first cell signal
    Xt_C,                      // Capacitive crosstalk signal
    dXt_C,                     // first derivative of the capacitive crosstalk signal
    Xt_L,                      // Inductive crosstalk signal
    dXt_L,                     // first derivative of the inductive crosstalk signal
    Noise,
    CompSignal,                // composite signal
    samp_Energy,
    samp_Signal ,              // Cell signal samples
    samp_dSignal,              //
    samp_Xt_C ,
    samp_dXt_C,
    samp_Xt_L ,
    samp_dXt_L,
    samp_Noise,
    samp_SigNoise ,
    samp_SigNoiseXt_C ,
    samp_SigNoiseXt_L ,
    samp_SigNoiseXt_CL,
    samp_SigXt_C ,
    samp_SigXt_L ,
    samp_SigXt_CL ,
    samp_Comp,
    EnergyPerCell_S1,          // energies on the S1 per cell
    EnergyPerCell_S2,          // energies on the S2 per cell
    EnergyPerCell_S3,          // energies on the S3 per cell
    sigmEnerCluster(3),
    E0_impact ,
    sumEnergyLayer(3),
    EnergyADC_S1,              // Energy in ADC counts
    EnergyADC_S2,              // Energy in ADC counts
    EnergyADC_S3;              // Energy in ADC counts
    TMatrixD EMshowerMatrix(5*5, 5*5 );
    
    
    TString pathOut = "/Users/mykola/Physics/CaloML/My_sim/" ;
    TString fileName = "calo_sim_output.root";
    TFile* calo_out = TFile::Open(pathOut+fileName,"recreate");
    TTree *SampleSignals = new TTree ( treename.c_str(), "Samples for each one cell on the each one Cluster") ;

//    TDatime *dateTime = new TDatime;
//    UInt_t iDate = dateTime->GetDate(),nClusters;
//    nClusters = 10;

//    TFile *fileSignSamp = TFile::Open( "/Users/mykola/Physics/CaloML/My_sim/hz.root", "RECREATE") ;
//    TTree *SampleSignals = new TTree ( "huy", "Samples for each one cell on the each one Cluster");
 //   CalcCellXtalkNoise( Signal, dSignal, Xt_C, dXt_C, Xt_L, dXt_L, Noise, CompSignal, iDate, pathOut );
    
    CaloClusterGenerator* new_cluster = new CaloClusterGenerator();
    
    //new_cluster->GenerateRandomDelays(5, 5, 4);
    double impact_energy = 50000;
    double tau_0 = 0.2;
    new_cluster->FillClusterEnergy(impact_energy);
    new_cluster->FillSignalSamples(tau_0);
    new_cluster->FillSignalSamples(tau_0);

    cout << "cells in l1: " <<  new_cluster->layers
    //new_cluster->FillRecoEnergy();
    //new_cluster->FillNTuples(TTree);
    //new_cluster->Reset();
    
    //new_cluster->GetClusterEnergy(0, E0_impact, sumEnergyLayer, EnergyPerCell_S1, EnergyADC_S1, EnergyPerCell_S2, EnergyADC_S2, EnergyPerCell_S3, EnergyADC_S3 );

    return 0;

}
