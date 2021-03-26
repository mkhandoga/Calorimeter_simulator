//
//  main.cpp
//  calo_sim
//
//  Created by Mykola Khandoga on 11/03/2021.
//  Copyright © 2021 Mykola Khandoga. All rights reserved.
//
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "CaloClusterGenerator.cpp"
#include "inv_functions.cpp"
using namespace std;

int start_sim_calo() {

    TString fileName = "calo_sim_output.root";
    TFile* calo_out = TFile::Open(fileName,"recreate");

    CaloClusterGenerator* new_cluster = new CaloClusterGenerator();

    double impact_energy = 50000;
    double tau_0 = 0.5;
    for (int i = 0; i <3; i++) {
        cout << "iteration: " << i << endl;
        new_cluster->FillClusterEnergy(impact_energy);
        new_cluster->FillSignalSamples(tau_0);
        new_cluster->FillRecoEnergyTime();
        new_cluster->FillTree();
        new_cluster->CleanUp();
    }
    
    new_cluster->ClusterTree->Write();    
    delete new_cluster;
    calo_out->Close();
    
    return 0;

}
