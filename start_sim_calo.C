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
#include <ctime>

using namespace std;

int start_sim_calo(int n_samp, string suffix, string output_dir) {

    time_t now = time(0);   
    tm *gmtm = localtime(&now);
    TString timestamp = to_string(1+gmtm->tm_mon) + to_string(gmtm->tm_mday) + to_string(gmtm->tm_hour) + to_string(gmtm->tm_min);
    
    double impact_energy_mean = 50000;
    double tau_0_mean = 0;
    double energu_std = 10000;
    double tau_std = 0.5;
    double tau_0 = 0;
    double impact_energy = 0;
    //double n_samples = 20000;
    int n_samples = n_samp;
    //TString output_dirname = "/eos/atlas/atlascerngroupdisk/phys-sm/LowMu2017WZ/histograms_mykola/Calo_ouput/";
    TString output_dirname  = output_dir;
    TString fileName = "calo_sim_" +TString(to_string(int(n_samples/1000)))+"k_"+ timestamp+ "_" +TString(suffix)+".root";
    TFile* calo_out = TFile::Open(output_dirname+fileName,"recreate");

    
    CaloClusterGenerator* new_cluster = new CaloClusterGenerator();

    TRandom3* RandomGen = new TRandom3();


    for (int i = 0; i <n_samples; i++) {
        if (i < 10 || (i%100 == 0) ) cout << "iteration: " << i << " out of " << n_samples << endl;
        impact_energy = abs(RandomGen->Gaus(impact_energy_mean,energu_std));
        new_cluster->FillClusterEnergy(impact_energy);
        tau_0 = RandomGen->Gaus(tau_0_mean,tau_std);
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
