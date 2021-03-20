//
//  CaloClusterGenerator.cpp
//  calo_sim
//
//  Created by Mykola Khandoga on 12/03/2021.
//  Copyright © 2021 Mykola Khandoga. All rights reserved.
//

#include "CaloClusterGenerator.hpp"

using namespace::ROOT;
#include <math.h>
#include <cmath>


CaloClusterGenerator:: CaloClusterGenerator(){
    
    RandomGen = new TRandom3();
    double impact_eta = 1.4*0.025;
    double impact_phi = 1.4*M_PI/128;

    int eta_size = 5; // cells in the second layer
    int phi_size = 5; // cells in the second layer
    
    
    InitCluster(impact_eta, impact_phi, eta_size, phi_size);

}


void CaloClusterGenerator::InitCluster( double impact_eta, double impact_phi, int eta_size, int phi_size){
    
    // cluster is formed defining the size of the second layer
    
    
    my_cluster->layers[1].layer_number = 2;
    my_cluster->layers[1].n_cells_eta = eta_size;
    my_cluster->layers[1].n_cells_phi = phi_size;

    my_cluster->layers[1].eta_cell_size_factor = 1;
    my_cluster->layers[1].phi_cell_size_factor = 1;

    
    // now we can compute the absolute dimensions, assuming that the particle hit the middle of the cluster

    my_cluster->layers[1].eta_max = impact_eta + ((eta_size%2)/2 + eta_size/2.)*g_DeltaEtaS2;
    my_cluster->layers[1].eta_min = impact_eta - ((eta_size%2)/2 + eta_size/2.)*g_DeltaEtaS2;

    my_cluster->layers[1].phi_max = impact_phi + ((phi_size%2)/2 + phi_size/2.)*g_DeltaPhiS2;
    my_cluster->layers[1].phi_min = impact_phi - ((phi_size%2)/2 + phi_size/2.)*g_DeltaPhiS2;
    
    
    //extrapolating cluster size to layer1
    //L1 cells are 8 times smaller in eta and 4 times larger in phi

    my_cluster->layers[0].layer_number = 1;
    my_cluster->layers[0].eta_cell_size_factor = 1/8.;
    my_cluster->layers[0].phi_cell_size_factor = 4;
    
    InitSideLayers(my_cluster->layers[1], my_cluster->layers[0]);

    
    //Extrapolating cluster size to layer3
    //L3 cells are 2 times larger in eta and same size in phi
    my_cluster->layers[2].layer_number = 3;
    my_cluster->layers[2].eta_cell_size_factor = 2;
    my_cluster->layers[2].phi_cell_size_factor = 1;

    InitSideLayers(my_cluster->layers[1],my_cluster->layers[2]);

    /*
    cout << "eta_max " << my_cluster->layers[0].eta_max << " layer1.eta_min " << my_cluster->layers[0].eta_min << " my_cluster->layers[0].n_cells_eta " << my_cluster->layers[0].n_cells_eta  << endl;
    cout << "phi_max" << my_cluster->layers[0].phi_max << " layer1.phi_min " << my_cluster->layers[0].phi_min << " my_cluster->layers[0].n_cells_phi " << my_cluster->layers[0].n_cells_phi  << endl;

    cout << "eta_max " << my_cluster->layers[1].eta_max << " layer2.eta_min " << my_cluster->layers[1].eta_min << " my_cluster->layers[1].n_cells_eta " << my_cluster->layers[1].n_cells_eta  << endl;
    cout << "phi_max" << my_cluster->layers[1].phi_max << " layer2.phi_min " << my_cluster->layers[1].phi_min << " my_cluster->layers[1].n_cells_phi " << my_cluster->layers[1].n_cells_phi  << endl;
    
    cout << "eta_max " << my_cluster->layers[2].eta_max << " layer3.eta_min " << my_cluster->layers[2].eta_min << " my_cluster->layers[2].n_cells_eta " << my_cluster->layers[2].n_cells_eta  << endl;
    cout << "phi_max" << my_cluster->layers[2].phi_max << " layer3.phi_min " << my_cluster->layers[2].phi_min << " my_cluster->layers[2].n_cells_phi " << my_cluster->layers[2].n_cells_phi  << endl;

    cout <<"layer: " << my_cluster->layers[0].layer_number << " layer.layer_cells.size()" << my_cluster->layers[0].layer_cells.size() << endl;
    cout <<"layer: " << my_cluster->layers[1].layer_number << " layer.layer_cells.size()" << my_cluster->layers[1].layer_cells.size() << endl;
    cout <<"layer: " << my_cluster->layers[2].layer_number << " layer.layer_cells.size()" << my_cluster->layers[2].layer_cells.size() << endl;
  */
    my_cluster->layers[0].r_min = g_R0;
    my_cluster->layers[0].r_max = g_R0+g_Rs1;
    
    my_cluster->layers[1].r_min = g_R0+g_Rs1;
    my_cluster->layers[1].r_max = g_R0+g_Rs1+g_Rs2;
    
    my_cluster->layers[2].r_min = g_R0+g_Rs1+g_Rs2;
    my_cluster->layers[2].r_max = g_R0+g_Rs1+g_Rs2+g_Rs3;
    
    
    my_cluster->eta_min = my_cluster->layers[1].eta_min;
    my_cluster->eta_max = my_cluster->layers[1].eta_max;
    my_cluster->phi_min = my_cluster->layers[1].phi_min;
    my_cluster->phi_max = my_cluster->layers[1].phi_max;
    my_cluster->r_min = g_R0;
    my_cluster->r_max = g_R0+g_Rs1+g_Rs2+g_Rs3;
    
    for (int l = 0; l < n_layers; l ++)  InitLayerCells(my_cluster->layers[l]);

}

void CaloClusterGenerator::InitSideLayers(ClusterLayer &central_layer, ClusterLayer &side_layer ){
    
    side_layer.n_cells_eta = ceil ((central_layer.eta_max - central_layer.eta_min) / (side_layer.eta_cell_size_factor*g_DeltaEtaS2));
    side_layer.n_cells_phi = ceil ( (central_layer.phi_max - central_layer.phi_min) / (side_layer.phi_cell_size_factor*g_DeltaEtaS2));
    side_layer.eta_min = central_layer.eta_min;
    side_layer.eta_max = central_layer.eta_min + side_layer.n_cells_eta*g_DeltaEtaS2*side_layer.eta_cell_size_factor;
    side_layer.phi_min = central_layer.phi_min;
    side_layer.phi_max = central_layer.phi_min + side_layer.n_cells_phi*g_DeltaPhiS2*side_layer.phi_cell_size_factor;

}

void CaloClusterGenerator::InitLayerCells(ClusterLayer &layer ){
    ClusterCell cur_cell;

    for (int e = 0; e < layer.n_cells_eta; e++) {
        for (int p = 0; p < layer.n_cells_phi; p++) {

            layer.layer_cells.push_back(cur_cell);
            
            layer.layer_cells.back().eta_in_cluster = e;
            layer.layer_cells.back().phi_in_cluster = p;
            
            layer.layer_cells.back().eta_min = layer.eta_min + e*layer.eta_cell_size_factor*g_DeltaEtaS2;
            layer.layer_cells.back().eta_max = layer.eta_min + (e+1)*layer.eta_cell_size_factor*g_DeltaEtaS2;
            
            layer.layer_cells.back().phi_min = layer.phi_min + p*layer.phi_cell_size_factor*g_DeltaPhiS2;;
            layer.layer_cells.back().phi_max = layer.phi_min + (p+1)*layer.phi_cell_size_factor*g_DeltaPhiS2;
            //cout <<"layer: " << layer.layer_number << " phi_min " << layer.layer_cells.back().phi_min << endl;
            //cout <<"layer: " << layer.layer_number << " phi_max " << layer.layer_cells.back().phi_max << endl;

            layer.layer_cells.back().r_min = layer.r_min;
            layer.layer_cells.back().r_max = layer.r_max;
            //cout <<"layer: " << layer.layer_number << " layer.layer_cells.back().r_min " << layer.layer_cells.back().r_min << endl;
            //cout <<"layer: " << layer.layer_number << " layer.layer_cells.back().r_max " << layer.layer_cells.back().r_max << endl;

        }
    }
    //cout <<"layer: " << layer.layer_number << " layer.layer_cells.size()" << layer.layer_cells.size() << endl;
}

void CaloClusterGenerator::FillClusterEnergy(double e_impact){
    double layer_energy_theoretical, layer_energy_sumcells,  cluster_energy_theoretical, cluster_energy_sumlayers;
    cluster_energy_sumlayers = 0;
    cout << "Starting filling the cluster, impact energy: " << e_impact << endl;
    for (int l = 0; l < n_layers; l++){
        cout << "Starting filling layer " << my_cluster->layers[l].layer_number << endl;
        layer_energy_sumcells = 0;
        for (int c = 0; c < my_cluster->layers[l].layer_cells.size(); c++){
            cout << "cell number: " << c << " eta_min: " << my_cluster->layers[l].layer_cells[c].eta_min << endl;
            cout << "cell number: " << c << " eta_max: " << my_cluster->layers[l].layer_cells[c].eta_max << endl;
            cout << "cell number: " << c << " phi_min: " << my_cluster->layers[l].layer_cells[c].phi_min << endl;
            cout << "cell number: " << c << " phi_max: " << my_cluster->layers[l].layer_cells[c].phi_max << endl;
             my_cluster->layers[l].layer_cells[c].E_truth  = e_impact*this->IntegEnerCell(my_cluster->layers[l].layer_cells[c].eta_min, my_cluster->layers[l].layer_cells[c].eta_max, my_cluster->layers[l].layer_cells[c].phi_min, my_cluster->layers[l].layer_cells[c].phi_max, my_cluster->layers[l].layer_cells[c].r_min, my_cluster->layers[l].layer_cells[c].r_max);
            layer_energy_sumcells+= my_cluster->layers[l].layer_cells[c].E_truth;
            cout << "cell number: " << c << " energy: " << my_cluster->layers[l].layer_cells[c].E_truth << endl;
            cout << "sum_cells in layer: " << my_cluster->layers[l].layer_number << " energy: " << layer_energy_sumcells << endl;

        }
        
        layer_energy_theoretical = e_impact*this->IntegEnerCell(my_cluster->layers[l].eta_min, my_cluster->layers[l].eta_max, my_cluster->layers[l].phi_min, my_cluster->layers[l].phi_max, my_cluster->layers[l].r_min, my_cluster->layers[l].r_max);
        cout << "energy sum_cells: " << layer_energy_sumcells << " energy theoretical: " << layer_energy_theoretical << endl;

        cluster_energy_sumlayers+=layer_energy_theoretical;
    }
    cluster_energy_theoretical = e_impact*this->IntegEnerCell(my_cluster->eta_min, my_cluster->eta_max, my_cluster->phi_min, my_cluster->phi_max, my_cluster->r_min, my_cluster->r_max);
    cout << " cluster energy sum layers: " << cluster_energy_sumlayers << " energy theoretical: " << cluster_energy_theoretical << endl;
};

void CaloClusterGenerator::GetClusterLimits( Double_t Lay, Double_t& dPhi, Double_t& dEta, Double_t& nPhi, Double_t& nEta, Double_t& Rmin, Double_t& Rmax ){
    if( Lay == 1){
        /* Layer = 1;    */
        dPhi  = 4*g_DeltaPhiS2 ;          // ~0.09817
        dEta  = g_DeltaEtaS2/8 ;          // Eta increment on layer
        nPhi  = ceil(g_nPhiS2*g_DeltaPhiS2/dPhi) ;
        nEta  = ceil(g_nEtaS2*g_DeltaEtaS2/dEta) ;
        Rmin  = g_R0 ;
        Rmax  = g_R0 + g_Rs1 ;
    }
    else if(Lay == 2){
        /* Layer = 2; */
        dPhi  = g_DeltaPhiS2;             // M_PI/128 ~ 0.02545
        dEta  = g_DeltaEtaS2;             // 0.025 ;
        nPhi  = g_nPhiS2;
        nEta  = g_nEtaS2;
        Rmin  = g_R0 + g_Rs1 ;
        Rmax  = Rmin + g_Rs2 ;
    }
    else if(Lay == 3){
        /* Layer = 3 ; */
        dPhi  = g_DeltaPhiS2 ;            // same size of S2
        dEta  = 2*g_DeltaEtaS2 ;
        nPhi  = ceil(g_nPhiS2*g_DeltaPhiS2/dPhi) ;
        nEta  = ceil(g_nEtaS2*g_DeltaEtaS2/dEta) ;
        Rmin  = g_R0 + g_Rs1 + g_Rs2 ;
        Rmax  = Rmin + g_Rs3 ;
    }
}

void CaloClusterGenerator::AdjustCoordinates( Double_t &dEta, Double_t &nEta, Double_t &EtaMin, Double_t &dPhi, Double_t &nPhi, Double_t &PhiMin ){
    Double_t eta0, phi0 ;
    
    if ( Int_t(nPhi)%2 != true ){
        phi0 = nPhi/2 ;                 // nPhi odd
    }else{
        phi0 = floor(nPhi/2 + 1) ;      // nPhi even
    }
    if ( Int_t(nEta)%2 != true ){
        eta0 = nEta/2 ;                 // nEta odd
    }else{
        eta0 = floor(nEta/2 + 1) ;      // nEta even
    }
    
    // Adjust of the PhiMin with the cordinates of impact to calculate the integral always on a full cell
    if ( g_PhiImpact/dPhi != floor(g_PhiImpact/dPhi) ){
        PhiMin = dPhi*( floor(g_PhiImpact/dPhi) - floor(nPhi/2.)  ) ;
        //PhiMin = dPhi*(floor( g_PhiImpact/dPhi ) - phi0 ) ;
    }
    else{
        // Adjust the coordinate to the first bottom cell of the cluster
        PhiMin = g_PhiImpact - dPhi*(floor(nPhi/2.)  ) ;
        //PhiMin = g_PhiImpact - dPhi*phi0 ;
    }
    // Adjust of the EtaMin with the coordinates of the impact to calculate the integral always on a full cell
    if ( g_EtaImpact/dEta != floor(g_EtaImpact/dEta) ){
        EtaMin = dEta*( floor(g_EtaImpact/dEta) - floor(nEta/2.)  ) ;
        //EtaMin = dEta*( floor(g_EtaImpact/dEta) - eta0 ) ;
    }
    else{
        // Adjust the coordinate to the first left cell of the cluster. **This implies that the maximum energy is on the center of the matrix
        EtaMin = g_EtaImpact - dEta*(floor(nEta/2.)  ) ;
        //EtaMin = g_EtaImpact - dEta*eta0 ;
    }
}// end of function

vector <double> CaloClusterGenerator::GenerateRandomDelays(Int_t n_etaL2, Int_t n_PhiL2, Int_t n_samples){
    vector <double> cell_delays_collection;
    for (Int_t i = 0; i < n_etaL2*n_PhiL2*n_samples; i++){
        if (i == 0){
            cell_delays_collection.push_back(0) ;
        }
        else{
            // The delay beetwen the source clock and each one cell 0.3 ns
            cell_delays_collection.push_back(abs(25 - RandomGen->Gaus(25, 0.3))) ;
            cout << cell_delays_collection.back() << endl;
        }
    }
    return cell_delays_collection;
}

// Function to calculatin Energie by profiles radial and longitudinal
/* Double_t IntegEnerCell( Double_t tMin , Double_t tMax, Double_t EtaMin , Double_t EtaMax , Double_t PhiMin, Double_t PhiMax ){ */
Double_t CaloClusterGenerator::IntegEnerCell( Double_t EtaMin ,
                       Double_t EtaMax ,
                       Double_t PhiMin,
                       Double_t PhiMax,
                       Double_t Rmin,
                       Double_t Rmax ){
  
    //TF3* IntegraCell = new TF3("IntegraCell", "cos([3]) ");
    //cout << "sdfsdfs"<<sech(4) << endl;
    cout << "EtaMin: " << EtaMin << " EtaMax " << EtaMax << "PhiMin: " << PhiMin << " PhiMax " << PhiMax << endl;
    cout << "Rmin: " << Rmin << " Rmax " << Rmax << endl;

    //EtaMin: -0.0375 EtaMax -0.034375PhiMin: -0.0981748 PhiMax 0;
    /*EtaMin = -0.0375;
    EtaMax = -0.034375;
    PhiMin = -0.0981748;
    PhiMax = 0;
    Rmin = 1500;
    Rmax = 1590;
     */
    
//EtaMin: 0.0125 EtaMax 0.0375PhiMin: 0.0122718 PhiMax 0.0368155
/*
    EtaMin = 0.0125;
    EtaMax = 0.0375;
     PhiMin = 0.0368155;
     PhiMax = 0.0613592;
    Rmin = 1590;
    Rmax = 1927;
 */
    TF3* IntegraCell = new TF3("IntegraCell", "(0.0262470945510689 * pow(0.11661276855717342,-1 + [0]) * pow([1],[0]) * abs((z * cosh([3]) * cosh(x)) / sqrt(cosh(2 * [3]) + cosh(2 * x) - 2 * cos([2] - y) * (cos([2] - y) + 2 * sinh([3]) * sinh(x)))) * pow(abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3])),-1 + [0]) * ((2 * (1 - 2.6150800000000003 * pow(exp(1),-pow(exp(1),(0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) + (0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) - 0.30939226519337015 * (0.348 - 0.4491923844822321 / pow(exp(1),pow(-1 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]),2)))) * pow(-0.04331491712707183 - 0.3463399226148051 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + 0.60338 * (pow(exp(1),-2.59 * (-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]))) + pow(exp(1),(-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) * (0.3585 + 0.0421 * TMath::Log([5])))),2) * z * sqrt(1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3]))) / pow(pow(-0.04331491712707183 - 0.3463399226148051 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + 0.60338 * (pow(exp(1),-2.59 * (-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]))) + pow(exp(1),(-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) * (0.3585 + 0.0421 * TMath::Log([5])))),2) + pow(z,2) * (1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3])),2) + (2 * (2.6150800000000003 * pow(exp(1),-pow(exp(1),(0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) + (0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) + 0.30939226519337015 * (0.348 - 0.4491923844822321 / pow(exp(1),pow(-1 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]),2)))) * z * pow(0.018819337016574587 + 0.027777161470318713 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + (0.012750673339578456 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]) + 0.00319 * TMath::Log([5]),2) * sqrt(1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3]))) / pow(pow(0.018819337016574587 + 0.027777161470318713 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + (0.012750673339578456 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]) + 0.00319 * TMath::Log([5]),2) + pow(z,2) * (1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3])),2))) / (pow(exp(1),0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) * TMath::Gamma([0]))") ;
    
    IntegraCell-> SetParameter(0, g_a) ;
    IntegraCell-> SetParameter(1, g_b) ;
    IntegraCell-> SetParameter(2, g_PhiImpact) ;
    IntegraCell-> SetParameter(3, g_EtaImpact) ;
    IntegraCell-> SetParameter(4, g_R0) ;
    IntegraCell-> SetParameter(5, g_E0) ;
    IntegraCell-> SetParameter(6, g_Rmoleff) ;
    
    //cout << "PhiImpact " << PhiImpact << " EtaImpact = " << EtaImpact << endl;
    // /* ********************************* */
    // zero as a limit is a singular value and the integral diverges and produce an inconsistent value
    if ( EtaMin < 0 & EtaMax == 0 ){
        EtaMax = -1e-6 ;
    }
    if ( EtaMin == 0 & EtaMax > 0 ){
        EtaMin = 1e-6 ;
    }
    if ( PhiMin == 0 & PhiMax > 0 ){
        PhiMin = 1e-6 ;
    }
    if ( PhiMin < 0 & PhiMax == 0 ){
        PhiMax = -1e-6 ;
    }
    // x = eta; y = phi; z = R;
  //  cout << " integral: " << IntegraCell -> Integral( EtaMin, EtaMax, PhiMin, PhiMax, Rmin, Rmax, 1e-5 ) << endl;
    return IntegraCell -> Integral( EtaMin, EtaMax, PhiMin, PhiMax, Rmin, Rmax, 1e-5 ) ;
} // end of function

Int_t CaloClusterGenerator::EnergyMeVtoADC(Double_t E_cell, Int_t Layer){
    Double_t F1 = 1. , // F1 = F(µA  to MeV) => MeV/µA
    F2 = 1. , // F2 = F(DAC to µA) => µA/DAC
    /// F1*F2 = MeV/DAC => 1/(F1*F2) => DAC/MeV
    E_dac = 0. ;
    
    if (Layer == 1){
        F1 = 374.14 ;
        F2 = 0.024  ;
    }
    else if (Layer == 2){
        F1 = 374.14 ;
        F2 = 0.074  ;
    }
    else if (Layer == 3){
        F1 = 316.61 ;
        F2 = 0.074  ;
    }
    E_dac = E_cell/(F1*F2) + g_ADCpedestal ;
    
    return E_dac ;
} // end of function

void CaloClusterGenerator::GetClusterEnergy( Bool_t fluctEnergy,
                   vector<Double_t> &E0_impact,
                   vector<Double_t> &sumEnergyLayer,
                   vector<Double_t> &EnergyPerCell_S1,
                   vector<Double_t> &EnergyADC_S1,
                   vector<Double_t> &EnergyPerCell_S2,
                   vector<Double_t> &EnergyADC_S2,
                   vector<Double_t> &EnergyPerCell_S3,
                   vector<Double_t> &EnergyADC_S3 ){
    Double_t    PhiMin   = 0. ,
    PhiMin0  = 0. ,
    PhiMax   = 0. ,
    dPhi     = 0. ,
    nPhi     = 0. ,
    EtaMin   = 0. ,
    EtaMax   = 0. ,
    dEta     = 0. ,
    nEta     = 0. ,
    tMin     = 0. ,
    tMax     = 0. ,
    Rmin     = 0. ,
    Rmax     = 0. ,
    sumS1    = 0. ,
    sumS2    = 0. ,
    sumS3    = 0. ,
    EnerCell = 0. ;
    TMatrixD S2_Energy(5,5);
    TMatrixD S3_Energy(5,5);
    sumEnergyLayer.clear() ;
    EnergyPerCell_S1.clear() ;
    EnergyPerCell_S2.clear() ;
    EnergyPerCell_S3.clear() ;
    EnergyADC_S1.clear() ;
    EnergyADC_S2.clear() ;
    EnergyADC_S3.clear() ;
    for(Int_t layer = 1 ; layer <= 3 ; ++layer){
        EnerCell = 0 ;
        this->GetClusterLimits( layer, dPhi, dEta, nPhi, nEta, Rmin, Rmax ) ;
        
        // Adjust of the coordinates of impact point to calculating energy on full size of the each one cells
        this->AdjustCoordinates( dEta, nEta, EtaMin, dPhi, nPhi, PhiMin ) ;
        PhiMin0 = PhiMin;
        ///* *************************** */
        // loop to cover eta
        //cout << "rows " << nEta << " " << EnergyPerCell_S1.GetNrows() << " Column " << nPhi << " " << EnergyPerCell_S1.GetNcols() << endl;
        
        for( Int_t iEta = 0 ; iEta < nEta ; ++iEta ){//
            EtaMax  = EtaMin + dEta ;
            // loop to cover phi
            for( Int_t iPhi = 0 ; iPhi < nPhi ; ++iPhi ){//
                PhiMax = PhiMin + dPhi ;
                
                //EnerCell = g_E0*IntegEnerCell( EtaMin, EtaMax, PhiMin, PhiMax, Rmin, Rmax ) ;
                // /************************************/
                // /* 0.1*sqrt(g_E0/1e3) - This value is the fluctuation on the impact energy */
                // /************************************/
                if (fluctEnergy == true ){
                    E0_impact.push_back(1e3*RandomGen->Gaus( g_E0/1e3, 0.1*sqrt(g_E0/1e3))) ;
                }else{
                    //cout << "No fluctuations on Energy" << endl;
                    E0_impact.push_back(g_E0) ;
                }
                EnerCell = E0_impact.back()*this->IntegEnerCell( EtaMin, EtaMax, PhiMin, PhiMax, Rmin, Rmax ) ;
                
                if (layer == 1){
                    EnergyPerCell_S1.push_back(EnerCell) ;
                    EnergyADC_S1.push_back(this->EnergyMeVtoADC(layer, EnerCell)) ;
                    sumS1 += EnerCell/1e3;
                }else if (layer == 2){
                    EnergyPerCell_S2.push_back(EnerCell) ;
                    EnergyADC_S2.push_back(this->EnergyMeVtoADC(layer, EnerCell)) ;
                    S2_Energy(iEta,iPhi) = EnerCell/1e3 ;
                    sumS2 += EnerCell/1e3;
                }else if (layer == 3){
                    EnergyPerCell_S3.push_back(EnerCell) ;
                    EnergyADC_S3.push_back(this->EnergyMeVtoADC(layer, EnerCell)) ;
                    S3_Energy(iEta,iPhi) = EnerCell/1e3 ;
                    sumS3 += EnerCell/1e3;
                }
                
                sumEnergyLayer[layer-1] += EnerCell ;
                
                PhiMin += dPhi ;
            }//end loop for phi
            PhiMin = PhiMin0;
            EtaMin += dEta ;
        } // end loop for eta
        
        //S2_Energy.Print() ;
        //S3_Energy.Print() ;
        //cout << "sum S1 = " << setprecision(4)<< sumS1 << " GeV  => " << 100*sumS1/(sumS1+sumS2+sumS3) << " % \n" <<
        //        "sum S2 = " << setprecision(4)<< sumS2 << " GeV  => " << 100*sumS2/(sumS1+sumS2+sumS3) << " % \n" <<
        //        "sum S3 = " << setprecision(4)<< sumS3 << " GeV  => " << 100*sumS3/(sumS1+sumS2+sumS3) << " % \n" << endl;
        //cin.get() ;
    } // end of loop for layer
    
    //cout << "Sum S1: " << sumEnergyLayer[0]/1e3 << " | Sum S2: "  << sumEnergyLayer[1]/1e3 << " | Sum S3: "  << sumEnergyLayer[2]/1e3 << endl;
    //cout << "S1:  "<< sumEnergyLayer[0]/500 << " % | S2: "  << sumEnergyLayer[1]/500 << " % | S3: "<< sumEnergyLayer[2]/500 << " %" << endl;
} // end of function

void  CaloClusterGenerator::GenerateCellSamples(vector<Double_t> EnergyPerCell_S2,
                 vector<Double_t> &samp_Energy ,
                 vector<Double_t> &samp_Signal ,
                 vector<Double_t> &samp_dSignal,
                 vector<Double_t> &samp_Xt_C  ,
                 vector<Double_t> &samp_dXt_C ,
                 vector<Double_t> &samp_Xt_L  ,
                 vector<Double_t> &samp_dXt_L ,
                 vector<Double_t> &samp_Noise ,
                 vector<Double_t> &samp_SigNoise ,
                 vector<Double_t> &samp_SigXt_C ,
                 vector<Double_t> &samp_SigXt_L ,
                 vector<Double_t> &samp_SigXt_CL ,
                 vector<Double_t> &samp_SigNoiseXt_C ,
                 vector<Double_t> &samp_SigNoiseXt_L ,
                 vector<Double_t> &samp_SigNoiseXt_CL ){
    
    samp_Energy.clear();
    samp_Signal.clear();
    samp_dSignal.clear();
    samp_Xt_C.clear();
    samp_dXt_C.clear();
    samp_Xt_L.clear();
    samp_dXt_L.clear();
    samp_Noise.clear();
    samp_SigNoise.clear();
    samp_SigXt_C.clear();
    samp_SigXt_L.clear();
    samp_SigXt_CL.clear();
    samp_SigNoiseXt_C.clear();
    samp_SigNoiseXt_L.clear();
    samp_SigNoiseXt_CL.clear();
    // 100 samples is because of we have 5x5 cells on the cluster with 4 sample each
    g_tDelayCell = GenerateRandomDelays(g_nEtaS2, g_nPhiS2, int(g_nSamp));

    for (int i = 0; i < g_tDelayCell.size(); i++)  cout << "i:" << i << " " << g_tDelayCell[i] << endl;
    for (UInt_t i = 0; i < g_nEtaS2*g_nPhiS2*g_nSamp; i++){
        // cin.get();

        samp_Noise.push_back(g_AmpNoise/g_ToNormNoise*RandomGen->Gaus (0, 2)) ;

        samp_Signal.push_back(GenerateSignalSample()-> Eval( (i%4 + 1 )*25 + g_tDelayCell[ i ], 0., 0.)) ;
        samp_dSignal.push_back(GenerateSignalSample()-> Derivative( (i%4 + 1)*25 + g_tDelayCell[ i ])) ;

        samp_Xt_C.push_back(g_AmpXt_C*GenerateXTalkSample()-> Eval( (i%4 + 1)*25 + g_tDelayCell[ i ], 0., 0.)) ;

        samp_dXt_C.push_back(g_AmpXt_C*GenerateXTalkSample()-> Derivative( (i%4 + 1)*25 + g_tDelayCell[ i ])) ;
        samp_Xt_L.push_back(g_AmpXt_L*GenerateXTalkSample()-> Derivative( (i%4 + 1)*25 + g_tDelayCell[ i ])) ;
        //samp_dXt_L.push_back(g_AmpXt_C*XTalk()-> Derivative2( (i%4 + 1)*25 + g_tDelayCell[ i ])) ;
        samp_Energy.push_back(EnergyPerCell_S2[ i/4 ]*samp_Signal[i]) ;
        //// Signal combinations
        
        /*         samp_SigNoise.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise.back() + samp_Signal.back()) ) ;
         samp_SigXt_C.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Signal.back() + samp_Xt_C.back()) ) ;
         
         samp_SigXt_L.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Signal.back() + samp_Xt_L.back()) ) ;
         samp_SigXt_CL.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Signal.back() + samp_Xt_L.back() + samp_Xt_C.back()) ) ;
         
         samp_SigNoiseXt_C.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise.back() + samp_Signal.back() + samp_Xt_C.back()) ) ;
         samp_SigNoiseXt_L.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise.back() + samp_Signal.back() + samp_Xt_L.back()) ) ;
         
         samp_SigNoiseXt_CL.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise.back() + samp_Signal.back() + samp_Xt_C.back() + samp_Xt_L.back()) ) ;
         
         */
        
        samp_SigNoise.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise[i] + samp_Signal[i]) ) ;
        samp_SigXt_C.push_back (EnergyPerCell_S2[ i/4 ]*(samp_Signal[i] + samp_Xt_C[i]) ) ;
        
        samp_SigXt_L.push_back (EnergyPerCell_S2[ i/4 ]*(samp_Signal[i] + samp_Xt_L[i]) ) ;
        samp_SigXt_CL.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Signal[i] + samp_Xt_L[i] + samp_Xt_C[i]) ) ;
        
        samp_SigNoiseXt_C.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise[i] + samp_Signal[i] + samp_Xt_C[i]) ) ;
        samp_SigNoiseXt_L.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise[i] + samp_Signal[i] + samp_Xt_L[i]) ) ;
        
        samp_SigNoiseXt_CL.push_back(EnergyPerCell_S2[ i/4 ]*(samp_Noise[i] + samp_Signal[i] + samp_Xt_C[i] + samp_Xt_L[i]) ) ;
        
    }
}

void CaloClusterGenerator::FillSignalSamples(double tau_0){
    double cell_delay = 0;
    for (int l = 0; l < n_layers; l++){
        cout << "Starting filling layer " << my_cluster->layers[l].layer_number << endl;
        for (int c = 0; c < my_cluster->layers[l].layer_cells.size(); c++){
            for (int s = 1; s <= n_samples; s++){
                cell_delay = RandomGen->Gaus (0, 0.3);
                my_cluster->layers[l].layer_cells[c].sampling_noise.push_back(g_AmpNoise/g_ToNormNoise*RandomGen->Gaus (0, 2));
                my_cluster->layers[l].layer_cells[c].samples_truth.push_back(GenerateSignalSample()-> Eval( s*g_tSamp +cell_delay, 0., 0.));
                my_cluster->layers[l].layer_cells[c].dsamples_truth.push_back(GenerateSignalSample()-> Derivative( s*g_tSamp + cell_delay, 0, 0.));

                my_cluster->layers[l].layer_cells[c].Xt_C.push_back(g_AmpXt_C*GenerateXTalkSample()-> Eval( s*g_tSamp + cell_delay));
                my_cluster->layers[l].layer_cells[c].dXt_C.push_back(g_AmpXt_C*GenerateXTalkSample()-> Derivative( s*g_tSamp + cell_delay));

                my_cluster->layers[l].layer_cells[c].Xt_L.push_back(g_AmpXt_L*GenerateXTalkSample()-> Eval( s*g_tSamp + cell_delay));
                
                cout << "cell: " << c << " sample: " << s << " samples truth: " <<  my_cluster->layers[l].layer_cells[c].samples_truth.back() << endl;
                cout << "cell: " << c << " sample: " << s << " xtc: " <<  my_cluster->layers[l].layer_cells[c].Xt_C.back() << endl;
                cout << "cell: " << c << " sample: " << s << " Xt_L: " <<  my_cluster->layers[l].layer_cells[c].Xt_L.back() << endl;

            }
        }
    }
}


TF1* CaloClusterGenerator::GenerateSignalSample(){
    
    TF1* CellM = new TF1("CellM","((exp(-x/[0])*x*x)/(2 *[0]*[0]*([0] - [1])) - (exp(-(x/[0]))*x*[1])/([0]*pow([0] - [1],2)) + exp(-(x/[0]))*[1]*[1]/pow([0] - [1],3) + (exp(-(x/[1]))*[1]*[1])/pow(-[0] + [1],3) + (1/(2*[2]*[0] *pow(([0] - [1]),3)))*exp(-x* (1/[0] + 1/[1]))* (-2 *exp(x *(1/[0] + 1/[1]))*[0] *pow(([0] - [1]),3) - 2 *exp(x/[0])*[0]*pow([1],3) + exp(x/[1]) *(x*x *pow(([0] - [1]),2) + 2*x*[0]*([0]*[0] - 3*[0]*[1] + 2*[1]*[1]) + 2*[0]*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))) + ((1 - (exp((-x + [2])/[0])*(x - [2])*([0] - 2*[1]))/pow(([0] - [1]),2) - (exp((-x + [2])/[0])*(x - [2])*(x- [2]))/(2*[0]*([0] - [1])) + (exp((-x + [2])/[1])*[1]*[1]*[1])/pow(([0] - [1]),3) - (exp((-x + [2])/[0])*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))/pow(([0] - [1]),3))* Heaviside(x -[2]))/[2])*[3]*[4]*[3]*[0]*[0]",0., g_tmax2);
    CellM ->SetParameter(0, g_taud);
    CellM ->SetParameter(1, g_taupa);
    CellM ->SetParameter(2, g_td);
    CellM ->SetParameter(3, g_Rf);
    CellM ->SetParameter(4, g_C1);
    return CellM ;
}// end of function

TF1* CaloClusterGenerator::GenerateXTalkSample(){
    TF1* XTalk = new TF1("Xtalk","(1/[6]*(([0]*[1]*[2]*(2*exp(x/[3])*pow([3],2)*(x*([3] - [4])*([4] + [5]) + [4]*(3*[3]*[4] + (2*[3] + [4])*[5]))- exp(x/[4]) * (pow(x,2)*pow([3] - [4],2)*([3] + [5]) - 2*x*[3]*([3] - [4])*(2*[3]*[4] + ([3] + [4])*[5]) +  2*pow([3],2)*[4]*(3*[3]*[4] + (2*[3] + [4])*[5])) + [3]*(-2*exp(x/[3] + [5]/[4])*[3]*[4]*(x*([3] - [4]) + 3*[3]*[4] + (-[3] + [4])*[5]) +  exp(x/[4] + [5]/[3]) * (pow(x,2)*pow([3] - [4],2) + 6*pow([3],2)*pow([4],2) +  4*[3]*([3] - [4])*[4]*[5] + pow([3] - [4],2)*pow([5],2) -  2*x*([3] - [4]) * (2*[3] *[4] + ([3] - [4])*[5])))*Heaviside(x - [5]))))/(2.*exp(x*(1/[3] + 1/[4]))*[3]*pow([3]- [4],4)*[5]))",0., g_tmax2);
    
    XTalk ->SetParameter(0,g_Cx);
    XTalk ->SetParameter(1,g_Rf);
    XTalk ->SetParameter(2,g_Rin);
    XTalk ->SetParameter(3,g_taud);
    XTalk ->SetParameter(4,g_taupa);
    XTalk ->SetParameter(5,g_td);
    XTalk ->SetParameter(6,g_ToNormXtC);   // This value normalize amplitude of the Xt_C to unit = 0.022206
    
    return XTalk ;
}

void CaloClusterGenerator::FillTree(TTree* calo_tree){
    
    
    

}
