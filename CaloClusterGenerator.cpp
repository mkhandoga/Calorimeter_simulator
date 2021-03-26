//
//  CaloClusterGenerator.cpp
//  calo_sim
//
//  Created by Mykola Khandoga on 12/03/2021.
//  Copyright © 2021 Mykola Khandoga. All rights reserved.
//

#include "CaloClusterGenerator.hpp"

using namespace::std;
#include <math.h>
#include <cmath>
#include "TVectorD.h"
#include "TMatrixD.h"


CaloClusterGenerator:: CaloClusterGenerator(){
    
    RandomGen = new TRandom3();
//    this->cluster_center_eta = 0.4*0.025;
//    this->cluster_center_phi = 2.15*M_PI/128;
    cluster_center_eta = 0.4*0.025;
    cluster_center_phi = 0.4*M_PI/128;
    
    
    int eta_size = 5; // cells in the second layer
    int phi_size = 5; // cells in the second layer
    
    
    InitCluster(cluster_center_eta, cluster_center_phi, eta_size, phi_size);
    
    string TreeName = "CaloSimTree";
    InitTree(TreeName);
    InitFunctions();
    InitOpFiltCoefficients();
    
    //for (int i = 0; i < 100; i++) cout << " i: " << i << " value: " << GenerateSignalSample()->Eval(i) << endl;
    
}

void CaloClusterGenerator::InitOpFiltCoefficients(){
    
    vector<double> gt, dgt, noise;
    for (int s = 0; s < n_samples; s++) {
        gt.push_back(GenerateSignalSample()->Eval(g_tSamp*(s+1)));
        dgt.push_back(GenerateSignalSample()->Derivative(g_tSamp*(s+1)));
        noise.push_back(g_AmpNoise/g_ToNormNoise*RandomGen->Gaus(0,2));
    }

    auto ofcoefs = OptimalFilteringInit(gt, dgt,noise);
    
    ai = ofcoefs.first;
    bi = ofcoefs.second;
    
   
    
    //cout << "ai2" << ai[2] << endl;

}


void CaloClusterGenerator::InitFunctions() {
    IntegraCell = new TF3("IntegraCell", "(0.0262470945510689 * pow(0.11661276855717342,-1 + [0]) * pow([1],[0]) * abs((z * cosh([3]) * cosh(x)) / sqrt(cosh(2 * [3]) + cosh(2 * x) - 2 * cos([2] - y) * (cos([2] - y) + 2 * sinh([3]) * sinh(x)))) * pow(abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3])),-1 + [0]) * ((2 * (1 - 2.6150800000000003 * pow(exp(1),-pow(exp(1),(0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) + (0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) - 0.30939226519337015 * (0.348 - 0.4491923844822321 / pow(exp(1),pow(-1 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]),2)))) * pow(-0.04331491712707183 - 0.3463399226148051 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + 0.60338 * (pow(exp(1),-2.59 * (-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]))) + pow(exp(1),(-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) * (0.3585 + 0.0421 * TMath::Log([5])))),2) * z * sqrt(1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3]))) / pow(pow(-0.04331491712707183 - 0.3463399226148051 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + 0.60338 * (pow(exp(1),-2.59 * (-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]))) + pow(exp(1),(-0.645 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) * (0.3585 + 0.0421 * TMath::Log([5])))),2) + pow(z,2) * (1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3])),2) + (2 * (2.6150800000000003 * pow(exp(1),-pow(exp(1),(0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) + (0.43466000000000005 - (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) / (1.313 - 0.0686 * TMath::Log([5]))) + 0.30939226519337015 * (0.348 - 0.4491923844822321 / pow(exp(1),pow(-1 + (0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]),2)))) * z * pow(0.018819337016574587 + 0.027777161470318713 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + (0.012750673339578456 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]) + 0.00319 * TMath::Log([5]),2) * sqrt(1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3]))) / pow(pow(0.018819337016574587 + 0.027777161470318713 / pow(exp(1),(0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0])) + (0.012750673339578456 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) / (-1 + [0]) + 0.00319 * TMath::Log([5]),2) + pow(z,2) * (1 - ((2 + cos(2 * [2] - 2 * y) - cosh(2 * x)) * pow(sech([3]),2)) / 2. - 2 * cos([2] - y) * sech([3]) * sinh(x) * tanh([3])),2))) / (pow(exp(1),0.11661276855717342 * [1] * abs(-([4] * cosh([3])) + z * cos([2] - y) * sech([3]) + z * sinh(x) * tanh([3]))) * TMath::Gamma([0]))") ;

    CellM = new TF1("CellM","((exp(-x/[0])*x*x)/(2 *[0]*[0]*([0] - [1])) - (exp(-(x/[0]))*x*[1])/([0]*pow([0] - [1],2)) + exp(-(x/[0]))*[1]*[1]/pow([0] - [1],3) + (exp(-(x/[1]))*[1]*[1])/pow(-[0] + [1],3) + (1/(2*[2]*[0] *pow(([0] - [1]),3)))*exp(-x* (1/[0] + 1/[1]))* (-2 *exp(x *(1/[0] + 1/[1]))*[0] *pow(([0] - [1]),3) - 2 *exp(x/[0])*[0]*pow([1],3) + exp(x/[1]) *(x*x *pow(([0] - [1]),2) + 2*x*[0]*([0]*[0] - 3*[0]*[1] + 2*[1]*[1]) + 2*[0]*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))) + ((1 - (exp((-x + [2])/[0])*(x - [2])*([0] - 2*[1]))/pow(([0] - [1]),2) - (exp((-x + [2])/[0])*(x - [2])*(x- [2]))/(2*[0]*([0] - [1])) + (exp((-x + [2])/[1])*[1]*[1]*[1])/pow(([0] - [1]),3) - (exp((-x + [2])/[0])*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))/pow(([0] - [1]),3))* Heaviside(x -[2]))/[2])*[3]*[4]*[3]*[0]*[0]",0., g_tmax2);
    
    XTalk = new TF1("Xtalk","(1/[6]*(([0]*[1]*[2]*(2*exp(x/[3])*pow([3],2)*(x*([3] - [4])*([4] + [5]) + [4]*(3*[3]*[4] + (2*[3] + [4])*[5]))- exp(x/[4]) * (pow(x,2)*pow([3] - [4],2)*([3] + [5]) - 2*x*[3]*([3] - [4])*(2*[3]*[4] + ([3] + [4])*[5]) +  2*pow([3],2)*[4]*(3*[3]*[4] + (2*[3] + [4])*[5])) + [3]*(-2*exp(x/[3] + [5]/[4])*[3]*[4]*(x*([3] - [4]) + 3*[3]*[4] + (-[3] + [4])*[5]) +  exp(x/[4] + [5]/[3]) * (pow(x,2)*pow([3] - [4],2) + 6*pow([3],2)*pow([4],2) +  4*[3]*([3] - [4])*[4]*[5] + pow([3] - [4],2)*pow([5],2) -  2*x*([3] - [4]) * (2*[3] *[4] + ([3] - [4])*[5])))*Heaviside(x - [5]))))/(2.*exp(x*(1/[3] + 1/[4]))*[3]*pow([3]- [4],4)*[5]))",0., g_tmax2);

}


void CaloClusterGenerator::InitCluster( double cluster_center_eta, double cluster_center_phi, int eta_size, int phi_size){
    
    // cluster is formed defining the size of the second layer
    
    
    my_cluster->layers[1].layer_number = 2;
    my_cluster->layers[1].n_cells_eta = eta_size;
    my_cluster->layers[1].n_cells_phi = phi_size;

    my_cluster->layers[1].eta_cell_size_factor = 1;
    my_cluster->layers[1].phi_cell_size_factor = 1;

    
    // now we can compute the absolute dimensions, assuming that the particle hits the middle cell of the cluster

    my_cluster->layers[1].eta_max = cluster_center_eta + ((eta_size%2)/2 + eta_size/2.)*g_DeltaEtaS2;
    my_cluster->layers[1].eta_min = cluster_center_eta - ((eta_size%2)/2 + eta_size/2.)*g_DeltaEtaS2;
    
    my_cluster->layers[1].phi_max = cluster_center_phi + ((phi_size%2)/2 + phi_size/2.)*g_DeltaPhiS2;
    
    my_cluster->layers[1].phi_min = cluster_center_phi - ((phi_size%2)/2 + phi_size/2.)*g_DeltaPhiS2;
    
    
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
    
    for (int l = 0; l < n_layers; l ++)  {
        InitLayerCells(my_cluster->layers[l]);
        my_cluster->layers[l].eta_cell_size = my_cluster->layers[l].eta_cell_size_factor*g_DeltaEtaS2;
        my_cluster->layers[l].phi_cell_size = my_cluster->layers[l].phi_cell_size_factor*g_DeltaPhiS2;
    }
    
    cout <<"layer: " << my_cluster->layers[0].layer_number << " layer.layer_cells.size()" << my_cluster->layers[0].layer_cells.size() << endl;
    /*
    int k = 2;
    cout <<"layer etamax: " << my_cluster->layers[k].eta_max << " cells etamax:" << my_cluster->layers[k].layer_cells.back().eta_max << endl;
    cout <<"layer etamin: " << my_cluster->layers[k].eta_min << " cells etamin:" << my_cluster->layers[k].layer_cells[0].eta_min << endl;
    cout <<"layer phimax: " << my_cluster->layers[k].phi_max << " cells phimax:" << my_cluster->layers[k].layer_cells.back().phi_max << endl;
    cout <<"layer phimin: " << my_cluster->layers[k].phi_min << " cells phimin:" << my_cluster->layers[k].layer_cells[0].phi_min << endl;
    //cout <<"layer: " << my_cluster->layers[1].layer_number << " layer.layer_cells.size()" << my_cluster->layers[1].layer_cells.size() << endl;
    //cout <<"layer: " << my_cluster->layers[2].layer_number << " layer.layer_cells.size()" << my_cluster->layers[2].layer_cells.size() << endl;
     */
    
}

void CaloClusterGenerator::InitSideLayers(ClusterLayer &central_layer, ClusterLayer &side_layer ){
    
    side_layer.n_cells_eta = ceil ((central_layer.eta_max - central_layer.eta_min) / (side_layer.eta_cell_size_factor*g_DeltaEtaS2));
    side_layer.n_cells_phi = ceil ( (central_layer.phi_max - central_layer.phi_min) / (side_layer.phi_cell_size_factor*g_DeltaPhiS2));
    side_layer.eta_min = central_layer.eta_min;
    side_layer.eta_max = central_layer.eta_min + side_layer.n_cells_eta*g_DeltaEtaS2*side_layer.eta_cell_size_factor;
    side_layer.phi_min = central_layer.phi_min;
    side_layer.phi_max = central_layer.phi_min + side_layer.n_cells_phi*g_DeltaPhiS2*side_layer.phi_cell_size_factor;
    //cout << "side_layer.phi_max: " << side_layer.phi_max << " central_layer.phi_max: " << central_layer.phi_max << endl;
}

void CaloClusterGenerator::InitLayerCells(ClusterLayer &layer ){
    ClusterCell cur_cell;
    pair<int,int> address;
    for (int e = 0; e < layer.n_cells_eta; e++) {
        for (int p = 0; p < layer.n_cells_phi; p++) {
            address = make_pair(e,p);
            //layer.layer_cells.push_back(cur_cell);
            layer.layer_cells_map[address] = cur_cell;

            layer.layer_cells_map[address].eta_in_cluster = e;
            layer.layer_cells_map[address].phi_in_cluster = p;
            
            layer.layer_cells_map[address].eta_min = layer.eta_min + e*layer.eta_cell_size_factor*g_DeltaEtaS2;
            layer.layer_cells_map[address].eta_max = layer.eta_min + (e+1)*layer.eta_cell_size_factor*g_DeltaEtaS2;
            
            layer.layer_cells_map[address].phi_min = layer.phi_min + p*layer.phi_cell_size_factor*g_DeltaPhiS2;;
            layer.layer_cells_map[address].phi_max = layer.phi_min + (p+1)*layer.phi_cell_size_factor*g_DeltaPhiS2;
            //cout <<"layer: " << layer.layer_number << " phi_min " << layer.layer_cells_map[address].phi_min << endl;
            //cout <<"layer: " << layer.layer_number << " phi_max " << layer.layer_cells_map[address].phi_max << endl;

            layer.layer_cells_map[address].r_min = layer.r_min;
            layer.layer_cells_map[address].r_max = layer.r_max;
            layer.layer_cells_map[address].layer = layer.layer_number;

            //cout <<"layer: " << layer.layer_number << " layer.layer_cells_map[address].r_min " << layer.layer_cells_map[address].r_min << endl;
            //cout <<"layer: " << layer.layer_number << " layer.layer_cells_map[address].r_max " << layer.layer_cells_map[address].r_max << endl;

        }
    }
    //cout <<"layer: " << layer.layer_number << " layer.layer_cells.size()" << layer.layer_cells.size() << endl;
}

void CaloClusterGenerator::FillClusterEnergy(double e_impact){
    //defining the impact spot - uniformly random within the hottest cell
    
    double rand_eta = RandomGen->Uniform();
    double rand_phi = RandomGen->Uniform();
    my_cluster->E_impact = e_impact;

    //rand_eta = 1/2;
    //rand_phi = 1/2;
    
    impact_eta = cluster_center_eta - g_DeltaEtaS2/2 +rand_eta*g_DeltaEtaS2;
    impact_phi = cluster_center_phi - g_DeltaPhiS2/2 +rand_phi*g_DeltaPhiS2;
    
    cout << "cluster_center_eta: " << cluster_center_eta <<  " impact_eta: " << impact_eta <<" rand_eta: "<<rand_eta << endl;
    cout << "cluster_center_phi: " << cluster_center_phi <<  " impact_phi: " << impact_phi << " rand_phi: "<<rand_phi <<  endl;
    
    double layer_energy_theoretical, cluster_energy_sumcells = 0.0, layer_energy_sumcells,  cluster_energy_theoretical, cluster_energy_sumlayers;
    cluster_energy_sumlayers = 0;
    cluster_energy_sumcells = 0;
    
    cout << "Starting filling the cluster, impact energy: " << e_impact << endl;
    for (int l = 0; l < n_layers; l++){
        cout << "Starting filling layer " << my_cluster->layers[l].layer_number << endl;
        layer_energy_sumcells = 0;
        map<pair<int,int>,ClusterCell>::iterator cell_iterator;
        for (cell_iterator = my_cluster->layers[l].layer_cells_map.begin(); cell_iterator != my_cluster->layers[l].layer_cells_map.end(); cell_iterator++){

            /*
            cout << "cell number: " << cell_iterator->second.eta_in_cluster << " eta_min: " << cell_iterator->second.eta_min << endl;
            cout << "cell number: " << cell_iterator->second.eta_in_cluster << " eta_max: " << cell_iterator->second.eta_max << endl;
            cout << "cell number: " << cell_iterator->second.phi_in_cluster << " phi_min: " << cell_iterator->second.phi_min << endl;
            cout << "cell number: " <<  cell_iterator->second.phi_in_cluster << " phi_max: " << cell_iterator->second.phi_max << endl;
            cout << "cell layer: " <<  cell_iterator->second.layer << " phi_max: " << cell_iterator->second.phi_max << endl;
            */
            GetListOfNeighbours(cell_iterator->second.layer, cell_iterator->second.eta_in_cluster, cell_iterator->second.phi_in_cluster, true);

            //auto k = *cell_iterator->second;
             cell_iterator->second.E_truth  = e_impact*this->IntegEnerCell(cell_iterator->second.eta_min, cell_iterator->second.eta_max, cell_iterator->second.phi_min, cell_iterator->second.phi_max, cell_iterator->second.r_min, cell_iterator->second.r_max);
            
            layer_energy_sumcells+= cell_iterator->second.E_truth;
            
            //cout << "cell number: " << c << " energy: " << cell_iterator->second.E_truth << endl;
            //cout << "sum_cells in layer: " << my_cluster->layers[l].layer_number << " energy: " << layer_energy_sumcells << endl;

        }
        my_cluster->layers[l].E_truth = layer_energy_sumcells;

        layer_energy_theoretical = e_impact*this->IntegEnerCell(my_cluster->layers[l].eta_min, my_cluster->layers[l].eta_max, my_cluster->layers[l].phi_min, my_cluster->layers[l].phi_max, my_cluster->layers[l].r_min, my_cluster->layers[l].r_max);
        cout << "energy sum_cells: " << layer_energy_sumcells << " energy theoretical: " << layer_energy_theoretical << endl;
        cluster_energy_sumcells+=layer_energy_sumcells;

        cluster_energy_sumlayers+=layer_energy_theoretical;

    }
     cluster_energy_theoretical = e_impact*this->IntegEnerCell(my_cluster->eta_min, my_cluster->eta_max, my_cluster->phi_min, my_cluster->phi_max, my_cluster->r_min, my_cluster->r_max);
    cout << " cluster energy sum layers: " << cluster_energy_sumlayers << " energy theoretical: " << cluster_energy_theoretical << endl;
    my_cluster->E_truth = cluster_energy_sumcells;
    FillXtalkAmplitudes(true, true);

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
  

    
    IntegraCell->SetNpx(1000);
    IntegraCell->SetNpy(1000);
    IntegraCell->SetNpz(100);
    IntegraCell-> SetParameter(0, g_a) ;
    IntegraCell-> SetParameter(1, g_b) ;
    IntegraCell-> SetParameter(2, impact_phi) ;
    IntegraCell-> SetParameter(3, impact_eta) ;
    IntegraCell-> SetParameter(4, g_R0) ;
    IntegraCell-> SetParameter(5, g_E0) ;
    IntegraCell-> SetParameter(6, g_Rmoleff) ;
    

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
    double integral = IntegraCell -> Integral( EtaMin, EtaMax, PhiMin, PhiMax, Rmin, Rmax, 1e-6 );
    
    //It happens sometimes that the integral can't converge with some specific values of eta/phi
    // a very simple regularization is proposed here
    
    double epsilon = 1e-6;
    int r = 0;
    
    if (isnan(integral)){
        cout << "Integral is NaN. Trying to regularize." << endl;
        while (isnan(integral) || r < 3) {
            integral = IntegraCell -> Integral( EtaMin, EtaMax-r*epsilon, PhiMin+r*epsilon, PhiMax, Rmin, Rmax, 1e-6 );
            r++;
            cout << "Regularization iteration " << r << " out of 3. Integral: " << integral << endl;
        }
    }

    return integral ;
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
    my_cluster->tau_0 = tau_0;

    double cell_delay = 0;
    
    for (int l = 0; l < n_layers; l++){
        cout << "Starting filling layer with signal samples with tau_0: "<< tau_0 << " at layer: "  << my_cluster->layers[l].layer_number << endl;
        map<pair<int,int>,ClusterCell>::iterator cell_iterator;
        for (cell_iterator = my_cluster->layers[l].layer_cells_map.begin(); cell_iterator != my_cluster->layers[l].layer_cells_map.end(); cell_iterator++){
            for (int s = 1; s <= n_samples; s++){
                cell_delay = (RandomGen->Gaus(0, 0.3));
                cell_iterator->second.sampling_delay.push_back(cell_delay);
                cell_iterator->second.sampling_noise.push_back(g_AmpNoise/g_ToNormNoise*RandomGen->Gaus (0, 2));
                cell_iterator->second.samples_truth.push_back(GenerateSignalSample()-> Eval( s*g_tSamp +cell_delay + tau_0, 0., 0.));
                cell_iterator->second.dsamples_truth.push_back(GenerateSignalSample()-> Derivative( s*g_tSamp + cell_delay+ tau_0, 0, 0.));

                cell_iterator->second.Xt_C.push_back(GenerateXTalkSample()-> Eval( s*g_tSamp + cell_delay+tau_0));
                cell_iterator->second.dXt_C.push_back(GenerateXTalkSample()-> Derivative( s*g_tSamp + cell_delay+tau_0));

                cell_iterator->second.Xt_L.push_back(GenerateXTalkSample()-> Eval( s*g_tSamp + cell_delay + tau_0));
            //    cell_iterator->second.samples_truthXtC.push_back(cell_iterator->second.samples_truth.back()+cell_iterator->second.Xt_C.back());
            //    cell_iterator->second.dsamples_truthXtC.push_back(cell_iterator->second.dsamples_truth.back()+cell_iterator->second.dXt_C.back());

               // cout << "cell: " << c << " sample: " << s << " samples truth: " <<  cell_iterator->second.samples_truth.back() << endl;
                //cout << "cell: " << c << " sample: " << s << " xtc: " <<  cell_iterator->second.Xt_C.back() << endl;
                //cout << "cell: " << c << " sample: " << s << " Xt_L: " <<  cell_iterator->second.Xt_L.back() << endl;

            }
        
         //   auto fake_etime = OptimalFiltering(cell_iterator->second.E_truth, cell_iterator->second.samples_truth, cell_iterator->second.dsamples_truth, cell_iterator->second.sampling_noise);
        //    auto reco_etime = GetRecoEnergyTime(cell_iterator->second.E_truth*cell_iterator->second.samples_truth+cell_iterator->second.sampling_noise);
            
      //      cout << " now with xtalk" << endl;
      //  OptimalFiltering(cell_iterator->second.samples_truthXtC, cell_iterator->second.dsamples_truthXtC, cell_iterator->second.sampling_noise);


        }
    }
}

void CaloClusterGenerator::FillRecoEnergyTime(){
    
    vector<double> signal;
    pair <double, double> energy_time_reco;
    for (int l = 0; l < n_layers; l++) {
        for (auto &cell : my_cluster->layers[l].layer_cells_map){

            for (int s = 0; s < cell.second.samples_truth.size(); s++){
            //    cout << " s: " << s << endl;
            //    cout << "cell.second.E_truth" << cell.second.E_truth << " cell.second.samples_truth[s] " << cell.second.samples_truth[s] <<  " cell.second.sampling_noise[s] " << cell.second.sampling_noise[s] << endl;

                signal.push_back(cell.second.E_truth*cell.second.samples_truth[s]+cell.second.sampling_noise[s]);
            }
            energy_time_reco = GetRecoEnergyTime(signal);
            
           // cout << "energy_time_reco.first: " << energy_time_reco.first << endl;
          //  cout << "energy_time_reco.second: " << energy_time_reco.second << endl;

            cell.second.E_reco_noise = energy_time_reco.first;
            cell.second.tau_reco_noise = energy_time_reco.second;
            
            energy_time_reco = OptimalFiltering(cell.second.E_truth, cell.second.samples_truth, cell.second.dsamples_truth, cell.second.sampling_noise);
            
            cell.second.E_reco_fake = energy_time_reco.first;
            cell.second.tau_reco_fake = energy_time_reco.second;

            signal.clear();
            for (int s = 0; s < cell.second.samples_truth.size(); s++){
                signal.push_back(cell.second.E_truth*cell.second.samples_truth[s]+cell.second.sampling_noise[s]+cell.second.Xt_C_amplitude*cell.second.Xt_C[s]+cell.second.Xt_L_amplitude*cell.second.Xt_L[s]);
            }
            
            energy_time_reco = GetRecoEnergyTime(signal);
            cell.second.E_reco_noiseXT = energy_time_reco.first;
            cell.second.tau_reco_noiseXT = energy_time_reco.second;
            signal.clear();
        }
    }
}



TF1* CaloClusterGenerator::GenerateSignalSample(){
    
    CellM ->SetParameter(0, g_taud);
    CellM ->SetParameter(1, g_taupa);
    CellM ->SetParameter(2, g_td);
    CellM ->SetParameter(3, g_Rf);
    CellM ->SetParameter(4, g_C1);
    return CellM ;
}// end of function

TF1* CaloClusterGenerator::GenerateXTalkSample(){
    
    XTalk ->SetParameter(0,g_Cx);
    XTalk ->SetParameter(1,g_Rf);
    XTalk ->SetParameter(2,g_Rin);
    XTalk ->SetParameter(3,g_taud);
    XTalk ->SetParameter(4,g_taupa);
    XTalk ->SetParameter(5,g_td);
    XTalk ->SetParameter(6,g_ToNormXtC);   // This value normalize amplitude of the Xt_C to unit = 0.022206
    
    return XTalk ;
}

void  CaloClusterGenerator::InitTree(string TreeName){
    TString ts_TreeName = TreeName;
    ClusterTree = new TTree ( ts_TreeName, "Calorimeter simulation data") ;
    ClusterTree->Branch("Cell_E_truth_L2", &cell_truth_energies_l2);
    ClusterTree->Branch("Cell_E_reco_noise_L2", &cell_reco_noise_energies_l2);
    ClusterTree->Branch("Cell_E_reco_noiseXT_L2", &cell_reco_noiseXT_energies_l2);
    ClusterTree->Branch("Cell_E_fake_L2", &cell_fake_energies_l2);
    
    ClusterTree->Branch("Cell_t_reco_noise_L2", &cell_reco_noise_tau_l2);
    ClusterTree->Branch("Cell_t_reco_noiseXT_L2", &cell_reco_noiseXT_tau_l2);
    ClusterTree->Branch("Cell_t_fake_L2", &cell_fake_tau_l2);
    ClusterTree->Branch("Cell_t_delay_L2", &cell_t_delay_l2);

    ClusterTree->Branch("Cell_signal_samples_L2", &cell_signal_samples_l2);
    ClusterTree->Branch("Cell_noise_samples_L2", &cell_noise_samples_l2);
    ClusterTree->Branch("Cell_XTc_samples_L2", &cell_xtc_samples_l2);
    ClusterTree->Branch("Cell_XTl_samples_L2", &cell_xtl_samples_l2);
    ClusterTree->Branch("Energy_L2", &layer2_energy);
    ClusterTree->Branch("Impact_Energy", &impact_energy);
    ClusterTree->Branch("tau_0", &tau_0);

    ClusterTree->Branch("xtc_amplitude_l2", &cell_xtc_amplitudes_l2);
    ClusterTree->Branch("xtl_amplitude_l2", &cell_xtl_amplitudes_l2);
    
    ClusterTree->Branch("ai", &ai);
    ClusterTree->Branch("bi", &bi);


}


void CaloClusterGenerator::FillTree(){
    impact_energy = my_cluster->E_truth;
    tau_0 = my_cluster->tau_0;
    layer2_energy = my_cluster->layers[1].E_truth;

    map<pair<int,int>,ClusterCell>::iterator cell_iterator;
    for (cell_iterator = my_cluster->layers[1].layer_cells_map.begin(); cell_iterator != my_cluster->layers[1].layer_cells_map.end(); cell_iterator++){
    
        cell_truth_energies_l2.push_back(cell_iterator->second.E_truth);
        cell_xtc_amplitudes_l2.push_back(cell_iterator->second.Xt_C_amplitude);
        cell_xtl_amplitudes_l2.push_back(cell_iterator->second.Xt_L_amplitude);
        
        cell_reco_noise_energies_l2.push_back(cell_iterator->second.E_reco_noise);
        cell_reco_noiseXT_energies_l2.push_back(cell_iterator->second.E_reco_noiseXT);
        cell_fake_energies_l2.push_back(cell_iterator->second.E_reco_fake);

        cell_reco_noise_tau_l2.push_back(cell_iterator->second.tau_reco_noise);
        cell_reco_noiseXT_tau_l2.push_back(cell_iterator->second.tau_reco_noiseXT);
        cell_fake_tau_l2.push_back(cell_iterator->second.tau_reco_fake);
        
        for (int s = 0; s < n_samples; s++) {
            cell_signal_samples_l2.push_back(cell_iterator->second.samples_truth[s]);
            cell_noise_samples_l2.push_back(cell_iterator->second.sampling_noise[s]);
            cell_xtc_samples_l2.push_back(cell_iterator->second.Xt_C[s]);
            cell_xtl_samples_l2.push_back(cell_iterator->second.Xt_L[s]);
            
            cell_t_delay_l2.push_back(cell_iterator->second.sampling_delay[s]);

        }
       // cout << "eta : " << cell_iterator->second.eta_in_cluster  <<" phi : " << cell_iterator->second.phi_in_cluster  << " etruth "<< cell_iterator->second.E_truth  << endl;
      //  cout << "xtc : " << cell_iterator->second.Xt_C_amplitude  <<" xtl: " << cell_iterator->second.Xt_L_amplitude << endl;

    }
    ClusterTree->Fill();

    cell_truth_energies_l2.clear();
    cell_signal_samples_l2.clear();
    cell_noise_samples_l2.clear();
    cell_xtc_samples_l2.clear();
    cell_xtl_samples_l2.clear();
    cell_xtc_amplitudes_l2.clear();
    cell_xtl_amplitudes_l2.clear();
    
    cell_t_delay_l2.clear();
    cell_fake_tau_l2.clear();
    cell_reco_noiseXT_tau_l2.clear();
    cell_reco_noise_tau_l2.clear();
    
    cell_reco_noise_energies_l2.clear();
    cell_reco_noiseXT_energies_l2.clear();
    cell_fake_energies_l2.clear();

}

// Create the Correlation Matrix
TMatrixD CorrMatrix(TMatrixD M){
    Int_t  nCols = M.GetNcols(),                 // samples
    nRows = M.GetNrows(),                 // each row it is a vector
    k     = 1 ;                           // lag on the vector
    vector<Double_t> meanM(nRows) ;
    
    for (Int_t i = 0; i < nRows; i++){
        for (Int_t j = 0; j < nCols; j++){
            meanM[i] += M(i,j)/nCols ;
        }
        for (Int_t j = 0; j < nCols; j++){
            M(i,j) -= meanM[i] ;
        }
    }
    
    TMatrixD M1 = M ,
    S  = M1.T()*M ,
    R  = S ,
    r = M*M1 ;
    //M.Print();
    //// Matrix of the autocorrelation for a vector with samples
    if (nRows < 2){
        R.Zero();
        for (Int_t k = 0; k <= nCols-1; k++){
            for (Int_t i = 0; i < nCols - k; i++){
                R(0, k) += M(0, i)*M(0, i+k)/r(0,0) ;
            }
        }
        for (Int_t i = 0; i < Int_t(R.GetNcols()); i++){
            for (Int_t j = 0; j < Int_t(R.GetNrows()); j++){
                for(Int_t idx = 1; Int_t(idx < R.GetNcols()); idx++){
                    if (abs(i-j) == idx){
                        R(i,j) = R(0,idx) ;
                    }
                }
                if ( i == j ){
                    R(i,j) = R(0,0) ;
                }
            }
        }
        /* R.Print(); */
        //// Matrix of the correlation for a banch of the row vectors
    }else{
        //cin.get();
        for (Int_t i = 0; i < nRows; i++){
            for (Int_t j = 0; j < nRows; j++){
                R(i,j) = S(i,j) / sqrt(S(i,i)*S(j,j));
            }
        }
    }
    return R ;
}// end of function



// Calculation coefficients for the Optimal Filter to reconstruct energy and time
void FilterCoefficients(TMatrixD gSig ,
                        TMatrixD gSigT ,
                        TMatrixD DgSig ,
                        TMatrixD DgSigT ,
                        TMatrixD NoiseMatrix
                         ){
    
    //cout << "CellEner " << CellEner.size() << " EnergySamples " << EnergySamples.size() << "  gSamples  "<< gSamples.size() << "  dgSamples " << dgSamples.size() << " sNoise "<< sNoise.size() << endl;
    //cin.get() ;
    // Alvaro thesis, Marco Delmastro thesis and 1994_Signal processing... source
    /* Double_t ns = sNoise.size(); */
    Int_t n = 1, inc;
    TMatrixD aCoef(4,1),
    bCoef(4,1),
    A(1,1) ,
    Atau(1,1) ;
    TMatrixD Parameters(1,5) ;
    // Autocorrelation matrix

    TMatrixD Rxx = CorrMatrix(NoiseMatrix) ;

    //NoiseMatrix.Print() ;
    //Rxx.Print();
    //cout << "\nCorrelation Noise Matrix order " << Rxx.GetNcols() << " x " << Rxx.GetNrows() << endl;
    
    // *****************************
    // Inverse of Noise Correlation Matrix
    Rxx.Invert() ;
    // *****************************
    
    //Rxx.Print();
    //cout << "Inverse Correlation Noise Matrix"  << endl;
    //cin.get();
    
    //     InvAutCorr.Print();
    // **************************************************
    //// Parameters
    Double_t Delta,
    mu,
    lambda,
    k,
    ro;
    //cout << " here 1 " << endl;

    TMatrixD Q1 = (gSigT  * Rxx) * gSig ;              // Q1  = g^T * R^{1} - g
    TMatrixD Q2 = (DgSigT * Rxx) * DgSig ;             // Q2  = g'^T * R^{1} - g'
    TMatrixD Q3 = (DgSigT * Rxx) * gSig ;              // Q3  = g'^T * R^{1} - g
   // cout << " here 2" << endl;

    Parameters(0,0) =  Delta  =  Q1(0,0) * Q2(0,0) - Q3(0,0) * Q3(0,0) ;
    Parameters(0,1) =  mu     =  Q3(0,0)/Delta ;
    Parameters(0,2) =  lambda =  Q2(0,0)/Delta ;
    Parameters(0,3) =  ro     = -Q1(0,0)/Delta ;
    Parameters(0,4) =  k      = -Q3(0,0)/Delta ;
    //cout << " here 3" << endl;

    /*         cout << "Q1     = " << Q1(0,0) << endl ;
     cout << "Q2     = " << Q2(0,0) << endl ;
     cout << "Q3     = " << Q3(0,0) << endl ;
     cout << "Delta  = " << Delta << endl ;
     cout << "mu     = " << mu << endl ;
     cout << "lambda = " << lambda << endl ;
     cout << "ro     = " << ro << endl ;
     cout << "k      = " << k << endl ;  */
    
    //cin.get();
    //cout << "gSig  Matrix order             " << gSig.GetNcols() << " x " << gSig.GetNrows() << endl ;
    //gSig.Print() ;
    
    //cout << "DgSig Matrix order             " << DgSig.GetNcols() << " x " << DgSig.GetNrows() << endl ;
    //DgSig.Print() ;
    
    //cin.get() ;
    
    aCoef = lambda * (Rxx * gSig) +  k * (Rxx * DgSig);
    bCoef =  mu    * (Rxx * gSig) + ro * (Rxx * DgSig);
    //cout << " here 4" << endl;

    /*         cout << "a coefficients " << endl;
     aCoef.Print();
     cout << "b coefficients " << endl;
     bCoef.Print(); */
    
    A    = aCoef.T()*gSig ;
    Atau = bCoef.T()*gSig ;
    
    //cout << "from martons:: E: " <<A(0,0)  << " E*tau: " << Atau(0,0) << endl;

    //cout << "aCoef Matrix order             " << aCoef.GetNcols() << " x " << aCoef.GetNrows() << endl ;
    //aCoef.Print() ;
    //cout << "Si Matrix order                " << Si.GetNcols() << " x " << Si.GetNrows() << endl ;
    //Si.Print() ;
    /*         cout << " ================================= " << endl;
     cout << " ========= Verifying ============= " << endl;
     cout << " Amplitude of the Energy on cell = " << CellEner[cell] << endl ;
     cout << " Amplitude reconstructed     A   = " << A(0,0) << endl;
     cout << " Time of fly tau                 = " << Atau(0,0)/A(0,0) << endl ;
     cout << " ================================= " << endl; */
   // Amp = A(0,0) ;
    //tau = Atau(0,0)/A(0,0) ;
    //cin.get() ;
    
} // end function FilterCoefficients

void CaloClusterGenerator::CleanUp(){
    for (int l = 0; l < n_layers; l++ ){
        for (auto &cell : my_cluster->layers[l].layer_cells_map){
            cell.second.sampling_noise.clear();
            cell.second.samples_truth.clear();
            cell.second.dsamples_truth.clear();
            cell.second.Xt_C.clear();
            cell.second.dXt_C.clear();
            cell.second.Xt_L.clear();
            cell.second.dXt_L.clear();
            cell.second.samples_reco.clear();
            cell.second.sampling_delay.clear();
        }
    }
}

vector<CaloClusterGenerator::ClusterCell*> CaloClusterGenerator::GetListOfNeighbours(int layer, int eta, int phi, bool includeDiagonal){
    //layer number should start from 0
    layer -=1;
    
    //cout << "cell in layer " << layer << " eta: " << eta << " phi: " << phi << endl;
    
    int increments [3] = {0, -1, 1};
    pair<int,int> potential_address;
    vector<pair<int,int>> possible_neighbours;
    vector<ClusterCell*> neighbours;

    if (includeDiagonal) {
        for (int e = 0; e < 3; e++) {
            for (int p = 0; p < 3; p++) {
                if (e == 0 && p == 0) continue;
                potential_address = make_pair(eta+increments[e], phi+increments[p]);
                possible_neighbours.push_back(potential_address);
            }
        }
    }
    else {
        for (int e = 1; e < 3; e++) {
            potential_address = make_pair(eta+increments[e], phi);
            possible_neighbours.push_back(potential_address);
        }
        for (int p = 1; p < 3; p++) {
            potential_address = make_pair(eta, phi+increments[p]);
            possible_neighbours.push_back(potential_address);
        }
    }
    
    for (auto &n : possible_neighbours ){
        auto perhaps_a_cell = my_cluster->layers[layer].layer_cells_map.find(n);
        if (perhaps_a_cell!=my_cluster->layers[layer].layer_cells_map.end()) {
         //   cout << "eta of real neighbour: " << perhaps_a_cell->first.first << endl;
         //   cout << "phi of real neighbour: " << perhaps_a_cell->first.second << endl;
            neighbours.push_back(&perhaps_a_cell->second);
        }
      //  else cout << "neighbout with eta " << n.first << " and phi " <<n.second << " does not exist " << endl;
    }
    
    
    return neighbours;
}
TMatrixD CaloClusterGenerator::GetCorrelationMatrix( TVectorD vect){
    int size = vect.GetNrows();
    TMatrixD covMatrix(size, size);
    double vect_mean = vect.Sum()/size;
    double vect_std_sqared = 0;
    for (int i = 0; i < size; i ++) vect_std_sqared += abs(vect(i) - vect_mean)/size;
    
    cout << "mean: " << vect_mean << " vect_std_sqared: " << vect_std_sqared << endl;

    
    
    for (int row = 0; row < size; row++ ) {
        for (int col = 0; col < size; col++ ) {
            covMatrix[row][col] = (vect[row] - vect_mean)*(vect[col] - vect_mean)/(size*vect_std_sqared);
        }
    }
    
    //covMatrix.Print();
    return covMatrix;
}

void CaloClusterGenerator::FillXtalkAmplitudes(bool ifxtc, bool ifxtl){
    double sum_neighbours_energy;
    vector <CaloClusterGenerator::CellCluster*> neighbours;
    for (int l = 0; l < n_layers; l++){
        for (auto &cell: my_cluster->layers[l].layer_cells_map){
            //Counting XtC amplitude from immediate neighbours
            sum_neighbours_energy = 0;
            for (auto &neigh : GetListOfNeighbours(cell.second.layer, cell.second.eta_in_cluster, cell.second.phi_in_cluster, false)) {
                sum_neighbours_energy+= neigh->E_truth;
            }
            cell.second.Xt_C_amplitude =sum_neighbours_energy*g_AmpXt_C;
            
            //Counting XtL amplitude from immediate+diagonal neighbours
            sum_neighbours_energy = 0;
            for (auto &neigh : GetListOfNeighbours(cell.second.layer, cell.second.eta_in_cluster, cell.second.phi_in_cluster, true)) {
                sum_neighbours_energy+= neigh->E_truth;
            }
            cell.second.Xt_L_amplitude =sum_neighbours_energy*g_AmpXt_L;
        }
    }
    
}


Double_t RecTimeNoNoise(TVectorD gSignal, TVectorD DgSignal, TVectorD SigSamples, Double_t EnoNoise){
    Double_t time ;
    time = (gSignal(3)*gSignal(2)*gSignal(1)*SigSamples(0) - gSignal(3)*gSignal(2)*gSignal(0)*SigSamples(1))/
    ((gSignal(3)*gSignal(2)*gSignal(1)*DgSignal(0) - gSignal(3)*gSignal(2)*gSignal(0)*DgSignal(1))*EnoNoise) ;
    return time ;
}

pair<double,double> CaloClusterGenerator::OptimalFiltering( double e_truth, vector<double> gt, vector<double> dgt, vector<double> noise){
    double reco_energy = 0;
    double reco_time = 0;
    int n_signals = gt.size();
    TMatrixD noise_matrix(1,4);
    TMatrixD Rxx(4,4);
    TMatrixD m_noise_matrix(1,4), gsig(4,1), gsigt(1,4), dgsig(4,1), dgsigt(1,4);

    TVectorD noise_v(n_signals), gt_v(n_signals), dgt_v(n_signals), signal_v(n_signals), s_truth(n_signals);
    
    for (int s = 0; s < gt.size(); s++) {
   //     cout <<  noise[s] << endl;
        noise_v[s] = noise[s];
        noise_matrix[0][s] =noise[s];
        s_truth[s]  = gt[s]*e_truth;
        signal_v[s] = gt[s]*e_truth + noise[s];
        
        m_noise_matrix[0][s] =noise[s];
        gsig[s][0] =gt[s];
        gsigt[0][s] =gt[s];
        dgsig[s][0] =dgt[s];
        dgsigt[0][s] =dgt[s];
        

        gt_v[s] = gt[s];
        dgt_v[s] = dgt[s];
    }
    
   // FilterCoefficients(gsig, gsigt,  dgsig, dgsigt, m_noise_matrix);
    
   // TMatrixD noiseCovMatrix = this->GetCorrelationMatrix(noise_v);
   // noiseCovMatrix.Print();
    Rxx = CorrMatrix(noise_matrix);
    Rxx.Invert() ;

    
    double Q1 = gt_v  * (Rxx * gt_v) ;              // Q1  = g^T * R^{1} - g
    double Q2 = dgt_v * (Rxx * dgt_v) ;             // Q2  = g'^T * R^{1} - g'
    double Q3 = dgt_v * (Rxx * gt_v) ;              // Q3  = g'^T * R^{1} - g
    
    double Delta  =  Q1 * Q2 - Q3 * Q3;
    double mu     =  Q3/Delta ;
    double lambda =  Q2/Delta ;
    double ro     = -Q1/Delta ;
    double k      = -Q3/Delta ;
    
    TVectorD aCoef = lambda * (Rxx * gt_v) +  k * (Rxx * dgt_v);
    TVectorD bCoef =  mu    * (Rxx * gt_v) + ro * (Rxx * dgt_v);
    
    double E_reco = aCoef*signal_v;
    double tau = bCoef* signal_v;
    
   // aCoef.Print();
   // bCoef.Print();
    
    TVectorD ai_precalc(n_samples);
    TVectorD bi_precalc(n_samples);
    
    for (int r = 0; r < ai.size(); r++){
        ai_precalc[r] = ai[r];
        bi_precalc[r] = bi[r];
        
    }
   // ai_precalc.Print();
   // bi_precalc.Print();
   // cout << "tau = 0: " << endl;
   // cout << "E_truth: " << e_truth << " E_reco: " <<E_reco << " E*tau: " << tau/E_reco << endl;
//    cout << "Erectimenonoise: " << RecTimeNoNoise( gt_v,  dgt_v,  s_truth,  e_truth) << endl;
   // cout << "tau = 0.5: " << endl;
    //cout << "E_truth_precalc: " << e_truth << " E_reco: " <<ai_precalc*signal_v << " E*tau precalc: " << bi_precalc*signal_v/(ai_precalc*signal_v) << endl;
    //cout << "Erectimenonoise precalc: " << RecTimeNoNoise( gt_v,  dgt_v,  s_truth,  e_truth) << endl;
    
    //reco_energy = ai_precalc*signal_v;
    //reco_time = bi_precalc*signal_v/reco_energy;
    
    reco_energy = E_reco;
    reco_time = tau/E_reco;
    
    return make_pair(reco_energy, reco_time);
}

pair<vector<double>,vector<double>> CaloClusterGenerator::OptimalFilteringInit( vector<double> gt, vector<double> dgt, vector<double> noise){

    int n_signals = gt.size();
    TMatrixD noise_matrix(1,4);
    TMatrixD Rxx(4,4);
    TMatrixD m_noise_matrix(1,4), gsig(4,1), gsigt(1,4), dgsig(4,1), dgsigt(1,4);
    
    TVectorD noise_v(n_signals), gt_v(n_signals), dgt_v(n_signals), signal_v(n_signals), s_truth(n_signals);
    
    for (int s = 0; s < gt.size(); s++) {
        //     cout <<  noise[s] << endl;
        noise_v[s] = noise[s];
        noise_matrix[0][s] =noise[s];

        
        m_noise_matrix[0][s] =noise[s];
        gsig[s][0] =gt[s];
        gsigt[0][s] =gt[s];
        dgsig[s][0] =dgt[s];
        dgsigt[0][s] =dgt[s];
        
        
        gt_v[s] = gt[s];
        dgt_v[s] = dgt[s];
    }
    
    // FilterCoefficients(gsig, gsigt,  dgsig, dgsigt, m_noise_matrix);
    
    // TMatrixD noiseCovMatrix = this->GetCorrelationMatrix(noise_v);
    // noiseCovMatrix.Print();
    Rxx = CorrMatrix(noise_matrix);
    Rxx.Invert() ;
    
    
    double Q1 = gt_v  * (Rxx * gt_v) ;              // Q1  = g^T * R^{1} - g
    double Q2 = dgt_v * (Rxx * dgt_v) ;             // Q2  = g'^T * R^{1} - g'
    double Q3 = dgt_v * (Rxx * gt_v) ;              // Q3  = g'^T * R^{1} - g
    
    double Delta  =  Q1 * Q2 - Q3 * Q3;
    double mu     =  Q3/Delta ;
    double lambda =  Q2/Delta ;
    double ro     = -Q1/Delta ;
    double k      = -Q3/Delta ;
    
    TVectorD aCoef = lambda * (Rxx * gt_v) +  k * (Rxx * dgt_v);
    TVectorD bCoef =  mu    * (Rxx * gt_v) + ro * (Rxx * dgt_v);

    
    //cout << "E_truth: " << e_truth << " E_reco: " <<E_reco << " E*tau: " << tau << endl;
    //cout << "Erectimenonoise: " << RecTimeNoNoise( gt_v,  dgt_v,  s_truth,  e_truth) << endl;
    
    vector <double> ai, bi;
    
    for (int r = 0; r < aCoef.GetNrows(); r++) {
        ai.push_back(aCoef[r]);
        bi.push_back(bCoef[r]);
    }
    
    return make_pair(ai, bi);
}

pair<double,double> CaloClusterGenerator::GetRecoEnergyTime( vector<double> signal){
    
   // cout << "signal size: " << signal.size() << endl;
    
    double E_reco = 0;
    double tau_reco = 0;
 

    for (int r = 0; r < signal.size(); r++)  {
        E_reco += signal[r]*ai[r];
   //     cout << " signal: " <<  signal[r]  << " ai: " << ai[r] << " bi: " << bi[r]<<endl;
    }
    for (int r = 0; r < signal.size(); r++)  tau_reco += signal[r]*bi[r]/E_reco;
 
   // cout << " E_reco: " << E_reco << " tau: " << tau_reco << endl;
    
    return make_pair(E_reco, tau_reco);
}
