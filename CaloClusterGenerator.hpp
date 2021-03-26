//
//  CaloClusterGenerator.hpp
//  calo_sim
//
//  Created by Mykola Khandoga on 12/03/2021.
//  Copyright © 2021 Mykola Khandoga. All rights reserved.
//

#ifndef CaloClusterGenerator_hpp
#define CaloClusterGenerator_hpp

#include <iostream>

#include <stdio.h>
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TF3.h"
#include "TF1.h"
#include <math.h>
#include "TMatrixD.h"

using namespace TMath ;
using namespace std;

class CaloClusterGenerator
{

    public:

    pair<vector<double>,vector<double>> OptimalFilteringInit(vector<double> gt, vector<double> dgt, vector<double> noise);
    pair<double,double> GetRecoEnergyTime( vector<double> signal);

    void FillXtalkAmplitudes(bool ifxtc, bool ifxtl);

    CaloClusterGenerator();
    TTree *ClusterTree;

    struct ClusterCell {
        int layer;
        
        int eta_in_cluster;
        int phi_in_cluster;
        
        double eta_min;
        double eta_max;
        double phi_min;
        double phi_max;
        double r_min;
        double r_max;
        
        std::vector<double> sampling_noise;
        std::vector<double> samples_truth;
        std::vector<double> dsamples_truth;
        std::vector<double> samples_truthXtC;
        std::vector<double> dsamples_truthXtC;
        std::vector<double> sampling_delay;

        std::vector<double> Xt_C;
        std::vector<double> dXt_C;
        std::vector<double> Xt_L;
        std::vector<double> dXt_L;

        std::vector<double> samples_reco;
        
        double Xt_C_amplitude;
        double Xt_L_amplitude;
        
        double E_truth;
        
        double E_reco_noise;
        double tau_reco_noise;

        double E_reco_noiseXT;
        double tau_reco_noiseXT;

        
        double E_reco_fake;
        double tau_reco_fake;
    };

    
    int const static n_layers = 3;
    struct ClusterLayer {
        vector<ClusterCell*> layer_cells;
        map<pair<int,int>,ClusterCell> layer_cells_map;

        double eta_cell_size_factor;
        double phi_cell_size_factor;

        double eta_cell_size;
        double phi_cell_size;
        
        int layer_number;
        int n_cells_eta;
        int n_cells_phi;
        
        double eta_min;
        double eta_max;
        double phi_min;
        double phi_max;
        double r_min;
        double r_max;
        
        double E_reco;
        double E_truth;
    };
    
    struct CellCluster {
        ClusterLayer layers[n_layers];
        
        double eta_min;
        double eta_max;
        double phi_min;
        double phi_max;
        double r_min;
        double r_max;
        
        double E_impact;
        double tau_0;
        double E_reco;
        double E_truth;
    };
    void InitOpFiltCoefficients();
    
    void FillClusterEnergy(double e_impact);
    
    Double_t IntegEnerCell( Double_t EtaMin ,Double_t EtaMax , Double_t PhiMin, Double_t PhiMax, Double_t Rmin, Double_t Rmax );
    
    void InitCluster(double cluster_center_eta, double cluster_center_phi, int eta_size = 5, int phi_size = 5);

    void InitSideLayers(ClusterLayer &central_layer, ClusterLayer &side_layer );

    void InitLayerCells(ClusterLayer &layer );

    void FillSignalSamples(double tau_0);

    void InitTree(string TreeName);
    
    void FillRecoEnergyTime();
    
    void InitFunctions();
    
    void  FillTree();
    TMatrixD GetCorrelationMatrix( TVectorD vect);
    
    vector<ClusterCell*> GetListOfNeighbours(int layer, int eta, int phi, bool includeDiagonal);
    
    void CleanUp();
    
    TF1* GenerateSignalSample();
    TF1* GenerateXTalkSample();
    CellCluster* my_cluster = new CellCluster;
    
    vector<double> ai;
    vector<double> bi;
    
    private:
    int n_samples = 4;

    TF3* IntegraCell;
    TF1* CellM;
    TF1* XTalk;
    
    
    double cluster_center_eta;
    double cluster_center_phi;
    
    double impact_eta;
    double impact_phi;
    
    TRandom3 *RandomGen;
    
    Double_t g_DeltaEtaS2 = 0.025 ;       // Delta Eta on S2
    Double_t g_DeltaPhiS2 = M_PI/128 ;    // Delta Phi on S2
    Int_t    g_nEtaS2     = 5 ;           // Cluster on layer 2 (Eta x Phi) - rows
    Int_t    g_nPhiS2     = 5 ;           // Cluster on layer 2 (Eta x Phi) - columns
    //
    // Impact Point
    Double_t g_EtaImpact = 0.025  ;
    Double_t g_PhiImpact = M_PI/128 ;
    // *********************************
    
    Double_t g_gainAdj     = 1    ;         ///ConstNormCluster[g_RowClust-1][g_ColClust-1] ;      // adjust for normalization of
    Double_t g_E0          = 5e4  ;         // Impact Energy (in MeV)
    Double_t NoiseAmp      =  50  ;         // amplitude of the noise on MeV
    Double_t g_ADCpedestal = 900 ;
    Double_t g_theta       = 2*atan(exp(-g_EtaImpact)) ;
    Double_t g_R0          = 1500 ; //1385.0  ;      // distance between the beam axes and the first layer of calorimeter. Start point of EM Calorimeter (mm)
    // /* *************************************** */
    // Lenghts of the sampling layers
    // g_Rs1 + g_Rs2 + g_Rs3 = 470
    Double_t g_Rs1         =   90 ;         // length of the S1 on mm
    Double_t g_Rs2         =  337 ;         // length of the S2 on mm
    Double_t g_Rs3         =   43 ;         // length of the S3 on mm
    // /* *************************************** */
    
    Double_t g_tmax2   = 1000.   ;
    Double_t g_dt      =    0.1  ;          // increment of time to generate graphical, in ns
    Double_t g_window  =  600.0  ;          // size of window in ns
    Double_t g_nPoints = g_window/g_dt ;    // number of points for a window of 600 ns
    Double_t g_tSamp   =   25.0  ;          // sample time in ns
    Double_t g_nSamp   =    4.0  ;          // number of samples
    
    vector<Double_t> g_tDelayCell ;
    
    Double_t g_ToNormXtC = 0.022206 ;       // This value normalize amplitude of the Xt_C to unit
    Double_t g_ToNormXtL = 0.0539463;       // This value normalize amplitude of the Xt_L to unit
    Double_t g_ToNormNoise = 4.4787   ;
    //Double_t g_AmpXt_C   = 7.0/100  ;       // XTalk amplitude on % values of the Energy
    Double_t g_AmpXt_C   = 4.0/100  ;       // XTalk amplitude on % values of the Energy
    Double_t g_AmpXt_L   = 2.3/100  ;       // XTalk amplitude on % values of the Energy
    //Double_t g_AmpNoise  = NoiseAmp/g_E0 ;  // Noise amplitude
    Double_t g_AmpNoise  = 1./100 ;  // Noise amplitude
    
    // /* *************************************** */
    /// parameters for the cell and Xtalk signals
    Double_t g_taud    =   15.82 ;
    Double_t g_taupa   =   17.31 ;
    Double_t g_td      =  420.00 ;
    Double_t g_Rf      =    0.078 ;
    Double_t g_C1      =   50.00 ;
    Double_t g_Rin     =    1.20 ;
    Double_t g_Cx      =   47.00 ;
    // /* *************************************** */
    
    Double_t g_Eth[89];                     // Values of the theorectical Energy
    Double_t g_Tauth[89];                   // Values of the theoretical Tau
    
    // /* *************************************** */
    // Data for calculating Moliere Radius for a sampling calorimeter
    // http://pdg.lbl.gov/2019/AtomicNuclearProperties/index.html
    Double_t g_da      =   4.00 ;           // thickness of the active media (Ar)
    Double_t g_dp      =   2.00 ;           // thickness of the passive media (Pb)
    
    Double_t g_RmLAr   =  90.43 ;           // Moliere Radius in mm
    Double_t g_Z_LAr   =  18.00 ;           // Atomic number
    Double_t g_X0_LAr  = 140.00 ;           // Radiation length in mm
    Double_t g_wLAr    =   0.36 ;           // weight for the LAr
    Double_t g_EcLAr   =  32.84 ;           // critical energy for e- on LAr in MeV, 31.91 MeV (for e+)
    
    Double_t g_RmLead  =  16.02 ;           // Moliere Radius in mm
    Double_t g_Z_Lead  =  82.00 ;           // Atomic number
    Double_t g_X0_Lead =   5.612 ;          // Radiation length in mm
    Double_t g_wLead   =   1 - g_wLAr;      // weight for the Lead
    Double_t g_EcLead  =   7.43 ;           // critical energy for the Lead in MeV (e-), 7.16 for e+
    
    Double_t g_Es      =  21.2  ;           // multiple scattering energy in MeV
    
    Double_t g_a       =   4.36  ;          // Adjusted accord to information of electron deposition on layers
    Double_t g_b       =   0.25 ;           // Same as a
    
    //Double_t g_a       =   4.00 ;
    //Double_t g_b       =   0.26 ;
    
    //LAr resolution terms
    Double_t g_SampTerm  = 10./100 ;         // 10% of the Energy
    Double_t g_ConstTerm = 0.7/100 ;         // Constant term equal to 0.7%
    Double_t g_NoiseTerm = 0.40 ;            // Energy on GeV. Noise term is equal to 400 MeV
    
    // /* *************************************** */
    // Effective values for a sampling calorimeter. Reference: ATL-COM-PHYS-2004-015
    
    // Moliere Radius for a sampling calorimeter
    
    Double_t g_Rmoleff = 1/(1/g_Es*(g_wLAr*g_EcLAr/g_X0_LAr + g_wLead*g_EcLead/g_X0_Lead)) ;
    Double_t g_X0eff   = 1/(g_wLAr/g_X0_LAr + g_wLead/g_X0_Lead) ;
    Double_t g_Eceff   = g_X0eff * ( (g_wLAr*g_EcLAr)/g_X0_LAr + (g_wLead*g_EcLead)/g_X0_Lead) ;
    Double_t g_Zeff    = g_wLAr*g_Z_LAr + g_wLead*g_Z_Lead ;
    
    Double_t e = 1/(1 + 0.007*(g_Z_Lead - g_Z_LAr));
    Double_t T  = (g_a - 1)/g_b ;
    
    /// Grindhammer parameters to fit shower
    
    /// Homogeneus media
    Double_t z1 = 0.0251 + 0.00319*TMath::Log(g_E0) ;
    Double_t z2 = 0.1162 - 0.000381*g_Z_LAr ;
    Double_t k1 = 0.659 - 0.00309*g_Z_LAr ;
    Double_t k2 = 0.645 ;
    Double_t k3 = -2.59 ;
    Double_t k4 = 0.3585 + 0.0421*TMath::Log(g_E0) ;
    Double_t p1 = 2.632 - 0.00094*g_Z_LAr ;
    Double_t p2 = 0.401 + 0.00187*g_Z_LAr ;
    Double_t p3 = 1.313 - 0.0686*TMath::Log(g_E0) ;
    Double_t y  = g_E0/g_Eceff ;
    Double_t t1hom = -0.59 ;
    Double_t t1sam = -0.59 ;
    Double_t Thom  = t1hom + TMath::Log(y) ;
    Double_t t2 = -0.53 ;
    
    /// Sampling media
    Double_t Fs = g_X0eff/(g_da + g_dp) ;
    Double_t Tsamp = (1 - e)*t2 + t1sam/Fs + Thom ;

    vector <double> cell_truth_energies_l2;
    vector <double> cell_reco_noise_energies_l2;
    vector <double> cell_reco_noiseXT_energies_l2;
    vector <double> cell_fake_energies_l2;

    vector <double> cell_reco_noise_tau_l2;
    vector <double> cell_reco_noiseXT_tau_l2;
    vector <double> cell_fake_tau_l2;
    vector <double> cell_t_delay_l2;
    
    vector <double> cell_signal_samples_l2;
    vector <double> cell_noise_samples_l2;
    vector <double> cell_xtc_samples_l2;
    vector <double> cell_xtl_samples_l2;
    vector <double> cell_xtc_amplitudes_l2;
    vector <double> cell_xtl_amplitudes_l2;
    
    double layer2_energy;
    double impact_energy;
    double tau_0;


};




#endif /* CaloClusterGenerator_hpp */
