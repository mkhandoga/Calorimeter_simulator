//
//  inv_functions.cpp
//  calo_sim
//
//  Created by Mykola Khandoga on 13/03/2021.
//  Copyright © 2021 Mykola Khandoga. All rights reserved.
//

#include "inv_functions.hpp"
#include <TH1.h>

// secant function
Double_t sec(Double_t x){
    return 1/cos(x);
}
// Cossecant function
Double_t csc(Double_t x){
    return 1/sin(x);
}
// Cotangent function
Double_t cotan(Double_t x){
    return 1/tan(x);
}
// Hiperbolic secant function
Double_t sech(Double_t x){
    return 1/cosh(x);
}
// Hiperbolic cossecant function
Double_t csch(Double_t x){
    return 1/sinh(x);
}
// Hiperbolic Cotangent function
Double_t cotanh(Double_t x){
    return 1/tanh(x);
}

Double_t Heaviside(Double_t x){
    if ( x < 0 ){
        return 0. ;
    }
    else{
        return 1. ;
    }
} 
