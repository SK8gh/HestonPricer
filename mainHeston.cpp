#include <iostream>

#include "Heston.h"
#include "Heston.cpp"

#include <vector>
#include <list>
#include <math.h>

//Export library
#include <fstream>

using namespace std;

void save_to_csv(std::string path, std::vector<double> vect_){
    std::ofstream File(path);
    //Test
    if (File){
        std::cout << "Access granted";}
    else{std::cout << "Couldn't access the file";}

    double length_vect = vect_.size();

    for (int i=0; i < length_vect; i++){
        if (i==length_vect-1){
            File << vect_[i];
        }
        else{
            File << vect_[i] << ", ";
        }
        
    }

    File.close();
}

void write_matrix_to_csv(std::string path, std::vector<std::vector<double> > matrix){
    double size = matrix.size();

    std::ofstream File(path);
        //Test
        if (File){
            std::cout << "Access granted";}
        else{std::cout << "Couldn't access the file";}
        std::cout << "" << std::endl;

        File << "Call values" << std::endl;

    for (int i=0; i<size; i++){

        double length_vect = matrix[0].size();

        for (int j=0; j<length_vect; j++){
            File << matrix[i][j];
            File << ", ";
        }
        File << std::endl;  
    }
}


int main(){

    std::vector<double> strike_list;
    strike_list.push_back(100);
    strike_list.push_back(105);

    std::vector<double> maturities_list;
    maturities_list.push_back(1.0);
    maturities_list.push_back(2.0);

    // Parameters
    double mu(0.002), kappa(1.22), theta(0.56), sigma(0.3), rho(-0.66), spot(165);
    int number_steps(20);

    Heston model(mu, kappa, theta, sigma, rho);

    std::vector<std::vector<double> > call_prices_matrix = model.ComputePrices(maturities_list, strike_list, 100, 0.2, 
    20, 1000);

    std::cout << "Strike, maturity = 100, 1" << "\n" << "Call price : " << call_prices_matrix[0][0];
}