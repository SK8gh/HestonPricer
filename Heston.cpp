#include "Heston.h"

#include "random_numbers.h"
#include "random_numbers.cpp"

#include <vector>
#include <math.h>
#include <iostream>
#include <list>

//Export library
#include <fstream>

using namespace std;


// The following functions will be computing constants, which will be very useful later for code readability
// QE scheme constants

double Heston::phi_inv(double& u, double& p, double& beta) const{
    if (0 <= u <= p){
        return(0);
    }
    else {
        return( (1 / beta) * log((1 - p) / (1 - u)));
    };
}

double Heston::compute_s2(const double& vol, const double& time_step) const{

    double exp_const = (1 - exp(-kappa_m * time_step));

    double term1 = vol * pow(sigma_m,2) * exp(-kappa_m * time_step);
    term1 = term1 * exp_const / kappa_m;

    double term2 = theta_m * pow(sigma_m,2) * pow(exp_const,2) / (2 * kappa_m);
    
    return(term1 + term2);
}

double Heston::compute_m(const double& vol, const double& time_step) const{
    return(theta_m + (vol - theta_m) * exp(- kappa_m * time_step));
}

double Heston::compute_b2(const double& phi) const{
    double b2 = (2 / phi) - 1 + sqrt(2 / phi) * sqrt(2 / phi - 1);
    return(b2);
 }

double Heston::compute_a(const double& m, const double& b_squared) const{
    return(m / (1 + b_squared));
}

double Heston::compute_p(const double& phi) const{
    return((phi - 1) / (phi + 1));
}

double Heston::compute_beta(const double& p, const double&m) const{
    return((1 - p) / m);
}

double compute_call_value(const double& S_T, const double& strike){
    if (S_T > strike){
        return(S_T - strike);
    }
    else{
        return(0);
    }
}

double compute_mean(const std::vector<double>& vector_arg){
    double mean = 0;
    double size = vector_arg.size();

    for (int i=0; i < size; i++){
        mean += vector_arg[i] / size;
    }
    return(mean);
}

// Andersen stock scheme constants (vect)
std::vector<double> Heston::Andersen_constants(const double& time_step) const{
    std:vector<double> constants;

    // The scheme given by Andersen uses the gamma constants equal to 0.5 we might come back and change these values later
    double gamma1(0.5), gamma2(0.5);

    double k0 = - time_step * (rho_m * kappa_m * theta_m / sigma_m);
    
    double k1 = gamma1 * time_step * (kappa_m * rho_m / sigma_m - 0.5) - (rho_m / sigma_m);
    double k2 = gamma2 * time_step * (kappa_m * rho_m / sigma_m - 0.5) + (rho_m / sigma_m);
    
    double k3 = gamma1 * time_step * (1 - pow(rho_m,2));
    double k4 = gamma2 * time_step * (1 - pow(rho_m, 2)); 

    constants.push_back(k0);
    constants.push_back(k1);
    constants.push_back(k2);
    constants.push_back(k3);
    constants.push_back(k4);

    return(constants);
}

// Implementing the QE scheme method
std::vector<double> Heston::QE(const double& maturity, const double& vol_spot, const int& number_steps) const{
    std::vector<double> stock_vol;
    double phi_c = 1.5;
    double time_step = maturity / number_steps;

    // Initialisation
    stock_vol.push_back(vol_spot);

    // Loop
    for (int i=0; i < number_steps-1; i++){
        double old_vol = stock_vol.back();
        double m = compute_m(old_vol, time_step);
        double phi = compute_s2(old_vol, time_step) / pow(m,2);

        // QE condition
        if (phi <= phi_c){
            // Drawing a normal variable
            double Z = random_number::normal();

            // Computng constants
            double m = compute_m(old_vol, time_step);
            double b2 = compute_b2(phi);
            double a = compute_a(m, b2);

            double new_vol = a * pow((sqrt(b2) + Z),2);

            //Appending the new volatility in the vol vector
            stock_vol.push_back(new_vol);
        }
        else{
            double uniform = random_number::uniform();
            double p = compute_p(phi);
            double beta = compute_beta(p, m);

            double new_vol = phi_inv(uniform, p, beta);

            // Append the new volatility in the vol vector
            stock_vol.push_back(new_vol);
        }
    }
    // returning the volatility process
    return(stock_vol);
}

// Implementing the Andersen scheme
std::vector<double> Heston::Andersen(const double& maturity, const double& spot, const int& number_steps,
    const std::vector<double>& volatility) const{
        
        // Initialisation
        std::vector<double> stock;
        stock.push_back(spot);
        double time_step = maturity / number_steps;

        // To get the constant K0, write K[0] (same for other K constants)
        std::vector<double> K = Andersen_constants(time_step);

        // Loop
        for (int i=0; i<number_steps-1; i++){
            double old_price = stock.back();
            double normal_var = random_number::normal();

            double new_log_price = log(old_price) + mu_m * time_step + K[0] + K[1] * volatility[i] + 
            K[2] * volatility[i+1] + sqrt(K[3] * volatility[i] + K[4] * volatility[i+1]) * normal_var;

            stock.push_back(exp(new_log_price));
        }

        // Returning the stock price vector
        return(stock);
    }

// Call Pricing
double Heston::CallPricing(const double& maturity, const double& strike, const double& spot, const double& vol_spot,
    const int& number_steps, const int& number_simulations) const{
        //Initialisation
        std::vector<double> call_values;

    for (int i=0; i < number_simulations; i++){
        std::vector<double> vol_process = QE(maturity, vol_spot, number_steps);
        std::vector<double> price_process = Andersen(maturity, spot, number_steps, vol_process);

        double maturity_price = price_process.back();
        double call_value_iteration = compute_call_value(maturity_price, strike);

        call_values.push_back(exp(- mu_m * maturity) * call_value_iteration);
    }

    // Returning the average call value
    double mean_call_value = compute_mean(call_values);
    return(mean_call_value);
}

// Compute prices in the model given strikes and maturities
std::vector<std::vector<double> > Heston::ComputePrices(const std::vector<double>& maturities, 
    const std::vector<double>& strikes, const double& spot, const double& vol_spot,
     const int& number_steps, const int& number_simulations) const{
         // Initialisation
         std::vector<std::vector<double> > CallPrices;
         double length_strikes = strikes.size();
         double length_maturities = maturities.size();

        // Iteration over the strikes
         for (int i=0; i<length_strikes; i++){
            // Iteration over the maturities
            double strike = strikes[i];

            // Fixing the strike, creating a list that will contain the call values for that strike val, and every maturity
            std::vector<double> calls_fixed_strike;

            for (int j=0; j<length_maturities; j++){
                double maturity = maturities[j];
                double call_fixed_strike = CallPricing(maturity, strike, spot, vol_spot, number_steps, number_simulations);

                // Appending the call value in the calls_fixed_strike
                calls_fixed_strike.push_back(call_fixed_strike);
            }

            // Appending the vector of call values (with a single strike) in the vector of vectors
            CallPrices.push_back(calls_fixed_strike);
         }

         // Returning the call prices
         return(CallPrices);
     }



double Calibration::get_error(const double& vol_spot, const double& kappa, const double& theta, 
const double& rho, const double& sigma) const{
    // mu_m is the drift, already known = 3M EURIBOR rate
    Heston model_c(mu_m, kappa, theta, sigma, rho);

    std::vector<std::vector<double> > heston_prices = model_c.ComputePrices(maturities_m, strikes_m, spot_m,
    vol_spot, number_steps_m, number_simulations_m);
    
    double size0 = heston_prices.size();
    double size1 = heston_prices[0].size();

    // Testing shapes :
    std::string error_message = "Market prices and prices generated by the model don't have the same shape";
    if (heston_prices.size() !=  market_prices_m.size()){
        throw error_message;
    }
    else{
        if (heston_prices[0].size() !=  market_prices_m[0].size()){
            throw error_message;
        }
    }

    double error = 0;
    
    // Computing the L2 error
    for (int i=0; i<size0; i++){
        for (int j=0; j<size1; j++){
            error += pow((heston_prices[i][j] - market_prices_m[i][j]), 2);
        }
    }
    
    // Returning the average L2 error
    return(error / (size0 * size1));

}

// Printing market prices
void Calibration::print_market_prices() const{
    double size0(market_prices_m.size()), size1(market_prices_m[0].size());
    for (int i=0; i<size0; i++){
        for (int j=0; j<size1; j++){
            std::cout << market_prices_m[i][j] << ", ";
        }
        std::cout << "" << std::endl;
    }
}

std::vector<double> Calibration::call() const{
    std::vector<double> init_parameters;

    // Loop, update parameters by gradient descent, return final parameters
    return(init_parameters);
}



