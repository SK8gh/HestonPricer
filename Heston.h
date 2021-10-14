#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

// Call value
double compute_call_value(const double& S_T, const double& strike);

// Arithmetic mean of a vector
double compute_mean(const double& vector_arg);


class Heston{
public:
    // Constructor
    Heston(const double& mu_arg, const double& kappa_arg, const double& theta_arg, const double& sigma_arg, const double& rho_arg)
     :mu_m(mu_arg), kappa_m(kappa_arg), theta_m(theta_arg), sigma_m(sigma_arg), rho_m(rho_arg){}

    // QE scheme
    std::vector<double> QE(const double& maturity, const double& vol_spot, const int& number_steps) const;

    // Andersen, stock price scheme together with the quadradic exponential scheme
    std::vector<double> Andersen(const double& maturity, const double& spot, const int& number_steps,
    const std::vector<double>& volatility) const;
    
    // Call pricing
    double CallPricing(const double& maturity, const double& strike, const double& spot, const double& vol_spot,
    const int& number_steps, const int& number_simulations) const;

    // Compute prices
    std::vector<std::vector<double> > ComputePrices(const std::vector<double>& maturities, 
    const std::vector<double>& strikes, const double& spot, const double& vol_spot,
     const int& number_steps, const int& number_simulations) const;

private:
    // Parameters of the Heston model
    const double mu_m;
    const double kappa_m; 
    const double theta_m;
    const double sigma_m;
    const double rho_m;


    // The following functions will be computing constants, which will be very useful later for code readability
    double phi_inv(double& u, double& p, double& beta) const;

    double compute_s2(const double& vol, const double& time_step) const;

    double compute_m(const double& vol, const double& time_step) const;

    double compute_b2(const double& phi) const;

    double compute_a(const double& m, const double& b) const;

    double compute_p(const double& phi) const;

    double compute_beta(const double& p, const double&m) const;

    // Andersen scheme constants
    std::vector<double> Andersen_constants(const double& time_step) const;
};


class Calibration{
public:
    // Constructor
    Calibration(const std::vector<std::vector<double> >& market_prices_arg, std::vector<double>& strikes_arg,
    std::vector<double>& maturities_arg, const double& spot_arg, const double& mu_arg, const int& number_steps_arg,
    const int& number_simulations_arg)
    :market_prices_m(market_prices_arg), strikes_m(strikes_arg), maturities_m(maturities_arg),
    spot_m(spot_arg), mu_m(mu_arg), number_steps_m(number_steps_arg),
    number_simulations_m(number_simulations_arg){}

    double get_error(const double& vol_spot, const double& kappa, const double& theta, const double& rho, 
    const double& sigma) const;

    // Get private market prices
    void print_market_prices() const;

    std::vector<double> call() const;
    
private:
    // Arguments for the calibration algorithm : market prices, strikes, maturities, spot, vol_spot, nb_steps, nb_simulations
    const std::vector<std::vector<double> > market_prices_m;
    const std::vector<double> strikes_m;
    const std::vector<double> maturities_m;

    const double spot_m; 
    const double mu_m; 
    const int number_steps_m; 
    const int number_simulations_m;
};





