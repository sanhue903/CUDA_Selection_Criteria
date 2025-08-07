#ifndef CRITERIA
#define CRITERIA

#include "sketch/sketch.h"
#include <cmath>

float sigma(int p)
{
	switch(p) {
		case 4:
			return 1.106 / sqrt(1 << p);
		case 5:
			return 1.07 / sqrt(1 << p);
		case 6:
			return 1.054 / sqrt(1 << p);
		case 7:
			return 1.046  / sqrt(1 << p);
	}
	return 1.039 / sqrt(1 << p);
}

double cota_n (size_t card_A, size_t card_B, double t_hat, int p, float Z=2.0, int order_n=1){
	// card_A <= card_B
	double gamma = (double)card_A / card_B;
	float sigma_p = sigma(p);
	double S = 0;
	double num = 1;
	for(int k=1; k<order_n+1; k++){
		num *= Z*sigma_p;
		S += num;
	}
	double minimo = std::min(1.0, (1.0+Z*sigma_p)*card_B / t_hat);
	return minimo*(1+gamma)*S;
}

double kota_mas(size_t card_A, size_t card_B, double t_hat, int p, float Z=2.0){
	// card_A <= card_B
	double gamma = (double)card_A / card_B;
	float sigma_p = sigma(p);
	double t_hat_mas = t_hat / (1.0 + Z*sigma_p);
	double K_mas = ((1.0+gamma)*card_B-t_hat_mas) / t_hat_mas;
	return K_mas;
}

inline bool CB(double tau, double card_A, double card_B){
	// card_A <= card_B
	double gamma = (double)card_A / card_B;
	return (gamma >= tau);
}


bool hll_an(double tau, size_t card_A, size_t card_B, std::shared_ptr<sketch::hll_t> S_A, std::shared_ptr<sketch::hll_t> S_B, int p, float Z=2.0, int order_n=1){
	// card_A <= card_B
	double t_hat = S_A->union_size(*S_B);
	double J_hat = (double)(card_A+card_B-t_hat)/t_hat;
	double C = cota_n (card_A, card_B, t_hat, p, Z, order_n);
	return ((J_hat + C) >= tau);
}

bool hll_a(double tau, size_t card_A, size_t card_B, std::shared_ptr<sketch::hll_t> S_A, std::shared_ptr<sketch::hll_t> S_B, int p, float Z=2.0){
	size_t t_hat =  S_A->union_size(*S_B);
	double K_mas = kota_mas (card_A, card_B, t_hat, p, Z);
	return (K_mas >= tau);
}

bool smh_a(std::vector<uint64_t> v1, std::vector<uint64_t> v2, uint n_rows, uint n_bands){
	if(n_rows*n_bands!=v1.size()){
		std::cerr << "ERROR: Number of bands and rows doesnt match the MinHash sketch size." << "\n";
		return 0;
	}
	for(uint band_id=0;band_id<n_bands;band_id++){
		bool check = std::equal(
				v1.begin() + band_id*n_rows, v1.begin() + (band_id+1)*n_rows,
				v2.begin() + band_id*n_rows, v2.begin() + (band_id+1)*n_rows
				);
		if (check){
			return 1;
		}
	}
	return 0;
}

bool CB_hll_an(double tau, size_t card_A, size_t card_B, std::shared_ptr<sketch::hll_t> S_A, std::shared_ptr<sketch::hll_t> S_B, int p, float Z=2.0, int order_n=1){
	if (!CB(tau, card_A, card_B)) return 0;
	return hll_an(tau, card_A, card_B, S_A, S_B, p, Z, order_n);
}

bool CB_hll_a(double tau, size_t card_A, size_t card_B, std::shared_ptr<sketch::hll_t> S_A, std::shared_ptr<sketch::hll_t> S_B, int p, float Z=2.0){
	if (!CB(tau,card_A,card_B)) return 0;
	return hll_a(tau, card_A, card_B, S_A, S_B, p, Z);
}

bool CB_smh_a(double tau, size_t card_A, size_t card_B, std::vector<uint64_t> v1, std::vector<uint64_t> v2, int n_rows, int n_bands){
	if (!CB(tau, card_A,card_B)) return 0;
	return smh_a(v1,v2,n_rows,n_bands);
}


#endif
