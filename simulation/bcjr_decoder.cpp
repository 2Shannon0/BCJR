#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <tuple>
#include <string>

using namespace std;
namespace py = pybind11;

vector<long double> decode(
    vector<vector<tuple<string, int, string>>> &edg,
    vector<vector<tuple<string, int, string>>> &edg_bpsk,
    const vector<double> &llr_in,
    double sigma2
) {
    double a_priori = 0.5;
    double constant_coef = a_priori * (1 / (2 * M_PI * sigma2));

    vector<vector<tuple<string, long double, string>>> gammas;

    for (const auto& layer : edg_bpsk) {
        vector<tuple<string, long double, string>> new_layer;
        for (const auto& row : layer) {
            string s1 = get<0>(row);
            int int_value = get<1>(row);
            string s2 = get<2>(row);

            long double double_value = static_cast<long double>(int_value);
            new_layer.emplace_back(s1, double_value, s2);
        }
        gammas.push_back(new_layer);
    }

    vector<unordered_map<string, long double>> alphas(gammas.size() + 1);
    tuple<string, long double, string> first_edge = gammas[0][0];
    alphas[0][get<0>(first_edge)] = 1.0;

    for (size_t i = 0; i < gammas.size(); i++) {
        for (size_t j = 0; j < gammas[i].size(); j++) {
            string prev_vex = get<0>(gammas[i][j]);
            int edge_value = get<1>(edg_bpsk[i][j]);
            string next_vex = get<2>(gammas[i][j]);

            long double diff = pow(llr_in[i] - edge_value, 2) / (2 * sigma2);
            long double cur_gamma = constant_coef * exp(-diff);

            gammas[i][j] = make_tuple(prev_vex, cur_gamma, next_vex);
            long double new_alpha = cur_gamma * alphas[i][prev_vex];
            alphas[i + 1][next_vex] += new_alpha;
        }

        long double sum_alpha = 0;
        for (const auto& p : alphas[i + 1]) sum_alpha += p.second;
        if (sum_alpha != 0) {
            for (auto& p : alphas[i + 1]) p.second /= sum_alpha;
        }
    }

    vector<unordered_map<string, long double>> betas(gammas.size() + 1);
    betas[gammas.size()][get<0>(gammas[0][0])] = 1;
    vector<long double> llr_out(llr_in.size(), 0);

    for (int i = gammas.size() - 1; i >= 0; i--) {
        long double up = 0, down = 0;
        for (size_t j = 0; j < gammas[i].size(); j++) {
            string next_vex = get<0>(gammas[i][j]);
            long double cur_gamma = get<1>(gammas[i][j]);
            string prev_vex = get<2>(gammas[i][j]);

            long double new_beta = cur_gamma * betas[i + 1][prev_vex];
            betas[i][next_vex] += new_beta;

            long double cur_alpha = alphas[i][next_vex];
            long double cur_beta = betas[i + 1][prev_vex];
            long double cur_sigma = cur_gamma * cur_alpha * cur_beta;

            if (get<1>(edg_bpsk[i][j]) == 1) {
                up += cur_sigma;
            } else {
                down += cur_sigma;
            }
        }

        llr_out[i] = (down == 0) ? 9999 : log(up / down);

        long double sum_beta = 0;
        for (const auto& p : betas[i]) sum_beta += p.second;
        if (sum_beta != 0) {
            for (auto& p : betas[i]) p.second /= sum_beta;
        }
    }

    return llr_out;
}

// Оборачиваем функцию для Pybind11
PYBIND11_MODULE(bcjr_decoder, m) {
    m.def("decode", &decode, "BCJR decoder function");
}
