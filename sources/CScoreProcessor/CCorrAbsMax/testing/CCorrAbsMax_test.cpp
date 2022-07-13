#include "CScoreProcessor/CCorrAbsMax/CCorrAbsMax.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matioCpp/matioCpp.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CCorrAbsMax works for high snr inputs (q: 64, N: 60)", "[corrabsmax][high]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    matioCpp::File parm_file("../data/parameters_20210903.mat", matioCpp::FileMode::ReadOnly);

    matioCpp::Vector<double> PN64 = parm_file.read("PN64").asVector<double>();
    vector<float>            pn(PN64.size());

    std::copy(PN64.begin(), PN64.end(), pn.begin());

    if (PN64.size() != q) {
        throw "PN64 and q don't match.";
    }

    if (unsigned(parm_file.read("n_frame").asElement<double>()) != N) {
        throw "Fetched N and local N don't match.";
    }

    mat_t * data_file = Mat_Open("../data/test_data.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "iter_fcts_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "iter_fcts_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const vector<float> cabs_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CCorrAbsMax<q> * proc = new CCorrAbsMax<q>(pn.data());

    vector<float> results(cabs_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(re_in[i], im_in[i]);
    }

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(cabs_out[i], 1e-9));
    }

    delete proc;
}

TEST_CASE("CCorrAbsMax works for low snr inputs (q: 64, N: 60)", "[corrabsmax][low]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    matioCpp::File parm_file("../data/parameters_20210903.mat", matioCpp::FileMode::ReadOnly);

    matioCpp::Vector<double> PN64 = parm_file.read("PN64").asVector<double>();
    vector<float>            pn(PN64.size());

    std::copy(PN64.begin(), PN64.end(), pn.begin());

    if (PN64.size() != q) {
        throw "PN64 and q don't match.";
    }

    if (unsigned(parm_file.read("n_frame").asElement<double>()) != N) {
        throw "Fetched N and local N don't match.";
    }

    mat_t * data_file = Mat_Open("../data/test_data.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "iter_fcts_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "iter_fcts_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const vector<float> cabs_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CCorrAbsMax<q> * proc = new CCorrAbsMax<q>(pn.data());

    vector<float> results(cabs_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(re_in[i], im_in[i]);
    }

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(cabs_out[i], 1e-4f));
    }

    delete proc;
}
