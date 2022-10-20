#include "./utility_functions.hpp"
#include <cstring>

using std::vector;

namespace QCSP {
namespace StandaloneDetector {
namespace Utilities {

/*********************************************/
/* Function Definitions *********************/
/*********************************************/

void write_score_fa(const vector<float> & score, mat_t * score_mat) {
    size_t dims[2] = {1U, score.size()};

    matvar_t * score_fa_var = Mat_VarCreate(
        "score_fa",
        MAT_C_SINGLE, MAT_T_SINGLE,
        2, dims,
        (void *) score.data(),
        0);
    if (!bool(score_fa_var)) {
        throw std::runtime_error("Failed to create a new variable score_fa");
    }

    const int result = Mat_VarWriteAppend(score_mat, score_fa_var, MAT_COMPRESSION_ZLIB, 2);
    Mat_VarFree(score_fa_var);

    if (result != 0) {
        throw std::runtime_error("Failed to append score_fa variable to the matfile.\nFailed with code " + std::to_string(result) + ".");
    }
}

void write_score_md(const vector<float> & score, mat_t * score_mat) {
    size_t dims[2] = {1U, score.size()};

    matvar_t * score_md_var = Mat_VarCreate(
        "score_md",
        MAT_C_SINGLE, MAT_T_SINGLE,
        2, dims,
        (void *) score.data(),
        0);
    if (!bool(score_md_var)) {
        throw std::runtime_error("Failed to create a new variable score_md");
    }

    const int result = Mat_VarWriteAppend(score_mat, score_md_var, MAT_COMPRESSION_ZLIB, 2);
    Mat_VarFree(score_md_var);

    if (result != 0) {
        throw std::runtime_error("Failed to append score_md variable to the matfile.\nFailed with code " + std::to_string(result) + ".");
    }
}

void write_full_score(const vector<vector<float>> & scores, mat_t * score_mat) {
    size_t dims[2] = {scores[0].size(), scores.size()};

    float * lin_scores = new float[dims[0] * dims[1]];
    for (int i = 0; i < int(scores.size()); i++) {
        memcpy(lin_scores + i * scores[0].size(), scores[i].data(), scores[0].size() * sizeof(float));
    }

    matvar_t * score_var = Mat_VarCreate(
        "score_full",
        MAT_C_SINGLE, MAT_T_SINGLE,
        2, dims,
        (void *) lin_scores,
        0);
    if (!bool(score_var)) {
        throw std::runtime_error("Failed to create a new variable score_full");
    }

    const int result = Mat_VarWriteAppend(score_mat, score_var, MAT_COMPRESSION_ZLIB, 2);
    Mat_VarFree(score_var);
    delete[] lin_scores;
    if (result != 0) {
        throw std::runtime_error("Failed to append score_full variable to the matfile.\nFailed with code " + std::to_string(result) + ".");
    }
}

void write_full_score(const vector<vector<float>> & scores, FILE * score_file) {
    for (const auto & score : scores) {
        for (const auto & value : score) {
            fprintf(score_file, "%-16f", value);
        }
        fprintf(score_file, "\n");
    }
}


void safe_write_score_fa(const std::vector<float> & score, mat_t * score_mat, omp_lock_t * lock) {
    omp_set_lock(lock);
    write_score_fa(score, score_mat);
    omp_unset_lock(lock);
}

void safe_write_score_md(const std::vector<float> & score, mat_t * score_mat, omp_lock_t * lock) {
    omp_set_lock(lock);
    write_score_md(score, score_mat);
    omp_unset_lock(lock);
}

void safe_write_full_score(const std::vector<std::vector<float>> & scores, mat_t * score_mat, omp_lock_t * lock) {
    omp_set_lock(lock);
    write_full_score(scores, score_mat);
    omp_unset_lock(lock);
}

void safe_write_full_score(const std::vector<std::vector<float>> & scores, FILE * score_file, omp_lock_t * lock) {
    omp_set_lock(lock);
    write_full_score(scores, score_file);
    omp_unset_lock(lock);
}

} // namespace Utilities
} // namespace StandaloneDetector
} // namespace QCSP
