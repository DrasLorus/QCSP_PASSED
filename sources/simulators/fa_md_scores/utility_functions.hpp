#ifndef _QCSP_PASSED_FA_MD_UTILITY_HPP_
#define _QCSP_PASSED_FA_MD_UTILITY_HPP_

#include <iostream>
#include <string>
#include <vector>

#include <matio.h>
#include <omp.h>

namespace QCSP {
namespace StandaloneDetector {
namespace Utilities {

/*********************************************/
/* Template Function Definitions *************/
/*********************************************/

template <typename T1, typename T2>
static void progressBar(T1 v, T2 max) {
    const double progress = double(v) / double(max);
    const size_t barWidth = 60;
    const size_t pos      = size_t(barWidth * progress);

    std::cout << "[";
    for (size_t i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << size_t(progress * 100.0) << " % \r";
    std::cout.flush();
}

/*********************************************/
/* Function Declarations *********************/
/*********************************************/

void write_score_fa(const std::vector<float> & score, mat_t * score_mat);

void write_score_md(const std::vector<float> & score, mat_t * score_mat);

void write_full_score(const std::vector<std::vector<float>> & scores, mat_t * score_mat);

void write_full_score(const std::vector<std::vector<float>> & scores, FILE * score_file);

void safe_write_score_fa(const std::vector<float> & score, mat_t * score_mat, omp_lock_t * lock);

void safe_write_score_md(const std::vector<float> & score, mat_t * score_mat, omp_lock_t * lock);

void safe_write_full_score(const std::vector<std::vector<float>> & scores, mat_t * score_mat, omp_lock_t * lock);

void safe_write_full_score(const std::vector<std::vector<float>> & scores, FILE * score_file, omp_lock_t * lock);



void write_score_fa(const std::vector<uint32_t> & score, mat_t * score_mat);

void write_score_md(const std::vector<uint32_t> & score, mat_t * score_mat);

void write_full_score(const std::vector<std::vector<uint32_t>> & scores, mat_t * score_mat);

void write_full_score(const std::vector<std::vector<uint32_t>> & scores, FILE * score_file);

void safe_write_score_fa(const std::vector<uint32_t> & score, mat_t * score_mat, omp_lock_t * lock);

void safe_write_score_md(const std::vector<uint32_t> & score, mat_t * score_mat, omp_lock_t * lock);

void safe_write_full_score(const std::vector<std::vector<uint32_t>> & scores, mat_t * score_mat, omp_lock_t * lock);

void safe_write_full_score(const std::vector<std::vector<uint32_t>> & scores, FILE * score_file, omp_lock_t * lock);

} // namespace Utilities
} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_FA_MD_UTILITY_HPP_
