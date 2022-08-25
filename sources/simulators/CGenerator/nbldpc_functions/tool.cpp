/*!
 * \file tools.cpp
 * \brief tools for NB decoder
 * \author C. Monière, C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015, modified the 25/08/2022
 */

#include <algorithm>
#include <cmath>
#include <cstring>

#include "./init.hpp"
#include "./struct.hpp"
#include "./tools.hpp"

/*!
 * \fn Bin2GF
 * \brief compute a GF(q) symbol corresponding to a frame of log_2(GF) bits
 * Parameters    :
 * Inputs        :
 * 	- vector<bool> gf_bin : array representing p bits
 * 	- int p               : size of the array gf_bin. p = log2 (q)
 * 	- int q               : order of the field
 * 	- table_t table       : mapping table
 * Outputs       :
 *      - index of the non-binary symbol
 */

int32_t Bin2GF(const vector<bool> & gf_bin, const table_t & table) {
    const auto iter = std::find(table.BINGF.begin(), table.BINGF.end(), gf_bin);
    return int32_t(iter - table.BINGF.begin());
}

/**
 * \fn GaussianElimination
 * \brief Perform a Gaussian elimination on the parity-check matrix.
 * 		   The procedure stops when the
 * 		   The redundancy symbols are initialized to 0.
 * Inputs
 * 	- code->mat   : parity-check matrix
 * 	- table       : lookup tables for computations in GF(q)
 * Outputs       :
 *      - code.encoding_mat : Upper triangular matrix used for encoding
 *      - code->Perm : Column permutation
 */
void GaussianElimination(const table_t & table, code_t & code) {
    const int N = code.N;
    const int M = code.M;

    code.encoding_mat.resize(M);
    for (int32_t m = 0; m < M; m++) {
        code.encoding_mat[m].resize(N);
    }

    code.perm.resize(N);
    for (int32_t n = 0; n < N; n++)
        code.perm[n] = n;

    for (int32_t m = 0; m < M; m++) {
        for (int32_t k = 0; k < code.rowDegrees[m]; k++) {
            code.encoding_mat[m][code.node_matrix[m][k]] = code.coefficients[m][k];
        }
    }

    for (int32_t m = 0; m < M; m++) {
        int32_t ind = m;
        while ((ind < N) && (code.encoding_mat[m][ind] == 0)) {
            ind++;
        }

        if (ind == N) {
            printf("The matrix is not full rank (%d,%d)\n", m, N);
            exit(EXIT_FAILURE);
        }

        {
            const int32_t buf = code.perm[ind];
            code.perm[ind]    = code.perm[m];
            code.perm[m]      = buf;
        }

        for (int32_t m1 = 0; m1 < M; m1++) {
            const int32_t buf          = code.encoding_mat[m1][m];
            code.encoding_mat[m1][m]   = code.encoding_mat[m1][ind];
            code.encoding_mat[m1][ind] = buf;
        }

        for (int32_t m1 = m + 1; m1 < M; m1++) {
            if (code.encoding_mat[m1][m] != 0) {
                const int32_t buf = code.encoding_mat[m1][m];
                for (int32_t n = m; n < N; n++) {
                    if (code.encoding_mat[m1][n] != 0)
                        code.encoding_mat[m1][n] = table.DIVGF[code.encoding_mat[m1][n]][buf];
                }
                for (int32_t n = m; n < N; n++) {
                    if (code.encoding_mat[m1][n] != 0)
                        code.encoding_mat[m1][n] = table.MULGF[code.encoding_mat[m1][n]][code.encoding_mat[m][m]];
                }
                for (int32_t n = m; n < N; n++) {
                    vector<bool> temp(code.p);
                    for (int32_t i = 0; i < code.p; i++) {
                        temp[i] = (table.BINGF[code.encoding_mat[m1][n]][i]) ^ (table.BINGF[code.encoding_mat[m][n]][i]);
                    }
                    code.encoding_mat[m1][n] = Bin2GF(temp, table);
                }
            }
        }
    }
}

/**
 * \fn Encoding
 * \brief Encode the information bits into a codeword.
 * 		   matUT beeing upper triangular, the backsubstitution method is used.
 * 		   The M first symbols in NSYMB are redundancy symbols (before deinterleaving)
 * Inputs
 * 	- KSYMB  ( KSYMB are information symbols)
 * Outputs
 *      - Codeword
 *      - NBIN : binary copy of the codeword
 */

int Encoding(const code_t & code, const table_t & table, const vector<int32_t> & message, vector<int32_t> & codeword) {
    const int N = code.N;
    const int M = code.M;
    const int p = code.p;

    vector<int> NSYMB(N);

    for (int32_t k = 0; k < N - M; k++) {
        NSYMB[M + k] = message[k];
        // printf(" %d ",KSYMB[k]);
    }
    // getchar();

    /* Backsubstitution */
    for (int32_t m = M - 1; m >= 0; m--) {
        int32_t buf = 0;
        for (int32_t n = m + 1; n < N; n++) {
            if (code.encoding_mat[m][n] != 0) {
                vector<bool> temp(p);
                for (int32_t i = 0; i < p; i++) {
                    temp[i] = (table.BINGF[buf][i]) ^ (table.BINGF[table.MULGF[code.encoding_mat[m][n]][NSYMB[n]]][i]);
                }
                buf = Bin2GF(temp, table);
                //  buf = table->ADDGF[buf][table->MULGF[code.encoding_mat[m][n]][NSYMB[n]]];
            }
        }
        /* Systematic codeword (interleaved) */
        NSYMB[m] = table.DIVGF[buf][code.encoding_mat[m][m]];
    }

    /* De-interleaving */
    for (int32_t n = 0; n < N; n++) {
        codeword[code.perm[n]] = NSYMB[n];
    }

    return 0;
}

/**
 * \fn Syndrom
 * \brief Compute the syndom of a message
 * Inputs
 * 	- code structure code_t
 * 	- table_t tableGF : lookup table
 * 	- message
 * Outputs
 * 	- synd is 0 iff. the decided message is a codeword
 * 	(the value of synd is not meaningful if synd != 0 )
 */
int Syndrom(const code_t & code, const table_t & tableGF, vector<int32_t> & decide) {
    int synd = 0;
    for (int32_t k = 0; k < code.M; k++) {
        for (int32_t l = 0; l < code.rowDegrees[k]; l++) {
            vector<bool> temp(code.p, 0);
            for (int32_t i = 0; i < code.p; i++) {
                temp[i] = (tableGF.BINGF[synd][i]) ^ (tableGF.BINGF[tableGF.MULGF[code.coefficients[k][l]][decide[code.node_matrix[k][l]]]][i]);
            }
            synd = Bin2GF(temp, tableGF);
            // synd = tableGF->ADDGF[synd][tableGF->MULGF[code->matValue[k][l]][decide[code->mat[k][l]]]];
        }

        if (synd != 0) {
            break;
        }
    }

    return synd;
}
