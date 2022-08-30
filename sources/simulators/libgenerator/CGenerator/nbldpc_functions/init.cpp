/**
 * \file init.cpp
 * \brief initialization of Non-binary LDPC decoder
 * \author C. Monière, C.Marchand, A. Al-Ghouwahel, Oussama Abassi, L. Conde-Canencia, A. abdmoulah, E. Boutillon
 * \copyright BSD copyright
 * \date 03/03/2015, modified the 25/08/2022
 *
 */

// initialization
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// #include "./struct.h"
#include "./init.hpp"
#include "./tools.hpp"

#define exit_if(condition)               \
    if (condition) {                     \
        perror("Failed in system call"); \
        exit(EXIT_FAILURE);              \
    }

/*!
 * \fn Table_Add_GF
 * \brief Compute the addition table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Add_GF(table_t & table, int32_t q, int32_t p) {
    for (int32_t j = 0; j < q; j++) {
        for (int32_t k = 0; k < q; k++) {
            vector<bool> temp(p);
            for (int32_t i = 0; i < p; i++) {
                temp[i] = (table.BINGF[j][i]) ^ (table.BINGF[k][i]);
            }
            table.ADDGF[j][k] = Bin2GF(temp, table);
        }
    }
}

// multiply GF values and output decimal
void Table_Mul_DEC(table_t & table, int32_t q) {
    for (int32_t i = 0; i < q; i++) {
        for (int32_t j = 0; j < q; j++) {
            table.MULDEC[i][j] = table.DECGF[table.MULGF[i][j]];
        }
    }
}

// divide dicimal by GF and output GF
void Table_Div_DEC(table_t & table, int q) {
    for (int32_t i = 0; i < q; i++) {
        for (int32_t j = 0; j < q; j++) {
            // table->DIVDEC[table->DECGF[i]][table->DECGF[j]]=table->DECGF[table->DIVGF[i][j]];
            table.DIVDEC[table.DECGF[i]][j] = table.DIVGF[i][j];
        }
    }
}

void Table_dec_GF(table_t & table, int q, int p) {

    // bin2dec
    for (int32_t j = 0; j < q; j++) {
        int32_t sum = 0;
        for (int32_t i = 0; i < p; i++) {
            const int32_t tmp = int32_t(table.BINGF[j][i]);
            // printf("%d",tmp);
            sum = sum + (tmp << i);
        }
        table.DECGF[j]   = sum;
        table.GFDEC[sum] = j;
        // printf(" \n bin2dec of GF %d is %d \n",j,sum);
    }
    // getchar();
}

/*!
 * \fn Table_Mul_GF
 * \brief compute the multiplication table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Mul_GF(table_t & table, int q) {
    for (int32_t i = 0; i < q; i++) {
        for (int32_t j = 0; j < q; j++) {
            if (i == 0 || j == 0)
                table.MULGF[i][j] = 0;
            else if (i == 1)
                table.MULGF[i][j] = j;
            else if (j == 1)
                table.MULGF[i][j] = i;
            else {
                const int32_t temp = i + j - 2;
                if (temp < q - 1)
                    table.MULGF[i][j] = temp + 1;
                else
                    table.MULGF[i][j] = (temp % (q - 1)) + 1;
            }
        }
    }
}

/*!
 * \fn Table_Div_GF
 * \brief compute the division table in GF(q)
 * Parameters    :
 * Inputs        :
 * 	- table  : structure containing an allocated pointer to the addition table
 * 	- int logGF : logGF = log2 (GF)
 * 	- int GF    : order of the field
 * Outputs       :
 */
void Table_Div_GF(table_t & table, int q) {
    int nb = q - 1;
    for (int32_t i = 0; i < q; i++) {
        for (int32_t j = 0; j < q; j++) {
            if (j == 0) {
                table.DIVGF[i][j] = 0;
            } else if (i == 0) {
                table.DIVGF[i][j] = 0;
            } else if (j == 1) {
                table.DIVGF[i][j] = i;
            } else {
                table.DIVGF[i][j] = nb--;
            };
            if (nb < 1) {
                nb = q - 1;
            }
        }
    }
}

/**
 * \fn LoadCode
 * \brief Open a parity-check matrix file. Allocate and initialize the code structures
 * Parameters    :
 * Inputs        :
 * 	- FileMatrix  : Name of the parity-check matrix file
 * 	- code_t code : Structure that describes the code
 * Outputs       :
 */
void LoadCode(const std::string & alist_file, code_t & code) {
    /*
     * Load the files corresponding to code (graph, size, GF)
     */
    FILE * f = fopen(alist_file.c_str(), "r");

    int32_t N, M, q;

    exit_if(fscanf(f, "%d", &N) != 1);
    exit_if(fscanf(f, "%d", &M) != 1);
    exit_if(fscanf(f, "%d", &q) != 1);
    code.N    = N;
    code.M    = M;
    code.q    = q;
    code.p    = (int32_t) log2f(float(q));
    code.K    = N - M;
    code.rate = float(N - M) / float(code.N);

    code.columnDegrees.resize(N);
    for (int32_t n = 0; n < N; n++) {
        exit_if(fscanf(f, "%d", code.columnDegrees.data() + n) != 1);
    }

    code.rowDegrees.resize(M);
    for (int32_t m = 0; m < M; m++) {
        exit_if(fscanf(f, "%d", code.rowDegrees.data() + m) != 1);
    }

    code.node_matrix.resize(M);
    code.coefficients.resize(M);
    for (int32_t m = 0; m < M; m++) {
        code.node_matrix[m].resize(code.rowDegrees[m]);
        code.coefficients[m].resize(code.rowDegrees[m]);
    }

    for (int32_t m = 0; m < M; m++) {
        for (int32_t k = 0; k < code.rowDegrees[m]; k++) {
            int32_t temp_node, temp_coef;
            exit_if(fscanf(f, "%d %d", &temp_node, &temp_coef) != 2);
            code.node_matrix[m][k]  = temp_node - 1;
            code.coefficients[m][k] = temp_coef + 1;
        }
    }

    fclose(f);

    code.nbBranch = 0;
    for (int32_t m = 0; m < M; m++) {
        code.nbBranch += code.rowDegrees[m];
    }

    // printf("LDPC code parameters: \n");
    // printf(" \t N \t:%d \n \t K \t:%d \n \t M\t:%d \n \t CR\t:%g \n \t GF \t:%d \n \t logGF \t:%d\n", N, N - M, M, code.rate, q, code.p);
    // fflush(stdout);
}

/**
 * \fn void LoadTables (table_t *table, int GF, int logGF)
 * \brief Memory allocation for the tables and Initialization of the tables.
 * Inputs        :
 * 	- table_t table : Structure that contains pointers to the tables
 * 	- int GF    : order of the field
 * 	- int logGF : logGF = log2(GF)
 * Outputs       :
 */
void LoadTables(table_t & table, int q, int p) {

    if (q != 16 && q != 32 && q != 64 && q != 256 && q != 4096) {
        printf("The binary image of GF(%d) is not available in this version of the program. Please try GF(64) or GF(256)\n", q);
        exit(EXIT_FAILURE);
    }

    int32_t nbRow, nbCol;

    /* BINGF [GF][log2GF] */
    nbRow       = q;
    nbCol       = p;
    table.BINGF = vector<vector<bool>>(nbRow, vector<bool>(nbCol));

    /* ADDGF [GF][GF] */
    nbRow       = q;
    nbCol       = q;
    table.ADDGF = vector<vector<int32_t>>(nbRow, vector<int32_t>(nbCol));

    /*DECGF [GF] */
    table.DECGF = vector<int32_t>(q);
    /*GFDEC [GF] */
    table.GFDEC = vector<int32_t>(q);

    /* MULGF [GF][GF] */
    table.MULGF = vector<vector<int32_t>>(nbRow, vector<int32_t>(nbCol));

    /* DIVGF [GF][GF] */
    table.DIVGF = vector<vector<int32_t>>(nbRow, vector<int32_t>(nbCol));

    /* MULDEC [GF][GF] */
    table.MULDEC = vector<vector<int32_t>>(nbRow, vector<int32_t>(nbCol));

    /* DIVDEC [GF][GF] */
    table.DIVDEC = vector<vector<int32_t>>(nbRow, vector<int32_t>(nbCol));

    if (q == 16) {
        for (int32_t g = 0; g < q; g++)
            for (int32_t l = 0; l < p; l++)
                table.BINGF[g][l] = BINGF_16[g][l];
        // printf("Loading of the binary image of GF(64): Success\n");
        // fflush(stdout);
    }

    if (q == 32) {
        for (int32_t g = 0; g < q; g++)
            for (int32_t l = 0; l < p; l++)
                table.BINGF[g][l] = BINGF_32[g][l];
        // printf("Loading of the binary image of GF(64): Success\n");
        // fflush(stdout);
    }

    if (q == 64) {
        for (int32_t g = 0; g < q; g++)
            for (int32_t l = 0; l < p; l++)
                table.BINGF[g][l] = BINGF_64[g][l];
        // printf("Loading of the binary image of GF(64): Success\n");
        // fflush(stdout);
    }

    if (q == 256) {
        for (int32_t g = 0; g < q; g++)
            for (int32_t l = 0; l < p; l++)
                table.BINGF[g][l] = BINGF_256[g][l];
        // printf("Loading of the binary image of GF(256): Success\n");
        // fflush(stdout);
    }

    if (q == 4096) {
        for (int32_t g = 0; g < q; g++)
            for (int32_t l = 0; l < p; l++)
                table.BINGF[g][l] = BINGF_4096[g][l];
        // printf("Loading of the binary image of GF(256): Success\n");
        // fflush(stdout);
    }

    /*
     * Build the addition, multiplication and division tables (corresponding to GF[q])
     */
    // Table_Add_GF(table,GF,logGF);
    Table_dec_GF(table, q, p);
    Table_Mul_GF(table, q);
    Table_Div_GF(table, q);
    Table_Mul_DEC(table, q);
    Table_Div_DEC(table, q);
}
