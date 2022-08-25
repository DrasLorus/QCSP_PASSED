#ifndef TOOLS_HPP_INCLUDED
#define TOOLS_HPP_INCLUDED

/*!
 *
 * \file tools.h
 *
 */

#include "./struct.hpp"

#include <cstdio>
#include <cstdlib>

int32_t Bin2GF(const std::vector<bool> & gf_bin, const table_t & table);

void GaussianElimination(const table_t & table, code_t & code);

int Encoding(const code_t & code, const table_t & table, const vector<int32_t> & message, vector<int32_t> & codeword);

int Syndrom(const code_t & code, const table_t & tableGF, vector<int32_t> & decide);

#endif
