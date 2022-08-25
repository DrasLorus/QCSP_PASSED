#ifndef INIT_HPP_INCLUDED
#define INIT_HPP_INCLUDED

/*!
 * \file init.h
 */

#include "./struct.hpp"

#include <string>

void Table_Add_GF(table_t & table, int q, int p);

void Table_Mul_GF(table_t & table, int q);

void Table_Div_GF(table_t & table, int q);

void Table_dec_GF(table_t & table, int q, int p);

void Table_Mul_DEC(table_t & table, int q);

void Table_Div_DEC(table_t & table, int q);

void LoadCode(const std::string & filename, code_t & code);

void LoadTables(table_t & table, int q, int p);

#endif
