/**
 * @file types.h
 * @brief This file contains various typedefs for the systematics code.
 * @author justin.mueller@colostate.edu
*/

#ifndef TYPES_H
#define TYPES_H

typedef std::tuple<Double_t, Double_t, Double_t, Double_t> index_t;
typedef std::tuple<Double_t, Double_t, Double_t> meta_t;
typedef std::map<std::string, size_t> systs_t;
typedef std::map<std::string, TH1*> weights_t;

#endif