#ifndef ONMI_HPP__
#define ONMI_HPP__

#include"state.hpp"
#include<vector>
long double calculate_oNMI(lvalue_input :: in< std::vector< std::vector<int64_t> > > ground_truth, lvalue_input :: in<State> st);

#endif
