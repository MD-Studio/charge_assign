
%module dp

%{
#define SWIG_FILE_WITH_INIT
#include "dp.h"
%}

%include "carrays.i"
%array_functions(double, doublea);
%array_functions(unsigned short, ushorta);

double solve_dp(double weights[], double profits[],
                unsigned short sets[], unsigned short num_items, unsigned short num_sets,
                unsigned short digits, short total_charge, double total_charge_diff,
                unsigned short solution[]);
