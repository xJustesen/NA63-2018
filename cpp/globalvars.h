#ifndef GLOBALVARS_H
#define GLOBALVARS_H

default_random_engine generator;
normal_distribution<double> xdistro;
normal_distribution<double> ydistro;
uniform_int_distribution<int> hits_per_event;
normal_distribution<double> x_entry_angles;
normal_distribution<double> y_entry_angles;
uniform_real_distribution<double> R;

#endif
