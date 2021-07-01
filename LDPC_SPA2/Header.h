#pragma once
#pragma warning(disable:4996)
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void log_SPA(int, int, double, double*);
short sgn(double);
double input(short, double);
void rUpdate(int, int, double*, double**, double**, int**);
double boxsum(double, double);
int end_condition_check(int, int*, short*, int**);
int bit_error_count(int, short*);
void freee(int, int*, int*, int**, double*, double**, double**, short*);


