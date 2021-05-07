#pragma once
#pragma warning(disable:4996)
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double input(short);
double phi(double);
double sgn(double);
int searchIndex(int, int, int, int**);
void rUpdate(int, int, int*, int**, int**, double**, double**);
void qUpdate(int, int, int*, int**, int**, double**, double**, double*, double*);
int end_condition_check(int, int*, short*, int**);
int repeat_result_check(int, short*, short*);
void freee(int, int, int*, int*, int**, int**,short*,double*, double*, double* , double**, double**);

