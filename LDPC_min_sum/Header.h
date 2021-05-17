#pragma once
#pragma warning(disable:4996)
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

short sgn(double x);
double input(short, double);
void rUpdate(int, int, double*, double**, double**, int**);
int end_condition_check(int, int*, short*, int**);
void freee(int, int*, int*, int**, double*, double**, double**, short*);


