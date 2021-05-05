#pragma once
#pragma warning(disable:4996)
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double input(short b);
double phi(double x);
double sgn(double x);
int searchIndex(int start, int dest, int weight, int** N);
void cUpdate(int vNum, int* vWeight, int** C);
void downUpdate(int cNum, int* cWeight, int* vWeight, int** V, int** down, int** up);

