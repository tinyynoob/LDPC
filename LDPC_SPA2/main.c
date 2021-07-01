#include"Header.h"
#include"rand.h"


int main()
{
	FILE* f;
	int i;
	double* er;
	int loopNum[] = {2000};
	int max_iteration[] = {16};
	double SNR_db[] = {1.3};
	rand_init();
	er = malloc(sizeof(double) * 2);	//error_rate	er[0]=BER, er[1]=FER
	f = fopen("./analysis.csv", "w");
	fprintf(f, "#loop,max_iteration,SNR(db),BER,FER\n");
	fclose(f);
	for (i = 0; i < sizeof(loopNum) / sizeof(int); i++)
	{
		f = fopen("./analysis.csv", "a");
		printf("processing SNR number %d\n", i + 1);
		log_SPA(max_iteration[i], loopNum[i], SNR_db[i], er);
		fprintf(f, "%d,%d,%lf,%lf,%lf\n", loopNum[i], max_iteration[i], SNR_db[i], er[0], er[1]);	//output data
		fclose(f);
	}
	free(er);
	system("pause");
	return 0;
}


void log_SPA(int max_iteration, int loopNum, double SNR_db, double* er)
{
	FILE* f;
	int i, j, n, count_iteration, count_loop;
	int vNum, cNum, * vWeight, * cWeight, ** V;		//to store the Tanner graph
	double* Q, ** r, ** prev_r;	//to store the values and messages
	short* binary;		//to store the estimated codeword
	double bit_error_rate, frame_error_rate;
	double sigma;

	rand_init();

	//printf("LDPC log-SPA\n\n");
	/*---------------------read alist and build the Tanner graph----------------------*/
	{
		f = fopen("./Gallager_3_6.txt", "r");
		if (!fscanf(f, "%d %d", &vNum, &cNum))
			printf("ERROR\n");
		//printf("%d %d\n", vNum, cNum);

		fscanf(f, "%d", &n);	//no use
		fscanf(f, "%d", &n);	//no use

		vWeight = malloc(sizeof(int) * vNum);	//vWeight[i] stores the number of edges from ith v-node 
		cWeight = malloc(sizeof(int) * cNum);	//cWeight[i] stores the number of edges from ith c-node 
		V = malloc(sizeof(int*) * cNum);	//V[i] stores the numbers of v-nodes connecting to ith c-node
		r = malloc(sizeof(double*) * cNum);
		prev_r = malloc(sizeof(double*) * cNum);
		//note that the index starts from 0

		for (i = 0; i < vNum; i++)
			if (!fscanf(f, "%d", &vWeight[i]))
				printf("ERROR\n");
		for (i = 0; i < cNum; i++)
			if (!fscanf(f, "%d", &cWeight[i]))
				printf("ERROR\n");

		for (i = 0; i < vNum; i++)	//discard the information of edges from V
			for (j = 0; j < vWeight[i]; j++)
				fscanf(f, "%d", &n);

		for (i = 0; i < cNum; i++)	//
		{
			V[i] = malloc(sizeof(int) * cWeight[i]);
			for (j = 0; j < cWeight[i]; j++)
			{
				if (!fscanf(f, "%d", V[i] + j))
					printf("ERROR\n");
				V[i][j]--;		//modify so that index starts from 0
			}
			r[i] = malloc(sizeof(double) * cWeight[i]);
			prev_r[i] = malloc(sizeof(double) * cWeight[i]);
		}
		fclose(f);
		//printf("%d\n", V[2][2]);
	}
	/*---------------------end building the Tanner graph----------------------*/
	Q = malloc(sizeof(double) * vNum);
	binary = malloc(sizeof(short) * vNum);

	bit_error_rate = 0;
	frame_error_rate = 0;
	sigma = sqrt(pow(10, -SNR_db / 10));	//compute sigma from SNR_db where we assume mean_of_signal_energy = 1

	for (count_loop = 0; count_loop < loopNum; count_loop++)
	{
		//printf("%dth loop\n", count_loop);
		/*---------------initialization step-------------------*/
		{
			for (i = 0; i < vNum; i++)
			{
				Q[i] = input(0, sigma);	//y_i	default codeword:00000~
				Q[i] = 2 * Q[i] / sigma / sigma;
				//printf("%.6lf\t", Q[i]);
			}
			printf("\n\n");
			for (i = 0; i < cNum; i++)
				for (j = 0; j < cWeight[i]; j++)
				{
					r[i][j] = 0;
					prev_r[i][j] = 0;
				}
		}
		/*---------------end initialization step-------------------*/

		/*----------------iteration step-------------------*/
		count_iteration = 0;
		while (1)
		{
			count_iteration++;
			for (i = 0; i < cNum; i++)
				rUpdate(i, cWeight[i], Q, r, prev_r, V);
			for (i = 0; i < cNum; i++)
				for (j = 0; j < cWeight[i]; j++)
					Q[V[i][j]] += r[i][j] - prev_r[i][j];
			for (i = 0; i < vNum; i++)
			{
				if (Q[i] < 0)
					binary[i] = 1;
				else
					binary[i] = 0;
			}
			//printf(" %dth iteration:\n", count_iteration);
			/*
			for (i = 0; i < vNum; i++)
				printf("%.6lf\t", Q[i]);
			printf("\n");
			for (i = 0; i < vNum; i++)
				printf("%d\t\t", binary[i]);
			printf("\n");
			*/

			/*--------------check the algorithm ending condition--------------*/
			if (end_condition_check(cNum, cWeight, binary, V))	//check if cH^{T}==0
				break;

			if (count_iteration == max_iteration)	//check if the maximum limit is achieved
				break;

			for (i = 0; i < cNum; i++)
				for (j = 0; j < cWeight[i]; j++)
					prev_r[i][j] = r[i][j];
		}
		/*---------------end iteration step-------------------*/

		//printf("--------------algorithm ends--------------\n\n");
		/*
		printf("the estimated codeword:\n");
		for (i = 0; i < vNum; i++)
			printf("%d  ", binary[i]);
		printf("\n\n");
		*/
		n = bit_error_count(vNum, binary);
		bit_error_rate += n;
		if (n)		//if n is not 0, then the estimated codeword is wrong
			frame_error_rate++;
	}
	bit_error_rate = bit_error_rate / loopNum / vNum;
	frame_error_rate = frame_error_rate / loopNum;
	er[0] = bit_error_rate;
	er[1] = frame_error_rate;

	//free the pointers and end
	freee(cNum, vWeight, cWeight, V, Q, r, prev_r, binary);
	return;
}


short sgn(double x)
{
	if (x > 0)
		return 1;
	else if (x)
		return -1;
	return 0;
}
double input(short b, double sigma)	//transform the binary bit b to modulated bit +-1 (inverse BPSK) and add noise
{
	double AWGN, y;
	if (b)	//mapping the bit
		y = -1;
	else
		y = 1;
	AWGN = sigma * gaussian();	//generate a number from X~N(0,sigma^2)
	y = y + AWGN;
	return y;
}

void rUpdate(int start, int weight, double* Q, double** r, double** prev_r, int** V)	// update r from the specified c-node
{
	int i;
	double** bs;
	if (weight == 1)	//consider again
		return;
	else if (weight == 2)
	{
		r[start][0] = Q[V[start][1]] - prev_r[start][1];
		r[start][1] = Q[V[start][0]] - prev_r[start][0];
		return;
	}
	bs = malloc(sizeof(double*) * 2);
	bs[0] = malloc(sizeof(double) * (weight - 1));	//boxsum --->
	bs[1] = malloc(sizeof(double) * (weight - 1));	//boxsum <---
	bs[0][0] = Q[V[start][0]] - prev_r[start][0];
	bs[1][0] = Q[V[start][weight - 1]] - prev_r[start][weight - 1];
	for (i = 1; i < weight - 1; i++)
		bs[0][i] = boxsum(bs[0][i - 1], Q[V[start][i]] - prev_r[start][i]);
	for (i = 1; i < weight - 1; i++)
		bs[1][i] = boxsum(bs[1][i - 1], Q[V[start][weight - 1 - i]] - prev_r[start][weight - 1 - i]);

	r[start][0] = bs[1][weight - 2];
	r[start][weight - 1] = bs[0][weight - 2];
	for (i = 1; i < weight - 1; i++)
		r[start][i] = boxsum(bs[0][i - 1], bs[1][weight - i - 2]);
	free(bs[0]);
	free(bs[1]);
	free(bs);
}

double boxsum(double x, double y)	//boxsum(x,y) = sgn(x)sgn(y)min(|x|,|y|)+log(1+e^{-|x+y|})-log(1+e^{-|x-y|})
{
	double a, b, c;
	a = fabs(x);
	if (fabs(y) < a)
		a = fabs(y);
	a *= sgn(x);
	a *= sgn(y);
	b = exp(-fabs(x + y));
	b = log(1 + b);
	c = exp(-fabs(x - y));
	c = log(1 + c);
	return a + b - c;
}

int end_condition_check(int cNum, int* cWeight, short* binary, int** V)		//if the ending codition is satisfied, return 1
{
	int i, j, sum;
	for (i = 0; i < cNum; i++)
	{
		sum = 0;
		for (j = 0; j < cWeight[i]; j++)
			sum = sum ^ binary[V[i][j]];	//xor operation
		if (sum)
			return 0;	//if a c-node doesn't sum to 0, return 0
	}
	return 1;	//no error is found
}

int bit_error_count(int vNum, short* binary)	//count how many bit is not 0 in a frame
{
	int i, count = 0;
	for (i = 0; i < vNum; i++)
		if (binary[i])
			count++;
	return count;
}

void freee(int cNum, int* vWeight, int* cWeight, int** V, double* Q, double** r, double** prev_r, short* binary)
{
	int i;
	for (i = 0; i < cNum; i++)
	{
		free(V[i]);
		free(r[i]);
		free(prev_r[i]);
	}
	free(vWeight);
	free(cWeight);
	free(V);
	free(Q);
	free(r);
	free(prev_r);
	free(binary);
}
