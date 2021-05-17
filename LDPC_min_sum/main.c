#include"Header.h"
#include"rand.h"

int main()
{
	int max_iteration = 16;
	double SNR_db=0.8;
	FILE* f;
	int i, j, n, count_iteration;
	int vNum, cNum, * vWeight, * cWeight, ** V;		//to store the Tanner graph
	double* Q, ** r, ** prev_r;	//to store the values and messages
	short* binary;		//to store the estimated codeword
	double sigma;
	rand_init();
	
	printf("LDPC min-sum\n\n");
	/*---------------------read alist and build the Tanner graph----------------------*/
	{
		f = fopen("./alist.txt", "r");
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

	sigma = sqrt(pow(10, -SNR_db / 10));	//compute sigma from SNR_db where we assume mean_of_signal_energy = 1

	/*---------------initialization step-------------------*/
	{
		for (i = 0; i < vNum; i++)
		{
			Q[i] = input(0, sigma);	//y_i	default codeword:00000~
			printf("%.6lf\t", Q[i]);
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
		printf(" %dth iteration:\n", count_iteration);
		
		for (i = 0; i < vNum; i++)
			printf("%.6lf\t", Q[i]);
		printf("\n");
		for (i = 0; i < vNum; i++)
			printf("%d\t\t", binary[i]);
		printf("\n");
		

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

	printf("--------------algorithm ends--------------\n\n");
	
	printf("the estimated codeword:\n");
	for (i = 0; i < vNum; i++)
		printf("%d  ", binary[i]);
	printf("\n\n");
	

	//free the pointers and end
	freee(cNum, vWeight, cWeight, V, Q, r, prev_r, binary);
	system("pause");
	return 0;
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
	short flag, sign = 1;	//the flag is used to indicate first loop
	int i, j;
	double min, q;

	for (i = 0; i < weight; i++)
	{
		flag = 1;
		for (j = 0; j < weight; j++)
		{
			if (j == i)
				continue;
			q = Q[V[start][i]] - prev_r[i][j];
			sign *= sgn(q);
			if (flag = 1)		//first
			{
				min = abs(q);
				flag = 0;
			}
			else if (abs(q) < min)
				min = abs(q);
		}
		r[start][i] = (double)sign * min;
	}
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
