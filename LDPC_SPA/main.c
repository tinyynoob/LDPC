#include"rand.h"
#include"Header.h"

int main()
{
	FILE* f;
	int i;
	double *er;
	int loopNum[] = { 300000,300000 };
	int max_iteration[] = {16,16};
	double SNR_db[] = {2.05,2.05};
	rand_init();
	er = malloc(sizeof(double) * 2);	//error_rate	er[0]=BER, er[1]=FER
	f = fopen("./analysis.csv", "w");
	fprintf(f, "#loop,max_iteration,SNR(db),BER,FER\n");
	fclose(f);
	for (i = 0; i < sizeof(loopNum)/sizeof(int);i++)
	{
		f = fopen("./analysis.csv", "a");
		printf("processing SNR number %d\n", i+1);
		log_SPA(max_iteration[i], loopNum[i], SNR_db[i], er);
		fprintf(f, "%d,%d,%lf,%lf,%lf\n", loopNum[i], max_iteration[i], SNR_db[i], er[0], er[1]);	//output data
		fclose(f);
	}
	
	free(er);
	system("pause");
	return 0;
}

void log_SPA(int max_iteration, int loopNum, double SNR_db,double *er)	//analysis the error rate and store at er
{
	FILE* f;
	int n, i, j, count_iteration, count_loop;
	int cNum, vNum, * vWeight, * cWeight, ** V, ** C;	//to store the Tanner graph
	double *variable,*Q,**q,**r;	//to store the values and messages
	short* binary;		//to store the estimated codeword
	double bit_error_rate, frame_error_rate; 
	double sigma;
	
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
		C = malloc(sizeof(int*) * vNum);	//C[i] stores the numbers of c-nodes connecting to ith v-node
		V = malloc(sizeof(int*) * cNum);	//V[i] stores the numbers of v-nodes connecting to ith c-node
		r = malloc(sizeof(double*) * vNum);
		q = malloc(sizeof(double*) * cNum);
		//note that the index starts from 0

		for (i = 0; i < vNum; i++)
			if (!fscanf(f, "%d", &vWeight[i]))
				printf("ERROR\n");
		for (i = 0; i < cNum; i++)
			if (!fscanf(f, "%d", &cWeight[i]))
				printf("ERROR\n");

		for (i = 0; i < vNum; i++)	//
		{
			C[i] = malloc(sizeof(int) * vWeight[i]);
			for (j = 0; j < vWeight[i]; j++)
			{
				if (!fscanf(f, "%d", C[i] + j))
					printf("ERROR\n");
				C[i][j]--;		//modify so that index starts from 0
			}
			r[i] = malloc(sizeof(double) * vWeight[i]);
		}
		for (i = 0; i < cNum; i++)	//
		{
			V[i] = malloc(sizeof(int) * cWeight[i]);
			for (j = 0; j < cWeight[i]; j++)
			{
				if (!fscanf(f, "%d", V[i] + j))
					printf("ERROR\n");
				V[i][j]--;		//modify so that index starts from 0
			}
			q[i] = malloc(sizeof(double) * cWeight[i]);
		}
		fclose(f);
		//printf("%d\n", V[2][2]);
	}
	/*---------------------end building the Tanner graph----------------------*/
	binary = malloc(sizeof(short) * vNum);
	variable = malloc(sizeof(double) * vNum);
	Q = malloc(sizeof(double) * vNum);
	
	bit_error_rate = 0;
	frame_error_rate = 0;
	sigma = sqrt(pow(10, -SNR_db / 10));	//compute sigma from SNR_db where we assume mean_of_signal_energy = 1

	for (count_loop = 0; count_loop < loopNum; count_loop++)
	{
		//printf("%dth loop\n", count_loop);
		/*---------------initialization step-------------------*/
		{
			//printf("initial y_i values:\n");
			for (i = 0; i < vNum; i++)
			{
				variable[i] = input(0, sigma);	//y_i	default codeword:00000~
				variable[i] = 2 * variable[i] / sigma / sigma;
				//printf("%.6lf\t", variable[i]);
			}
			//printf("\n\n");

			for (i = 0; i < vNum; i++)
				for (j = 0; j < vWeight[i]; j++)
				{
					n = searchIndex(i, C[i][j], cWeight[C[i][j]], V);
					q[C[i][j]][n] = variable[i];		//update q
				}
		}
		/*---------------end initialization step-------------------*/


		/*---------------iteration step-------------------*/
		count_iteration = 0;
		while (1)
		{
			count_iteration++;
			for (i = 0; i < cNum; i++)
				rUpdate(i, cWeight[i], vWeight, V, C, q, r);
			for (i = 0; i < vNum; i++)
			{
				qUpdate(i, vWeight[i], cWeight, V, C, q, r, Q, variable);
				if (Q[i] < 0)		//set the binary value of each v-node
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
	freee(vNum, cNum, vWeight, cWeight, V, C, binary, variable, Q,  q, r);
	return;
}


double phi(double x){
	return -log(tanh(0.5*x));
}
double sgn(double x)
{
	if (!x)
		return 0;
	else if (x > 1)
		return 1;
	return -1;
}
double input(short b,double sigma)	//transform the binary bit b to modulated bit +-1 (inverse BPSK) and add noise
{
	double AWGN,y;
	if (b)	//mapping the bit
		y = -1;
	else
		y = 1;
	AWGN = sigma*gaussian();	//generate a number from X~N(0,sigma^2)
	y = y + AWGN;
	return y;
}

int searchIndex(int start,int dest,int weight,int **N)	//given start and dest, find the number of start in dest's list
{														//thus we need weight of dest and the connecting node list of the dest
	int i;
	for (i = 0; i < weight; i++)
	{
		if (N[dest][i] == start)
			return i;
	}
	return -1;	//if ERROR
}

void rUpdate(int start ,int weight, int *vWeight, int**V, int **C, double**q, double**r)
{	/*L(r)=2atanh(product of tanh( 0.5*(L(q)) ) )*/
	int i, j, n;
	double ans;
	for (i = 0; i < weight; i++)	//for each edge from start to V[start][i]
	{
		ans = 1;
		for (j = 0; j < weight; j++)
		{
			if (j == i)
				continue;
			ans = ans * tanh(0.5 * q[start][j]);
		}
		ans = 2 * atanh(ans);
		n = searchIndex(start,V[start][i],vWeight[V[start][i]],C);
		r[V[start][i]][n] = ans;
	}
}

void qUpdate(int start, int weight, int* cWeight, int** V, int** C, double**q, double**r, double *Q, double *variable)
{
	int i, n;
	double ans;
	ans = variable[start];
	for (i = 0; i < weight; i++)	//update Q[start]
		ans += r[start][i];
	Q[start] = ans;

	for (i = 0; i < weight; i++)	//for each edge from start to C[start][i]
	{
		ans = Q[start];
		n = searchIndex(start,C[start][i],cWeight[C[start][i]],V);
		ans -= r[start][i];
		q[C[start][i]][n] = ans;
	}
}

int end_condition_check(int cNum, int *cWeight, short *binary, int **V)		//if the ending codition is satisfied, return 1
{
	int i,j, sum;
	for (i = 0; i < cNum;i++)
	{
		sum = 0;
		for (j=0;j<cWeight[i];j++)
			sum = sum ^ binary[V[i][j]];	//xor operation
		if (sum)
			return 0;	//if a c-node doesn't sum to 0, return 0
	}
	return 1;	//no error is found
}

int bit_error_count(int vNum,short *binary)	//count how many bit is not 0 in a frame
{
	int i, count = 0;
	for (i = 0; i < vNum; i++)
		if (binary[i])
			count++;
	return count;
}

void freee(int vNum,int cNum,int *vWeight,int *cWeight,int **V, int**C, short*binary, double*variable, double*Q, double**q, double**r)
{
	int i;
	for (i = 0; i < vNum; i++)
	{
		free(C[i]);
		free(r[i]);
	}
	for (i = 0; i < cNum; i++)
	{
		free(V[i]);
		free(q[i]);
	}
	free(vWeight);
	free(cWeight);
	free(V);
	free(C);
	free(binary);
	free(variable);
	free(Q);
	free(q);
	free(r);
}
