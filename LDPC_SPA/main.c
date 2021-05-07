#include"rand.h"
#include"Header.h"


int main()
{
	FILE* f;
	int n, i, j,count_iteration,count_loop;
	int cNum, vNum, * vWeight, * cWeight, ** V, ** C;	//to store the Tanner graph
	double *variable,*Q,**q,**r;	//to store the values and messages
	short* binary,*last_binary;		//to store the estimated codewords
	double sigma;
	rand_init();
	sigma = 1;

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
	last_binary = malloc(sizeof(short) * vNum);
	variable = malloc(sizeof(double) * vNum);
	Q = malloc(sizeof(double) * vNum);
	f = fopen("./analysis.csv", "w");
	for (count_loop = 0; count_loop < 10000; count_loop++)
	{
		/*---------------initialization step-------------------*/
		{
			printf("LDPC log-SPA\n\n initial y_i values:\n");
			for (i = 0; i < vNum; i++)
			{
				variable[i] = input(0);	//y_i	default codeword: 00000~
				variable[i] = 2 * variable[i] / sigma / sigma;
				printf("%.6lf\t", variable[i]);
				last_binary[i] = 5;	// arbitrarily set but cant be 0 or 1
			}
			printf("\n\n");

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

			if (repeat_result_check(vNum, last_binary, binary))	//check if the result repeats
				break;

			for (i = 0; i < vNum; i++)		//update last_binary
				last_binary[i] = binary[i];
		}
		/*---------------end iteration step-------------------*/

		printf("--------------algorithm ends--------------\n\n");
		printf("the estimated codeword:\n");
		for (i = 0; i < vNum; i++)
			printf("%d  ", binary[i]);
		printf("\n\n");

		for (i = 0; i < vNum; i++)
			fprintf(f,"%d,",binary[i]);
		fprintf(f,"\n");
	}
	fclose(f);
	//free the pointers and end
	freee(vNum, cNum, vWeight, cWeight, V, C, binary, last_binary, variable, Q,  q, r);
	system("pause");
	return 0;
}

double input(short b)	//transform the binary bit b to modulated bit +-1 (BPSK) and add noise
{
	double AWGN,y;
	double mean_of_signal_energy = 1;  //in the case inverse BPSK (1 -> -1, 0 -> 1)
	double SNR = 0.5;   //SNR = mean_of_signal_energy / mean_of_WGN_energy ,
	double mean_of_AWGN_energy = mean_of_signal_energy / SNR;   // mean_of_signal_energy / mean_of_AWGN_energy = SNR   
	double AWGN_SNR;
	if (b)	//mapping the bit
		y = -1;
	else
		y = 1;
	AWGN = gaussian();	//get a number from X~N(0,1)
	AWGN_SNR = AWGN * sqrt(mean_of_AWGN_energy);
	y = y + AWGN_SNR;
	
	return y;
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

void rUpdate(int start ,int weight, int *vWeight,int**V,int **C,double**q,double**r)
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

void qUpdate(int start, int weight, int* cWeight, int** V, int** C, double**q,double**r,double *Q,double *variable)
{
	int i, n;
	double ans;
	ans = variable[start];
	for (i = 0; i < weight; i++)	//update each Q[i]
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

int end_condition_check(int cNum,int *cWeight,short *binary,int **V)		//if the ending codition is satisfied, return 1
{
	int i,j, sum;
	for (i = 0; i < cNum;i++)
	{
		sum = 0;
		for (j=0;j<cWeight[i];j++)
		{
			sum = sum ^ binary[V[i][j]];	//xor operation
			if (sum)
				return 0;	//if a c-node doesn't sum to 0, return 0
		}
	}
	return 1;	//no error is found
}

int repeat_result_check(int vNum,short* last_binary,short*binary)
{
	int i;
	for (i = 0; i < vNum; i++)	//if we have identical estimated codeword as the last iteration, then end the algorithm
		if (last_binary[i] - binary[i] != 0)
			return 0;
	return 1;
}

void freee(int vNum,int cNum,int *vWeight,int *cWeight,int **V,int**C,short*binary, short* last_binary,double*variable,double*Q,double**q,double**r)
{
	int i, j;
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
	free(last_binary);
	free(q);
	free(r);
}
