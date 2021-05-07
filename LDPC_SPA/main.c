#include"Header.h"
#include"rand.h"

int main()
{
	int n, i, j;
	int cNum, vNum, * vWeight, * cWeight, ** V, ** C,max_vWeight, max_cWeight;	//to store the Tanner graph
	double * check,*variable,**q,**r;	//to store the values and messages
	double sigma,value;
	rand_init();
	sigma = 1;

	/*---------------------read alist and build the Tanner graph----------------------*/
	{	
		FILE* f;
		f = fopen("./alist.txt", "r");
		if (!fscanf(f, "%d %d", &vNum, &cNum))
			printf("ERROR\n");
		//printf("%d %d\n", vNum, cNum);

		if(!fscanf(f, "%d", &max_vWeight))
			printf("ERROR\n"); 
		if (!fscanf(f, "%d", &max_cWeight))
			printf("ERROR\n");

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
				C[i][j]--;		//modify so that index from 0
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
				V[i][j]--;		//modify so that index from 0
			}
			q[i] = malloc(sizeof(double) * cWeight[i]);
		}
		fclose(f);
		//printf("%d\n", V[2][2]);
	}
	/*---------------------end reading the Tanner graph----------------------*/

	check = malloc(sizeof(double) * cNum);
	variable = malloc(sizeof(double) * vNum);
	/*---------------initialization step-------------------*/
	{
		for (i = 0; i < vNum; i++)
		{
			variable[i] = input(0);	//y_i
			variable[i] = 2 * variable[i] / sigma / sigma;
			//printf("%.20lf\n", variable[i]);
		}
		for (i = 0; i < cNum; i++)
			check[i] = 0;

		for (i = 0; i < vNum; i++)
			for (j = 0; j < vWeight[i]; j++)
			{
				n = searchIndex(i, C[i][j],cWeight[C[i][j]],V);
				q[C[i][j]][n] = variable[i];		//update q
			}
	}
	/*---------------end initialization step-------------------*/

	//while (1)
	{
		for (i = 0; i < vNum; i++)
		{
			/*L(r)=2atanh(product of tanh( 0.5*(L(r)) ) )*/
		}
	}

	
	//freee();
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
int searchIndex(int start,int dest,int weight,int **N)
{
	int i;
	for (i = 0; i < weight; i++)
	{
		if (N[start][i] == dest)
			return i;
	}
	return -1;	//if ERROR
}

int freee() {

}
