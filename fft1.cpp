//C语言版FFT算法编写
#include<stdio.h>
#include<math.h>
#define pi 3.1415926535

int *reverse(int x[]);
float cal_real(int N,int n);
float cal_imag(int N,int n);

int main()
{
	int xn[8] = {0,1,0,1,1,1,0,0};
//	int xn[8]={0,1,2,3,4,5,6,7};
	int N,M,i;
	N = 8;
	M = 3;
	printf("初始的x[n]:");
	for(i=0;i<N;i++)
	{
		printf("%d",xn[i]);
	}
	printf("\n");

	int *x = reverse(xn);
	printf("倒位序之后得到的x[n]:");
    for(i=0;i<N;i++)
	{
		printf("%d",xn[i]);
	}
	printf("\n");
	
	float X_r[8],X_i[8]={0},T_r,T_i;//因为C语言本身没有复数运算，所以这里将其实部与虚部分开保存
	for(i=0;i<N;i++)
	{
		X_r[i] = xn[i];
	}


	int m,B,IndexJ,P,k;
	for(m=1;m<=M;m++)
	{
		B = pow(2,m-1);
		for(IndexJ=0;IndexJ<B;IndexJ++)
		{
			P = pow(2,M-m) * IndexJ;
			for(k=IndexJ;k<N;k=k+pow(2,m))
			{
				T_r = X_r[k+B] * cal_real(N,P) - X_i[k+B] * cal_imag(N,P);
				T_i = X_r[k+B] * cal_imag(N,P) + X_i[k+B] * cal_real(N,P);
				X_r[k+B] = X_r[k] - T_r;
				X_i[k+B] = X_i[k] - T_i;
				X_r[k] = X_r[k] + T_r;
				X_i[k] = X_i[k] + T_i;
			}
		}
	}

printf("最终的X[k]为:\n");
for(i=0;i<N;i++)
{
	printf("%f %fi\n",X_r[i],X_i[i]);
}

	return 0;
}

//计算WNp
float cal_real(int N,int n)
{
	float r;
	r = cos(-2*pi*n/N);
	return r;

}
float cal_imag(int N,int n)
{
	float i;
	i = sin(-2*pi*n/N);
	return i;
}

//倒位序算法
int *reverse(int x[])
{
	int N,LH,N1,IndexJ;
	N = 2*sizeof(x);
//	printf("%d",N);
	LH = N / 2;
	IndexJ = LH;
	N1 = N - 2;
	int i,T,k;
	for(i=1;i<=N1;i++)
	{
		if(i<IndexJ)
		{
			T = x[i];
			x[i] = x[IndexJ];
			x[IndexJ] = T;
		}
		k = LH;
		while(IndexJ>=k)
		{
			IndexJ = IndexJ - k;
			k = k/2;
		}
		IndexJ = IndexJ + k;
	}
	return x;
}