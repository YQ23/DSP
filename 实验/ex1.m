clear;clc;
L=4;N=8;                                    %方波参数
x1=[ones(1,L),zeros(1,N-L)];                %方波x(n)
n = 0:1:N-1;
x2 = cos(n.*pi./4);
x3 = sin(n.*pi./8);

X1 = dft(x1,N);
X2 = dft(x2,N);
X3 = dft(x3,N);
set(0,'defaultfigurecolor','w')%使面板颜色为白色
figure(1);
stem(n,abs(X1));
title('x1[n]=R4[n]的DFT(N=8)');
figure(2);
stem(n,x1);
title('x1[n]=R4[n]矩形序列');


figure(3);
stem(n,abs(X2));
title('x2[n]=cos(0.25*pi*n)的DFT(N=8)');
figure(4);
stem(n,x2);
title('x2[n]=cos(0.25*pi*n)余弦序列');

figure(5);
stem(n,abs(X3));
title('x3[n]=sin(0.125*pi*n)的DFT(N=8)');
figure(6);
stem(n,x3);
title('x3[n]=sin(0.125*pi*n)正弦序列');

n2 = 0:1:15;N=16;
x1=[ones(1,L),zeros(1,N-L)];                %方波x(n)
x2 = cos(n2.*pi./4);
x3 = sin(n2.*pi./8);
X1 = dft(x1,N);
X2 = dft(x2,N);
X3 = dft(x3,N);
figure(7);
stem(n2,abs(X1));
title('x1[n]=R4[n]的DFT(N=16)');
figure(8);
stem(n2,x1);
title('x1[n]=R4[n]矩形序列');

figure(9);
stem(n2,abs(X2));
title('x2[n]=cos(0.25*pi*n)的DFT(N=16)');
figure(10);
stem(n2,x2);
title('x2[n]=cos(0.25*pi*n)余弦序列');

figure(11);
stem(n2,abs(X3));
title('x3[n]=sin(0.125*pi*n)的DFT(N=16)');
figure(12);
stem(n2,x3);
title('x3[n]=sin(0.125*pi*n)正弦序列');

fs = 64;%采样频率
T = 1./fs;
N = 16;
n = 0:1:N-1;
t = n.*T;

x4 = sin(8*pi.*t)+cos(16*pi.*t)+cos(20*pi.*t);%合成信号
X4 = dft(x4,N);%DFT变换
figure(20);
stem(n,abs(X4));
title('x4(t)=sin(8πt)+cos(16πt)+cos(20πt)的DFT的幅度谱');
figure(21);
%stem(n,x4);
plot(n,x4);
title('x4(t)=sin(8πt)+cos(16πt)+cos(20πt)合成序列');


N = 32;
n = 0:1:N-1;
t = n.*T;

x4 = sin(8*pi.*t)+cos(16*pi.*t)+cos(20*pi.*t);
X4 = dft(x4,N);
figure(22);
stem(n,abs(X4));
title('x4(t)=sin(8πt)+cos(16πt)+cos(20πt)的DFT的幅度谱');
figure(23);
%stem(n,x4);
plot(n,x4);
title('x4(t)=sin(8πt)+cos(16πt)+cos(20πt)合成序列');


N = 64;
n = 0:1:N-1;
t = n.*T;

x4 = sin(8*pi.*t)+cos(16*pi.*t)+cos(20*pi.*t);
X4 = dft(x4,N);
figure(24);
stem(n,abs(X4));
title('x4(t)=sin(8πt)+cos(16πt)+cos(20πt)的DFT的幅度谱');
figure(25);
%stem(n,x4);
plot(n,x4);
title('x4(t)=sin(8πt)+cos(16πt)+cos(20πt)合成序列');


N = 8;
n = 0:1:N-1;
x2 = cos(n.*pi./4);
x3 = sin(n.*pi./8);
x5 = x2+1i.*x3;
X5 = fft(x5);%FFT变换
figure(51);
stem(n,abs(X5));
title('X(k)(N=8)');

%利用共轭对称性进行变换
X5 = [X5,X5(1)];
X5_r = fliplr(X5);%反转

X5_1 = 0.5.*(X5+conj(X5_r));
X5_2 = 0.5.*(X5-conj(X5_r));
X5_1 = X5_1(1:N);
X5_2 = X5_2(1:N);
figure(30);
stem(n,abs(X5_1));
title('共轭对称性得到的余弦序列的谱(N=8)');
figure(31);
stem(n,abs(X5_2));
title('共轭对称性得到的正弦序列的谱(N=8)');

figure(40);
X51 = fft(x2);
stem(n,abs(X51));
title('余弦序列直接DFT(N=8)');

figure(41);
X52 = fft(x3);
stem(n,abs(X52));
title('正弦序列直接DFT(N=8)');

N = 16;
n = 0:1:N-1;
x2 = cos(n.*pi./4);
x3 = sin(n.*pi./8);
x5 = x2+1i.*x3;
X5 = fft(x5);
figure(52);
stem(n,abs(X5));
title('X(k)(N=16)');

X5 = [X5,X5(1)];
X5_r = fliplr(X5);

X5_1 = 0.5.*(X5+conj(X5_r));
X5_2 = 0.5.*(X5-conj(X5_r));
X5_1 = X5_1(1:N);
X5_2 = X5_2(1:N);

figure(32);

stem(n,abs(X5_1));
title('共轭对称性得到的余弦序列的谱(N=16)');
figure(33);

stem(n,abs(X5_2));
title('共轭对称性得到的正弦序列的谱(N=16)');


figure(42);
X51 = fft(x2);
stem(n,abs(X51));
title('余弦序列直接DFT(N=16)');

figure(43);
X52 = fft(x3);
stem(n,abs(X52));
title('正弦序列直接DFT(N=16)');

function Xk=dft(xn,N)
n=0:1:N-1;
k=0:1:N-1;
WN=exp(-1i*2*pi/N);
nk=n'*k;
WNnk=WN.^(nk);
Xk=xn*WNnk;
end

function [xn]=idft(Xk,N)
n=[0:1:N-1];
k=[0:1:N-1];
WN=exp(-j*2*pi/N);
nk=n'*k;
WNnk=WN.^(-nk);
end
