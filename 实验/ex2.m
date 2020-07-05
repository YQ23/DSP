%清除无关变量
clear;
close all;
%数字滤波器指标：
wp=0.2*pi;ws=0.3*pi;Rp=1;As=15;
%转换成为模拟指标：
Fs=1;T=1/Fs;
Omgp=(2/T)*tan(wp/2);
Omgs=(2/T)*tan(ws/2);
%如果是脉冲响应不变法,则Omgp=wc*Fs;Omgs=ws*FS;
%模拟原型滤波器计算：
[n,Omgc]=buttord(Omgp,Omgs,Rp,As,'s');%计算阶数和截止频率
%双线性变换法计算数字滤波器系数
[B,A] = butter(n,Omgc,'s');
[bd,ad] = bilinear(B,A,Fs);
tf(bd,ad,1)%输出传递函数
%求数字系统的频率特性:
[H1,w1]=freqz(bd,ad);
%心电信号：
xn = [-4,-2,0,-4,-6,-4,-2,-4,-6,-6,-4,-4,-6,-6,-2,6,12,8,0,-16,-38,-60,-84,-90,-66,...
  -32,-4,2,-4,8,12,12,10,6,6,6,4,0,0,0,0,0,-2,-4,0,0,0,-2,-2,0,0,-2,-2,-2,-2,0,...
  -2,-4,-2,0,-2,-4,-4,2,0,0,-2,-4,-2,0,0,-2,-4,-2,0,0,-4,-4,-2,-2,-4,-6,-6,-4,-4,8,-10,-8,...
  -6,-6,-8,-12,-10,-8,-8,-10,-12,-10,-8,-8,-10,-10,-8,-6,-6,-8,-8,-4,-2,-4,-4,-4,...
  0,0,-2,-4,-2,-2,0,-4];

yn=filter(bd,ad,xn);%进行滤波处理

%窗函数法设计FIR滤波器
wp=0.2*pi;ws=0.3*pi;
deltaw=ws-wp;%过渡带宽
wc=(ws+wp)/2;

N0 = ceil(1.8*pi/deltaw);%使用矩形窗
N = N0+mod(N0+1,2);%确保N为奇数
windows = (boxcar(N))';

% N0 = ceil(6.1*pi/deltaw);%使用三角形窗
% N = N0+mod(N0+1,2);%确保N为奇数
% windows = (bartlett(N))';

% N0 = ceil(6.6*pi/deltaw);%使用海明窗
% N = N0+mod(N0+1,2);%确保N为奇数
% windows = (hamming(N))';
% 
% N0 = ceil(6.2*pi/deltaw);%使用汉宁窗
% N = N0+mod(N0+1,2);%确保N为奇数
% windows = (hanning(N))';

% N0 = ceil(11*pi/deltaw);%使用布莱克曼窗
% N = N0+mod(N0+1,2);%确保N为奇数
% windows = (blackman(N))';

% %使用凯撒窗
% N = N0+mod(N0+1,2);%确保N为奇数
% windows = (kaiser(N))';

n = 0:ceil(N-1);
m=n-(N-1)/2+eps;%+eps转换成浮点数
hd = sin(wc*m)./(pi*m);
%hh = fir1(N-1,wc/pi,'low',boxcar(N));
%加窗施加在hn上
b=hd.*windows;

N = 13;
fc = 0.23;
b = fir1(N-1,fc,rectwin(N));
[H2,w2]=freqz(b,1);
hb = tf(b,1,1)
%y = filter(b,1,xn);%窗函数滤波处理
y = filter(b,1,xn);%窗函数滤波处理
%yy = fftfilt(hh,xn,length(xn));
set(0,'defaultfigurecolor','w')%使面板颜色为白色
figure(10);
xt = 0:(length(xn)-1);
plot(xt,xn,'blue');
hold on;
plot(xt,yn,'green');
plot(xt,y,'red');
title('不同方法下的时域波形');
legend('原始信号','双线性变换法滤波后的信号','窗函数滤波后的信号');

figure(11);
plot(xt,abs(fft(xn)),'blue');
hold on;
plot(xt,abs(fft(yn)),'green');
plot(xt,abs(fft(y)),'red');
title('不同方法下的频域波形');
legend('原始频谱','双线性变换法滤波后的频谱','窗函数滤波后的频谱');

figure(12);
subplot(211);

dbH1=20*log10((abs(H1)+eps)/max(abs(H1)));%化为分贝值
plot(w1/pi,dbH1);
axis([0 1 -400 200]);
title('幅频和相频特性曲线(双线性变换法)');
xlabel('Normalized Frequency(\pi rad/sample)');
ylabel('Magnitude(dB)');
grid
subplot(212);
plot(w1/pi,phase(H1)*180/pi);
xlabel('Normalized Frequency(\pi rad/sample)');
ylabel('Phase(degrees)');
grid

figure(13);
subplot(211);
dbH2=20*log10((abs(H2)+eps)/max(abs(H2)));%化为分贝值
plot(w2/pi,dbH2);
axis([0 1 -100 50]);
title('幅频和相频特性曲线(窗函数法)');
xlabel('Normalized Frequency(\pi rad/sample)');
ylabel('Magnitude(dB)');
grid
subplot(212);
h2 = abs(phase(H2)*180/pi);
h3 = mod(h2,360);%将超过360度的部分转换为360度以内
plot(w2/pi,-h3);
xlabel('Normalized Frequency(\pi rad/sample)');
ylabel('Phase(degrees)');
grid
