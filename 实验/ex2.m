%����޹ر���
clear;
close all;
%�����˲���ָ�꣺
wp=0.2*pi;ws=0.3*pi;Rp=1;As=15;
%ת����Ϊģ��ָ�꣺
Fs=1;T=1/Fs;
Omgp=(2/T)*tan(wp/2);
Omgs=(2/T)*tan(ws/2);
%�����������Ӧ���䷨,��Omgp=wc*Fs;Omgs=ws*FS;
%ģ��ԭ���˲������㣺
[n,Omgc]=buttord(Omgp,Omgs,Rp,As,'s');%��������ͽ�ֹƵ��
%˫���Ա任�����������˲���ϵ��
[B,A] = butter(n,Omgc,'s');
[bd,ad] = bilinear(B,A,Fs);
tf(bd,ad,1)%������ݺ���
%������ϵͳ��Ƶ������:
[H1,w1]=freqz(bd,ad);
%�ĵ��źţ�
xn = [-4,-2,0,-4,-6,-4,-2,-4,-6,-6,-4,-4,-6,-6,-2,6,12,8,0,-16,-38,-60,-84,-90,-66,...
  -32,-4,2,-4,8,12,12,10,6,6,6,4,0,0,0,0,0,-2,-4,0,0,0,-2,-2,0,0,-2,-2,-2,-2,0,...
  -2,-4,-2,0,-2,-4,-4,2,0,0,-2,-4,-2,0,0,-2,-4,-2,0,0,-4,-4,-2,-2,-4,-6,-6,-4,-4,8,-10,-8,...
  -6,-6,-8,-12,-10,-8,-8,-10,-12,-10,-8,-8,-10,-10,-8,-6,-6,-8,-8,-4,-2,-4,-4,-4,...
  0,0,-2,-4,-2,-2,0,-4];

yn=filter(bd,ad,xn);%�����˲�����

%�����������FIR�˲���
wp=0.2*pi;ws=0.3*pi;
deltaw=ws-wp;%���ɴ���
wc=(ws+wp)/2;

N0 = ceil(1.8*pi/deltaw);%ʹ�þ��δ�
N = N0+mod(N0+1,2);%ȷ��NΪ����
windows = (boxcar(N))';

% N0 = ceil(6.1*pi/deltaw);%ʹ�������δ�
% N = N0+mod(N0+1,2);%ȷ��NΪ����
% windows = (bartlett(N))';

% N0 = ceil(6.6*pi/deltaw);%ʹ�ú�����
% N = N0+mod(N0+1,2);%ȷ��NΪ����
% windows = (hamming(N))';
% 
% N0 = ceil(6.2*pi/deltaw);%ʹ�ú�����
% N = N0+mod(N0+1,2);%ȷ��NΪ����
% windows = (hanning(N))';

% N0 = ceil(11*pi/deltaw);%ʹ�ò���������
% N = N0+mod(N0+1,2);%ȷ��NΪ����
% windows = (blackman(N))';

% %ʹ�ÿ�����
% N = N0+mod(N0+1,2);%ȷ��NΪ����
% windows = (kaiser(N))';

n = 0:ceil(N-1);
m=n-(N-1)/2+eps;%+epsת���ɸ�����
hd = sin(wc*m)./(pi*m);
%hh = fir1(N-1,wc/pi,'low',boxcar(N));
%�Ӵ�ʩ����hn��
b=hd.*windows;

N = 13;
fc = 0.23;
b = fir1(N-1,fc,rectwin(N));
[H2,w2]=freqz(b,1);
hb = tf(b,1,1)
%y = filter(b,1,xn);%�������˲�����
y = filter(b,1,xn);%�������˲�����
%yy = fftfilt(hh,xn,length(xn));
set(0,'defaultfigurecolor','w')%ʹ�����ɫΪ��ɫ
figure(10);
xt = 0:(length(xn)-1);
plot(xt,xn,'blue');
hold on;
plot(xt,yn,'green');
plot(xt,y,'red');
title('��ͬ�����µ�ʱ����');
legend('ԭʼ�ź�','˫���Ա任���˲�����ź�','�������˲�����ź�');

figure(11);
plot(xt,abs(fft(xn)),'blue');
hold on;
plot(xt,abs(fft(yn)),'green');
plot(xt,abs(fft(y)),'red');
title('��ͬ�����µ�Ƶ����');
legend('ԭʼƵ��','˫���Ա任���˲����Ƶ��','�������˲����Ƶ��');

figure(12);
subplot(211);

dbH1=20*log10((abs(H1)+eps)/max(abs(H1)));%��Ϊ�ֱ�ֵ
plot(w1/pi,dbH1);
axis([0 1 -400 200]);
title('��Ƶ����Ƶ��������(˫���Ա任��)');
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
dbH2=20*log10((abs(H2)+eps)/max(abs(H2)));%��Ϊ�ֱ�ֵ
plot(w2/pi,dbH2);
axis([0 1 -100 50]);
title('��Ƶ����Ƶ��������(��������)');
xlabel('Normalized Frequency(\pi rad/sample)');
ylabel('Magnitude(dB)');
grid
subplot(212);
h2 = abs(phase(H2)*180/pi);
h3 = mod(h2,360);%������360�ȵĲ���ת��Ϊ360������
plot(w2/pi,-h3);
xlabel('Normalized Frequency(\pi rad/sample)');
ylabel('Phase(degrees)');
grid
