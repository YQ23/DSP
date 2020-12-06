function varargout = gui1(varargin)
% GUI1 MATLAB code for gui1.fig
%      GUI1, by itself, creates a new GUI1 or raises the existing
%      singleton*.
%
%      H = GUI1 returns the handle to a new GUI1 or the handle to
%      the existing singleton*.
%
%      GUI1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI1.M with the given input arguments.
%
%      GUI1('Property','Value',...) creates a new GUI1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui1

% Last Modified by GUIDE v2.5 03-Jul-2020 11:22:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui1_OpeningFcn, ...
                   'gui_OutputFcn',  @gui1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui1 is made visible.
function gui1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui1 (see VARARGIN)

% Choose default command line output for gui1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%按下录音按钮
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global flag;%全局变量,其他函数也需要用到
global fpath;%文件路径

flag = 0;           %录音还是打开音频文件,如果是录音flag=0,打开音频文件flag=1
fs=8000;           %采样频率
a = get(handles.edit3,'String');%获取录音时间
duration = str2double(a);%将录音时间转换为数字
%fprintf('duration:%d\n',duration);

fprintf('开始录音...\n');
%y=audiorecorder(duration*fs,fs);
recObj = audiorecorder(fs,16,1);%16位单通道,采样频率为8000Hz
disp('结束录音...');
recordblocking(recObj, duration);

%回放录音 
%fprintf('Playback recording...\n');
play(recObj); 
%sound(recObj);

%获取录音数据 
y = getaudiodata(recObj);
%fprintf('Press any key to play the recording...\n');
%sound(y);

%频域信息
yf = abs(fft(y));%进行快速傅里叶变换
fpath = './test.wav';%将录音存储为wav文件
audiowrite(fpath,y,fs);%进行存储
fprintf('完成存储...\n');

axes(findobj('tag','axes1'));%在第一个图绘制时域波形
plot(y);

axes(findobj('tag','axes2'));%在第二个图绘制频域波形
plot(yf);
set(handles.edit1,'string','原始');%设置编辑框文字为原始

%按下变女人声音按钮
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global fpath; %录音文件存储路径
    [data,fs]=audioread(fpath);           % 载入存储的语音  
    data = data + 1e-20;  %加一个很小的数,避免有的地方数据为0导致后面运算出现NaN报错的情况
    data = data/max(data);	%归一化处理
    L = length(data);          % 读入语音长度
    Framelength = 80;           % 帧长
    Windowlength = 240;         % 窗长
    N = 10;                     % LPC法运算得到的预测系数个数
    FrameNum = floor(L/Framelength)-2;     % 计算帧数
	%预测滤波器
    en_pre = zeros(L,1);       % 激励信号
    st_pre = zeros(N,1);    % 预测滤波器的状态
    %重建滤波器
    data_re = zeros(L,1);     % 重建语音
    st_re = zeros(N,1);      %重建滤波器的状态
%     % 合成滤波器
%     en_syn = zeros(L,1);    % 合成的激励信号（脉冲串）
%     data_syn = zeros(L,1);  % 合成语音
% 	last_syn = 0;           %存储上一个（或多个）段的最后一个脉冲的下标
% 	st_syn = zeros(N,1);    % 合成滤波器的状态
	% 变调不变速滤波器
    en_syn_t = zeros(L,1);   % 合成的激励信号（脉冲串）
    data_syn_t = zeros(L,1);     % 合成语音
	last_syn_t = 0;         %存储上一个（或多个）段的最后一个脉冲的下标
	st_syn_t = zeros(N,1);   % 合成滤波器的状态(变调不变速)
    hw = hamming(Windowlength);       % 海明窗
    %依次处理每帧语音
    for n = 3:FrameNum
        % 计算预测系数
        data_w = data(n*Framelength-Windowlength+1:n*Framelength).*hw;    %使用汉明窗得到加权后的语音
        [A,E] = lpc(data_w, N);   %用线性预测法计算P个预测系数,A是预测系数,E会被用来计算合成激励的能量                           
        data_f = data((n-1)*Framelength+1:n*Framelength);       % 对本帧语音进行处理
        %用filter函数对data_f计算激励,滤波器状态在不断更新
		[en1,st_pre] = filter(A,1,data_f,st_pre);
        en_pre((n-1)*Framelength+1:n*Framelength) = en1; %计算得到的激励
        %用filter函数和en1重建语音,滤波器状态在不断更新
		[data_re1,st_re] = filter(1,A,en1,st_re);
        data_rec((n-1)*Framelength+1:n*Framelength) = data_re1; %计算得到的重建语音
        %只有在得到en_pre后才会计算正确
        data_Pitch = en_pre(n*Framelength-222:n*Framelength);
        PT = fp(data_Pitch);    % 计算基音周期PT
        G = sqrt(E*PT);           % 计算合成激励的能量G
        %将基音周期减小一半,并将共振峰频率增加100Hz,重新合成语音
        delta = 100*2*pi/fs; %共振峰频率变化量
		new_PT =floor(PT/2);   %减小基音周期
        poles = roots(A);   %求出预测系数的特征根
		
		for p=1:N   %增加共振峰频率,将实轴上方的极点逆时针旋转,下方的极点顺时针旋转
			if imag(poles(p))>0 
                poles(p) = poles(p)*exp(1i*delta);
            elseif imag(poles(p))<0 
                poles(p) = poles(p)*exp(-1i*delta);
			end
		end
		A1=poly(poles);%计算出新的系数
        temp_syn_t =[1:n*Framelength-last_syn_t]';
		en_syn1_t = zeros(length(temp_syn_t),1);
		en_syn1_t(mod(temp_syn_t,new_PT)==0) = G; %某一段算出的脉冲
		en_syn1_t = en_syn1_t((n-1)*Framelength-last_syn_t+1:n*Framelength-last_syn_t);
        last_syn_t = last_syn_t+new_PT*floor((n*Framelength-last_syn_t)/new_PT);
		[data_syn1_t,st_syn_t] = filter(1,A1,en_syn1_t,st_syn_t);  %使用filter函数进行处理
		en_syn_t((n-1)*Framelength+1:n*Framelength) =  en_syn1_t;   %计算得到的合成激励
		data_syn_t((n-1)*Framelength+1:n*Framelength) = data_syn1_t;   %计算得到的合成语音
    end
     

axes(findobj('tag','axes1'));
plot(data_syn_t);%时域波形

axes(findobj('tag','axes2'));
plot(abs(fft(data_syn_t)));%频域波形
set(handles.edit1,'string','女人');
sound(data_syn_t);%播放变声后的音频



%按下变老人声按钮
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global fpath;%原始音频路径
    [data,fs]=audioread(fpath);         % 载入语音
    data = data + 1e-20;      %增加一个很小的数避免后面出现NaN
	data = data/max(data);	  %进行归一化处理
    L = length(data);         % 读入语音长度
    Framelength = 80;                % 帧长
    Windowlength = 240;               % 窗长
    N = 10;                 % lpc法的预测系数个数
    FrameNum = floor(L/Framelength)-2;     % 计算帧数
	% 预测滤波器
    en_pre = zeros(L,1);       % 激励信号
    st_pre = zeros(N,1);    % 预测滤波器的状态
    %重建滤波器
    data_re = zeros(L,1);     % 重建语音
    st_re = zeros(N,1);       %重建滤波器的状态
% 	% 合成滤波器
%     exc_syn = zeros(L,1);   % 合成的激励信号（脉冲串）
%     x1_syn = zeros(L,1);     % 合成语音
% 	last_syn = 0;   %存储上一个（或多个）段的最后一个脉冲的下标
% 	zi_syn = zeros(N,1);   % 合成滤波器的状态
    % 变速不变调滤波器（速度减慢一倍）
	v=0.5;
    en_syn_v = zeros(v\L,1);   % 合成的激励信号（脉冲串）
    data_syn_v = zeros(v\L,1);     % 合成语音
	last_syn_v = 0;   %存储上一个（或多个）段的最后一个脉冲的下标
	st_syn_v = zeros(N,1);   % 合成滤波器的状态
    hw = hamming(Windowlength);       % 汉明窗
    % 依次处理每帧语音
    for n = 3:FrameNum
        data_w = data(n*Framelength-Windowlength+1:n*Framelength).*hw;    %汉明窗加权后的语音
        [A,E] = lpc(data_w, N);            %用线性预测法计算P个预测系数,A是预测系数,E会被用来计算合成激励的能量
        data_f = data((n-1)*Framelength+1:n*Framelength);       % 对本帧语音做处理
        %用filter函数data_f计算激励,滤波器状态不断更新
		[en1,st_pre] = filter(A,1,data_f,st_pre);
        en_pre((n-1)*Framelength+1:n*Framelength) = en1; %计算得到的激励
        %用filter函数和en重建语音,滤波器状态不断更新
		[data_re1,st_re] = filter(1,A,en1,st_re);
        data_re((n-1)*Framelength+1:n*Framelength) = data_re1; %计算得到的重建语音
        %只有在得到en后才会计算正确
        data_Pitch = en_pre(n*Framelength-222:n*Framelength);
        PT = fp(data_Pitch);    % 计算基音周期PT
        G = sqrt(E*PT);           % 计算合成激励的能量G
        % 这里不改变基音周期和预测系数，将合成激励的长度增加一倍，再作为filter
        % 的输入得到新的合成语音，从而使速度变慢了，但音调没有变。
		Framelength_v = floor(Framelength/v);%改变速度
		temp_syn_v = [1:n*Framelength_v-last_syn_v]';
		en_syn1_v = zeros(length(temp_syn_v),1);
		en_syn1_v(mod(temp_syn_v,PT)==0) = G; %某一段算出的脉冲
		en_syn1_v = en_syn1_v((n-1)*Framelength_v-last_syn_v+1:n*Framelength_v-last_syn_v);
		[data_syn1_v,st_syn_v] = filter(1,A,en_syn1_v,st_syn_v);%使用filter函数进行处理		
     	last_syn_v = last_syn_v+PT*floor((n*Framelength_v-last_syn_v)/PT);%计算最后一个脉冲下标   
        en_syn_v((n-1)*Framelength_v+1:n*Framelength_v) =en_syn1_v;  %计算得到的加长合成激励
        data_syn_v((n-1)*Framelength_v+1:n*Framelength_v) = data_syn1_v;   %计算得到的加长合成语音 
    end

axes(findobj('tag','axes1'));
plot(data_syn_v);%时域波形

axes(findobj('tag','axes2'));
plot(abs(fft(data_syn_v)));%快速傅里叶变换得到频域波形
set(handles.edit1,'string','老人');
sound(data_syn_v);%播放处理后的音频


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global fpath;%原始音频路径
    [data,fs] = audioread(fpath);             % 载入需要处理的语音
    data = data + 1e-20;    %加一个很小的数避免后面运算出现NaN
	data = data/max(data);	%对数据进行归一化处理
    L = length(data);          % 读入语音长度
    Framelength = 80;                % 帧长
    Windowlength = 240;               % 窗长
    N = 10;                       % lpc法得到的预测系数个数
    FrameNum = floor(L/Framelength)-2;     % 计算帧数
	% 预测滤波器
    en = zeros(L,1);       % 激励信号
    st_pre = zeros(N,1);    % 预测滤波器的状态
    % 重建滤波器
    data_re = zeros(L,1);     % 重建语音
    st_re = zeros(N,1);       % 重建滤波器的状态
% 	% 合成滤波器
%     exc_syn = zeros(L,1);   % 合成的激励信号（脉冲串）
%     x1_syn = zeros(L,1);     % 合成语音
% 	last_syn = 0;   %存储上一个（或多个）段的最后一个脉冲的下标
% 	zi_syn = zeros(N,1);   % 合成滤波器的状态
	% 变调不变速滤波器
    en_syn_t = zeros(L,1);   % 合成的激励信号（脉冲串）
    data_syn_t = zeros(L,1);     % 合成语音
	last_syn_t = 0;   %存储上一个（或多个）段的最后一个脉冲的下标
	st_syn_t = zeros(N,1);   % 合成滤波器的状态
    hw = hamming(Windowlength);       %汉明窗
    % 依次处理每帧语音
    for n = 3:FrameNum
        % 计算预测系数
        data_w =data(n*Framelength-Windowlength+1:n*Framelength).*hw;    %汉明窗加权后的语音
        [A,E] = lpc(data_w, N);  %用线性预测法计算P个预测系数,A是预测系数，E会被用来计算合成激励的能量
        data_f = data((n-1)*Framelength+1:n*Framelength);       % 对本帧语音做处理
        % 用filter函数data_f计算激励,滤波器状态不断更新
		[en1,st_pre] = filter(A,1,data_f,st_pre);
        en((n-1)*Framelength+1:n*Framelength) = en1; %计算得到的激励
        %  用filter函数和en1重建语音,滤波器状态不断更新
		[data_re1,st_re] = filter(1,A,en1,st_re);
        data_re((n-1)*Framelength+1:n*Framelength) = data_re1; %计算得到的重建语音
        % 注意下面只有在得到en后才会计算正确
        data_Pitch = en(n*Framelength-222:n*Framelength);
        PT = fp(data_Pitch);    % 计算基音周期PT
        G = sqrt(E*PT);                % 计算合成激励的能量G
        % 将基音周期减小一半,并将共振峰频率增加700Hz,再重新合成语音
        delta = 700*2*pi/fs;
		new_PT =floor(PT/2);   %减小基音周期
        poles = roots(A);  %计算预测系数的特征根
		for p=1:N   %增加共振峰频率，将实轴上方的极点逆时针旋转，下方极点顺时针旋转
			if imag(poles(p))>0 
                poles(p) = poles(p)*exp(1i*delta);
			elseif imag(poles(p))<0 
                poles(p) = poles(p)*exp(-1i*delta);
			end
		end
		A1=poly(poles);%新的预测系数的特征根
        temp_syn_t = [1:n*Framelength-last_syn_t]';
		en_syn1_t = zeros(length(temp_syn_t),1);
		en_syn1_t(mod(temp_syn_t,new_PT)==0) = G; %某一段算出的脉冲
		en_syn1_t = en_syn1_t((n-1)*Framelength-last_syn_t+1:n*Framelength-last_syn_t);
        last_syn_t = last_syn_t+new_PT*floor((n*Framelength-last_syn_t)/new_PT);
		[data_syn1_t,st_syn_t] = filter(1,A1,en_syn1_t,st_syn_t);%使用filter函数进行处理
		en_syn_t((n-1)*Framelength+1:n*Framelength) =  en_syn1_t;   %计算得到的合成激励
		data_syn_t((n-1)*Framelength+1:n*Framelength) = data_syn1_t;   %计算得到的合成语音
    end

axes(findobj('tag','axes1'));
plot(data_syn_t);%时域波形

axes(findobj('tag','axes2'));
plot(abs(fft(data_syn_t)));%快速傅里叶变换得到频域波形
set(handles.edit1,'string','小孩');
sound(data_syn_t);%播放处理后的音频



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%按下打开音频文件按钮
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('.wav','选择音频文件');
global fpath;
global flag;
flag = 1;
fprintf('pathname:%s\n',pathname);
fprintf('filename:%s\n',filename);
fpath = strcat(pathname,filename);%得到文件路径
[y,fs]=audioread(fpath); %读取文件
fprintf('fs:%d\n',fs);
yf = abs(fft(y));%进行快速傅里叶变换
axes(findobj('tag','axes1'));
sound(y);%播放音频文件
plot(y);

axes(findobj('tag','axes2'));
plot(yf);
set(handles.edit1,'string','原始');
 



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fpath;
global flag;
%[y,fs]=audioread(fpath); 
if flag == 0 %如果是录音,则播放录音文件
    [y,fs] = audioread('./test.wav');
else     %如果是打开的音频文件,则播放对应文件
    [y,fs] = audioread(fpath);
end

axes(findobj('tag','axes1'));
plot(y);%时域波形

axes(findobj('tag','axes2'));
plot(abs(fft(y)));%快速傅里叶变换得到频域波形
set(handles.edit1,'string','原始');
sound(y,fs);


 %计算基音周期
 function PT = fp(s)%s为223*1的向量
 fs = 8000;
 N = 5;
 fc = 700;
[B, A] = butter(N, fc/(0.5*fs));
s = filter(B,A,s);
R = zeros(143,1);
for k=1:143
    R(k) = s(144:223)'*s(144-k:223-k);%计算自相关函数
end
[R1,T1] = max(R(80:143));
T1 = T1 + 79;
R1 = R1/(norm(s(144-T1:223-T1))+1);
[R2,T2] = max(R(40:79));
T2 = T2 + 39;
R2 = R2/(norm(s(144-T2:223-T2))+1);
[R3,T3] = max(R(20:39));
T3 = T3 + 19;
R3 = R3/(norm(s(144-T3:223-T3))+1);
Top = T1;
Rop = R1;
if R2 >= 0.85*Rop %设置阈值为0.85
    Rop = R2;
    Top = T2;
end
if R3 > 0.85*Rop
    Rop = R3;
    Top = T3;
end
PT = Top;


