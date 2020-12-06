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


%����¼����ť
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global flag;%ȫ�ֱ���,��������Ҳ��Ҫ�õ�
global fpath;%�ļ�·��

flag = 0;           %¼�����Ǵ���Ƶ�ļ�,�����¼��flag=0,����Ƶ�ļ�flag=1
fs=8000;           %����Ƶ��
a = get(handles.edit3,'String');%��ȡ¼��ʱ��
duration = str2double(a);%��¼��ʱ��ת��Ϊ����
%fprintf('duration:%d\n',duration);

fprintf('��ʼ¼��...\n');
%y=audiorecorder(duration*fs,fs);
recObj = audiorecorder(fs,16,1);%16λ��ͨ��,����Ƶ��Ϊ8000Hz
disp('����¼��...');
recordblocking(recObj, duration);

%�ط�¼�� 
%fprintf('Playback recording...\n');
play(recObj); 
%sound(recObj);

%��ȡ¼������ 
y = getaudiodata(recObj);
%fprintf('Press any key to play the recording...\n');
%sound(y);

%Ƶ����Ϣ
yf = abs(fft(y));%���п��ٸ���Ҷ�任
fpath = './test.wav';%��¼���洢Ϊwav�ļ�
audiowrite(fpath,y,fs);%���д洢
fprintf('��ɴ洢...\n');

axes(findobj('tag','axes1'));%�ڵ�һ��ͼ����ʱ����
plot(y);

axes(findobj('tag','axes2'));%�ڵڶ���ͼ����Ƶ����
plot(yf);
set(handles.edit1,'string','ԭʼ');%���ñ༭������Ϊԭʼ

%���±�Ů��������ť
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global fpath; %¼���ļ��洢·��
    [data,fs]=audioread(fpath);           % ����洢������  
    data = data + 1e-20;  %��һ����С����,�����еĵط�����Ϊ0���º����������NaN��������
    data = data/max(data);	%��һ������
    L = length(data);          % ������������
    Framelength = 80;           % ֡��
    Windowlength = 240;         % ����
    N = 10;                     % LPC������õ���Ԥ��ϵ������
    FrameNum = floor(L/Framelength)-2;     % ����֡��
	%Ԥ���˲���
    en_pre = zeros(L,1);       % �����ź�
    st_pre = zeros(N,1);    % Ԥ���˲�����״̬
    %�ؽ��˲���
    data_re = zeros(L,1);     % �ؽ�����
    st_re = zeros(N,1);      %�ؽ��˲�����״̬
%     % �ϳ��˲���
%     en_syn = zeros(L,1);    % �ϳɵļ����źţ����崮��
%     data_syn = zeros(L,1);  % �ϳ�����
% 	last_syn = 0;           %�洢��һ�����������ε����һ��������±�
% 	st_syn = zeros(N,1);    % �ϳ��˲�����״̬
	% ����������˲���
    en_syn_t = zeros(L,1);   % �ϳɵļ����źţ����崮��
    data_syn_t = zeros(L,1);     % �ϳ�����
	last_syn_t = 0;         %�洢��һ�����������ε����һ��������±�
	st_syn_t = zeros(N,1);   % �ϳ��˲�����״̬(���������)
    hw = hamming(Windowlength);       % ������
    %���δ���ÿ֡����
    for n = 3:FrameNum
        % ����Ԥ��ϵ��
        data_w = data(n*Framelength-Windowlength+1:n*Framelength).*hw;    %ʹ�ú������õ���Ȩ�������
        [A,E] = lpc(data_w, N);   %������Ԥ�ⷨ����P��Ԥ��ϵ��,A��Ԥ��ϵ��,E�ᱻ��������ϳɼ���������                           
        data_f = data((n-1)*Framelength+1:n*Framelength);       % �Ա�֡�������д���
        %��filter������data_f���㼤��,�˲���״̬�ڲ��ϸ���
		[en1,st_pre] = filter(A,1,data_f,st_pre);
        en_pre((n-1)*Framelength+1:n*Framelength) = en1; %����õ��ļ���
        %��filter������en1�ؽ�����,�˲���״̬�ڲ��ϸ���
		[data_re1,st_re] = filter(1,A,en1,st_re);
        data_rec((n-1)*Framelength+1:n*Framelength) = data_re1; %����õ����ؽ�����
        %ֻ���ڵõ�en_pre��Ż������ȷ
        data_Pitch = en_pre(n*Framelength-222:n*Framelength);
        PT = fp(data_Pitch);    % �����������PT
        G = sqrt(E*PT);           % ����ϳɼ���������G
        %���������ڼ�Сһ��,���������Ƶ������100Hz,���ºϳ�����
        delta = 100*2*pi/fs; %�����Ƶ�ʱ仯��
		new_PT =floor(PT/2);   %��С��������
        poles = roots(A);   %���Ԥ��ϵ����������
		
		for p=1:N   %���ӹ����Ƶ��,��ʵ���Ϸ��ļ�����ʱ����ת,�·��ļ���˳ʱ����ת
			if imag(poles(p))>0 
                poles(p) = poles(p)*exp(1i*delta);
            elseif imag(poles(p))<0 
                poles(p) = poles(p)*exp(-1i*delta);
			end
		end
		A1=poly(poles);%������µ�ϵ��
        temp_syn_t =[1:n*Framelength-last_syn_t]';
		en_syn1_t = zeros(length(temp_syn_t),1);
		en_syn1_t(mod(temp_syn_t,new_PT)==0) = G; %ĳһ�����������
		en_syn1_t = en_syn1_t((n-1)*Framelength-last_syn_t+1:n*Framelength-last_syn_t);
        last_syn_t = last_syn_t+new_PT*floor((n*Framelength-last_syn_t)/new_PT);
		[data_syn1_t,st_syn_t] = filter(1,A1,en_syn1_t,st_syn_t);  %ʹ��filter�������д���
		en_syn_t((n-1)*Framelength+1:n*Framelength) =  en_syn1_t;   %����õ��ĺϳɼ���
		data_syn_t((n-1)*Framelength+1:n*Framelength) = data_syn1_t;   %����õ��ĺϳ�����
    end
     

axes(findobj('tag','axes1'));
plot(data_syn_t);%ʱ����

axes(findobj('tag','axes2'));
plot(abs(fft(data_syn_t)));%Ƶ����
set(handles.edit1,'string','Ů��');
sound(data_syn_t);%���ű��������Ƶ



%���±���������ť
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global fpath;%ԭʼ��Ƶ·��
    [data,fs]=audioread(fpath);         % ��������
    data = data + 1e-20;      %����һ����С��������������NaN
	data = data/max(data);	  %���й�һ������
    L = length(data);         % ������������
    Framelength = 80;                % ֡��
    Windowlength = 240;               % ����
    N = 10;                 % lpc����Ԥ��ϵ������
    FrameNum = floor(L/Framelength)-2;     % ����֡��
	% Ԥ���˲���
    en_pre = zeros(L,1);       % �����ź�
    st_pre = zeros(N,1);    % Ԥ���˲�����״̬
    %�ؽ��˲���
    data_re = zeros(L,1);     % �ؽ�����
    st_re = zeros(N,1);       %�ؽ��˲�����״̬
% 	% �ϳ��˲���
%     exc_syn = zeros(L,1);   % �ϳɵļ����źţ����崮��
%     x1_syn = zeros(L,1);     % �ϳ�����
% 	last_syn = 0;   %�洢��һ�����������ε����һ��������±�
% 	zi_syn = zeros(N,1);   % �ϳ��˲�����״̬
    % ���ٲ�����˲������ٶȼ���һ����
	v=0.5;
    en_syn_v = zeros(v\L,1);   % �ϳɵļ����źţ����崮��
    data_syn_v = zeros(v\L,1);     % �ϳ�����
	last_syn_v = 0;   %�洢��һ�����������ε����һ��������±�
	st_syn_v = zeros(N,1);   % �ϳ��˲�����״̬
    hw = hamming(Windowlength);       % ������
    % ���δ���ÿ֡����
    for n = 3:FrameNum
        data_w = data(n*Framelength-Windowlength+1:n*Framelength).*hw;    %��������Ȩ�������
        [A,E] = lpc(data_w, N);            %������Ԥ�ⷨ����P��Ԥ��ϵ��,A��Ԥ��ϵ��,E�ᱻ��������ϳɼ���������
        data_f = data((n-1)*Framelength+1:n*Framelength);       % �Ա�֡����������
        %��filter����data_f���㼤��,�˲���״̬���ϸ���
		[en1,st_pre] = filter(A,1,data_f,st_pre);
        en_pre((n-1)*Framelength+1:n*Framelength) = en1; %����õ��ļ���
        %��filter������en�ؽ�����,�˲���״̬���ϸ���
		[data_re1,st_re] = filter(1,A,en1,st_re);
        data_re((n-1)*Framelength+1:n*Framelength) = data_re1; %����õ����ؽ�����
        %ֻ���ڵõ�en��Ż������ȷ
        data_Pitch = en_pre(n*Framelength-222:n*Framelength);
        PT = fp(data_Pitch);    % �����������PT
        G = sqrt(E*PT);           % ����ϳɼ���������G
        % ���ﲻ�ı�������ں�Ԥ��ϵ�������ϳɼ����ĳ�������һ��������Ϊfilter
        % ������õ��µĺϳ��������Ӷ�ʹ�ٶȱ����ˣ�������û�б䡣
		Framelength_v = floor(Framelength/v);%�ı��ٶ�
		temp_syn_v = [1:n*Framelength_v-last_syn_v]';
		en_syn1_v = zeros(length(temp_syn_v),1);
		en_syn1_v(mod(temp_syn_v,PT)==0) = G; %ĳһ�����������
		en_syn1_v = en_syn1_v((n-1)*Framelength_v-last_syn_v+1:n*Framelength_v-last_syn_v);
		[data_syn1_v,st_syn_v] = filter(1,A,en_syn1_v,st_syn_v);%ʹ��filter�������д���		
     	last_syn_v = last_syn_v+PT*floor((n*Framelength_v-last_syn_v)/PT);%�������һ�������±�   
        en_syn_v((n-1)*Framelength_v+1:n*Framelength_v) =en_syn1_v;  %����õ��ļӳ��ϳɼ���
        data_syn_v((n-1)*Framelength_v+1:n*Framelength_v) = data_syn1_v;   %����õ��ļӳ��ϳ����� 
    end

axes(findobj('tag','axes1'));
plot(data_syn_v);%ʱ����

axes(findobj('tag','axes2'));
plot(abs(fft(data_syn_v)));%���ٸ���Ҷ�任�õ�Ƶ����
set(handles.edit1,'string','����');
sound(data_syn_v);%���Ŵ�������Ƶ


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global fpath;%ԭʼ��Ƶ·��
    [data,fs] = audioread(fpath);             % ������Ҫ���������
    data = data + 1e-20;    %��һ����С������������������NaN
	data = data/max(data);	%�����ݽ��й�һ������
    L = length(data);          % ������������
    Framelength = 80;                % ֡��
    Windowlength = 240;               % ����
    N = 10;                       % lpc���õ���Ԥ��ϵ������
    FrameNum = floor(L/Framelength)-2;     % ����֡��
	% Ԥ���˲���
    en = zeros(L,1);       % �����ź�
    st_pre = zeros(N,1);    % Ԥ���˲�����״̬
    % �ؽ��˲���
    data_re = zeros(L,1);     % �ؽ�����
    st_re = zeros(N,1);       % �ؽ��˲�����״̬
% 	% �ϳ��˲���
%     exc_syn = zeros(L,1);   % �ϳɵļ����źţ����崮��
%     x1_syn = zeros(L,1);     % �ϳ�����
% 	last_syn = 0;   %�洢��һ�����������ε����һ��������±�
% 	zi_syn = zeros(N,1);   % �ϳ��˲�����״̬
	% ����������˲���
    en_syn_t = zeros(L,1);   % �ϳɵļ����źţ����崮��
    data_syn_t = zeros(L,1);     % �ϳ�����
	last_syn_t = 0;   %�洢��һ�����������ε����һ��������±�
	st_syn_t = zeros(N,1);   % �ϳ��˲�����״̬
    hw = hamming(Windowlength);       %������
    % ���δ���ÿ֡����
    for n = 3:FrameNum
        % ����Ԥ��ϵ��
        data_w =data(n*Framelength-Windowlength+1:n*Framelength).*hw;    %��������Ȩ�������
        [A,E] = lpc(data_w, N);  %������Ԥ�ⷨ����P��Ԥ��ϵ��,A��Ԥ��ϵ����E�ᱻ��������ϳɼ���������
        data_f = data((n-1)*Framelength+1:n*Framelength);       % �Ա�֡����������
        % ��filter����data_f���㼤��,�˲���״̬���ϸ���
		[en1,st_pre] = filter(A,1,data_f,st_pre);
        en((n-1)*Framelength+1:n*Framelength) = en1; %����õ��ļ���
        %  ��filter������en1�ؽ�����,�˲���״̬���ϸ���
		[data_re1,st_re] = filter(1,A,en1,st_re);
        data_re((n-1)*Framelength+1:n*Framelength) = data_re1; %����õ����ؽ�����
        % ע������ֻ���ڵõ�en��Ż������ȷ
        data_Pitch = en(n*Framelength-222:n*Framelength);
        PT = fp(data_Pitch);    % �����������PT
        G = sqrt(E*PT);                % ����ϳɼ���������G
        % ���������ڼ�Сһ��,���������Ƶ������700Hz,�����ºϳ�����
        delta = 700*2*pi/fs;
		new_PT =floor(PT/2);   %��С��������
        poles = roots(A);  %����Ԥ��ϵ����������
		for p=1:N   %���ӹ����Ƶ�ʣ���ʵ���Ϸ��ļ�����ʱ����ת���·�����˳ʱ����ת
			if imag(poles(p))>0 
                poles(p) = poles(p)*exp(1i*delta);
			elseif imag(poles(p))<0 
                poles(p) = poles(p)*exp(-1i*delta);
			end
		end
		A1=poly(poles);%�µ�Ԥ��ϵ����������
        temp_syn_t = [1:n*Framelength-last_syn_t]';
		en_syn1_t = zeros(length(temp_syn_t),1);
		en_syn1_t(mod(temp_syn_t,new_PT)==0) = G; %ĳһ�����������
		en_syn1_t = en_syn1_t((n-1)*Framelength-last_syn_t+1:n*Framelength-last_syn_t);
        last_syn_t = last_syn_t+new_PT*floor((n*Framelength-last_syn_t)/new_PT);
		[data_syn1_t,st_syn_t] = filter(1,A1,en_syn1_t,st_syn_t);%ʹ��filter�������д���
		en_syn_t((n-1)*Framelength+1:n*Framelength) =  en_syn1_t;   %����õ��ĺϳɼ���
		data_syn_t((n-1)*Framelength+1:n*Framelength) = data_syn1_t;   %����õ��ĺϳ�����
    end

axes(findobj('tag','axes1'));
plot(data_syn_t);%ʱ����

axes(findobj('tag','axes2'));
plot(abs(fft(data_syn_t)));%���ٸ���Ҷ�任�õ�Ƶ����
set(handles.edit1,'string','С��');
sound(data_syn_t);%���Ŵ�������Ƶ



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

%���´���Ƶ�ļ���ť
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('.wav','ѡ����Ƶ�ļ�');
global fpath;
global flag;
flag = 1;
fprintf('pathname:%s\n',pathname);
fprintf('filename:%s\n',filename);
fpath = strcat(pathname,filename);%�õ��ļ�·��
[y,fs]=audioread(fpath); %��ȡ�ļ�
fprintf('fs:%d\n',fs);
yf = abs(fft(y));%���п��ٸ���Ҷ�任
axes(findobj('tag','axes1'));
sound(y);%������Ƶ�ļ�
plot(y);

axes(findobj('tag','axes2'));
plot(yf);
set(handles.edit1,'string','ԭʼ');
 



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
if flag == 0 %�����¼��,�򲥷�¼���ļ�
    [y,fs] = audioread('./test.wav');
else     %����Ǵ򿪵���Ƶ�ļ�,�򲥷Ŷ�Ӧ�ļ�
    [y,fs] = audioread(fpath);
end

axes(findobj('tag','axes1'));
plot(y);%ʱ����

axes(findobj('tag','axes2'));
plot(abs(fft(y)));%���ٸ���Ҷ�任�õ�Ƶ����
set(handles.edit1,'string','ԭʼ');
sound(y,fs);


 %�����������
 function PT = fp(s)%sΪ223*1������
 fs = 8000;
 N = 5;
 fc = 700;
[B, A] = butter(N, fc/(0.5*fs));
s = filter(B,A,s);
R = zeros(143,1);
for k=1:143
    R(k) = s(144:223)'*s(144-k:223-k);%��������غ���
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
if R2 >= 0.85*Rop %������ֵΪ0.85
    Rop = R2;
    Top = T2;
end
if R3 > 0.85*Rop
    Rop = R3;
    Top = T3;
end
PT = Top;


