clc;
clear;
close all;
% @time: Dec 2021
% @author: yanhao
% @software: MATLAB

%待传输数据96个十进制数
data=0:95;
%生成比特流
bits_stream=dec2bin(data,8)-48;%每行对应一个字节
bits_stream=bits_stream';
%QPSK映射
bits_stream=bits_stream*2-1;%1表示数值的1，0表示数值的-1
qpsk_symbol=(1i*bits_stream(1:2:end)+bits_stream(2:2:end))/sqrt(2);%MSB
%%绘制星座图
figure(1)
hold on,grid on
plot(real(qpsk_symbol),imag(qpsk_symbol),'bo');
xlabel("Real");
ylabel("Image");
title("Binary Sequences");
%串并交换
qpsk_symbol_s2p=reshape(qpsk_symbol,48,8);
%OFDM符号数据子载波分配（8个OFDM数据符号/帧）
ofdm_DataCarrier_allocator = zeros(64,8);
ofdm_DataCarrier_allocator(33+[-21 -7 7],:)=1;%导频子载波波段前三个为1最后一个为-1
ofdm_DataCarrier_allocator(33+21,:)=-1;
ofdm_DataCarrier_allocator(33+(-26:-22),:)=qpsk_symbol_s2p(1:5,:);
ofdm_DataCarrier_allocator(33+(-20:-8),:)=qpsk_symbol_s2p(6:18,:);
ofdm_DataCarrier_allocator(33+(-6:-1),:)=qpsk_symbol_s2p(19:24,:);
ofdm_DataCarrier_allocator(33+(1:6),:)=qpsk_symbol_s2p(25:30,:);
ofdm_DataCarrier_allocator(33+(8:20),:)=qpsk_symbol_s2p(31:43,:);
ofdm_DataCarrier_allocator(33+(22:26),:)=qpsk_symbol_s2p(44:48,:);
%%绘制星座图
figure(2)
hold on,grid on
plot(real(ofdm_DataCarrier_allocator),imag(ofdm_DataCarrier_allocator),'bo');
xlabel("Real");
ylabel("Image");
title("Binary Sequences with SubCarrier");

load sync_word_data.mat
sync_word=[sync_word1;sync_word2];
ofdm_freq=[sync_word',ofdm_DataCarrier_allocator];%在信号前加入同步字
ofdm_time=ifft(ifftshift(ofdm_freq,1),64,1);%IFFT
ofdm_time_addCP=[ofdm_time(49:64,:);ofdm_time];%加循环前缀
ofdm_time_addCP=ofdm_time_addCP(:);%并串变换
%比较第k+1个OFDM信号(待发送的一个OFDM信号)
figure(3)
hold on
k=3;
a1=64*fftshift(fft(ofdm_time_addCP(80*k+17:80*k+80)));
stem(real(a1),'r*-.');
stem(imag(a1),'bo-.');
xlabel("Time");
ylabel("Amplitude");
title("One of OFDM Symbol of Transmition Signal");
%比较整个时域信号(整个待发送的OFDM信号)
figure(4)
hold on
plot(64*real(ofdm_time_addCP),'r-*');
xlabel("Time");
ylabel("Amplitude");
title("Whole Transmition Signal");

%PlutoSDR发信号
Tdata=ofdm_time_addCP;
txPluto=sdrtx('Pluto');
txPluto.CenterFrequency=2500000000;%设置中心频率2.5GHz
txPluto.BasebandSampleRate=1e6;%设置采样率1MHz
txPluto.Gain=0;
%PlutoSDR反复发送数据
txPluto.transmitRepeat(Tdata);
%PlutoSDR收信号
rxPluto=sdrrx('Pluto','OutputDataType','double','GainSource','Manual');
rxPluto.CenterFrequency=2500000000;%设置中心频率2.5GHz
rxPluto.BasebandSampleRate=1e6;%设置采样率1MHz
rxPluto.SamplesPerFrame=1e5;%设置采样点数100000
rxPluto.Gain=40;
%PlutoSDR接收数据
for k=1:3
   rxPluto(); 
end
Rdata=rxPluto();
release(txPluto);
%save OFDM_RX_Pluto.mat Rdata

clc;
clear;
close all;
load sync_word_data.mat;
load OFDM_RX_Pluto.mat;
Rx_Sig=Rdata(end-2000:end);
N=0;%画图编号
%% 接收信号
N=N+1;
figure(N);
hold on
plot(real(Rx_Sig),'r');
plot(imag(Rx_Sig),'b');
legend("Real","Imag");
title("Received Signal")
%% 帧同步
%同步码转为时域
sync_word1_time=ifft(ifftshift(sync_word1));
sync_word2_time=ifft(ifftshift(sync_word2,1));
%时域同步字加CP
sync_word1_time_cp=[sync_word1_time(end-15:end),sync_word1_time];
sync_word2_time_cp=[sync_word2_time(end-15:end),sync_word2_time];
N=N+1;
figure(N);
subplot(2,1,1);
hold on,grid on
plot(real(sync_word1_time_cp(17:end)),'b');
set(gca,'XTickMode','manual','XTick',[1,17,33,49,64]);
title("Real of Sync_word1_time with Circuit Prefix");
subplot(2,1,2);
hold on,grid on
plot(imag(sync_word1_time_cp(17:end)),'r');
set(gca,'XTickMode','manual','XTick',[1,17,33,49,64]);
title("Imag of Sync_word1_time with Circuit Prefix");
%% 帧同步
%时域同步字1加CP的对称特性体现在：
%出现plateau，即17个连续的最大值；发现第18个点也较高，但仍认为是17个，原因是只有一个点没有匹配，其他点都配对了
sdata=[zeros(1,100),sync_word1_time_cp,zeros(1,100)];
for d=1:length(sdata)-63
    data_corr(d)=sum(conj(sdata(d:d+31)).*sdata(d+32:d+63));
end
N=N+1;
figure(N);
hold on
plot(abs(data_corr),'b');
scatter([101,117],abs(data_corr([101,117])),'filled');
title("Correlation of Synchronic Words 1 with Circuit Prefix");
%% 帧同步
%利用同步字1加CP实现时间同步
for d=1:length(Rx_Sig)-63
    data_corr(d)=sum(conj(Rx_Sig(d:d+31)).*Rx_Sig(d+32:d+63));
end
N=N+1;
figure(N);
hold on
plot(abs(data_corr),'b');
title("Correlation of Received Signal");
index=find(abs(data_corr)>max(abs(data_corr))*9/10);
start_data=index(1)%数据的起始点
scatter(start_data,abs(data_corr(start_data)),'filled');
Rx_Sig_frame=Rx_Sig(start_data:start_data+799);%读取一帧数据
%% 频偏估计
x_ofdm1=Rx_Sig_frame(17:80);
phase=angle(conj(x_ofdm1(1:32)).*x_ofdm1(33:64));
delta_f=angle(sum(conj(x_ofdm1(1:32)).*x_ofdm1(33:64)));
if delta_f>pi/2
    delta_f=delta_f-pi;
elseif delta_f<-pi/2
    delta_f=delta_f+pi;
else 
end
for k=1:length(phase)
    if phase(k)>pi/2
        phase(k)=phase(k)-pi;
    elseif phase(k)<-pi/2
        phase(k)=phase(k)+pi;
    else
    end
end
N=N+1;
figure(N);
hold on
plot(phase*180/pi,'b');
title("Frequency Offset Estimation");
%% 频偏补偿
%Rx_Sig_frame=Rx_Sig_frame.*exp(-2j*delta_f*(0:799)'/64);
Rx_Sig_frame=Rx_Sig_frame.*exp(-j*delta_f/2);
%FFT
N=N+1;
figure(N);
grid on
for d=0
    x_frame_s2p=reshape(Rx_Sig_frame,80,10);
    x_frame_s2p(:,1:2)=[];%去掉同步字共两个OFDM信号
    z_frame=fftshift(fft(x_frame_s2p(17-d:80-d,:),64,1),1);
    %绘制星座图
    N=N+1;
    figure(N);
    grid on
    plot(real(z_frame),imag(z_frame),'bo');
    title("Constellation Diagram");
    pause(1);
end
%% 信道估计
data_pilot=[1;1;1;-1];
h_est=ones(64,8);
index_pilot=[-21,-7,7,21];
h_est_pilot=z_frame(33+index_pilot,:)./(data_pilot*ones(1,8));
h_est(((-26:26)+33),:)=interp1(index_pilot,h_est_pilot,(-26:26),'linear','extrap');
%Channel Equalization
z_frame=z_frame./h_est;
%绘制星座图
N=N+1;
figure(N);
hold on,grid on
plot(real(z_frame),imag(z_frame),'bo');
title("Constellation Diagram");
%De-mapping
z_frame(33+[-32:-27,-21,-7,0,7,21,27:31],:)=[];%去掉非数据子载波
data_complex=z_frame(:).';

pause(1);
data_bin=[imag(data_complex)>0;real(data_complex)>0];
c=dec2bin(0:95,8)-48;
c=c.';
c=c(:);
N=N+1;
figure(N);
hold on
plot(c-data_bin(:),'ro-');
hold on;
xlabel("Time");
ylabel("Error");
title("Error Conditions");
ylim([-1.1 1.1])
BER=1-mean(data_bin(:)==c)

%绘制星座图
N=N+1;
figure(N);
hold on,grid on
plot(real(data_complex),imag(data_complex),'bo');
error_index=find(c-data_bin(:))
scatter(real(data_complex(error_index)),imag(data_complex(error_index)),'r','filled');
title("Constellation Diagram");