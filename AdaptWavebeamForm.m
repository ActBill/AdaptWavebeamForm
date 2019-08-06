%% hello word!!!非自适应波束形成
clc;
clear all;
close all;
%% 波形参数
RF = 10e9; % 载频10GHz
c = 3e8;
lambda = c/RF;

%% 阵列参数
N = 16; % 阵元个数
d = lambda/2; % 阵元间隔
theta0 = 0; % 天线导向
win = hamming(N); % 窗矢量
k_theta0 = 2*pi*d*sind(theta0)/lambda; 
a_theta = exp(-1j*(0:1:N-1).*k_theta0).'; % 空域导向矢量
h = win.*conj(a_theta);

%% 空域信号参数
Amp = 10; % 信号幅度
Pnoise =  Amp^2; % 噪声功率
% 干扰源1
J_theta1 = 18; % 到达角
JSR1 = 10^(50/10); % 干扰源信号比
k_theta1 = 2*pi*d*sind(J_theta1)/lambda; % 空域频率
% 干扰源2
J_theta2 = -2;
JSR2 = 10^(30/10); 
k_theta2 = 2*pi*d*sind(J_theta2)/lambda;

%% 匹配滤波器参数
kappa = Pnoise;
a_J1 = exp(-1j*(0:1:N-1)*k_theta1).'; % 干扰源1导向矢量
a_J2 = exp(-1j*(0:1:N-1)*k_theta2).'; % 干扰源2导向矢量
S_I1 = Amp^2*JSR1*conj(a_J1)*a_J1.'; % 干扰源1协方差矩阵
S_I2 = Amp^2*JSR1*conj(a_J2)*a_J2.'; % 干扰源1协方差矩阵
S_I = Pnoise*eye(N) + S_I1 + S_I2; % 干扰协方差矩阵
h = kappa*inv(S_I)*conj(a_theta); % 匹配滤波器矢量

%% 波束形成 
theta = -90:0.5:90;% 到达角
y = Amp^2*exp(-1j*2*pi*(0:1:N-1).'*d*sind(theta)/lambda); % 阵列快拍
z_theta = h.'*y; % 波束形成输出
z_theta_nowin = conj(a_theta).'*y;

%% display
figure;
plot(theta, db(abs(z_theta)./(Amp^2*N)), 'linewidth', 1.5)
hold on
plot([J_theta1 J_theta1], [-60 0], 'g--', 'linewidth', 1)
hold on
plot([J_theta2 J_theta2], [-60 0], 'k--', 'linewidth', 1)
xlabel('到达角(\circ)')
ylabel('归一化阵列响应(dB)')
grid on
axis([-90 90 -60 0])
figure;
plot(theta, db(abs(z_theta_nowin)./(Amp^2*N)), 'linewidth', 1.5)
hold on
plot([J_theta1 J_theta1], [-60 0], 'g--', 'linewidth', 1)
hold on
plot([J_theta2 J_theta2], [-60 0], 'k--', 'linewidth', 1)
xlabel('到达角(\circ)')
ylabel('归一化阵列响应(dB)')
grid on
axis([-90 90 -60 0])

