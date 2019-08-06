%% 非自适应波束形成
clc;
clear all;
close all;
%% 波形参数
RF = 10e9; % 载频10GHz
c = 3e8;
lambda = c/RF;

%% 阵列参数
N = 11; % 阵元个数
d = lambda/2; % 阵元间隔
theta0 = 30; % 天线导向
win = hamming(N); % 窗矢量
k_theta0 = 2*pi*d*sind(theta0)/lambda; 
a_theta = exp(-1j*(0:1:N-1).*k_theta0).'; % 空域导向矢量
h = win.*conj(a_theta);

%% 波束形成
Amp = 10; 
theta = linspace(-90, 90, 360); % 到达角
y = Amp*exp(-1j*2*pi*(0:1:N-1).'*d*sind(theta)/lambda); % 阵列快拍

z_theta = h.'*y; % 波束形成输出
z_theta_nowin = conj(a_theta).'*y;

%% display
figure;
plot(theta, db(abs(z_theta)./(Amp*N)), 'r', 'linewidth', 1.5)
hold on 
plot(theta, db(abs(z_theta_nowin)./(Amp*N)), 'linewidth', 1.5)
xlabel('到达角(\circ)')
ylabel('归一化阵列响应(dB)')
grid on
axis([-90 90 -60 0])
legend('加窗','不加窗', 'Location', 'northwest')

