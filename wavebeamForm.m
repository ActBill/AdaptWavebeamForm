%% ������Ӧ�����γ�
clc;
clear all;
close all;
%% ���β���
RF = 10e9; % ��Ƶ10GHz
c = 3e8;
lambda = c/RF;

%% ���в���
N = 11; % ��Ԫ����
d = lambda/2; % ��Ԫ���
theta0 = 30; % ���ߵ���
win = hamming(N); % ��ʸ��
k_theta0 = 2*pi*d*sind(theta0)/lambda; 
a_theta = exp(-1j*(0:1:N-1).*k_theta0).'; % ������ʸ��
h = win.*conj(a_theta);

%% �����γ�
Amp = 10; 
theta = linspace(-90, 90, 360); % �����
y = Amp*exp(-1j*2*pi*(0:1:N-1).'*d*sind(theta)/lambda); % ���п���

z_theta = h.'*y; % �����γ����
z_theta_nowin = conj(a_theta).'*y;

%% display
figure;
plot(theta, db(abs(z_theta)./(Amp*N)), 'r', 'linewidth', 1.5)
hold on 
plot(theta, db(abs(z_theta_nowin)./(Amp*N)), 'linewidth', 1.5)
xlabel('�����(\circ)')
ylabel('��һ��������Ӧ(dB)')
grid on
axis([-90 90 -60 0])
legend('�Ӵ�','���Ӵ�', 'Location', 'northwest')

