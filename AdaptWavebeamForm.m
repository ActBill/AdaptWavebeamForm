%% hello word!!!������Ӧ�����γ�
clc;
clear all;
close all;
%% ���β���
RF = 10e9; % ��Ƶ10GHz
c = 3e8;
lambda = c/RF;

%% ���в���
N = 16; % ��Ԫ����
d = lambda/2; % ��Ԫ���
theta0 = 0; % ���ߵ���
win = hamming(N); % ��ʸ��
k_theta0 = 2*pi*d*sind(theta0)/lambda; 
a_theta = exp(-1j*(0:1:N-1).*k_theta0).'; % ������ʸ��
h = win.*conj(a_theta);

%% �����źŲ���
Amp = 10; % �źŷ���
Pnoise =  Amp^2; % ��������
% ����Դ1
J_theta1 = 18; % �����
JSR1 = 10^(50/10); % ����Դ�źű�
k_theta1 = 2*pi*d*sind(J_theta1)/lambda; % ����Ƶ��
% ����Դ2
J_theta2 = -2;
JSR2 = 10^(30/10); 
k_theta2 = 2*pi*d*sind(J_theta2)/lambda;

%% ƥ���˲�������
kappa = Pnoise;
a_J1 = exp(-1j*(0:1:N-1)*k_theta1).'; % ����Դ1����ʸ��
a_J2 = exp(-1j*(0:1:N-1)*k_theta2).'; % ����Դ2����ʸ��
S_I1 = Amp^2*JSR1*conj(a_J1)*a_J1.'; % ����Դ1Э�������
S_I2 = Amp^2*JSR1*conj(a_J2)*a_J2.'; % ����Դ1Э�������
S_I = Pnoise*eye(N) + S_I1 + S_I2; % ����Э�������
h = kappa*inv(S_I)*conj(a_theta); % ƥ���˲���ʸ��

%% �����γ� 
theta = -90:0.5:90;% �����
y = Amp^2*exp(-1j*2*pi*(0:1:N-1).'*d*sind(theta)/lambda); % ���п���
z_theta = h.'*y; % �����γ����
z_theta_nowin = conj(a_theta).'*y;

%% display
figure;
plot(theta, db(abs(z_theta)./(Amp^2*N)), 'linewidth', 1.5)
hold on
plot([J_theta1 J_theta1], [-60 0], 'g--', 'linewidth', 1)
hold on
plot([J_theta2 J_theta2], [-60 0], 'k--', 'linewidth', 1)
xlabel('�����(\circ)')
ylabel('��һ��������Ӧ(dB)')
grid on
axis([-90 90 -60 0])
figure;
plot(theta, db(abs(z_theta_nowin)./(Amp^2*N)), 'linewidth', 1.5)
hold on
plot([J_theta1 J_theta1], [-60 0], 'g--', 'linewidth', 1)
hold on
plot([J_theta2 J_theta2], [-60 0], 'k--', 'linewidth', 1)
xlabel('�����(\circ)')
ylabel('��һ��������Ӧ(dB)')
grid on
axis([-90 90 -60 0])

