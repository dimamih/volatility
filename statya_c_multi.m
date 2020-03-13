%% stock simulation
%clc;
clear;
multiflag = 1;
NP = 4; %����� �������� (�� 1 �� 7)
rng('default');
%rng(59947);
rng(357);
if multiflag == 1
    T = 6; %months
else
    T = 3;
end
%sigma = 0.1;
sig = @(t) 0.1*sqrt(t);
%L = 9; Ns = 900; hs = 0.01;
K = 0.9; %strike
r = 0.05; 
ht = 0.0025;
% t0 = linspace(0,T,Nt+1);
t0 = 0:ht:T;
% Nt = 1200; 
Nt = size(t0,2)-1;
s = ones(1, Nt+1);
xi = rand(1, Nt);

for k = 2: Nt+1
    d = exp(-sig(s(k-1))*sqrt(ht));
    u = exp(sig(s(k-1))*sqrt(ht));
    p = (u-1-r*ht)/(u-d);
    if xi(k-1) < p
        s(k) = d * s(k-1);
    else
       s(k) = u * s(k-1); 
    end
%     if xi(k-1) > 1 - p
%         s(k) = u * s(k-1);
%     end
end
%%
figure;
hold on;
grid on;
plot(t0 ,s,'g', 'LineWidth', 1);
xlabel('Time');
ylabel('Stock price');
title('s=\phi(t) curve');
hold off;

% figure;
% hold on;
% grid on;
% histogram(s, 'Normalization', 'count')
% title('s=\phi(t) histogram');
% hold off;
%% --- ����� �� ������� � ���� �����
L = 1.5; %9 �.�.
hs = 0.01; %��� ����� � ����� ��������
T = 3; %3 ������
%x = linspace(0,L,90+1);
x = 0:hs:L; %����� �� ���� �����
%x = [0:10*hs:0.7,0.8:hs:1.1,1.2:10*hs:L]; %�����, ����������� ����� s = 1
t = linspace(0,T,1200+1); %����� �� �������
r_true = 0.05;  
K = 0.9; %strike �������
u0 = max(K-x,0); %��������� �������

%% ---��������� ������---
tic
truesig = @(x) 0.1*sqrt(x); %��������� �������
truevsig = truesig(x);
trueu = u_c(x,t,truevsig,r_true*ones(1,size(x,2)),u0,K); %��������� u (������� �� ��������� sig(s))
phi = s; %phi(t) �� ������� ����� - �������� ���� �����
t1 = 3:-ht:0;
g = zeros(NP, 1201);
for k = 1:NP
    g(k,:) = interp2(x,t,trueu,phi(1+(k-1)*200:1:1201+(k-1)*200),t1); %�������� ���� ������� �����
end
disp(['����� ����������: ',num2str(toc)])
%% ---������������ ��������� ������


for k= 1:NP
    %������� ���� ����� � �������
    figure; 
    subplot(2,1,1);
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200), K*ones(length(t)), 'r--');
    hold on;
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200) ,phi(1+(k-1)*200:1:1201+(k-1)*200),'g', 'LineWidth', 1);
    xlabel('t - ����� � �������');
    ylabel('s - ���� �����');
    title('s=\phi(t) ���� �����');
    subplot(2,1,2);
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200), zeros(length(t)), 'k--');
    hold on;
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200) ,g(k,:),'b', 'LineWidth', 1);
    xlabel('t - ����� � �������');
    ylabel('u - ���� �������');
    title(['g_{',num2str(k),'}(t) ���� �������']);
end



%% ---����---


rflag = 0; %��� �������� = 1 ��������� ������ �� ������ ������������ r
r0 = 0.06; 

a = find(x==0.75);
b = find(x==1.04);

eps = x(2)/5000; %������� � ����������� ������-������� 5000
delta = @(x,phi) 1/(2*sqrt(pi*eps))*exp(-(x-phi).^2/(4*eps)); %����������� ������-�������

G = zeros(1201,151,NP);
Q1 = zeros(1201,151,NP);
for k = 1:NP
    [x1,phi1] = meshgrid(x,phi(1201+(k-1)*200:-1:1+(k-1)*200));
    G(:,:,k) = repmat(g(k,end:-1:1)',1,size(x,2)); %G(T-t1)
    Q1(:,:,k) = 2*delta(x1,phi1); %Q=Q1.*(u-G), Q1=2*delta(s-phi(T-t1))
end
N = 5000; %����� �������� �����

J = zeros(1,N+1); %�������� �����������
acc = zeros(1,N+1); %��������
sig = @(x) 0.1*sqrt(x)/sqrt(3); %��������� �����������
vsig = zeros(N+1,size(x,2));
vsig(1,:) = sig(x); %�������� ����������� ������� sigma(s) �� �����

alpha = (0.0015/NP)*ones(1,N); 
%alpha = 0.05./sqrt(1:1:N);
disp('----------');
disp(['�������� ���������������: alpha=',num2str(alpha(1)),'; eps=',num2str(eps)]);

if rflag == 1
   vr = zeros(1,N+1);
   vr(1) = r0;
   S = repmat(x,size(t,2),1);
   acc_r = zeros(1,N+1); %��������
   beta = 1.5;
   disp(['beta=',num2str(beta)]);
   r = vr(1);
else
   r = r_true;
end


tic
for i = 1:N
    u = u_c(x,t,vsig(i,:),r*ones(1,size(x,2)),u0,K); %�������
    g_i = zeros(NP, 1201);
    for k = 1:NP
        g_i(k,:) = interp2(x,t,u,phi(1+(k-1)*200:1:1201+(k-1)*200),t1); %�������� ���� ������� �����
        J(i) = J(i) + t(2)*trapz((g_i(k,:)-g(k,:)).^2); %�������� �����������
    end   
    
    acc(i) = x(2)*trapz((vsig(i,a:1:b)-truevsig(a:1:b)).^2); %�������� �� �������
    Q = zeros(1201,151);
    for k = 1:NP
        Q = Q + Q1(:,:,k).*(u-G(:,:,k)); %Q
    end
    psi = psi_c(x,t,vsig(i,:),r*ones(1,size(x,2)),Q);
    psi = psi(end:-1:1,:); %����������� �������

    u_ss = ddif(u,x(2)); 
    vsig(i+1,:) = vsig(i,:) + alpha(i)/2*t(2)*trapz(u_ss.*psi).*x.^2; %���������� ������� sigma(s)
    if rflag == 1  
        acc_r(i) = abs(vr(i)-r_true);       
        u_s = dif(u,x(2));
        vr(i+1) = vr(i) + beta*x(2)*t(2)*trapz(trapz((S.*u_s-u).*psi)); 
        r = vr(i+1);
    end
end

u = u_c(x,t,vsig(N+1,:),r*ones(1,size(x,2)),u0,K); %�������
for k = 1:NP
    g_i(k,:) = interp2(x,t,u,phi(1+(k-1)*200:1:1201+(k-1)*200),t1); %�������� ���� ������� �����
    J(N+1) = J(N+1) + t(2)*trapz((g_i(k,:)-g(k,:)).^2); %�������� �����������
end
acc(N+1) = x(2)*trapz((vsig(N+1,a:1:b)-truevsig(a:1:b)).^2); %��������
if rflag == 1  
    acc_r(N+1) = abs(vr(N+1)-r_true); 
    [~,minacc_r] = min(acc_r);
end

tau = toc;  %����� ���������� �����  
disp(['����� ���������� ����� �� ',num2str(N),' ��������: ',num2str(tau)]);

[~,minacc] = min(acc);
[~,minj] = min(J);
        


%%
if rflag ~= 1
    figure; %��������� ������������
    subplot(2,1,1);
    plot(0:N ,J,'g', 'LineWidth', 1);
    hold on;
    xlabel('����� ��������');
    title('�������� ����������� J');
    subplot(2,1,2);
    plot(0:N ,acc,'b', 'LineWidth', 1);
    ylim([min(acc) max(acc)]);
    hold on;
    xlabel('����� ��������');
    title('��������');
end
if rflag == 1  
    figure; %��������� ������������
    subplot(3,1,1);
    plot(0:N ,J,'g')%, 'LineWidth', 1);
    hold on;
    xlabel('����� ��������');
    title('�������� ����������� J');
    subplot(3,1,2);
    plot(0:N ,acc,'b')%, 'LineWidth', 1);
    ylim([min(acc) max(acc)]);
    hold on;
    xlabel('����� ��������');
    title('�������� ��������������� sigma(s)');
    subplot(3,1,3);
    plot(0:N ,acc_r,'b')%, 'LineWidth', 1);
    %ylim([min(acc_r) max(acc_r)]);
    hold on;
    xlabel('����� ��������');
    title('�������� ��������������� r(s)');
end;


%%
figure; %��������� ������������
semilogy(0:N ,J,'g', 'LineWidth', 1);
hold on;
xlabel('����� ��������');
title('�������� ����������� J');

disp(['�������� ����������� J �� ������� ���� ',num2str(J(1),'%10.2e')]);
disp(['�������� ����������� J �� ',num2str(N),'-� ���� ',num2str(J(N+1),'%10.2e')]);

%% ��������� ��������� ������� � ������������-������ �� acc

a = find(x==0.75); %0.77
b = find(x==1.04); %1.02
x0 = x(a:1:b);
vminacc = vsig(minacc,a:1:b);
figure;
plot(x0,truesig(x0),'g-','LineWidth', 1);
hold on;
plot(x0,vminacc,'r--','LineWidth', 1);
xlabel('s - ���� �����');
xlim([x(a) x(b)]);
title(['��������� ������������� ������� \sigma_{',num2str(minacc),'}(s) � \sigma(s)'])
legend('sigma(s)',['sigma_{',num2str(minacc),'}(s)'],'Location','southeast');
hold off;
disp(['���, �� ������� ����� ����� ������������ ������� sigma: ',num2str(minacc)]);
disp(['�������� ����������� acc �� ',num2str(minacc),'-� ���� ',num2str(acc(minacc),'%10.2e')]);
disp(['�������� ����������� J �� ',num2str(minacc),'-� ���� ',num2str(J(minacc),'%10.2e')]);

%% ��������� ��������� ������� � ������������-������ �� J

a = find(x==0.77);     %0.77
b = find(x==1.04);   %1.04
x0 = x(a:1:b);
vminj = vsig(minj,a:1:b);
figure;
plot(x0,truesig(x0),'g-','LineWidth', 1);
hold on;
plot(x0,vminj,'r--','LineWidth', 1);
xlabel('s - ���� �����');
xlim([x(a) x(b)]);
title(['��������� ������������� ������� \sigma_{',num2str(minj),'}(s) � \sigma(s)'])
legend('sigma(s)',['sigma_{',num2str(minj),'}(s)'],'Location','southeast');
hold off;
disp(['���, �� ������� ����������� ���������� �������� ����������� J:',num2str(minj)]);
disp(['�������� ����������� acc �� ',num2str(minj),'-� ���� ',num2str(acc(minj),'%10.2e')]);
disp(['�������� ����������� J �� ',num2str(minj),'-� ���� ',num2str(J(minj),'%10.2e')]);

if rflag == 1
figure;
plot(x,r_true*ones(1,size(x,2)),'g-','LineWidth', 1);
hold on;
plot(x,vr(minj)*ones(1,size(x,2)),'r--','LineWidth', 1);
xlabel('s - ���� �����');
%xlim([x(a) x(b)]);
title(['��������� ������������� ������� r_{',num2str(minj),'}(s) � r(s)'])
legend('r(s)',['r_{',num2str(minj),'}(s)'],'Location','southeast');
hold off;
disp(['���, �� ������� ����������� ���������� �������� ����������� J:',num2str(minj)]);
disp(['�������� ����������� acc �� ',num2str(minj),'-� ���� ',num2str(acc(minj),'%10.2e')]);
disp(['�������� ����������� J �� ',num2str(minj),'-� ���� ',num2str(J(minj),'%10.2e')]);
end
%% ��������� ��������� ������� � ������������-������ �� acc_r

if rflag == 1
    a = find(x==0.75); %0.77
    b = find(x==1.04); %1.02
    x0 = x(a:1:b);
    vminacc_r = vr(minacc_r,a:1:b);
    figure;
    plot(x0,r_true(a:1:b),'g-','LineWidth', 1);
    hold on;
    plot(x0,vminacc_r,'r--','LineWidth', 1);
    xlabel('s - ���� �����');
    xlim([x(a) x(b)]);
    title(['��������� ������������� ������� r_{',num2str(minacc_r),'}(s) � r(s)'])
    legend('r(s)',['r_{',num2str(minacc_r),'}(s)'],'Location','southeast');
    hold off;
    disp(['���, �� ������� ����� ����� ������������ ������� r: ',num2str(minacc_r)]);
    disp(['�������� ����������� acc �� ',num2str(minacc_r),'-� ���� ',num2str(acc(minacc_r),'%10.2e')]);
    disp(['�������� ����������� J �� ',num2str(minacc_r),'-� ���� ',num2str(J(minacc_r),'%10.2e')]);
    disp(['�������� ����������� acc_r �� ',num2str(minacc_r),'-� ���� ',num2str(acc_r(minacc_r),'%10.2e')]);
end
%% �������� (��������� ������������� ������� � ���������)
figure;
vid = VideoWriter('movie_multi.mp4','MPEG-4');
vid.FrameRate = 50;
open(vid);
mov(1:N+1) = struct('cdata', [], 'colormap', []);
for i=1:N+1
    plot(x,truesig(x),'g-','LineWidth', 1);
    hold on;
    plot(x,vsig(i,:),'r--','LineWidth', 1);
    xlabel('s - ���� �����');
    title(['��������� ������������� ������� \sigma_{',num2str(i-1),'}(s) � \sigma(s)'])
    hold off;
    %pause(0.1);
    mov(i) = getframe(gcf);
    writeVideo(vid,mov(i));
    
end;
close(vid);

