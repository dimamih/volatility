%% stock simulation
%clc;
clear;
multiflag = 1;
NP = 4; %число опционов (от 1 до 7)
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
%% --- сетки по времени и цене акции
L = 1.5; %9 у.е.
hs = 0.01; %шаг сетки в месте сгущения
T = 3; %3 месяца
%x = linspace(0,L,90+1);
x = 0:hs:L; %сетка по цене акции
%x = [0:10*hs:0.7,0.8:hs:1.1,1.2:10*hs:L]; %сетка, сгущающаяся возле s = 1
t = linspace(0,T,1200+1); %сетка по времени
r_true = 0.05;  
K = 0.9; %strike опциона
u0 = max(K-x,0); %начальное условие

%% ---настоящие данные---
tic
truesig = @(x) 0.1*sqrt(x); %настоящая функция
truevsig = truesig(x);
trueu = u_c(x,t,truevsig,r_true*ones(1,size(x,2)),u0,K); %настоящая u (создана по настоящей sig(s))
phi = s; %phi(t) на крупной сетке - реальная цена акции
t1 = 3:-ht:0;
g = zeros(NP, 1201);
for k = 1:NP
    g(k,:) = interp2(x,t,trueu,phi(1+(k-1)*200:1:1201+(k-1)*200),t1); %реальная цена опциона акции
end
disp(['Время выполнения: ',num2str(toc)])
%% ---Визуализация настоящих данных


for k= 1:NP
    %графики цены акции и опциона
    figure; 
    subplot(2,1,1);
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200), K*ones(length(t)), 'r--');
    hold on;
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200) ,phi(1+(k-1)*200:1:1201+(k-1)*200),'g', 'LineWidth', 1);
    xlabel('t - время в месяцах');
    ylabel('s - цена акции');
    title('s=\phi(t) цена акции');
    subplot(2,1,2);
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200), zeros(length(t)), 'k--');
    hold on;
    plot(t0(1+(k-1)*200:1:1201+(k-1)*200) ,g(k,:),'b', 'LineWidth', 1);
    xlabel('t - время в месяцах');
    ylabel('u - цена опциона');
    title(['g_{',num2str(k),'}(t) цена опциона']);
end



%% ---цикл---


rflag = 0; %при значении = 1 обобщение задачи на случай неизвестного r
r0 = 0.06; 

a = find(x==0.75);
b = find(x==1.04);

eps = x(2)/5000; %епсилон в приближении дельта-функции 5000
delta = @(x,phi) 1/(2*sqrt(pi*eps))*exp(-(x-phi).^2/(4*eps)); %приближение дельта-функции

G = zeros(1201,151,NP);
Q1 = zeros(1201,151,NP);
for k = 1:NP
    [x1,phi1] = meshgrid(x,phi(1201+(k-1)*200:-1:1+(k-1)*200));
    G(:,:,k) = repmat(g(k,end:-1:1)',1,size(x,2)); %G(T-t1)
    Q1(:,:,k) = 2*delta(x1,phi1); %Q=Q1.*(u-G), Q1=2*delta(s-phi(T-t1))
end
N = 5000; %число итераций цикла

J = zeros(1,N+1); %значение функционала
acc = zeros(1,N+1); %точность
sig = @(x) 0.1*sqrt(x)/sqrt(3); %начальное приближение
vsig = zeros(N+1,size(x,2));
vsig(1,:) = sig(x); %значения приближения функции sigma(s) на сетке

alpha = (0.0015/NP)*ones(1,N); 
%alpha = 0.05./sqrt(1:1:N);
disp('----------');
disp(['Значения гиперпараметров: alpha=',num2str(alpha(1)),'; eps=',num2str(eps)]);

if rflag == 1
   vr = zeros(1,N+1);
   vr(1) = r0;
   S = repmat(x,size(t,2),1);
   acc_r = zeros(1,N+1); %точность
   beta = 1.5;
   disp(['beta=',num2str(beta)]);
   r = vr(1);
else
   r = r_true;
end


tic
for i = 1:N
    u = u_c(x,t,vsig(i,:),r*ones(1,size(x,2)),u0,K); %прогноз
    g_i = zeros(NP, 1201);
    for k = 1:NP
        g_i(k,:) = interp2(x,t,u,phi(1+(k-1)*200:1:1201+(k-1)*200),t1); %реальная цена опциона акции
        J(i) = J(i) + t(2)*trapz((g_i(k,:)-g(k,:)).^2); %значение функционала
    end   
    
    acc(i) = x(2)*trapz((vsig(i,a:1:b)-truevsig(a:1:b)).^2); %точность на отрезке
    Q = zeros(1201,151);
    for k = 1:NP
        Q = Q + Q1(:,:,k).*(u-G(:,:,k)); %Q
    end
    psi = psi_c(x,t,vsig(i,:),r*ones(1,size(x,2)),Q);
    psi = psi(end:-1:1,:); %сопряженная система

    u_ss = ddif(u,x(2)); 
    vsig(i+1,:) = vsig(i,:) + alpha(i)/2*t(2)*trapz(u_ss.*psi).*x.^2; %обновление функции sigma(s)
    if rflag == 1  
        acc_r(i) = abs(vr(i)-r_true);       
        u_s = dif(u,x(2));
        vr(i+1) = vr(i) + beta*x(2)*t(2)*trapz(trapz((S.*u_s-u).*psi)); 
        r = vr(i+1);
    end
end

u = u_c(x,t,vsig(N+1,:),r*ones(1,size(x,2)),u0,K); %прогноз
for k = 1:NP
    g_i(k,:) = interp2(x,t,u,phi(1+(k-1)*200:1:1201+(k-1)*200),t1); %реальная цена опциона акции
    J(N+1) = J(N+1) + t(2)*trapz((g_i(k,:)-g(k,:)).^2); %значение функционала
end
acc(N+1) = x(2)*trapz((vsig(N+1,a:1:b)-truevsig(a:1:b)).^2); %точность
if rflag == 1  
    acc_r(N+1) = abs(vr(N+1)-r_true); 
    [~,minacc_r] = min(acc_r);
end

tau = toc;  %Время выполнения цикла  
disp(['Время выполнения цикла на ',num2str(N),' итераций: ',num2str(tau)]);

[~,minacc] = min(acc);
[~,minj] = min(J);
        


%%
if rflag ~= 1
    figure; %изменение функционалов
    subplot(2,1,1);
    plot(0:N ,J,'g', 'LineWidth', 1);
    hold on;
    xlabel('Номер итерации');
    title('Значение функционала J');
    subplot(2,1,2);
    plot(0:N ,acc,'b', 'LineWidth', 1);
    ylim([min(acc) max(acc)]);
    hold on;
    xlabel('Номер итерации');
    title('Точность');
end
if rflag == 1  
    figure; %изменение функционалов
    subplot(3,1,1);
    plot(0:N ,J,'g')%, 'LineWidth', 1);
    hold on;
    xlabel('Номер итерации');
    title('Значение функционала J');
    subplot(3,1,2);
    plot(0:N ,acc,'b')%, 'LineWidth', 1);
    ylim([min(acc) max(acc)]);
    hold on;
    xlabel('Номер итерации');
    title('Точность прогнозирования sigma(s)');
    subplot(3,1,3);
    plot(0:N ,acc_r,'b')%, 'LineWidth', 1);
    %ylim([min(acc_r) max(acc_r)]);
    hold on;
    xlabel('Номер итерации');
    title('Точность прогнозирования r(s)');
end;


%%
figure; %изменение функционалов
semilogy(0:N ,J,'g', 'LineWidth', 1);
hold on;
xlabel('Номер итерации');
title('Значение функционала J');

disp(['Значение функционала J на нулевом шаге ',num2str(J(1),'%10.2e')]);
disp(['Значение функционала J на ',num2str(N),'-м шаге ',num2str(J(N+1),'%10.2e')]);

%% сравнение настоящей функции с приближенной-лучшей по acc

a = find(x==0.75); %0.77
b = find(x==1.04); %1.02
x0 = x(a:1:b);
vminacc = vsig(minacc,a:1:b);
figure;
plot(x0,truesig(x0),'g-','LineWidth', 1);
hold on;
plot(x0,vminacc,'r--','LineWidth', 1);
xlabel('s - Цена акции');
xlim([x(a) x(b)]);
title(['Сравнение предсказанной функции \sigma_{',num2str(minacc),'}(s) с \sigma(s)'])
legend('sigma(s)',['sigma_{',num2str(minacc),'}(s)'],'Location','southeast');
hold off;
disp(['Шаг, на котором лучше всего приближается функция sigma: ',num2str(minacc)]);
disp(['Значение функционала acc на ',num2str(minacc),'-м шаге ',num2str(acc(minacc),'%10.2e')]);
disp(['Значение функционала J на ',num2str(minacc),'-м шаге ',num2str(J(minacc),'%10.2e')]);

%% сравнение настоящей функции с приближенной-лучшей по J

a = find(x==0.77);     %0.77
b = find(x==1.04);   %1.04
x0 = x(a:1:b);
vminj = vsig(minj,a:1:b);
figure;
plot(x0,truesig(x0),'g-','LineWidth', 1);
hold on;
plot(x0,vminj,'r--','LineWidth', 1);
xlabel('s - Цена акции');
xlim([x(a) x(b)]);
title(['Сравнение предсказанной функции \sigma_{',num2str(minj),'}(s) с \sigma(s)'])
legend('sigma(s)',['sigma_{',num2str(minj),'}(s)'],'Location','southeast');
hold off;
disp(['Шаг, на котором достигается наименьшее значение функционала J:',num2str(minj)]);
disp(['Значение функционала acc на ',num2str(minj),'-м шаге ',num2str(acc(minj),'%10.2e')]);
disp(['Значение функционала J на ',num2str(minj),'-м шаге ',num2str(J(minj),'%10.2e')]);

if rflag == 1
figure;
plot(x,r_true*ones(1,size(x,2)),'g-','LineWidth', 1);
hold on;
plot(x,vr(minj)*ones(1,size(x,2)),'r--','LineWidth', 1);
xlabel('s - Цена акции');
%xlim([x(a) x(b)]);
title(['Сравнение предсказанной функции r_{',num2str(minj),'}(s) с r(s)'])
legend('r(s)',['r_{',num2str(minj),'}(s)'],'Location','southeast');
hold off;
disp(['Шаг, на котором достигается наименьшее значение функционала J:',num2str(minj)]);
disp(['Значение функционала acc на ',num2str(minj),'-м шаге ',num2str(acc(minj),'%10.2e')]);
disp(['Значение функционала J на ',num2str(minj),'-м шаге ',num2str(J(minj),'%10.2e')]);
end
%% сравнение настоящей функции с приближенной-лучшей по acc_r

if rflag == 1
    a = find(x==0.75); %0.77
    b = find(x==1.04); %1.02
    x0 = x(a:1:b);
    vminacc_r = vr(minacc_r,a:1:b);
    figure;
    plot(x0,r_true(a:1:b),'g-','LineWidth', 1);
    hold on;
    plot(x0,vminacc_r,'r--','LineWidth', 1);
    xlabel('s - Цена акции');
    xlim([x(a) x(b)]);
    title(['Сравнение предсказанной функции r_{',num2str(minacc_r),'}(s) с r(s)'])
    legend('r(s)',['r_{',num2str(minacc_r),'}(s)'],'Location','southeast');
    hold off;
    disp(['Шаг, на котором лучше всего приближается функция r: ',num2str(minacc_r)]);
    disp(['Значение функционала acc на ',num2str(minacc_r),'-м шаге ',num2str(acc(minacc_r),'%10.2e')]);
    disp(['Значение функционала J на ',num2str(minacc_r),'-м шаге ',num2str(J(minacc_r),'%10.2e')]);
    disp(['Значение функционала acc_r на ',num2str(minacc_r),'-м шаге ',num2str(acc_r(minacc_r),'%10.2e')]);
end
%% анимация (сравнение предсказанных функций с настоящей)
figure;
vid = VideoWriter('movie_multi.mp4','MPEG-4');
vid.FrameRate = 50;
open(vid);
mov(1:N+1) = struct('cdata', [], 'colormap', []);
for i=1:N+1
    plot(x,truesig(x),'g-','LineWidth', 1);
    hold on;
    plot(x,vsig(i,:),'r--','LineWidth', 1);
    xlabel('s - Цена акции');
    title(['Сравнение предсказанной функции \sigma_{',num2str(i-1),'}(s) с \sigma(s)'])
    hold off;
    %pause(0.1);
    mov(i) = getframe(gcf);
    writeVideo(vid,mov(i));
    
end;
close(vid);

