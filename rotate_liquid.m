clc;
clear;

alpha = 1;  % = 2*omega(зависит от угловой скорости вращения жидкости) - переменная
%alpha = 0:0.1:100;

k = 5;      % переменная
%k = 0:0.01:20;

R_m = inf;  % магнитное число Рейнольдса - переменная

H0 = 3;   % толщина слоя в состоянии покоя   - переменная

g = 4;    % ускорение свободного падения - переменная
%g = 0.01:0.01:5;

mu = 2;   % магнитная проницаемость - константа
ro = 3;   % плотность - константа

%b_0 = 2;% однородный фон магнитной индукции(x-компонента) - переменная
b_0 = 0:0.01:50;

con1 = k.^2;        
con2 = b_0*k; 

d = (g.*H0.*con1-alpha.^2-(2.*con2.^2)/(mu*ro)).^2 + (4.*g.*H0.*con1.*con2.^2)/(mu*ro) - (4.*con2.^4)/(mu*ro).^2;

sigma_2plus = 1/2*(alpha.^2 + (2.*con2.^2)/(mu*ro)-g.*H0.*con1+sqrt(d));
sigma_2minus = 1/2*(alpha.^2 + (2.*con2.^2)/(mu*ro)-g.*H0.*con1-sqrt(d));
sigma1 = sqrt(sigma_2plus);
sigma2 = sqrt(sigma_2minus);
sigma3 = -sqrt(sigma_2plus);
sigma4 = -sqrt(sigma_2minus);

%subplot(1,2,1)
%plot(alpha, real(sigma1));
%xlabel('alpha');
%ylabel('Re(sigma)');
%subplot(1,2,2)
%plot(alpha, imag(sigma1));
%xlabel('alpha');
%ylabel('Im(sigma)');

%subplot(1,2,1)
%plot(k, real(sigma1));
%xlabel('k');
%ylabel('Re(sigma)');
%subplot(1,2,2)
%plot(k, imag(sigma1));
%xlabel('k');
%ylabel('Im(sigma)');

%subplot(1,2,1)
%plot(g, real(sigma4));
%xlabel('g');
%ylabel('Re(sigma)');
%subplot(1,2,2)
%plot(g, imag(sigma4));
%xlabel('g');
%ylabel('Im(sigma)');

subplot(1,2,1)
plot(b_0, real(sigma4));
xlabel('b_0');
ylabel('Re(sigma)');
subplot(1,2,2)
plot(b_0, imag(sigma4));
xlabel('b_0');
ylabel('Im(sigma)');
