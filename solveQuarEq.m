clc;
clear;

% const
H0 = 3;
mu = 2;
ro = 3;

%vars
alpha = 2;
k = 5;
g = 4;
b0 = 2;
R_m = 5;

var1 = 0.1:0.1:20;
%vars-end

sigmaSize = length(var1);

sigma1 = zeros(1,sigmaSize);
sigma2 = zeros(1,sigmaSize);
sigma3 = zeros(1,sigmaSize);
sigma4 = zeros(1,sigmaSize);

count = 0;
for i=1:sigmaSize
    %coefficients
    b0 = var1(i);
    
    A = 1;
    B = 1j*((2*k^2)/(R_m));
    C = (g*H0*k^2) - (alpha^2) - 2*(b0^2*k^2)/(mu*ro) - (2*k^4)/(R_m^2);
    D = 1j*((g*H0*k^4 - 2*(alpha^2)*(k^2))/(R_m) + (2*(b0^2)*(k^4))/(mu*ro*R_m));
    E = ((alpha^2)*(k^4))/(R_m^2) - (g*H0*(b0^2)*(k^4))/(mu*ro) + (b0*k)^4/(mu*ro)^2;

    b_div_a = B/A;
    c_div_a = C/A;
    d_div_a = D/A;

    b_div_a_sqr = b_div_a^2;
    b_div_a_cube = b_div_a^3;

    P = c_div_a - 3/8 * b_div_a_sqr;    
    Q = d_div_a + 1/8 * b_div_a_cube - 1/2 * (B * C)/(A^2);
    R = 16/256 * b_div_a_sqr * c_div_a - 64/256 * b_div_a * d_div_a - 3/256 * (b_div_a_sqr^2) + E/A;    

    % solving of cubic equation with coefs
    A_3 = 1;
    B_3 = P / 2;
    C_3 = (P^2 - 4 * R) / 16;
    D_3 = -(Q^2) / 64;

    sqr_a_3x = 3 * (A_3^2);    
    P_3 = C_3/A_3 - (B_3^2)/sqr_a_3x;
    Q_3 = 2/27 * (B_3/A_3)^3 - (C_3*B_3)/sqr_a_3x + D_3/A_3;

    if (abs(Q_3) >= 0.0001)
        P_3_1 = (P_3 / Q_3)^2 * P_3 / 27 + 1/4;
        sqrt_Q_3 = Q_3 * sqrt(P_3_1);
    else
        sqrt_Q_3 = (P_3 / 3)^(3/2);
    end;
    
    ALPHA = (-Q_3/2 + sqrt_Q_3)^(1/3);

    if (abs(ALPHA) < 0.0001)
        BETA = 0;
    else
        BETA = -P_3/(3 * ALPHA);
    end;
        
    tmp1 = (ALPHA - BETA)/2 * 1j * sqrt(3);
    tmp2 = (ALPHA + BETA)/2;
    substr_const = B_3/(3*A_3);

    roots3 = zeros(1, 3);
    
    roots3(1) = tmp2 * 2 - substr_const;
    roots3(2) = -tmp2 + tmp1 - substr_const;
    roots3(3) = -tmp2 - tmp1 - substr_const;

    %control
    %x = roots3(3);
    %eval = A_3*x.^3+B_3*x.^2+C_3*x+D_3;
    %wrong = find(abs(eval) > 0.0001)
    %if (not(isempty(wrong)))
    %    i
    %    roots(wrong)
    %end;
    % solving cubic equation end
    
    %if (abs(roots3(1)) < 0.001)
    %    roots3(1) = 0;
    %else
    %    roots3(1) = sqrt(roots3(1));
    %end;
    
    %if (abs(roots3(2)) < 0.001)
    %    roots3(2) = 0;
    %else
    %    roots3(2) = sqrt(roots3(2));
    %end;
    
    %if (abs(roots3(3)) < 0.001)
    %    roots3(3) = 0;
    %else
    %    roots3(3) = sqrt(roots3(3));
    %end;
    
    roots3 = sqrt(roots3);

    P = -B / (4 * A);
    Q = -Q / 8;
    
    if (abs(roots3(1) * roots3(2) * roots3(3) - Q) < 0.001)
        sigma1(i) = (roots3(1) + roots3(2) + roots3(3)) + P;
        sigma2(i) = (-roots3(1) - roots3(2) + roots3(3)) + P;
        sigma3(i) = (roots3(1) - roots3(2) - roots3(3)) + P;
        sigma4(i) = (-roots3(1) + roots3(2) - roots3(3)) + P;
    else
        sigma1(i) = (roots3(1) + roots3(2) - roots3(3)) + P;
        sigma2(i) = (-roots3(1) - roots3(2) - roots3(3)) + P;
        sigma3(i) = (roots3(1) - roots3(2) + roots3(3)) + P;
        sigma4(i) = (-roots3(1) + roots3(2) + roots3(3)) + P;
    end;
    
    %control
    x = sigma1(i);
    eval = abs(A*x.^4+B*x.^3+C*x.^2+D*x+E);
    if (abs(eval) > 0.0001)
        i
        eval
        count = count + 1;
    end;
    x = sigma2(i);
    eval = abs(A*x.^4+B*x.^3+C*x.^2+D*x+E);
    if (abs(eval) > 0.0001)
        i
        eval
        count = count + 1;
    end;
    x = sigma3(i);
    eval = abs(A*x.^4+B*x.^3+C*x.^2+D*x+E);
    if (abs(eval) > 0.0001)
        i
        eval
        count = count + 1;
    end;
    x = sigma4(i);
    eval = abs(A*x.^4+B*x.^3+C*x.^2+D*x+E);
    if (abs(eval) > 0.0001)
        i
        eval
        count = count + 1;
    end;
end;

%var1(127)
count

% smooth zeros functions
sum1 = sum(real(sigma1));
if (abs(sum1) < 10^-4)
    sigma1 = imag(sigma1).*1j;
end;
sum1 = sum(imag(sigma1));
if (abs(sum1) < 10^-4)
    sigma1 = real(sigma1);
end;
sum1 = sum(real(sigma2));
if (abs(sum1) < 10^-4)
    sigma2 = imag(sigma2).*1j;
end;
sum1 = sum(imag(sigma2));
if (abs(sum1) < 10^-4)
    sigma2 = real(sigma2);
end;
sum1 = sum(real(sigma3));
if (abs(sum1) < 10^-4)
    sigma3 = imag(sigma3).*1j;
end;
sum1 = sum(imag(sigma3));
if (abs(sum1) < 10^-4)
    sigma3 = real(sigma3);
end;
sum1 = sum(real(sigma4));
if (abs(sum1) < 10^-4)
    sigma4 = imag(sigma4).*1j;
end;
sum1 = sum(imag(sigma4));
if (abs(sum1) < 10^-4)
    sigma4 = real(sigma4);
end;

% sorting of decisions
sigma = [sigma1 sigma2 sigma3 sigma4];
var2 = [var1 var1 var1 var1];

sigma_table = zeros(length(var1), 5);
for i=1:length(var1)
    sigma_table(i,1) = var1(i);
    sigma_table(i,2) = sigma1(i);
    sigma_table(i,3) = sigma2(i);
    sigma_table(i,4) = sigma3(i);
    sigma_table(i,5) = sigma4(i);
end;

sigma_table
%

subplot(1,2,1)
plot(var2, real(sigma), '.');
xlabel('b_0');
ylabel('Re(sigma)');
subplot(1,2,2)
plot(var2, imag(sigma), '.');
xlabel('b_0');
ylabel('Im(sigma)');
