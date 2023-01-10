clear all
close all
clc

% choose the number of Chebyshev polinomials and the discretization step:

N=40;
step=0.01;

% x=a^2

x=0.7:step:40; 

% determination of the numerical curve of marginal stability:

for i=1:length(x)

    [A,B]=fullmatrices_B(N,x(i));
    R=eig(A,B,'qz');

    R_pos_r=zeros(length(R),1);

% selection of real positive eigenvalues:
    
    for j=1:length(R)
        if R(j)>0 && abs(imag(R(j)))<10^(-10)
            R_pos_r(j)=R(j);
        end
    end

% ordering of positive eigenvalues:

    R_pos_r_ord=sort(R_pos_r);

% selection of the first positive eigenvalue and its square:

    flag=0;
    for j=1:length(R_pos_r_ord) 
        if R_pos_r_ord(j)>0 && flag==0
            R_x(i)=R_pos_r_ord(j);
            R_x_sqr(i)=R_x(i)^2;
            flag=1;
        end
    end

end

% plot of the numerical marginal stability curve:

figure
plot(x,R_x_sqr,'k--','LineWidth', 2)
title('Benard Problem (Chebyshev-tau)')
xlabel('a^2');
ylabel('R_{Chebyshev}^2');
legend('R_{Chebyshev}^2','Location',"north");

% determination of the critical Rayleigh number:

[M,i]=min(R_x_sqr);
wave_c=double(x(i))
R_c=R_x_sqr(i)

% comparison with the analytical curve:

syms aa
F(aa)=((aa+pi^2)^3)/aa;
Ra = matlabFunction(F);

figure(2)
plot(x,Ra(x),x,R_x_sqr,'--k','LineWidth', 2)
title('Benard Problem (camparison)')
xlabel('a^2');
ylabel('R^2');
legend('R_{Analytic}^2', 'R_{Chebyshev}^2','Location',"north");