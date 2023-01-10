clear all
close all
clc

% choose the number of Chebyshev polynomials, and TU interval:

N=30;
TU=4:0.1:8;

% initialize critical Rayleigh number vector:

RcritL=zeros(1,length(TU));

% x=a^2

x = 0.4:0.01:20; 

for k=1:length(TU)

for i=1:length(x)

    [A,B]=fullmatrices_PClear(N,x(i),TU(k));
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

% determination of the critical Rayleigh number:

[M,i]=min(R_x_sqr);
R_c=R_x_sqr(i);

RcritL(k)=R_c;

end

% plot of the behaviour curve:

figure
plot(TU,RcritL,'--k','LineWidth', 2)
title('Behaviour of R_c with respect to T_U')
xlabel('T_U');
ylabel('R_c');
legend('R_{c,L}','Location',"north");
