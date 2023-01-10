% Theta(z) in correspondence with upper temperature TU = TU(k):
funzTk = 0;
for i=1:N
   fun = Cheby_base(i)*Matrix_Eig(k,i+N)+funzTk;
   funzTk = fun;
end
Tk=matlabFunction(funzTk);