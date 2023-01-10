function[coeff]=CoeffExtraction(t,order)
%%%% t is a symbolic polynomial function
%%%% ORDER is the degree of polynomial
syms z
for i = 1:order
    coeff(i)=subs(t,z,0);
    syms z
    t(z)=simplify((t-coeff(i))/z);
end
end