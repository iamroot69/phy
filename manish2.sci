clc;clear;clf();
function y = z(e, g,T)
    k = 8.617333262D-5;
    b = 1./ (k * T);
    y = sum(exp(-b*e)*g',2);  // Vectorized operation
endfunction

e = [0, 1,2,3,4];
g=[1, 2,3,4,4];
T = (0:1:100000)';
k = 8.617333262D-5;
b = 1./ (k * T);

subplot(2,2,1)
pl = 1./ z(e,g, T);
pu = exp(-b)./ z(e,g, T);
plot(T, pl, T, pu);

subplot(2,2,2)
E = exp(-b)./ z(e,g, T);
plot(k * T, E);

subplot(2,2,3)
s = log(z(e,g, T)) + b.*E;
plot(k * T, s);
