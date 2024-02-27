clear all;
clc;

lb = [40 20 0 0]';
ub = [60 40 60 40]';

%lb = [20 20 5 5]';
%ub = [60 40 55 35]';

function r = h(x)
  r = [ x(1)-x(3);
        x(2)-x(4) ];
endfunction

function obj1 = phi(x)
  d_ad = x(1);
  L_bc = x(2);
  d_ae = x(3);
  d_bf = x(4);
  d_v = (d_ad + L_bc)/2;
  N = (11*(d_ad)^2) - (16*d_ad*d_ae) + (12*(d_ae)^2) + (6*d_ad*L_bc) + ...
      (8*(L_bc)^2) - (12*L_bc*d_bf) + (12*(d_bf)^2);
  K = (N^2) - (8*((d_ad+L_bc)^2)*((4*(d_ad^2)) + (3*(L_bc^2))));
  sigma1 = sqrt(N-sqrt(K))/2*d_v*sqrt(2);
  sigma3 = sqrt(N+sqrt(K))/2*d_v*sqrt(2);
  obj1 = sigma3/sigma1;
endfunction

x0 = [60 30 10 20]';
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, [], @h, lb, ub)