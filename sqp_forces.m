clear all;
clc;

lb = [20 20 5 5]';
ub = [60 40 55 35]';

%lb = [40 20 0 0]';
%ub = [60 40 60 40]';

function r = g(x)
  r = [ x(3) - (2/3)*x(1);
        x(4) - x(2)/2 ];
endfunction

function r = h(x)
  r = [ x(1)-x(3);
        x(2)-x(4) ];
endfunction

function obj2 = phi(x)
  d_ad = x(1);
  L_bc = x(2);
  d_ae = x(3);
  d_bf = x(4);
  d_v = (d_ad + L_bc)/2;
  J_star = [1 0 1 0 1 0;
            0 -1 0 -1 0 -1;
            -(d_bf-(L_bc/2))/d_v -d_ae/d_v -d_bf/d_v (d_ad-d_ae)/d_v (L_bc-d_bf)/d_v (d_ad-d_ae)/d_v];
  q_c = [0 100000 0]';
  %q_c = [15000 100000 20000]';
  f_c = (J_star'*inv(J_star*(J_star)'))*q_c;
  fA = sqrt(f_c(1)^2 + f_c(2)^2);
  phiA = atan2(f_c(1),f_c(2));
  fB = sqrt(f_c(3)^2 + f_c(4)^2);
  phiB = atan2(f_c(3),f_c(4));
  fC = sqrt(f_c(5)^2 + f_c(6)^2);
  phiC = atan2(f_c(5),f_c(6));
  
  obj1 = fA + fB + fC;
  obj2 = norm(f_c);
endfunction

x0 = [30 30 20 20]';

[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [], lb, ub)

x(3)/x(1)
x(4)/x(2)