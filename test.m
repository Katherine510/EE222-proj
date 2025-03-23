sympref('FloatingPointOutput',true)
g = 9.81;
rg = 0.0254;
L = 0.4255;
K = 1.5;
t = 0.025;


%% get jacobian
syms x1 x2 x3 x4;
nonLeqn = [ x2; (5/7)*(rg/L)*g*sin(x3) - (5/7)*(rg/L)^2 * x4^4 * (cos(x3))^2; x4; -(x4/t)];
A  = jacobian(nonLeqn, [x1, x2, x3, x4]);
B = [0; 0; 0; K/t];
disp(A)
disp(B)


%% Calculate LQR matrix

A = double(subs(A, [x1, x2, x3, x4], [1, 1, 1, 1]));


Q = [1,0,0,0;
     0,0,0,0;
     0,0,1,0;
     0,0,0,0];
R = 1;


[K,S,P] = lqr(A,B,Q,R);