%% matan 1

syms x y z
u = y*x/(z^2*(2*y-3*z));
p = [3, -1, -1];
h = [3, -2, -2];

grad_u = gradient(u, [x, y, z]);
grad_u_at_p = subs(grad_u, [x, y, z], p);
dpf_h = dot(grad_u_at_p, h);
disp(dpf_h);

%% matan 2

% Определение переменных
syms x y z;
u = (x/y^3) - 3*y/z^3;
p = [-1, 1, 2];
h = [1, -2, 2];

% Вычисление первых производных
u_x = diff(u, x);
u_y = diff(u, y);
u_z = diff(u, z);

% Вычисление вторых производных
u_xx = diff(u_x, x);
u_yy = diff(u_y, y);
u_zz = diff(u_z, z);
u_xy = diff(u_x, y);
u_xz = diff(u_x, z);
u_yz = diff(u_y, z);

% Создание матрицы Гессе
Hessian = [u_xx u_xy u_xz; u_xy u_yy u_yz; u_xz u_yz u_zz];

% Подстановка значений p в матрицу Гессе
Hessian_p = subs(Hessian, [x y z], p);

% Вычисление d²pf(h)
d2pf_h = h * Hessian_p * h';

% Вывод результата
disp(d2pf_h);

%% matan 3
syms x y;
u = (y^4)/(x^4);
p = [1, -2];x
h = [1, -3];

% Вычисляем частные производные
du_dx = diff(u, x);
du_dy = diff(u, y);

% Вычисляем вторые производные
d2u_dx2 = diff(du_dx, x);
d2u_dy2 = diff(du_dy, y);
d2u_dxdy = diff(du_dx, y);

% Вычисляем третьи производные
d3u_dx3 = diff(d2u_dx2, x);
d3u_dy3 = diff(d2u_dy2, y);
d3u_dx2dy = diff(d2u_dx2, y);
d3u_dxdy2 = diff(d2u_dxdy, y);
d3u_dy2dx = diff(d2u_dy2, x);
d3u_dxdydx = diff(d2u_dxdy, x);

% Вычисляем значения производных в точке p
d3u_dx3_val = subs(d3u_dx3, {x, y}, p);
d3u_dy3_val = subs(d3u_dy3, {x, y}, p);
d3u_dx2dy_val = subs(d3u_dx2dy, {x, y}, p);
d3u_dxdy2_val = subs(d3u_dxdy2, {x, y}, p);
d3u_dy2dx_val = subs(d3u_dy2dx, {x, y}, p);
d3u_dxdydx_val = subs(d3u_dxdydx, {x, y}, p);

% Вычисляем d3pf(h)
d3pf_h = h(1)^3 * d3u_dx3_val + 3 * h(1)^2 * h(2) * d3u_dx2dy_val + 3 * h(1) * h(2)^2 * d3u_dxdy2_val + h(2)^3 * d3u_dy3_val;

disp(d3pf_h);
