function [t, v, y, h] = bungee_project_RK4(a, b, n)
% bungee_project_RK4 Bungee Project Model Using Fourth Order Runge-Kutta
% Method
% [t, v, y, h] = bungee_project_RK4(f, a, b, alpha, n) computes the bungee
% jumping using RK4 taking n steps from t = a to t = b with v(a) = alpha.
% let y(1) = 0 and v(1) = 0 which are the distance and velocity of the
% jumper at the platform (y0 & v0). 

% declare the constants
c = 0.9; % drag coefficient in kg/m
m = 80; % mass of jumper in kg
k = 90; % spring constant of bungee rope in N/m (originally 90)
L = 25; % length of bungee rope in m (originally 25)
g = 9.8; % gravitational acceleration in m/s^2
C = c/m;
K = k/m;
f1 = @(t,v) v;
f2 = @(t,y,v) g - C*abs(v)*v - max(0, K*(y - L));

% calculate h
h = (b-a)/n;

% create t array
t = a:h:b;

% initialise v & y array
v = zeros(size(t));
v(1) = 0;
y = zeros(size(t));
y(1) = 0;

% perform the iterations
for i = 1:n
    k1_v = h * f2(t(i), y(i), v(i)); 
    k1_y = h * f1(t(i), v(i));
    
    k2_v = h * f2(t(i) + (h/2), y(i) + (k1_y/2), v(i) + (k1_v/2));
    k2_y = h * f1(t(i) + (h/2), v(i) + (k1_v/2));
    
    k3_v = h * f2(t(i) + (h/2), y(i) + (k2_y/2), v(i) + (k2_v/2));
    k3_y = h * f1(t(i) + (h/2), v(i) + (k2_v/2));
    
    k4_v = h * f2(t(i) + h, y(i) + k3_y, v(i) + k3_v);
    k4_y = h * f1(t(i) + h, v(i) + k3_v);
    
    v(i+1) = v(i) + (1/6 * (k1_v + (2*k2_v) + (2*k3_v) + k4_v));  
    y(i+1) = y(i) + (1/6 * (k1_y + (2*k2_y) + (2*k3_y) + k4_y)); 
end

end