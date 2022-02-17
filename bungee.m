format compact;

a = 0; b= 60; n = 1000;

% RUN RK4 BUNGEE JUMPING MODEL
[t, v, y, h] = bungee_project_RK4(a, b, n);

% PLOT POSITION VS TIME GRAPH
figure(1);
plot(t, y);
title('Jumper''s Position Against Time');
xlabel('time (s)');
ylabel('position (m)');

% 4.1 TIMING AND BOUNCES
bounces = 0; % create new variable to store number of bounces

% start for loop to count bounces which starts from i = 3
for i = 3:n
    % compare y(i) with y(i-1) and y(i-2) to check for bounces
    if(y(i - 2) < y(i - 1) && y(i - 1) > y(i))
        bounces = bounces + 1; % add 1 to bounces if a new bounce is found
    end
end

% print the results
fprintf('There are %d "bounces" in %d seconds\n', bounces, b);

% 4.2 MAXIMUM SPEED EXPERIENCED BY JUMPER
% plot velocity vs time graph
figure(2);
plot(t, v);
title('Jumper''s Velocity Against Time');
xlabel('time (s)');
ylabel('velocity (m/s)');

max_speed = 0; % create new variable to store the maximum speed

% start for loop to find what the maximum speed is
for i = 1:n
    if(abs(v(i)) > max_speed)
        max_speed = abs(v(i)); % set current velocity as max_speed if its larger
    end
end

max_speed_index = find(v == max_speed); % find the index of the max speed
max_speed_time = t(max_speed_index); % find the time the max speed occurs

% print the results
fprintf('The maximum speed experienced by the jumper is %.4f m/s and it occurs at %.4f seconds\n', max_speed, max_speed_time);

% 4.3 MAXIMUM ACCELERATION EXPERIENCED BY THE JUMPER
h = (b-a)/n; % the size between points

% calculate the acceleration for each point using first order forward difference method
for i = 1:n
    acc1(i) = (v(i + 1) - v(i)) ./ h; % finding acceleration for acc1(i)
end

% fill the last point because we do not have n+1 points
acc1(n+1) = acc1(n);

% plot acceleration vs time graph
figure(3);
plot(t, acc1);
xlabel('time (s)');
ylabel('acceleration (m/s^2)');
title('Jumper''s Acceleration Against Time');
hold on;

% find coordinates of the largest point
[max_acceleration position_y] = max(abs(acc1));

% finding time using the y position
time = t(position_y);

% plot the maximum acceleration
hold on;
plot(time, acc1(position_y), 'ro');
legend('acceleration vs time', 'max acceleration');

% print the results
fprintf('The maximum acceleration is %.4f meters per second square and occurs at %.2f seconds\n',max_acceleration,time);

% 4.4 DISTANCE TRAVELLED BY THE JUMPER
S = 0; % initialising S

% Trapezoidal Rule, summing the values from 2 to n
for j = 2:n
    S = S + abs(v(j)); % add velocity at j (v(j)) to the sum
end

% calculate the integration
I = h/2 * (abs(v(1)) + 2*S + abs(v(n+1)));

% print the results
fprintf('The distance travelled by the jumper is %.4f meters\n', I);

% 4.5 AUTOMATED CAMERA SYSTEM
H = 74; D = 31; % declare the height of jump point and deck height

H_to_D = H - D; % calculate H-D

% set yi and yi+1 to be the min y and y1+2 and y1+3 to be the max y
y_1 = min(y); y_2 = min(y); y_3 = max(y); y_4 = max(y);

% use for loop to find the four closest value to H-D
for i = 1:n
    if round(y(i)) < H_to_D % to find 2 values which are less than H-D
        if y(i) > y_2
            y_1 = y_2;
            y_2 = y(i);
        end
    else
        if round(y(i)) > H_to_D % to find 2 values which are more than H-D
            if y(i) < y_4
                y_4 = y_3;
                y_3 = y(i);
            end
        end
    end
end

% find the index of y_1 to y_4
index_1 = find(y == y_1); index_2 = find(y == y_2); index_3 = find(y == y_3); index_4 = find(y == y_4);

% find the time when y_1, y_2, y_3 and y_4 occurs
t_1 = t(index_1); t_2 = t(index_2); t_3 = t(index_3); t_4 = t(index_4);

Y = [y_1 y_2 y_3 y_4]; % create Y array to store the four y values
T = [t_1 t_2 t_3 t_4]; % create T array to store the corresponding t for Y

% use the bisection method to find an approximation for t
for i = 1:n
    p = (a+b)/2; % compute the midpoint value
    if sign(divided_difference_q5(T, Y, a) - H_to_D) == sign(divided_difference_q5(T, Y, p) - H_to_D)
        a = p;
    else
        b = p;
    end
end

% print the results
fprintf('The camera should trigger to capture an image of the jumper at %.4f seconds\n', p);

% 4.6 WATER TOUCH OPTION
a = 0; b = 60; % reset the values of a and b
fprintf('\nBungee jumping model with modified spring constant & rope length:\n');
% run the RK4 bungee jumping model where the spring constant and rope
% length has been modified
k = 90; % spring constant of bungee rope in N/m
L = 25; % length of bungee rope in m
g = 9.8; % gravitational acceleration in m/s^2
target_height = 74 - 1.75;
max_height = max(abs(y));
accel_limit = 2 * g;
max_accel = max(abs(acc1));

%set incrementing factor
inc_factor = 0.1;
run_count = 0;
loop_limit = 10000;

k_initial = 80;
L_initial = 0;

L_max = 74;

tol = 0.01;

while run_count < loop_limit
    if L < L_max
        L = L + inc_factor;
    end
    %reset k value for testing
    k = k_initial;
    
    % Recalculate
    [t, v, y, h, acc1] = bungee_project_RK4_modified(a, b, n, k, L);
    max_height = max(abs(y));
    max_accel = max(abs(acc1));
    
    
   	if max_height < (target_height + tol) && max_height > (target_height - tol) && max_accel < accel_limit
        break
    end    
      
    for counter = 1:loop_limit
        k = k + inc_factor;
        % Recalculate
        [t, v, y, h, acc1] = bungee_project_RK4_modified(a, b, n, k, L);
        max_height = max(abs(y));
        max_accel = max(abs(acc1));
        
        if max_height < (target_height + tol) && max_height > (target_height - tol) || max_accel > accel_limit
            break
        end
    end
    
    if max_height < (target_height + tol) && max_height > (target_height - tol) && max_accel < accel_limit
        break
    end
    
   run_count = run_count + 1;
   if run_count == loop_limit
       fprintf('loop limit reached\n');
   end
       
end

fprintf('Current max height is %.2f\n', max_height);
fprintf('Current max acceleration is %.4f\n', max_accel);
fprintf('Rope Length L = %.2f\n', L);
fprintf('Spring constant k = %.2f\n', k);

%plotting
figure(4);
plot(t,acc1);
xlabel('time (s)');
ylabel('acceleration (m/s^2)');
title({'Jumper''s Acceleration Against Time using Forward Difference Method'; '"Water Touch" Option'});


% finding the coordinates of the largest point
[max_acceleration, position_y] = max(abs(acc1));

% finding time using the y position
time = t(position_y);
fprintf('The maximum acceleration is %.4f meters per second and occurs at %.2f seconds\n',max_acceleration,time)

% Check to see if result should be discounted because g-force threshold
% exceeded

% Display Water Touch Scenario
figure(5);
plot(t, y);
yline(72.25, '--', 'Required point for Water Touch - 72.25m');
title("Height Displacement over Time");
xlabel("Time");
ylabel("Height Displacement");