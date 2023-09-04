%% Excercise 1:
%% question
% 
% <<Exercise_1.png>>
% 

%% Part A

clc
close all
clear

Taw = 10;  %in ms
I1 = 1;     I2 = 2;
T = 100;  %total simulation time
dt = 0.01; %time step

alpha = 0;      betha = 0;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = betha = 0')
legend('activity of r1', 'activity of r2')

%% Part B_a

alpha = 0.5;      betha = 0.5;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = betha = 0.5')
legend('activity of r1', 'activity of r2')

%% Part B_b

alpha = 0.5;      betha = 1;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = 0.5 & betha = 1')
legend('activity of r1', 'activity of r2')

alpha = 1;      betha = 1;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = betha = 1')
legend('activity of r1', 'activity of r2')

alpha = 0.5;      betha = -0.5;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = 0.5 & betha = -0.5')
legend('activity of r1', 'activity of r2')

alpha = 1;      betha = -1;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = 1 & betha = -1')
legend('activity of r1', 'activity of r2')

alpha = 5;      betha = -5;
[r1, r2 , t] = neuron_act(alpha, betha, Taw, I1, I2, T, dt);
figure
plot(t,r1),     hold on,        plot(t,r2),     xlabel('time (ms)'),        title('alpha = 5 & betha = -5')
legend('activity of r1', 'activity of r2')

 %% Result:
 disp('')

%%
% *When recurrent is present, different behaviors can be seen: ocilation, divergance, increasing linearly and stability.*
%
% *Behavior is varient due to different eigenvalues (lambda1 & llambda2)These values are determined by alpha and betha*
%
% *Lamda1_2 = alpha +- betha*

 %%
function [r1, r2, time] = neuron_act(alpha, betha, Taw, I1, I2, Simulation_time, time_step)
      total_data_points = Simulation_time / time_step;

      time = (0:total_data_points)*time_step;

      r1 = zeros(1, total_data_points);
      r2 = zeros(1, total_data_points);
      
      for i = 1:total_data_points
            r1(i+1) = r1(i) + (time_step/Taw) * (-r1(i) + alpha*r2(i) + I1);
            r2(i+1) = r2(i) + (time_step/Taw) * (-r2(i) + betha*r1(i) + I2);
      end
end

