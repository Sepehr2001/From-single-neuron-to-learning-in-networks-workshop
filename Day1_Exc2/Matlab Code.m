%% Excercise 2:
%% question
%%
% 
% <<Excercise2.png>>
% 


%% Part a

close all
clear
clc

Taw_m = 10; %in milli seconds
V_th = -54;  %in milli volts
V_reset = -80;  %in milli volts

tot_data_points = 3000000;
dt = 0.1; %time step in ms

sigma_v = 0:0.5:10;
std_of_V = zeros(1,length(sigma_v));
V = zeros(1, tot_data_points);
V(1) = V_reset;
for j = 1:length(sigma_v)
      Eff = -56 + sigma_v(j)*sqrt(2*Taw_m/dt)*randn(1, tot_data_points);
      for i = 1:tot_data_points
            V(i+1) = V(i) + (dt/Taw_m)*(-V(i) + Eff(i));
      end
      std_of_V(j) = std(V);
      V(1) = V_reset;
end
figure
plot(sigma_v, std_of_V),    xlabel('sigma of V'),    ylabel('std of V')
%%
% *we can see that without implementing the threshold, the response of the neuron is linear* 

%% Part b
close all
clear
clc

Taw_m = 10; %in milli seconds
V_th = -54;  %in milli volts
V_reset = -80;  %in milli volts

tot_data_points = 3000000;
dt = 0.1; %time step in ms

sigma_v = 0:0.5:10;
std_of_V = zeros(1,length(sigma_v));
V = zeros(1, tot_data_points);
V(1) = V_reset;
a=0;  %this represents the number of spikes

for j = 1:length(sigma_v)
      Eff = -56 + sigma_v(j)*sqrt(2*Taw_m/dt)*randn(1, tot_data_points);
      for i = 1:tot_data_points
            V(i+1) = V(i) + (dt/Taw_m)*(-V(i) + Eff(i));
            if V(i+1) > V_th
                V(i+1) = V_reset;
                a = a+1;
            end
      end
      std_of_V(j) = a/(tot_data_points * dt/1000);
      % std_of_V(j) = std(V);
      a=0;
      V(1) = V_reset;
end
figure
plot(sigma_v, std_of_V),       xlabel('sigma of V'),    ylabel('Average of firing rate')
hold on
plot(sigma_v, sigma_v)

%%
% *But, with adding the threshold and reseting to our model, the response beacme nonlinear.*
 
