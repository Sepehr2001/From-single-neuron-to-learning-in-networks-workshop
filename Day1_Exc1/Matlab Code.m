%% Excercise 1
%% question
% 
% <<Excercise1.png>>
% 

clc
clear
close all

E_L = -70; % mV
E_x = 0;
g_L = 1; % \micro S / mm^2
c_m = 10; % \nano F / mm^2
V_th = -54; % mV
V_reset = -80; % mV

tau_exc = 10; % ms  
Delta_g_exc = 0.5; % \micro S / mm^2

tot_data_points = 100000;    %total number of datapoints
dt = 0.01; %time step in ms
t = (0:tot_data_points) * dt;         %time in milli seconds = number of data points * step time


Presyn_SpikeTimes= [100 200 230 300 320 400 410] / dt;   %time/step time = data point number

V = zeros (1, tot_data_points) ; % membrane potential
g_exc = zeros (1, tot_data_points); 
X = zeros (1, tot_data_points); % Spikes
I_exc = zeros (1, tot_data_points); % Excitatory currents

V(1) = V_reset;
    
for i = 1 : tot_data_points
           
    V(i + 1) = V( i ) - dt/c_m * ( g_L * ( V( i ) - E_L ) + g_exc ( i ) * ( V( i ) - E_x));  %Euler method
 
    g_exc(i + 1) = g_exc ( i ) - (dt / tau_exc)  * g_exc ( i );                                           %Euler method
    
    if ismember(i+1, Presyn_SpikeTimes)
          g_exc (i+1) = g_exc (i+1) + Delta_g_exc ;
    end
    
    I_exc (i + 1) = g_exc (i + 1) * ( V (i + 1) - E_x );
    
    if V ( i + 1 ) >= V_th
         X(i + 1) = 1 ;
         V(i + 1) = V_reset;
    end 

end

n_spikes = sum(X);
spike_times = find(X==1)*dt;  % in milli seconds
figure
subplot (211)
plot (t, V)
title ( ' membrane potential versus time ' )
xlabel ( ' time in milli seconds ' )
ylabel ( ' V ' )

subplot(212)
plot( t , I_exc )
title ( ' Excitatory current ' )
xlabel ( ' time in milli seconds ' )
ylabel ( ' I_e_x_c ' )

fprintf('Number of spikes \t Time of spiking(in milli seconds)\n')

for i =1:length(spike_times)
        fprintf('\t%d\t\t\t   %.2f\n',i,spike_times(i))
end
%% 
% * The spike time was 412ms. because there are 2 EPSPs occuring at 400 & 410  which are so close together and this caused a spike! 
% * This is a good presentation of the integration in the LIF model.
% * We can also change the presynaptic times to generate more spikes as well.