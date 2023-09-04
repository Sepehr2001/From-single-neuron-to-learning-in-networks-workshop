%% Excercise 2:
%% question
% 
% <<Exercise_2.png>>
% 

%% Answer
clc
clear
close all

N_samples = 1000;
Reward_y = normrnd(0.5, 0.25, 1, N_samples);
Reward_b = normrnd(1, 0.25, 1, N_samples);


N_actions = 100;

beta = 1;  %this is a trade-off between exploration and exploitation (something like a learning rate!)
% very large beta values may lead to the wrong decision (choosing yellow)
epsilon = 0.1;

my_exp = zeros(1,N_actions);
mb_exp = zeros(1,N_actions);
ch = zeros(1,N_actions);

ch(1) = round (rand);
for i = 1:N_actions-1

    if ch(i) == 0
        my_exp(i+1) = my_exp(i) + epsilon * (Reward_y(i) - my_exp(i));
        mb_exp(i+1) = mb_exp(i);
    elseif ch(i) == 1
        mb_exp(i+1) = mb_exp(i) + epsilon * (Reward_b(i) - mb_exp(i));
        my_exp(i+1) = my_exp(i);
    end


    P_dens = exp(beta*my_exp(i+1)) + exp(beta*mb_exp(i+1));
    py = exp(beta * my_exp(i+1)) / P_dens;
    pb = exp(beta * mb_exp(i+1)) / P_dens;

    xy = py * rand(1);
    xb = pb * rand (1);
    if xy > xb
        ch(i+1) = 0;
    elseif xb > xy
        ch(i+1) = 1;
    end

end

% - plotting

bins = -2:0.1:3;

figure
subplot(2,2,1); hold on
hist(Reward_y, bins)
xlabel('Reward value (Yellow)')
ylabel('Probability')

subplot(2,2,2); hold on
hist(Reward_b, bins)
xlabel('Reward value (Blue)')
ylabel('Probability')

subplot(2,2,3); hold on
plot(cumsum(ch==0), 'LineWidth',2, 'Color','y')
plot(cumsum(ch==1), 'LineWidth',2, 'Color','b')
xlabel('Action history')
ylabel('Sum action')


subplot(2,2,4); hold on
plot(my_exp, 'LineWidth',2, 'color', 'y')
plot(mb_exp, 'LineWidth',2, 'color', 'b')
xlabel('Action history')
ylabel('Expected reward')