% Recurrent connection of two neurons
clearvars

tau1 = 10;
tau2 = 10;
dt = 0.1;
time = 20000;

W_RNN = [ 0.2 -0.7; -0.7 0.2];

r1 = zeros ( time , 1);
r2 = zeros ( time , 1);

tau_eff1 = tau1 / (1-W_RNN(1,1));
tau_eff2 = tau2 / (1-W_RNN(2,2));

delt1 = dt / tau_eff1;
delt2 = dt / tau_eff2;

[Eigenvects,Eigenvals] = eig(W_RNN);
AmplifyFactors = 1./(1 - diag(Eigenvals));
eigvect1 = Eigenvects(:,1);
eigvect2 = Eigenvects(:,2);

I = [63 57]';

r1( 1 ) = I(1);
r2( 1 ) = I(2);

for i=1:time

   % if ( i > 10000 )
   %    I(2) = 0;
   %    I(1) = 0;
   % end
   
r1(i+1) = r1(i)* ( 1 - delt1 ) + delt1 * ( r2 (i) * W_RNN(1,2)/(1-W_RNN(1,1)) + I(1)/ (1-W_RNN(1,1)));
r2(i+1) = r2(i)* ( 1 - delt2 ) + delt2 * ( r1 (i) * W_RNN(2,1)/(1-W_RNN(2,2)) + I(2)/ (1-W_RNN(2,2)));

end

time_ms = 0.1:0.1:real((time+1)/10);


figure
subplot (2,1,1)
plot (time_ms , r1, 'b','LineWidth',2)
hold on;
plot( time_ms , r2, 'r','LineWidth',2)
title ( ' Firing rate over time ' )
grid on
legend('r_1' , 'r_2')
xlabel ( ' time ' )
ylabel ( ' Firing rate ' )
ylim([0 100])
xlim([-10 2000])


subplot (2,1,2)
plot ( r1 , r2 , 'r', 'LineWidth' ,2)
xlabel ( ' r1 ' )
ylabel ( ' r2 ' )
hold on;
g(1)=line([0 -100*eigvect1(1)], [0 -100*eigvect1(2)],'linewidth',1,'Color',[0 0 1])% eigenvector1
g(2)=line([0 -100*eigvect2(1)], [0 -100*eigvect2(2)],'linewidth',1,'Color',[0 1 0])% eigenvector2
legend('trajectory r1,r2 ' , 'e1','e2','Location','SouthWest')



