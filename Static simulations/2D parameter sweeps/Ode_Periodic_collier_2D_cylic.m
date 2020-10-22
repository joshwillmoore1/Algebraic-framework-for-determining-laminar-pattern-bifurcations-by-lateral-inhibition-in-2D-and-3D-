function dydt = Ode_Periodic_collier_2D_cylic(t,y,W_MAT,w1,w2,n1,n2);
dydt = zeros(12,1);

%paramter values from collier et al. 1996

a = 0.01;
b = 100;

k = 2;
h = 1;

mu1 = 1;
mu2 = 1;


Y = transpose([y(2),y(4),y(6),y(8),y(10),y(12)]);
d = (W_MAT*Y)./(n1*w1 + n2*w2);


%% The model %%


%Her the geometry is a grid and we assume period BCs - No stroma

%%Cell i-1,j+1%%

d0 = d(1);
dydt(1) = (d0)^k/( a + d0^k ) - mu1*y(1);
dydt(2) = 1/(1+ b*y(1)^h) - mu2*y(2); 


%%Cell i,j+1%%

d1 = d(2);

dydt(3) = (d1)^k/( a + d1^k ) - mu1*y(3);
dydt(4) = 1/(1+ b*y(3)^h) - mu2*y(4); 


%%Cell i+1,j+1%%

d2 = d(3);

dydt(5) = (d2)^k/( a + d2^k ) - mu1*y(5);
dydt(6) = 1/(1+ b*y(5)^h) - mu2*y(6); 


%%Cell i-1,j%%

d3 = d(4);

dydt(7) = (d3)^k/( a + d3^k ) - mu1*y(7);
dydt(8) = 1/(1+ b*y(7)^h) - mu2*y(8); 


%%Cell i,j%%

d4 = d(5);

dydt(9) = (d4)^k/( a + d4^k ) - mu1*y(9);
dydt(10) = 1/(1+ b*y(9)^h) - mu2*y(10); 


%%Cell i+1,j%%

d5 = d(6);

dydt(11) = (d5)^k/( a + d5^k ) - mu1*y(11);
dydt(12) = 1/(1+ b*y(11)^h) - mu2*y(12); 
end