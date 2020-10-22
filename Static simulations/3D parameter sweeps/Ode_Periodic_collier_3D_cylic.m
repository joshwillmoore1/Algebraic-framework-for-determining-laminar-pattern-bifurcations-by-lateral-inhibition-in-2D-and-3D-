function dydt = Ode_Periodic_collier_3D_cylic(t,y,W_MAT,w1,w2,n1,n2);
dydt = zeros(36,1);

%paramter values from collier et al. 1996

a = 0.01;
b = 100;

k = 2;
h = 1;

mu1 = 1;
mu2 = 1;


Y = transpose([y(2),y(4),y(6),y(8),y(10),y(12),y(14),y(16),y(18),y(20),y(22),y(24),y(26),y(28),y(30),y(32),y(34),y(36)]);
d = (W_MAT*Y)./(n1*w1 + n2*w2);


%% The model %%


%Her the geometry is a grid and we assume period BCs - No stroma

%basal layer

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


d6 = d(7);

dydt(13) = (d6)^k/( a + d6^k ) - mu1*y(13);
dydt(14) = 1/(1+ b*y(13)^h) - mu2*y(14);


d7 = d(8);

dydt(15) = (d7)^k/( a + d7^k ) - mu1*y(15);
dydt(16) = 1/(1+ b*y(15)^h) - mu2*y(16);


d8 = d(9);

dydt(17) = (d8)^k/( a + d8^k ) - mu1*y(17);
dydt(18) = 1/(1+ b*y(17)^h) - mu2*y(18);




%luminal layer

d9 = d(10);

dydt(19) = (d9)^k/( a + d9^k ) - mu1*y(19);
dydt(20) = 1/(1+ b*y(19)^h) - mu2*y(20);


d10 = d(11);

dydt(21) = (d10)^k/( a + d10^k ) - mu1*y(21);
dydt(22) = 1/(1+ b*y(21)^h) - mu2*y(22);


d11 = d(12);

dydt(23) = (d11)^k/( a + d11^k ) - mu1*y(23);
dydt(24) = 1/(1+ b*y(23)^h) - mu2*y(24);


d12 = d(13);

dydt(25) = (d12)^k/( a + d12^k ) - mu1*y(25);
dydt(26) = 1/(1+ b*y(25)^h) - mu2*y(26);


d13 = d(14);

dydt(27) = (d13)^k/( a + d13^k ) - mu1*y(27);
dydt(28) = 1/(1+ b*y(27)^h) - mu2*y(28);


d14 = d(15);

dydt(29) = (d14)^k/( a + d14^k ) - mu1*y(29);
dydt(30) = 1/(1+ b*y(29)^h) - mu2*y(30);


d15 = d(16);

dydt(31) = (d15)^k/( a + d15^k ) - mu1*y(31);
dydt(32) = 1/(1+ b*y(31)^h) - mu2*y(32);


d16 = d(17);

dydt(33) = (d16)^k/( a + d16^k ) - mu1*y(33);
dydt(34) = 1/(1+ b*y(33)^h) - mu2*y(34);

d17 = d(18);

dydt(35) = (d17)^k/( a + d17^k ) - mu1*y(35);
dydt(36) = 1/(1+ b*y(35)^h) - mu2*y(36);

end