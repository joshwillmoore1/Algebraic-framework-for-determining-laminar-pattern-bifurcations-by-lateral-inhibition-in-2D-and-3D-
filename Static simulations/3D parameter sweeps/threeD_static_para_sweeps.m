close all
clear
clc
plot_setup %remove on users machine

%ODEs setup
tspan = [0 100];

IC_str_basal_small = [ 1-0.1,  1+0.1,  1-0.1,  1+0.1, 1-0.1, 1+0.1,1-0.1, 1+0.1,1-0.1, 1+0.1,1-0.1, 1+0.1,1-0.1, 1+0.1,1-0.1, 1+0.1,1-0.1, 1+0.1]; %% if we were to start polarised
IC_str_lum_small = [1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1];
IC_POLAR_SMALL = [IC_str_basal_small,IC_str_lum_small];

%reduce the mesh size to increase efficiency
trials = 150;
Length = 9; %3 cell periodic
Height = 2; %bilayer
pattern_tol = 0.1;

%paramter space for w1 and w2
w1=linspace(1e-4,0.25,trials);
w2=linspace(1e-4,1,trials);

[W1,W2] = meshgrid(w1,w2);


%Polar initial conditions
avg_delta_diff_3DVN_POLAR = zeros(trials,trials);
avg_delta_diff_3DT_POLAR = zeros(trials,trials);
avg_delta_diff_3DOS1_POLAR = zeros(trials,trials);
avg_delta_diff_3DOS2_POLAR = zeros(trials,trials);
avg_delta_diff_3DM_POLAR = zeros(trials,trials);

delta_3DVN_POLAR = [];
delta_3DT_POLAR = [];
delta_3DOS1_POLAR = [];
delta_3DOS2_POLAR = [];
delta_3DM_POLAR = [];



for u = 1:numel(w1)
    for v = 1:numel(w2)
        
        [t_3DVN,y_3DVN] = ode45(@(t,y) Ode_Periodic_collier_3D_cylic(t,y,W_mat_3DVN(w1(u),w2(v)),w1(u),w2(v),4,1), tspan, IC_POLAR_SMALL);
        [t_3DT,y_3DT] = ode45(@(t,y) Ode_Periodic_collier_3D_cylic(t,y,W_mat_3DT(w1(u),w2(v)),w1(u),w2(v),6,3), tspan, IC_POLAR_SMALL);
        [t_3DOS1,y_3DOS1] = ode45(@(t,y) Ode_Periodic_collier_3D_cylic(t,y,W_mat_3DOS1(w1(u),w2(v)),w1(u),w2(v),4,4), tspan, IC_POLAR_SMALL);
        [t_3DOS2,y_3DOS2] = ode45(@(t,y) Ode_Periodic_collier_3D_cylic(t,y,W_mat_3DOS1(w1(u),w2(v)),w1(u),w2(v),8,4), tspan, IC_POLAR_SMALL);
        [t_3DM,y_3DM] = ode45(@(t,y) Ode_Periodic_collier_3D_cylic(t,y,W_mat_3DM(w1(u),w2(v)),w1(u),w2(v),8,5), tspan, IC_POLAR_SMALL);
        
        
        for i = 1:Length
            for j = 1:Height
                
                if j == 1
                    index =  mod(i,9);
                    if index == 0
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,20);
                        delta_3DT_POLAR(i,j) = y_3DT(end,20);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,20);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,20);
                        delta_3DM_POLAR(i,j) = y_3DM(end,20);
                    end
                    
                    if index == 1
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,22);
                        delta_3DT_POLAR(i,j) = y_3DT(end,22);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,22);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,22);
                        delta_3DM_POLAR(i,j) = y_3DM(end,22);
                    end
                    if index == 2
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,24);
                        delta_3DT_POLAR(i,j) = y_3DT(end,24);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,24);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,24);
                        delta_3DM_POLAR(i,j) = y_3DM(end,24);
                    end
                    if index == 3
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,26);
                        delta_3DT_POLAR(i,j) = y_3DT(end,26);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,26);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,26);
                        delta_3DM_POLAR(i,j) = y_3DM(end,26);
                    end
                    if index == 4
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,28);
                        delta_3DT_POLAR(i,j) = y_3DT(end,28);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,28);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,28);
                        delta_3DM_POLAR(i,j) = y_3DM(end,28);
                    end
                    if index == 5
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,30);
                        delta_3DT_POLAR(i,j) = y_3DT(end,30);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,30);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,30);
                        delta_3DM_POLAR(i,j) = y_3DM(end,30);
                    end
                    if index == 6
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,32);
                        delta_3DT_POLAR(i,j) = y_3DT(end,32);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,32);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,32);
                        delta_3DM_POLAR(i,j) = y_3DM(end,32);
                    end
                    if index == 7
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,34);
                        delta_3DT_POLAR(i,j) = y_3DT(end,34);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,34);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,34);
                        delta_3DM_POLAR(i,j) = y_3DM(end,34);
                    end
                    if index == 8
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,36);
                        delta_3DT_POLAR(i,j) = y_3DT(end,36);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,36);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,36);
                        delta_3DM_POLAR(i,j) = y_3DM(end,36);
                    end
                end
                
                
                if j == 2
                    index =  mod(i,9);
                    
                    if index == 0
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,2);
                        delta_3DT_POLAR(i,j) = y_3DT(end,2);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,2);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,2);
                        delta_3DM_POLAR(i,j) = y_3DM(end,2);
                        
                    end
                    
                    if index == 1
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,4);
                        delta_3DT_POLAR(i,j) = y_3DT(end,4);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,4);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,4);
                        delta_3DM_POLAR(i,j) = y_3DM(end,4);
                        
                    end
                    
                    if index == 2
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,6);
                        delta_3DT_POLAR(i,j) = y_3DT(end,6);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,6);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,6);
                        delta_3DM_POLAR(i,j) = y_3DM(end,6);
                        
                    end
                    
                    if index == 3
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,8);
                        delta_3DT_POLAR(i,j) = y_3DT(end,8);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,8);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,8);
                        delta_3DM_POLAR(i,j) = y_3DM(end,8);
                        
                    end
                    
                    if index == 4
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,10);
                        delta_3DT_POLAR(i,j) = y_3DT(end,10);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,10);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,10);
                        delta_3DM_POLAR(i,j) = y_3DM(end,10);
                        
                    end
                    
                    if index == 5
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,12);
                        delta_3DT_POLAR(i,j) = y_3DT(end,12);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,12);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,12);
                        delta_3DM_POLAR(i,j) = y_3DM(end,12);
                        
                    end
                    
                    if index == 6
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,14);
                        delta_3DT_POLAR(i,j) = y_3DT(end,14);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,14);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,14);
                        delta_3DM_POLAR(i,j) = y_3DM(end,14);
                        
                    end
                    
                    if index == 7
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,16);
                        delta_3DT_POLAR(i,j) = y_3DT(end,16);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,16);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,16);
                        delta_3DM_POLAR(i,j) = y_3DM(end,16);
                        
                    end
                    
                    if index == 8
                        
                        delta_3DVN_POLAR(i,j) = y_3DVN(end,18);
                        delta_3DT_POLAR(i,j) = y_3DT(end,18);
                        delta_3DOS1_POLAR(i,j) = y_3DOS1(end,18);
                        delta_3DOS2_POLAR(i,j) = y_3DOS2(end,18);
                        delta_3DM_POLAR(i,j) = y_3DM(end,18);
                        
                    end
                end
                
            end
        end
        
        
        temp_del_b_3DVN=sum(delta_3DVN_POLAR(:,2));
        temp_del_l_3DVN=sum(delta_3DVN_POLAR(:,1));
        
        
        temp_del_b_3DT=sum(delta_3DT_POLAR(:,2));
        temp_del_l_3DT=sum(delta_3DT_POLAR(:,1));
        
        
        temp_del_b_3DOS1=sum(delta_3DOS1_POLAR(:,2));
        temp_del_l_3DOS1=sum(delta_3DOS1_POLAR(:,1));
        
        temp_del_b_3DOS2=sum(delta_3DOS2_POLAR(:,2));
        temp_del_l_3DOS2=sum(delta_3DOS2_POLAR(:,1));
        
        temp_del_b_3DM=sum(delta_3DM_POLAR(:,2));
        temp_del_l_3DM=sum(delta_3DM_POLAR(:,1));
        
        
        avg_delta_diff_3DVN_POLAR(u,v) = (temp_del_b_3DVN-temp_del_l_3DVN)./Length;
        avg_delta_diff_3DT_POLAR(u,v) = (temp_del_b_3DT-temp_del_l_3DT)./Length;
        avg_delta_diff_3DOS1_POLAR(u,v) = (temp_del_b_3DOS1-temp_del_l_3DOS1)./Length;
        avg_delta_diff_3DOS2_POLAR(u,v) = (temp_del_b_3DOS1-temp_del_l_3DOS1)./Length;
        avg_delta_diff_3DM_POLAR(u,v) = (temp_del_b_3DM-temp_del_l_3DM)./Length;
        
        
        
        
    end
    disp(['processing... ',num2str(u/numel(w1)*100),'%']);
end



figure();

subplot(1,5,1);

[M,c] = contourf(W2,W1,avg_delta_diff_3DVN_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(2,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('3DVN','Interpreter','latex');
box off
ylim([0,0.25])


subplot(1,5,2);
[M,c] = contourf(W2,W1,avg_delta_diff_3DT_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(1,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('3DT','Interpreter','latex');
box off
ylim([0,0.25])




subplot(1,5,3);
[M,c] = contourf(W2,W1,avg_delta_diff_3DOS1_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(2/3,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('3DOS1','Interpreter','latex');
box off
ylim([0,0.25])


subplot(1,5,4);
[M,c] = contourf(W2,W1,avg_delta_diff_3DOS2_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(2/3,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('3DOS2','Interpreter','latex');
box off
ylim([0,0.25])



subplot(1,5,5);
[M,c] = contourf(W2,W1,avg_delta_diff_3DM_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(2/3,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('3DM','Interpreter','latex');
box off
ylim([0,0.25])



%existence bound from Theorem 4.1 - Moore et al. 2020
function out = bound_w1(R,w2)
out  = 0.2106.*w2./R;
end

%full adjacency matrices for the 3D geometries
function out = W_mat_3DVN(w1,w2)

M1 = [ [0,w1,w1,w1,0,0,w1,0,0];[w1,0,w1,0,w1,0,0,w1,0]; [w1,w1,0,0,0,w1,0,0,w1]; [w1,0,0,0,w1,w1,w1,0,0]; [0,w1,0,w1,0,w1,0,w1,0];[0,0,w1,w1,w1,0,0,0,w1]; [w1,0,0,w1,0,0,0,w1,w1] ; [0,w1,0,0,w1,0,w1,0,w1];[0,0,w1,0,0,w1,w1,w1,0] ];
M2 = w2.*eye(9,9);

out = [[M1,M2];[M2,M1]];
end

function out = W_mat_3DT(w1,w2)

M2 = [ [w2,0,w2,0,0,0,w2,0,0]; [w2,w2,0,0,0,0,0,w2,0]; [0,w2,w2,0,0,0,0,0,w2]; [w2,0,0,w2,0,w2,0,0,0]; [0,w2,0,w2,w2,0,0,0,0]; [0,0,w2,0,w2,w2,0,0,0]; [0,0,0,w2,0,0,w2,0,w2];  [0,0,0,0,w2,0,w2,w2,0]; [0,0,0,0,0,w2,0,w2,w2]   ];
M1 = [[0,w1,w1,w1,0,w1,w1,w1,0 ];  [w1,0,w1,w1,w1,0,0,w1,w1]; [w1,w1,0,0,w1,w1,w1,0,w1]; [w1,w1,0,0,w1,w1,w1,0,w1]; [0,w1,w1,w1,0,w1,w1,w1,0]; [w1,0,w1,w1,w1,0,0,w1,w1]; [w1,0,w1,w1,w1,0,0,w1,w1]; [w1,w1,0,0,w1,w1,w1,0,w1]; [0,w1,w1,w1,0,w1,w1,w1,0]];

out = [[M1,M2];[M2,M1]];
end

function out = W_mat_3DOS1(w1,w2)

M1 = [[0,w1,w1,w1,0,0,w1,0,0]; [w1,0,w1,0,w1,0,0,w1,0]; [w1,w1,0,0,0,w1,0,0,w1]; [w1,0,0,0,w1,w1,w1,0,0]; [0,w1,0,w1,0,w1,0,w1,0]; [0,0,w1,w1,w1,0,0,0,w1];[w1,0,0,w1,0,0,0,w1,w1]; [0,w1,0,0,w1,0,w1,0,w1];[0,0,w1,0,0,w1,w1,w1,0]];
M2= [[w2,0,w2,0,0,0,w2,0,w2]; [w2,w2,0,0,0,0,w2,w2,0]; [0,w2,w2,0,0,0,0,w2,w2]; [w2,0,w2,w2,0,w2,0,0,0]; [w2,w2,0,w2,w2,0,0,0,0]; [0,w2,w2,0,w2,w2,0,0,0]; [0,0,0,w2,0,w2,w2,0,w2]; [0,0,0,w2,w2,0,w2,w2,0]; [0,0,0,0,w2,w2,0,w2,w2]];

out = [[M1,M2];[M2,M1]];
end


function out = W_mat_3DOS2(w1,w2)

M1 = w1.*ones(9,9);
M1(1,1) = 0;
M1(2,2) = 0;
M1(3,3) = 0;
M1(4,4) = 0;
M1(5,5) = 0;
M1(6,6) = 0;
M1(7,7) = 0;
M1(8,8) = 0;
M1(9,9) = 0;
M2= [[w2,0,w2,0,0,0,w2,0,w2]; [w2,w2,0,0,0,0,w2,w2,0]; [0,w2,w2,0,0,0,0,w2,w2]; [w2,0,w2,w2,0,w2,0,0,0]; [w2,w2,0,w2,w2,0,0,0,0]; [0,w2,w2,0,w2,w2,0,0,0]; [0,0,0,w2,0,w2,w2,0,w2]; [0,0,0,w2,w2,0,w2,w2,0]; [0,0,0,0,w2,w2,0,w2,w2]];

out = [[M1,M2];[M2,M1]];

end


function out = W_mat_3DM(w1,w2)

M2 = [ [w2,w2,w2,w2,0,0,w2,0,0];[w2,w2,w2,0,w2,0,0,w2,0]; [w2,w2,w2,0,0,w2,0,0,w2]; [w2,0,0,w2,w2,w2,w2,0,0]; [0,w2,0,w2,w2,w2,0,w2,0];[0,0,w2,w2,w2,w2,0,0,w2]; [w2,0,0,w2,0,0,w2,w2,w2] ; [0,w2,0,0,w2,0,w2,w2,w2];[0,0,w2,0,0,w2,w2,w2,w2] ];

M1 = w1.*ones(9,9);
M1(1,1) = 0;
M1(2,2) = 0;
M1(3,3) = 0;
M1(4,4) = 0;
M1(5,5) = 0;
M1(6,6) = 0;
M1(7,7) = 0;
M1(8,8) = 0;
M1(9,9) = 0;


out = [[M1,M2];[M2,M1]];

end

