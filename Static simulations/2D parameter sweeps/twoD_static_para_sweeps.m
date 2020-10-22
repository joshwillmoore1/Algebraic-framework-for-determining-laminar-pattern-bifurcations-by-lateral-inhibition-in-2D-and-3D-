close all
clear 
clc
plot_setup %remove on users machine

%ODEs setup
tspan = [0 100];

IC_str_basal_small = [ 1-0.1,  1+0.1,  1-0.1,  1+0.1, 1-0.1, 1+0.1]; %% if we were to start polarised
IC_str_lum_small = [1+0.1,1-0.1,1+0.1,1-0.1,1+0.1,1-0.1];
IC_POLAR_SMALL = [IC_str_basal_small,IC_str_lum_small];

%reduce the mesh size to increase efficiency
trials = 150;
Length = 3; %3 cell periodic
Height = 2; %bilayer
pattern_tol = 0.1;

%paramter space for w1 and w2
w1=linspace(1e-4,0.25,trials);
w2=linspace(1e-4,1,trials);

[W1,W2] = meshgrid(w1,w2);


%Polar initial conditions
avg_delta_diff_5_POLAR = zeros(trials,trials);
avg_delta_diff_4_POLAR = zeros(trials,trials);
avg_delta_diff_3_POLAR = zeros(trials,trials);

delta_5Moore_POLAR = [];
delta_hex_POLAR = [];
delta_3Moore_POLAR = [];

area_count_3 = 0;
area_count_4 = 0;
area_count_5 = 0;
    
   for u = 1:numel(w1)
        for v = 1:numel(w2)
            
            [t_5Moore,y_5Moore] = ode45(@(t,y) Ode_Periodic_collier_2D_cylic(t,y,W_mat_M(w1(u),w2(v)),w1(u),w2(v),2,3), tspan, IC_POLAR_SMALL);
            [t_hexagonal,y_hexagonal] = ode45(@(t,y) Ode_Periodic_collier_2D_cylic(t,y,W_mat_T(w1(u),w2(v)),w1(u),w2(v),2,2), tspan, IC_POLAR_SMALL);
            [t_3Moore,y_3Moore] = ode45(@(t,y) Ode_Periodic_collier_2D_cylic(t,y,W_mat_VN(w1(u),w2(v)),w1(u),w2(v),2,1), tspan, IC_POLAR_SMALL);
            
            
            for i = 1:Length
                for j = 1:Height
                    
                    if j == 1
                        index =  mod(i,3);
                        if index == 0
                            
                            delta_5Moore_POLAR(i,j) = y_5Moore(end,8);
                            delta_hex_POLAR(i,j) = y_hexagonal(end,8);
                            delta_3Moore_POLAR(i,j) = y_3Moore(end,8);
                        end
                        if index == 1
                            
                            delta_5Moore_POLAR(i,j) = y_5Moore(end,10);
                            delta_hex_POLAR(i,j) = y_hexagonal(end,10);
                            delta_3Moore_POLAR(i,j) = y_3Moore(end,10);
                        end
                        if index == 2
                            
                            delta_5Moore_POLAR(i,j) = y_5Moore(end,12);
                            delta_hex_POLAR(i,j) = y_hexagonal(end,12);
                            delta_3Moore_POLAR(i,j) = y_3Moore(end,12);
                        end
                    end
                    
                    
                    if j == 2
                        index =  mod(i,3);
                        
                        if index == 0
                            
                            delta_5Moore_POLAR(i,j) = y_5Moore(end,2);
                            delta_hex_POLAR(i,j) = y_hexagonal(end,2);
                            delta_3Moore_POLAR(i,j) = y_3Moore(end,2);
                            
                        end
                        
                        if index == 1
                            
                            delta_5Moore_POLAR(i,j) = y_5Moore(end,4);
                            delta_hex_POLAR(i,j) = y_hexagonal(end,4);
                            delta_3Moore_POLAR(i,j) = y_3Moore(end,4);
                            
                        end
                        
                        if index == 2
                            
                            delta_5Moore_POLAR(i,j) = y_5Moore(end,6);
                            delta_hex_POLAR(i,j) = y_hexagonal(end,6);
                            delta_3Moore_POLAR(i,j) = y_3Moore(end,6);
                            
                        end
                    end
                    
                end
            end
            
            %5 connectivity
            temp_del_b_5=sum(delta_5Moore_POLAR(:,2));
            temp_del_l_5=sum(delta_5Moore_POLAR(:,1));
            
            %4 connectivity
            temp_del_b_4=sum(delta_hex_POLAR(:,2));
            temp_del_l_4=sum(delta_hex_POLAR(:,1));
            
            %3 connectivity
            temp_del_b_3=sum(delta_3Moore_POLAR(:,2));
            temp_del_l_3=sum(delta_3Moore_POLAR(:,1));
            
            
            avg_delta_diff_5_POLAR(u,v) = (temp_del_b_5-temp_del_l_5)./Length;
            avg_delta_diff_4_POLAR(u,v) = (temp_del_b_4-temp_del_l_4)./Length;
            avg_delta_diff_3_POLAR(u,v) = (temp_del_b_3-temp_del_l_3)./Length;
            
            if avg_delta_diff_3_POLAR(u,v) >= pattern_tol
                area_count_3 = area_count_3 + 1;
            end
            
            if avg_delta_diff_4_POLAR(u,v) >= pattern_tol
                area_count_4 = area_count_4 + 1;
            end
            
            if avg_delta_diff_5_POLAR(u,v) >= pattern_tol
                area_count_5 = area_count_5 + 1;
            end
            
            
            
        end
       disp(['processing... ',num2str(u/numel(w1)*100),'%']);
    end
    
  

%area calculation
polarised_area_3 = (area_count_3/(trials^2)) * 0.25
polarised_area_4 = (area_count_4/(trials^2)) * 0.25
polarised_area_5 = (area_count_5/(trials^2)) * 0.25


figure();

subplot(1,3,1);

[M,c] = contourf(W2,W1,avg_delta_diff_3_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(2,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('2D Neumann','Interpreter','latex');
box off
ylim([0,0.25])


subplot(1,3,2);
[M,c] = contourf(W2,W1,avg_delta_diff_4_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(1,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('2D Triangular','Interpreter','latex');
box off
ylim([0,0.25])




subplot(1,3,3);
[M,c] = contourf(W2,W1,avg_delta_diff_5_POLAR',[pattern_tol pattern_tol],'--k','linewidth',3,'fill','on');
c.FaceColor = [0.75 0.75 0.75];
hold on
plot3(w2,bound_w1(2/3,w2),20.*ones(1,numel(w2)),'-k')
xlabel("$w_{2}$",'Interpreter','latex');
ylabel({'$w_{1}$'},'Interpreter','latex');
shading interp
view(0,90);
title('2D Moore','Interpreter','latex');
box off
ylim([0,0.25])



%existence bound from Theorem 4.1 - Moore et al. 2020
function out = bound_w1(R,w2)
out  = 0.2106.*w2./R;
end

%full adjacency matrices for the 2D geometries
function out = W_mat_VN(w1,w2)
    out = [ [0,w1,w1,w2,0,0]; [w1,0,w1,0,w2,0]; [w1,w1,0,0,0,w2]; [w2,0,0,0,w1,w1];[0,w2,0,w1,0,w1];[0,0,w2,w1,w1,0]   ];
end

function out = W_mat_T(w1,w2)
    out = [ [0,w1,w1,w2,w2,0]; [w1,0,w1,0,w2,w2]; [w1,w1,0,w2,0,w2]; [w2,0,w2,0,w1,w1];[w2,w2,0,w1,0,w1];[0,w2,w2,w1,w1,0]   ];
end

function out = W_mat_M(w1,w2)
    out = [ [0,w1,w1,w2,w2,w2]; [w1,0,w1,w2,w2,w2]; [w1,w1,0,w2,w2,w2]; [w2,w2,w2,0,w1,w1];[w2,w2,w2,w1,0,w1];[w2,w2,w2,w1,w1,0]   ];
end
