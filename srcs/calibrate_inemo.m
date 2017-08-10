clear all;
close all;
clc;
experim=1;
format longg
% Importing data
switch (experim)
    case 0
        % David Tedaldi's setting for Xsens
        IMU0x2Dalpha = importdata('IMU0x2Dalpha.mat');
        IMU0x2Domega = importdata('IMU0x2Domega.mat');
        % gravity magnitude at the experiment site
        magnitude = 9.81744;
        % the accelerometer calibration parameters initial values
        theta_pr = [0, 0, 0, 0, 0, 0, 0, 0, 0]; 
        % gyro calibration parameter settings
        theta_pr_gyro = [1/6258,0,0,0,1/6258,0,0,0,1/6258]; 
        % values from Datasheet for Xsens
        DS_a_scale = diag([415, 413, 415]);
        DS_a_misal = [1.00, 0.00, -0.01; 0.01, 1.00, 0.01; 0.02, 0.01, 1.00];
        DS_g_scale = diag([4778, 4758, 4766]);
        DS_g_misal = [1.00, -0.01, -0.02; 0.00, 1.00, 0.04; -0.01, 0.01, 1.00];
        % offset used by Xsens
        offset_acc_x = 33123;
        offset_acc_y = 33276;
        offset_acc_z = 32360;
        offset_gyro_x = 32768;
        offset_gyro_y = 32466;
        offset_gyro_z = 32485;
        % device characteristics of Xsens
        freq=100;
        Tinit=30;
    case 1        
        % Huai{ iNEMO data collected on June 19 2014 has a Tinit 112 sec. twait is
        % approximately 3 secs. 112 sec is reasonable beacuse in Tedaldi's
        % experiment, they found 50 s is good for 100HZ IMU. INEMO is 50HZ, so I
        % choose some value around 100 sec as Tinit. By Allan
        % variance analysis as done in Tedaldi's paper, 2014, from 50 to
        % 120 second the allan variance doesnot change much
        chunk=load('H:\OpenShoe-Matlab-Implementation\PoseKinIMU\keepme19062014111706.tsv');
        chunk(:,1)=chunk(:,1)/1000;
        chunk(:,6)=-chunk(:,6); % wrong y sign
        IMU0x2Dalpha = chunk(:,1:4);
        IMU0x2Domega = chunk(:,[1,5:7]);
        % gravity magnitude at the experiment site
        magnitude = 9.8008;
        % the accelerometer calibration parameters initial values
        theta_pr = [0, 0, 0, 9.8/1000, 9.8/1000, 9.8/1000,  -0.2617,  0.3032, 1.1519];
%         bound=0.2; % artifically 0.2 gives better results
%         lb=[-.5, -.5, -.5, 9.8/2000, 9.8/2000, 9.8/2000, -bound, -bound, -bound]; % lower bound
%         ub=[.5, .5, .5, 9.8/500, 9.8/500, 9.8/500, bound, bound, bound]; % upper bound
        % gyro calibration parameter settings
        theta_pr_gyro = [pi/180,0,0,0,pi/180,0,0,0,pi/180];
        
        % my supposed value for iNemo
        DS_a_scale = diag([1000/9.8, 1000/9.8, 1000/9.8]);
        DS_a_misal = [1.00, 0.00, 0.0; 0.0, 1.00, 0.0; 0.0, 0.0, 1.00];
        DS_g_scale = diag([180/pi, 180/pi, 180/pi]);
        DS_g_misal = [1.00, 0.0, 0.0; 0.00, 1.00, 0.0; 0.0, 0.0, 1.00];
        
        % offset for iNemo
        offset_acc_x = 0;
        offset_acc_y = 0;
        offset_acc_z = 0;
        offset_gyro_x = 0;
        offset_gyro_y = 0;
        offset_gyro_z = 0;
        % characteristics of iNemo
        freq=50; % how many samples we have in one second
        Tinit= 110; % how long is the initial static period
end
% } Huai

time = IMU0x2Domega(:,1)';

intervals = zeros(1, length(time));
intervals(1) = ((time(2) - time(1))/2) + time(1)/2;
intervals(2:length(intervals)) = time(2:length(time))-time(1:length(time)-1);

total_sample = length(IMU0x2Domega(:,1));

a_xp = IMU0x2Dalpha(:,2)' - offset_acc_x*ones(1,total_sample);
a_yp = IMU0x2Dalpha(:,3)' - offset_acc_y*ones(1,total_sample);
a_zp = IMU0x2Dalpha(:,4)' - offset_acc_z*ones(1,total_sample);
omega_x = IMU0x2Domega(:,2)' - offset_gyro_x*ones(1,total_sample);
omega_y = IMU0x2Domega(:,3)' - offset_gyro_y*ones(1,total_sample);
omega_z = IMU0x2Domega(:,4)' - offset_gyro_z*ones(1,total_sample);

%% Static State Statistical Filter
range=freq*Tinit; % what is the largest index for a sample belong to the initial static period
var_3D = (var(a_xp(1:range))^2 + var(a_yp(1:range))^2 + var(a_zp(1:range))^2);  
w_d= freq+1; % Window's dimension
normal_x = zeros(1, total_sample);
normal_y = zeros(1, total_sample);
normal_z = zeros(1, total_sample);

half_w_d = floor(w_d/2);

% inizialize integral
integral_x = sum(a_xp(1:w_d)/freq);
integral_y = sum(a_yp(1:w_d)/freq);
integral_z = sum(a_zp(1:w_d)/freq);

for i = (half_w_d + 1):total_sample - (half_w_d + 1)
    
   normal_x(i) = var(a_xp(i - half_w_d:i + half_w_d));
   normal_y(i) = var(a_yp(i - half_w_d:i + half_w_d));
   normal_z(i) = var(a_zp(i - half_w_d:i + half_w_d));
   
end

s_square = (normal_x.^2 + normal_y.^2 + normal_z.^2);
figure
plot(s_square);
title('s^2, variance in a window');

s_filter = zeros(1, total_sample);

%% Cycle used to individuate the optimal threshold

max_times_the_var = 20;

res_norm_vector = zeros(9 + 1 + 1,max_times_the_var);

for times_the_var = 1:max_times_the_var  
    for i = half_w_d:total_sample - (half_w_d + 1)
        
        if s_square(i) < times_the_var*var_3D            
            s_filter(i) = 1;            
        end        
    end    
    [QS_time_interval_info_matrix, s_filter]=compute_time_interval_info_matrix(s_filter, freq);
    % data selection - accelerometer
    qsTime = 1;
    sample_period =1/freq;
    num_samples = qsTime/sample_period;
    
    signal = [a_xp;a_yp; a_zp];
    
    selected_data = zeros(3, 1);
    l = 1;
    
    % select static samples with a step of selection_step
    for j = 1:length(QS_time_interval_info_matrix(:,1))
        
        if QS_time_interval_info_matrix(j,3) < num_samples
        % do nothing when there is a mis detection of static period           
        else            
            selection_step = floor(QS_time_interval_info_matrix(j,3)/num_samples); % interval between static samples           
            for i = 1:(num_samples - 1)                
                selected_data(1:3, l) = signal(1:3, QS_time_interval_info_matrix(j, 1) + (i - 1)*selection_step);
                l = l + 1;                
            end            
            selected_data(1:3, l) = signal(1:3, QS_time_interval_info_matrix(j, 2));
            l = l + 1;
            
        end
        
    end
    % minimization
    selectedAccData = selected_data;
    % Huai {
    ObjectiveFunction = @(theta_pr) accCostFunctLSQNONLIN(theta_pr, selectedAccData, magnitude);
    %     options = optimset('MaxFunEvals', 150000, 'MaxIter', 6000,
    %     'TolFun', 10^(-10)); % setting by Tedaldi
    options = optimset('MaxFunEvals', 2000, 'MaxIter', 600, 'TolFun', 1e-6); % these parameters typically do not make much difference
    if(exist('lb', 'var'))
        [theta_pr rsnorm] = lsqnonlin(ObjectiveFunction, theta_pr, lb,ub, options); 
    else
        [theta_pr rsnorm] = lsqnonlin(ObjectiveFunction, theta_pr, [],[], options);
    end
    % Huai }    
    res_norm_vector(:,times_the_var) = [theta_pr';rsnorm;times_the_var*var_3D];    
end

vec = res_norm_vector(10,:);  
z=find(vec==min(vec));
% Huai {
z=15; % from my observation, 8 is the best value to distinguish static session in my June 6 2014 iNemo dataset 
% } Huai
threshold_opt = res_norm_vector(11,z);

theta_pr_opt = res_norm_vector(1:9,z)';

estimated_misalignmentMatrix = [1, -theta_pr_opt(1), theta_pr_opt(2); 0, 1, -theta_pr_opt(3); 0, 0, 1];
estimated_scalingMatrix = diag([theta_pr_opt(4), theta_pr_opt(5), theta_pr_opt(6)]);
estimated_biasVector = [theta_pr_opt(7); theta_pr_opt(8); theta_pr_opt(9)];

s_filter = zeros(1, total_sample);
for i = half_w_d:total_sample - (half_w_d + 1)    
    if s_square(i) < threshold_opt        
        s_filter(i) = 1;        
    end    
end
[QS_time_interval_info_matrix, s_filter]=compute_time_interval_info_matrix(s_filter, freq);
figure
plot(a_xp)
hold on
plot(a_yp, 'red')
plot(a_zp, 'green')
plot(floor(max(a_zp)/2)*s_filter, 'black')
title('specific force in x, y, z and the dynamic session with value 0');
legend('ax', 'ay', 'az', 'dynamic');

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     GYROSCOPE BIAS REMUVAL                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QS filter for the first static region individuation
init_long_qs_interval_start = 0;
init_long_qs_interval_end = 0;
flag_is_first_long_static_interval = 1;

for i=1:total_sample
    
    if s_filter(i) == 0 && flag_is_first_long_static_interval == 1
        0;
    elseif s_filter(i) == 1 && flag_is_first_long_static_interval == 1
        init_long_qs_interval_start = i;
        flag_is_first_long_static_interval = 2;
    elseif s_filter(i) == 1 && flag_is_first_long_static_interval == 2
    elseif s_filter(i) == 0 && flag_is_first_long_static_interval == 2
        init_long_qs_interval_end = i;
        break
    end
end

figure
plot(omega_x);
hold on
plot(omega_y, 'red');
plot(omega_z, 'green');
title('raw angular rates in x, y, z');
legend('x', 'y', 'z');
estimate_bias_x = mean(omega_x(init_long_qs_interval_start:init_long_qs_interval_end));
estimate_bias_y = mean(omega_y(init_long_qs_interval_start:init_long_qs_interval_end));
estimate_bias_z = mean(omega_z(init_long_qs_interval_start:init_long_qs_interval_end));

omega_x = omega_x - estimate_bias_x;
omega_y = omega_y - estimate_bias_y;
omega_z = omega_z - estimate_bias_z;

figure
plot(omega_x);
hold on
plot(omega_y, 'red');
plot(omega_z, 'green');
title('calibrated angular rates in x, y, z');
legend('x', 'y', 'z');
% Huai {
disp(estimated_biasVector); % accelerator bias vector
disp(estimated_scalingMatrix); % acc
disp(estimated_misalignmentMatrix); % acc
disp([estimate_bias_x, estimate_bias_y, estimate_bias_z]); % gyro
% Huai}

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     GYROSCOPE MINIMIZATION                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimated_acc_misalignmentMatrix = estimated_misalignmentMatrix;
estimated_acc_scalingMatrix = estimated_scalingMatrix;
estimated_acc_biasVector = estimated_biasVector;
% Calibrating acceleromter data
calib_acc = estimated_acc_misalignmentMatrix*estimated_acc_scalingMatrix*[a_xp;a_yp;a_zp] - (diag(estimated_acc_biasVector)*ones(3, length(a_xp)));
a_xp_calib = calib_acc(1,:);
a_yp_calib = calib_acc(2,:);
a_zp_calib = calib_acc(3,:);

num_samples = qsTime/sample_period;

signal = [a_xp_calib;a_yp_calib; a_zp_calib];

selected_data = zeros(3, 1);                %% 3xN vector, where N is the number of selected data
l = 1;
QS_time_interval_info_matrix=[QS_time_interval_info_matrix'; zeros(size(QS_time_interval_info_matrix, 2), size(QS_time_interval_info_matrix, 1))];
for g = 1:length(QS_time_interval_info_matrix(1,:))    
    selected_acc_data = zeros(3,1);
    selection_step = floor(QS_time_interval_info_matrix(3,g)/num_samples);
    
    for i = 1:(num_samples - 1)        
        selected_acc_data(1:3, i) = signal(1:3, QS_time_interval_info_matrix(1, g) + (i - 1)*selection_step);        
    end    
    selected_data(1:3, num_samples) = signal(1:3, QS_time_interval_info_matrix(2, g));
    QS_time_interval_info_matrix(4:6, l) = mean(selected_acc_data, 2); % the mean gravity vector in this static session
    l = l + 1;    
end
QS_time_interval_calib_info_matrix = QS_time_interval_info_matrix;
% Minimizing LSQNONLIN

option = optimset('TolX', 10^-7, 'TolFun' , 10^-6, 'MaxFunEvals', 400);
theta_pr_gyro = lsqnonlin(@(theta_pr_gyro) gyroCostFunctLSQNONLIN(theta_pr_gyro, QS_time_interval_calib_info_matrix, omega_x, omega_y, omega_z, 1/freq), theta_pr_gyro, [], [], option);

misal_matrix = [1, theta_pr_gyro(2), theta_pr_gyro(3); theta_pr_gyro(4), 1, theta_pr_gyro(6); theta_pr_gyro(7), theta_pr_gyro(8), 1];
scale_matrix = [theta_pr_gyro(1), 0, 0; 0, theta_pr_gyro(5), 0; 0, 0, theta_pr_gyro(9)];

disp(scale_matrix); % gyro
disp(misal_matrix); % gyro

disp(['Accelerometer''','s Scaling Matrix taken from datasheet:']);
disp(DS_a_scale);
disp(['Accelerometer''','s Misalignment Matrix taken from datasheet:']);
format bank
disp(DS_a_misal);
disp(['Gyroscope''','s Scaling Matrix taken from datasheet:']);
format short
disp(DS_g_scale);
disp(['Gyroscope''','s Misalignment Matrix taken from datasheet:']);
format bank
disp(DS_g_misal);
disp('-----------------------------------------------------------------')

format short
[comp_a_scale, comp_a_misal, comp_g_scale, comp_g_misal] = ...
    obtainComparableMatrix(estimated_scalingMatrix, estimated_misalignmentMatrix, scale_matrix, misal_matrix);

disp(['Accelerometer''','s Estimated Scaling Matrix:']);
disp(comp_a_scale);
disp(['Accelerometer''','s Estimated Misalignment Matrix:']);
disp(comp_a_misal);
disp(['Gyroscope''','s Estimated Scaling Matrix:']);
disp(comp_g_scale);
disp(['Gyroscope''','s Estimated Misalignment Matrix:']);
disp(comp_g_misal);