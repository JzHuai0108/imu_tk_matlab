function [ comp_a_scale_matrix, comp_a_misal_matrix,  ...
    comp_g_scale_matrix, comp_g_misal_matrix ] = obtainComparableMatrix(...
    acc_scale_matrix, acc_misal_matrix, gyro_scale_matrix, ...
    gyro_misal_matrix )

% acc misal parameter taken from datasheet
alpha_xz_6 = 0.01;
alpha_xy_6 = -0.02;
alpha_yx_6 = 0.01;

theta = - alpha_xz_6;
w = [0 0 1];

q = [cos(theta/2) sin(theta/2)*w];

r_11 = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
r_12 = 2*(q(2)*q(3) - q(1)*q(4));
r_13 = 2*(q(2)*q(4) + q(1)*q(3));
r_21 = 2*(q(2)*q(3) + q(1)*q(4));
r_22 = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
r_23 = 2*(q(3)*q(4) - q(1)*q(2));
r_31 = 2*(q(2)*q(4) - q(1)*q(3));
r_32 = 2*(q(3)*q(4) + q(1)*q(2));
r_33 = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;

R_xz = [r_11, r_12, r_13; r_21, r_22, r_23; r_31, r_32, r_33];

theta = - alpha_xy_6;
w = [0 1 0];

q = [cos(theta/2) sin(theta/2)*w];

r_11 = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
r_12 = 2*(q(2)*q(3) - q(1)*q(4));
r_13 = 2*(q(2)*q(4) + q(1)*q(3));
r_21 = 2*(q(2)*q(3) + q(1)*q(4));
r_22 = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
r_23 = 2*(q(3)*q(4) - q(1)*q(2));
r_31 = 2*(q(2)*q(4) - q(1)*q(3));
r_32 = 2*(q(3)*q(4) + q(1)*q(2));
r_33 = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;

R_xy = [r_11, r_12, r_13; r_21, r_22, r_23; r_31, r_32, r_33];

theta = - alpha_yx_6;
w = [1 0 0];

q = [cos(theta/2) sin(theta/2)*w];

r_11 = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
r_12 = 2*(q(2)*q(3) - q(1)*q(4));
r_13 = 2*(q(2)*q(4) + q(1)*q(3));
r_21 = 2*(q(2)*q(3) + q(1)*q(4));
r_22 = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
r_23 = 2*(q(3)*q(4) - q(1)*q(2));
r_31 = 2*(q(2)*q(4) - q(1)*q(3));
r_32 = 2*(q(3)*q(4) + q(1)*q(2));
r_33 = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;

R_yx = [r_11, r_12, r_13; r_21, r_22, r_23; r_31, r_32, r_33];
 
comp_a_scale_matrix = inv(acc_scale_matrix); 
comp_a_misal_matrix = inv(R_xz*R_xy*R_yx*acc_misal_matrix);
comp_g_scale_matrix = inv(gyro_scale_matrix); 
comp_g_misal_matrix = inv(R_xz*R_xy*R_yx*gyro_misal_matrix);



end

