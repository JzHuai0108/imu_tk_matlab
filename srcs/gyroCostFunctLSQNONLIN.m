function [res_vector] = gyroCostFunctLSQNONLINUsingOnlyTheFilter(E, QS_time_interval_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat)

omega_hat = [omega_x_hat; omega_y_hat; omega_z_hat];

misalignmentMatrix = [1, E(2), E(3); E(4), 1, E(6); E(7), E(8), 1];
scalingMatrix = diag([E(1), E(5), E(9)]);

omega_bar = misalignmentMatrix*scalingMatrix*omega_hat;

omega_x = omega_bar(1,:);
omega_y = omega_bar(2,:);
omega_z = omega_bar(3,:);

vector = zeros(3,5);

for pr = 1:size(QS_time_interval_info_matrix, 2) - 1
    
    vector((pr-1)*3 + 1:(pr)*3, 1) = QS_time_interval_info_matrix(4:6,pr);
    vector((pr-1)*3 + 1:(pr)*3, 5) = QS_time_interval_info_matrix(4:6,pr + 1);
    gyroUnbiasUncalibratedValues = [omega_x(QS_time_interval_info_matrix(2,pr) + 1:QS_time_interval_info_matrix(1,pr + 1) - 1); omega_y(QS_time_interval_info_matrix(2,pr) + 1:QS_time_interval_info_matrix(1,pr + 1) - 1); omega_z(QS_time_interval_info_matrix(2,pr) + 1:QS_time_interval_info_matrix(1,pr + 1) - 1)];
    R = rotationRK4(gyroUnbiasUncalibratedValues);

    vector((pr-1)*3 + 1:(pr)*3, 2:4) = R;
    
end

residuals = zeros(length(vector(:,1))/3, 1);

for i = 1:length(vector(:,1))/3
    
    v = vector((i-1)*3 + 1:(i)*3, 5)/(vector((i-1)*3 + 1, 5)^2 + vector((i-1)*3 + 2, 5)^2 + vector((i)*3, 5)^2)^(1/2) - vector((i-1)*3 + 1:(i)*3, 2:4)*vector((i-1)*3 + 1:(i)*3, 1)/(vector((i-1)*3 + 1, 5)^2 + vector((i-1)*3 + 2, 5)^2 + vector((i)*3, 5)^2)^(1/2);
    v = v';
    residuals(i,1) = (v(1)^2 + v(2)^2 + v(3)^2)^(1/2);
    
end

res_vector = residuals;

end

