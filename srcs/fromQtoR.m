%% from rotation versor(quaternion) to raotation matrix

function [R] = fromQtoR(q)

r_11 = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
r_12 = 2*(q(2)*q(3) - q(1)*q(4));
r_13 = 2*(q(2)*q(4) + q(1)*q(3));
r_21 = 2*(q(2)*q(3) + q(1)*q(4));
r_22 = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
r_23 = 2*(q(3)*q(4) - q(1)*q(2));
r_31 = 2*(q(2)*q(4) - q(1)*q(3));
r_32 = 2*(q(3)*q(4) + q(1)*q(2));
r_33 = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;

R = [r_11, r_12, r_13; r_21, r_22, r_23; r_31, r_32, r_33];

end