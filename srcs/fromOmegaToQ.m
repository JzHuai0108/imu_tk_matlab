%% [q, angularRotation, direction] = fromOmegaToQ(omega, time)

function [ q, angularRotation, direction ] = fromOmegaToQ( omega, intervals )

s = size(omega);
angularRotation = zeros(1,1);
direction = zeros(3,1);
q = zeros(1,4);

for i = 1:s(2)
    
    angularRotation(i) = (omega(1,i)^2 + omega(2,i)^2 + omega(3,i)^2)^(1/2)*intervals(i);
    direction(:,i) = omega(:,i)*(intervals(i)/angularRotation(i));
    dir = direction(:,i)';
    q(i, :) = [cos(angularRotation(i)/2), sin(angularRotation(i)/2)*dir];
    
end

end