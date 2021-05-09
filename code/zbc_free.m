function a = zbc_free (a)
% This function implements z-direction free boundary condition to a given field a

a(:, :, 1) = a(:, :, 2);
a(:, :, end) = a(:, :, end-1);

