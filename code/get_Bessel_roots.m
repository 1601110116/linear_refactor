% This script generates the roots of Bessel functions.
%  The (m+1)-th line for the Bessel function of order m
%  The l-th colomn for the l-th positive zero


% Up to the order of m_max
m_max = 20;
% First l_max positive zeros
l_max = 30;
% The accuracy of the roots
accuracy = 1e-5;
% The maximum root should be less than root_max
root_max = 200;
% The output file
output = 'alpha.mat'

% The output, zeros of the Bessel functions
alpha = zeros(m_max+1, l_max);

x = accuracy: accuracy: root_max;
for m = 0: m_max
	y = besselj(m, x);
	sgn = sign(y);
	sgn_diff = sgn(1:end-1) - sgn(2:end);
	zero_ind = find(sgn_diff~=0);
	for l = 1: l_max
		alpha(m+1, l) = x(zero_ind(l));
	end
end
save(output, 'alpha');

