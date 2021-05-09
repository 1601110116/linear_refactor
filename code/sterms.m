function sterms()
% This function calculates terms that will be further used for time advancing.
%  They are calculated here because they are shared by more than 1 euquations.
%  We don't have boundary grids for these terms because they are never needed
global ddz_ve ve
%ddz_ve = ddz(ve);
