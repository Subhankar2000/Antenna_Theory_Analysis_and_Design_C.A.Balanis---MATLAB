%***********************************************************************
%	[w, error] = LMS(w,Mu,x,d)
%***********************************************************************
%	LMS is a MATLAB function that computes the weight coefficients
%   using the Least Mean Square algorithm
%
%	Input Parameters Description
%	----------------------------
%	- w          initial weight coefficients (row or column vector)
%   - Mu         convergence factor
%	- x          input data (column vector)
%	- d          sample of desired or reference signal
%
%	Output Parameters Description
%	-----------------------------
%	- w          updated weight coefficients (row or column vector)
%	- error      d - w' * x
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function [w, error] = LMS(w,Mu,x,d)

error = d - w' * x;
w = w + 2 * Mu * x * conj(error);
