%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the (inverse) Fractional Fast Fourier Transform.
%
% Reference: David H. Bailey and Paul N. Swarztrauber 1991.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performs the Fractional FFT.
% param vd_series: Series to be transformed. Required to be of length 2^N.
% param d_beta:    Beta parameter of FRFT.     
%
function out_transform = FRFT(vd_series, d_beta)

  i_N        = numel(vd_series);
  vi_indices = 0:(i_N-1);
  
  % create 2N sized vectors
  vd_y = [(vd_series .* exp((1i*pi*d_beta)*(vi_indices.^2))) zeros(1, i_N)];
  vd_z = [exp((-1i*pi*d_beta)*(vi_indices.^2)) ...
          exp((-1i*pi*d_beta)*(vi_indices - i_N).^2)];

  % perform one inverse and two ordinary DFTs via FFT algorithm
  vd_tmp = 2/i_N * FFT((i_N^2 * (FFT(vd_y, true) .* FFT(vd_z, true))), false);
  
  % return first half of transformation.
  out_transform = exp((1i*pi*d_beta)*(vi_indices.^2)) .* vd_tmp(1:i_N);

end