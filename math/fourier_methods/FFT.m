%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the (inverse) Fast Fourier Transform.
%
% Algorithm: Radix-2 by Cooley & Tucker.
% x^(k) = sum_(j=0)^(N-1)x(j)*omega_N^(kj) , 0 <= k < N
% See e.g. Neubauer 2012.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Performs the FFT.
% param vd_series: Series to be transformed. Required to be of length 2^N.
% param b_inverse: If TRUE calculate inverse Fourier transform.      
%
function vd_transform = FFT(vd_series, b_inverse)

  i_N = numel(vd_series);
  
  if i_N == 1
  
    vd_transform = vd_series;
  
  else

    vd_indices  = 1:i_N;
    vd_ind_even = vd_indices(2:2:i_N);
    vd_ind_odd  = vd_indices(1:2:i_N);
    
    vd_even = FFT(vd_series(vd_ind_even), b_inverse);
    vd_odd  = FFT(vd_series(vd_ind_odd), b_inverse);

    if b_inverse
      d_W = exp(-2*pi*1i/i_N);
    else
      d_W = exp(2*pi*1i/i_N);
    end
    
    i_mid    = i_N / 2;
    vi_first  = 1:i_mid;
    vi_second = (i_mid+1):i_N;

    vd_transform = zeros(1, i_N);
    vd_transform(vi_first)  = vd_odd + vd_even .* d_W.^(vi_first-1);
    vd_transform(vi_second) = vd_odd - vd_even .* d_W.^(vi_first-1);
    if b_inverse
      vd_transform = vd_transform / 2;
    end
    
  end

end