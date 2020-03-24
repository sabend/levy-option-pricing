%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the (inverse) Discrete Fourier Transform.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates the DFT.
% param vd_x:      Arbitrary row vector.
% param b_inverse: If TRUE calculate inverse Fourier transform.
%
function out_transform = DFT(vd_series, b_inverse)

  i_N = length(vd_series);

  % create coefficient matrix
  if b_inverse
    d_w = exp(-2*pi*1i/i_N);
  else
    d_w = exp(2*pi*1i/i_N);
  end
  
  % fill matrix element wise
  m_W = ones(i_N, i_N);
  for j = 2 : i_N
    for k = 2 : i_N
    
      m_W(j, k) = d_w^((j-1)*(k-1));
 
    end
  end
  
  % return transformed vector as coefficient matrix times x
  out_transform = transpose(m_W * transpose(vd_series));
  if b_inverse
    out_transform = out_transform / i_N;
   end

end