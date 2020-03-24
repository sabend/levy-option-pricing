%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extracting interpolated values from vector via specified
% interpolation method.
%
%####################################################################

% Interpolates elements.
% param vd_which:    Vector of X grid for which to determine interpolated values.
% param vd_values_x: Available grid of X values.
% param vd_values_y: Available grid of Y values.
% param s_method:    Specifies how to perform interpolate.
%
function out_interpolated = extractInterpolatedElements(vd_which, vd_values_x, vd_values_y, s_method)

  if numel(vd_values_x) ~= numel(vd_values_y)
    error('[MATH] size mismatch between interpolation X and Y values!');
  end

  vd_interpolated           = zeros(1, numel(vd_which));
  [vd_values_x, vi_indices] = sort(vd_values_x);
  vd_values_y               = vd_values_y(vi_indices);

  switch s_method
      
    case 'loglinear'

      for i = 1 : numel(vd_which)
   
        d_which = vd_which(i);
        i_after = find(d_which <= vd_values_x, 1);

        if isempty(i_after)
          error(['[MATH] element ' num2str(d_which) ' lies outside available X values!']);
        end
        
        if i_after == 1
          if d_which ~= vd_values_x(1)
            error(['[MATH] element ' num2str(d_which) ' lies outside available X values!']);
          else
            vd_interpolated(i) = d_which;
            continue
          end
        end
        
        d_value_x_after  = vd_values_x(i_after);
        d_value_x_before = vd_values_x(i_after - 1);
        d_value_y_after  = vd_values_y(i_after);
        d_value_y_before = vd_values_y(i_after - 1);
        d_diff           = d_value_x_after - d_value_x_before;
        
        vd_interpolated(i) = d_value_y_after^((d_which-d_value_x_before)/d_diff)*...
                             d_value_y_before^((d_value_x_after-d_which)/d_diff);
      
      end
    
    otherwise
      error(['[MATH] specified interpolation method ' num2str(s_method) ' not implemented!']);
    
  end
  
  out_interpolated = vd_interpolated;

end