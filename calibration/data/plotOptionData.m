% Plots observed option prices on strikes.
%
function plotOptionData(c_option_data, d_lower, d_upper)

  figure;
  m_colors     = jet(numel(c_option_data));
  c_maturities = cell();
  
  for i = 1:numel(c_option_data)
  
    m_tmp        = c_option_data{i, 1};
    vd_prices    = m_tmp(:, 3);
    vd_strikes   = m_tmp(:, 2);
    c_maturities = [c_maturities num2str(unique(m_tmp(:, 1)))];
    
    vi_sub = (vd_strikes >= d_lower) & (vd_strikes <= d_upper);
    
    hold on;
    scatter(vd_strikes(vi_sub), vd_prices(vi_sub), [], m_colors(i, :));
    
  end
  
  title('Option prices on strikes');
  h_legend = legend(c_maturities);
  set(h_legend, 'title', 'maturities');
  set(h_legend, 'location', 'east');

end