%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for building option surface.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reads relevant data from file for single maturity
% param s_path:       Path to option data.
% param s_asof:       As-of date in format 'YYYYMMDD'.
% param s_maturity:   Month of expiration in format 'YYYYMM'.
% param s_min_strike: Minimal strike value to be added to option table.
% param s_max_strike: Maximal strike value to be added to option table.
%
function opt_table = txt2surface(s_path, s_asof, s_maturity, s_min_strike, s_max_strike)

  % Determine time to maturity as year fraction (252 days).
  % For Eurex options the final settlement day is the third
  % friday within the month of expiration. For details, see
  % http://www.eurexchange.com/exchange-en/products/idx/stx/blc/EURO-STOXX-50--Index-Options/19066
  % Fridays correspond to the Octave weekday 6.
  i_asof = datenum(s_asof, 'yyyymmdd');
  i_mat  = datenum([s_maturity '01'], 'yyyymmdd');
  i_weekday = weekday(i_mat);
  i_datediff = 0;
  if (i_weekday == 6)
    i_datediff = 14;
  elseif (i_weekday == 7)
    i_datediff = 14 + 6;
  else
     i_datediff = 14 + 6 - i_weekday;
  end
  i_mat      = i_mat + i_datediff;
  vi_dates    = i_asof:i_mat;
  i_days_to_mat = sum((weekday(vi_dates) ~= 1) & (weekday(vi_dates) ~= 7));
  d_maturity = i_days_to_mat / 252;

  % read file row-wise and extract strike and price
  f_data = fopen(s_path);
  
  s_line = fgetl(f_data);
  i = 1;
  m_tmp = [0 0 0];
  while ischar(s_line)
  
      s_line = fgetl(f_data);
      if (s_line == -1) break; end;
      s_line = strrep(s_line, ',', '');
      c_line = strsplit(s_line, '\t');
      
      d_strike = str2double(c_line{1, 1});
      if (d_strike >= s_min_strike) && (d_strike <= s_max_strike)
        m_tmp(i, 1) = d_maturity;                % maturity
        m_tmp(i, 2) = d_strike;                  % strike
        m_tmp(i, 3) = str2double(c_line{1, 11}); % price
        
        i = i + 1;
      end
      
  end

  fclose(f_data);

  % return
  disp(['successfully loaded data for maturity ' datestr(i_mat)]);
  opt_table = m_tmp;

end