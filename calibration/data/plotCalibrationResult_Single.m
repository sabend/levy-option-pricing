% Plots observed and model option prices on strikes.
%
function [] = plotCalibrationResult_Single(c_market_data, c_calibrated_surface, b_plot_iters)

    i_count = numel(c_market_data);

    % prepare
    figure
    if ~b_plot_iters
        m_colors = jet(i_count);
    else
        i_n = 3;
        vd_b = transpose([(0:1:i_n-1)/i_n, ones(1, i_n)]);
        vd_g = zeros(2*i_n, 1);
        vd_r = transpose([ones(1, i_n), (i_n-1:-1:0)/i_n]);
        m_colors = flipud([vd_r, vd_g, vd_b]);
    end
    
    % gather maturities
    c_maturities = cell(i_count, 1);
    for i = 1 : i_count
        
        c_maturities{i, 1} = num2str(c_market_data{i, 1}(1, 1));
        
    end

    % add market prices
    for i = 1 : i_count

        m_tmp = c_market_data{i, 1};
      
        hold on
        scatter(m_tmp(:, 2), m_tmp(:, 3), [], m_colors(i, :));
        hold off

    end

    xlabel('strikes')
    ylabel('prices')
    h_legend = legend(c_maturities);
    set(h_legend, 'location', 'northeast')

    % add model prices
    i_start = numel(c_calibrated_surface);
    if b_plot_iters i_start = 1; end
    
    i_iters   = numel(c_calibrated_surface);
    vi_good   = [1:5 i_iters];
    i_no_good = numel(vi_good);
    vs_align  = {'right', 'left', 'right', 'right', 'right', 'right'};
    
    k = 0;
    for j = i_iters : -1 : i_start
        
        if isempty(find(j==vi_good, 1)) continue; end
 
        c_tmp_surface = c_calibrated_surface{j};
        for i = 1 : numel(c_tmp_surface)

            m_tmp = c_tmp_surface{i, 1};

            hold on
            %d_shadow = 1 - 0.7*k/i_no_good; 
            d_shadow = 1;
            vd_x = m_tmp(:, 2);
            vd_y = m_tmp(:, 3);
            if ~b_plot_iters
                scatter(vd_x, vd_y, [], d_shadow*m_colors(i, :), '+');
            else
                scatter(vd_x, vd_y, [], d_shadow*m_colors(k+1, :), '+');
            end
            hold off
            
            if b_plot_iters
                i_mid = ceil(numel(vd_x)/2);
                text(vd_x(i_mid), vd_y(i_mid), [num2str(j) ' '], ...
                    'HorizontalAlignment', vs_align{k+1}, 'FontSize', 10)
            end

        end
        k = k + 1;
        
    end

end

