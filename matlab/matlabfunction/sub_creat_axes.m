function h_axes=sub_creat_axes(nrow,nclm,margin,gap)
    % creat multiple axes
    % by zhong, 2014/6/19
    % modified by guo, 15/9/12
    % output: h_axes(k,l), handle of each axis
    h_N = nrow; % row number
    w_N = nclm; % column number
    h_gap = gap(1); % gap between rows
    w_gap = gap(2);  % gap between columns
    h_marg = margin(1:2); % top and bottom margins
    w_marg = margin(3:4); % left and right margins
    h_ax = (1 - sum(h_marg) - (h_N-1)*h_gap)./h_N;
    w_ax = (1 - sum(w_marg) - (w_N-1)*w_gap)./w_N;
    h_axes = zeros(h_N, w_N);
    y_pos = 1 - h_ax - h_marg(1);
    for h_i = 1:h_N
        x_pos = w_marg(1);
        for w_i = 1:w_N
            h_axes(h_i, w_i) = axes('Units','normalized', ...
                'Position',[x_pos y_pos w_ax h_ax], ...
                'XTickLabel','');
            x_pos = x_pos+w_ax+w_gap;
        end
        y_pos = y_pos - h_ax - h_gap; 
    end
end
