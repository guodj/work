function h=get_plot_imf_mdays(lrdays,str_imf)
    % Draw a plot of imf variation during lrdays
    imf=get_imf_mdays(lrdays);
    x=imf(:,1);
    switch str_imf
        case 'bx'
            b=imf(:,2);
        case 'bye'
            b=imf(:,3);
        case 'bze'
            b=imf(:,4);
        case 'bym'
            b=imf(:,5);
        case 'bzm'
            b=imf(:,6);
        otherwise
            error('wrong input of str_imf!');
    end
    h=plot(x,b);
end
