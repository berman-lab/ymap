%--------------------------------------------------------------------------
% alters Xscale labelling so we can have histogram with resolution<1.
%--------------------------------------------------------------------------
function [] = setXaxis(Xstart,Xend)
    range       = Xend-Xstart;
    original    = gca;
    position    = get(original, 'Position');
    temporary   = axes('Position', position, 'Visible', 'off');
    xlim_       = get(original, 'XLim');
    set(temporary, 'XLim', xlim_);
    xtick       = get(temporary, 'XTick');
    xticklabel  = num2str(str2num(get(temporary, 'XTickLabel'))+Xstart);
    set(original, 'XTick', xtick+0.5, 'XTickLabel', xticklabel);
end