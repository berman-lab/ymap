%--------------------------------------------------------------------------
% alters Xscale labelling so we can have histogram with resolution<1.
%--------------------------------------------------------------------------
function [] = scaleXaxis(multiplier)
    original    = gca;
    position    = get(original, 'Position');
    temporary   = axes('Position', position, 'Visible', 'off');
    xlim_       = get(original, 'XLim')/multiplier;
    set(temporary, 'XLim', xlim_);
    xtick       = get(temporary, 'XTick');
    xticklabel  = get(temporary, 'XTickLabel');
    set(original, 'XTick', xtick*multiplier+0.5, 'XTickLabel', xticklabel);
end