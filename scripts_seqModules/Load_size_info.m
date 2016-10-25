function [linear_fig_height,linear_fig_width,linear_left_padding,linear_chr_gap,linear_chr_max_width,linear_height...
            ,linear_base,rotate,linear_chr_font_size,linear_axis_font_size,linear_gca_font_size,stacked_fig_height,stacked_fig_width,...
            stacked_chr_font_size,stacked_title_size,stacked_axis_font_size,gca_stacked_font_size,stacked_copy_font_size,max_chrom_label_size] = Load_size_info(chr_in_use,num_chrs,chr_label,chr_size)
fprintf('\n---------------------------------Load_size_info.m started---------------------------------------------------\n');

%% calculating data used to determine size 

% calculate the number of chrs used
num_chrs_used = 0;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
        num_chrs_used = num_chrs_used + 1;
	end;
end;
fprintf('num_chrs_used - %d\n',num_chrs_used);

% calculating the maximum length of chromosome label size for the height of linear figure
max_chrom_label_size = 1;
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        if (numel(chr_label{chr}) > numel(max_chrom_label_size))
            max_chrom_label_size = numel(chr_label{chr});
        end;
    end;
end;
fprintf('max_chrom_label_size - %d\n',max_chrom_label_size);

%% global definitions
% the default dpi used by the system, in matlab no screen mode before 2014 it's 72dpi, and in later versions 150dpi
if verLessThan('matlab','8.4')
    % -- Code to run in MATLAB R2014a and earlier here --
    system_dpi = 72; 
else
    % -- Code to run in MATLAB R2014b and later here --
    system_dpi = 150; 
end
fprintf('using %d dpi\n',system_dpi);

%% linear figure
fprintf('-----------------------------------Linear Figure -----------------------------------------------\n');
% setting size for linear figure
linear_fig_plot_height = 130; % the height of each chromosom in the figure in px (including y-axis)
linear_fig_charc_height = 20; % the size to be allocated in px for each character in the label 

linear_fig_height_px = linear_fig_plot_height + linear_fig_charc_height*max_chrom_label_size; % the total height of the linear figure in px
% normalize height according to dpi
linear_fig_height = linear_fig_height_px / system_dpi;

linear_fig_width_px = 2400;
linear_fig_width = linear_fig_width_px / system_dpi;

% base value in matlab (the scaling here is from 0 to 1 and represent
% relative position
linear_left_padding    = 0.02;               % left margin
linear_right_padding   = 0.02;               % right margin
linear_total_gap     = 0.07;               % the size in precentage for total gap accross all figure
linear_cartoon_height_px = 111;            % the size of the chromosome cartoon in pixles (without y-axis) for proper scaling
linear_chr_gap  = linear_total_gap/(num_chrs-1);  % gaps between chr subfigures.
linear_chr_max_width = 1 - linear_total_gap - linear_left_padding - linear_right_padding;  % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
linear_height        = linear_cartoon_height_px/ linear_fig_height_px; % the size 
linear_base          = 0.1;

% setting rotation
% if there are more than5 charcters in label 90 degrees rotation, else
% rotating accroding to lowest chromosome size
% 45 degrees boundries 
lower_boundary = 0.10;
upper_boundary = 0.25;
% removing zero enteries in chromosom sizes
chr_size_cleaned = chr_size(chr_size ~= 0);
% calculate ration between smallest chromosome size to highest
ratio = min(chr_size_cleaned)/max(chr_size_cleaned);
rotate = 0;
if (max_chrom_label_size > 5)
    rotate = 90;
elseif (lower_boundary <= ratio && ratio <= upper_boundary)
    rotate = 45;
elseif(ratio < lower_boundary)
    rotate = 90;
end;

% font definitions
linear_axis_font_size = 10;
linear_gca_font_size = 12;
linear_chr_font_size = 12;

fprintf('linear figure paramters:\n');
fprintf('height:%d px, width:%d px, rotate:%d\n',linear_fig_height_px,linear_fig_width_px,rotate);
fprintf('chr font size:%d, axis font size:%d px, gca font size:%d\n',linear_chr_font_size,linear_axis_font_size,linear_gca_font_size);

%% stacked figure %%
fprintf('-----------------------------------Stacked Figure -----------------------------------------------\n');

stacked_plot_height = 205; % size in pixels of each cartoon height including gap
stacked_fig_height_px = stacked_plot_height*num_chrs_used; % the total height of the linear figure in px
% normalize height according to dpi
stacked_fig_height = stacked_fig_height_px / system_dpi;

stacked_fig_width_px = 2400;
stacked_fig_width = stacked_fig_width_px / system_dpi;

% the size of the stacked title normalized to 8 chromosomes (that appear
% right), 8 is what appeared correct in candida albicans, just increasing size makes figure look not good
stacked_title_size = 18*(8/num_chrs_used);
% the size of axis values in stacked figure
stacked_axis_font_size = 10*(8/num_chrs_used);
% the size of gca font in stacked mode
gca_stacked_font_size = 12*(8/num_chrs_used);
% the size of chromosome text in stacked figure
stacked_chr_font_size = 16*(8/num_chrs_used);
% the size of the copy beside figure in stacked
stacked_copy_font_size = 20*(8/num_chrs_used);

fprintf('stacked figure paramters:\n');
fprintf('height:%d px, width:%d px, title size:%d\n',stacked_fig_height_px,stacked_fig_width_px,stacked_title_size);
fprintf('chr font size:%d, axis font size:%d px, gca font size:%d\n',stacked_chr_font_size,stacked_axis_font_size,gca_stacked_font_size);



fprintf('\n---------------------------------Load_size_info.m ended---------------------------------------------------\n');


end
