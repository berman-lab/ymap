function [] = ChARM_v1a(projectName, workingDir,figureDir,displayBREAKS)
%% ========================================================================
% Generate graphs involved in analysis of CNV probe data from 1' Agilent
% slide.
%==========================================================================

fprintf(['\nGenerating ChARM figure from project "' projectName '" data.\n']);

% Centromere_format : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
Centromere_format = 0;
HistPlot          = true;
ChrNum            = true;
Chr_max_width     = 0.8;
show_annotations  = true;
Run_algorithms_1  = true;
   temp_figures_1 = true;
Run_algorithms_2  = true;
   temp_figures_2 = true;
Run_algorithms_3  = false;

%%=========================================================================
% Load common_CNV file for project : 'CNVplot2', 'genome_CNV'.
%--------------------------------------------------------------------------
dataFile = [workingDir 'matlab_dir/' projectName '.common_CNV.mat'];
fprintf(['\nLoading common_CNV file for "' projectName '" : ' dataFile '\n']);
load(dataFile);

vars = who('-file',dataFile)

genome   = genome_CNV;

%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines annotation locations in bp.
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,figureDir,genome);            

for i = 1:length(chr_sizes)
    chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
end;
for i = 1:length(centromeres)
    cen_start(centromeres(i).chr) = centromeres(i).start;
    cen_end(centromeres(i).chr)   = centromeres(i).end;
end;
if (length(annotations) > 0)
    fprintf(['\nAnnotations for ' genome '.\n']);
    for i = 1:length(annotations)
        annotation_chr(i)       = annotations(i).chr;
        annotation_type{i}      = annotations(i).type;
        annotation_start(i)     = annotations(i).start;
        annotation_end(i)       = annotations(i).end;
        annotation_fillcolor{i} = annotations(i).fillcolor;
        annotation_edgecolor{i} = annotations(i).edgecolor;
        annotation_size(i)      = annotations(i).size;
        fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' num2str(annotations(i).end) ':' annotations(i).fillcolor ':' annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
    end;
end;
for i = 1:length(figure_details)
    if (figure_details(i).chr == 0)
        key_posX   = figure_details(i).posX;
        key_posY   = figure_details(i).posY;
        key_width  = figure_details(i).width;
        key_height = figure_details(i).height;
    else
        chr_id    (figure_details(i).chr) = figure_details(i).chr;
        chr_label {figure_details(i).chr} = figure_details(i).label;
        chr_name  {figure_details(i).chr} = figure_details(i).name;
        chr_posX  (figure_details(i).chr) = figure_details(i).posX;
        chr_posY  (figure_details(i).chr) = figure_details(i).posY;
        chr_width (figure_details(i).chr) = figure_details(i).width;
        chr_height(figure_details(i).chr) = figure_details(i).height;
    end;
end;
num_chr = length(chr_size);
%bases_per_bin			= 5000;
bases_per_bin                   = max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;


%% ###################################################################################################
%===================================================================================================== 
% Implementation of ChARM algorithms.
%-----------------------------------------------------------------------------------------------------
% default values from ChARMRunner.java.   Values in brackets represent values chosen for Seq data.
median_window_width = 5;         % Median filter window size.        [9]
smooth_window_width = 5;         % Smoothing filter window size.     [9]

adjacent_dist       = 2;         % The minimum index (in genes) between any two initial edge estimates from filtering stage.
init_snr            = 0.5;       % Edge removal SNR criteria BEFORE edge convergence.
max_snr             = 0.75;      % Edge removal SNR criteria AFTER edge convergence.
percent_edges       = 0.2;       % The number of initial edge estimates reported from filtering stage: num_estimates = percent_edges*chromosome_length.
num_permutations    = 200;       % Number of permutations required for mean permutation test (with Gaussian approximation).
merge_snr_thresh    = 1.5;       % (?)

% default values from Edge.java
max_roi             = 20;        % The maximum allowable radius of influence (i.e. the number of genes on either side of an edge that influence edge placement).
percent_window_size = 0.5;       % The radius of influence is determined as follows: min(percent_window_size*window_size,max_ROI).
start_likelihood    = -1000000;  % Starting likelihood value for EM iterations.
min_window          = 2;         % Minimum allowable window size (windows smaller than this are discarded during EM iterations).



%% ###################################################################################################
%===================================================================================================== 
% Implementation of ChARM algorithm #1, "an edge detection filter that identifies points on chromosomes
%	where potential aneuploidies start or end".
%-----------------------------------------------------------------------------------------------------
if (Run_algorithms_1 == true)
    %% =================================================================================================== 
    % Median filter.
    %-----------------------------------------------------------------------------------------------------
    window_width = median_window_width;
    window_halfwidth = (window_width-1)/2;
    fprintf('\nMedian Filter');
    for chr = 1:num_chr
        fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(CNVplot2{chr})) ]);
        for data = 1:length(CNVplot2{chr})
            window_start   = max(data-window_halfwidth, 1);
            window_end     = min(data+window_halfwidth, length(CNVplot2{chr}));
            window         = CNVplot2{chr}(window_start:window_end);
            if (window_start == 1)
                if (length(window) < window_width)
                    for jj = 1:(window_width - length(window))
                         window = [CNVplot2{chr}(1) window];
                    end;
                end;
            elseif (window_end == length(CNVplot2{chr}))
                if (length(window) < window_width)
                    for jj = 1:(window_width - length(window))
                        window = [window CNVplot2{chr}(end)];
                    end;
                end;
            end;
            CNV_median{chr,1}(data) = median(window);
            CNV_median{chr,2}(data) = bases_per_bin;
        end;
	% fprintf(['\n' num2str(CNV_median{chr,1}) ':' num2str(CNV_median{chr,2}) '\n']);
    end;

    %% =================================================================================================== 
    % Smoothing filter.
    %-----------------------------------------------------------------------------------------------------
    fprintf('\nSmoothing Filter');
    window_width = smooth_window_width;
    window_halfwidth = (window_width-1)/2;
    for chr = 1:num_chr
        fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(CNV_median{chr,1})) ]);
        for data = 1:length(CNV_median{chr,1})
            window_start   = max(data-window_halfwidth, 1);
            window_end     = min(data+window_halfwidth, length(CNVplot2{chr}));
            window         = CNV_median{chr,1}(window_start:window_end);
            if (window_start == 1)
                if (length(window) < window_width)
                    for jj = 1:(window_width - length(window))
                        window = [CNV_median{chr,1}(1) window];
                    end;
                end;
            elseif (window_end == length(CNV_median{chr,1}))
                if (length(window) < window_width)
                    for jj = 1:(window_width - length(window))
                        window = [window CNV_median{chr,1}(end)];
                    end;
                end;
            end;
            CNV_smoothing{chr,1}(data) = sum(window)/window_width;
	    CNV_smoothing{chr,2}(data) = bases_per_bin;
        end;
    end;

    %% =================================================================================================== 
    % Differentiation filter.
    %-----------------------------------------------------------------------------------------------------
    fprintf('\nDifferentiation Filter');
    for chr = 1:num_chr
        fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(CNV_smoothing{chr,1})) ]);
        for data = 1:length(CNV_smoothing{chr,1})
            window_start   = max(data-1, 1);
            window_end     = min(data+1, length(CNV_smoothing{chr,1}));
	    if (window_start == 1)
	        window(1) = CNV_smoothing{chr,1}(data  );
	        window(2) = CNV_smoothing{chr,1}(data+1);
	    elseif (window_end == length(CNV_smoothing{chr,1}))
	        window(1) = CNV_smoothing{chr,1}(data-1);
	        window(2) = CNV_smoothing{chr,1}(data  );
	    else
	        window(1) = CNV_smoothing{chr,1}(data-1);
	        window(2) = CNV_smoothing{chr,1}(data+1);
	    end;
            CNV_differentiation{chr,1}(data) = (window(2)-window(1))/2*8;
	    CNV_differentiation{chr,2}(data) = bases_per_bin;
        end;
    end;

    %% ===================================================================================================
    % Find local maxima/minima.
    %-----------------------------------------------------------------------------------------------------
    fprintf('\nFinding local maxima & minima');
    for chr = 1:num_chr
        [pks1{chr},locs1{chr}] = findpeaks( CNV_differentiation{chr,1});
	[pks2{chr},locs2{chr}] = findpeaks(-CNV_differentiation{chr,1});
        fprintf(['Peak positions on chr : "' num2str(chr) '"\n']);
        fprintf('\t[');
        for edge = 1:(length(locs1{chr})-1)
            fprintf([num2str(locs1{chr}(edge)) ', ']);
	    if (mod(edge,30) == 0);    fprintf('\n\t');    end;
        end;
        edge = length(locs1{chr});
        fprintf([num2str(locs1{chr}(edge)) '] (' num2str(length(locs1{chr})) ')\n']);
	locs{chr} = sort([1 locs1{chr} locs2{chr} length(CNVplot2{chr,1})]);
    end;
end;

%% Assigning raw and filtered data to convenient names for later use.
data1 = CNVplot2;
for chr = 1:num_chr
    data2{chr} = CNV_median{chr,1};                 % ./CNV_median{chr,2};
    data3{chr} = CNV_smoothing{chr,1};              % ./CNV_smoothing{chr,2};
    data4{chr} = CNV_differentiation{chr,1} + 1;    % ./CNV_differentiation{chr,2};
end;

if (temp_figures_1 == true)
    maxY = 2;
    fig = figure(1);    dataShow = data1;
    set(gcf, 'Position', [0 70 1024 600]);
    for chr = 1:num_chr
	left   = chr_posX(chr);    bottom = chr_posY(chr);
	width  = chr_width(chr);   height = chr_height(chr);
	subplot('Position',[left bottom width height]);
	hold on;
	c_ = [0 0 0];
	fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
	for i = 1:length(dataShow{chr});
	    x_ = [i i i-1 i-1];
	    if (dataShow{chr}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr}(i);    end;
	    startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);
	    set(f,'linestyle','none');
	end;
	x2 = chr_size(chr)*chr_length_scale_multiplier;
	plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
	hold off;
	xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
	set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''});
    end;
    saveas(fig,[figureDir projectName '.ChARM_test.1.eps'], 'epsc');
    delete(fig);

    fig = figure(1);    dataShow = data2;
    set(gcf, 'Position', [0 70 1024 600]);
    for chr = 1:num_chr
        left   = chr_posX(chr);    bottom = chr_posY(chr);
        width  = chr_width(chr);   height = chr_height(chr);
        subplot('Position',[left bottom width height]);
        hold on;
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
        for i = 1:length(dataShow{chr});
            x_ = [i i i-1 i-1];
            if (dataShow{chr}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr}(i);    end;
            startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);   
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
        hold off;
        xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
        set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''}); 
    end;
    saveas(fig,[figureDir projectName '.ChARM_test.2.eps'], 'epsc');
    delete(fig);

    fig = figure(1);    dataShow = data3;
    set(gcf, 'Position', [0 70 1024 600]);
    for chr = 1:num_chr
        left   = chr_posX(chr);    bottom = chr_posY(chr);
        width  = chr_width(chr);   height = chr_height(chr);
        subplot('Position',[left bottom width height]);
        hold on;
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
        for i = 1:length(dataShow{chr});
            x_ = [i i i-1 i-1];
            if (dataShow{chr}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr}(i);    end;
            startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);   
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
        hold off;
        xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
        set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''}); 
    end;
    saveas(fig,[figureDir projectName '.ChARM_test.3.eps'], 'epsc');
    delete(fig);

    fig = figure(1);    dataShow = data4;
    set(gcf, 'Position', [0 70 1024 600]);
    for chr = 1:num_chr
        left   = chr_posX(chr);    bottom = chr_posY(chr);
        width  = chr_width(chr);   height = chr_height(chr);
        subplot('Position',[left bottom width height]);
        hold on;
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
        for i = 1:length(dataShow{chr});
            x_ = [i i i-1 i-1];
            if (dataShow{chr}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr}(i);    end;
            startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);   
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
        hold off;
        xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
        set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''}); 
    end;
    saveas(fig,[figureDir projectName '.ChARM_test.4.eps'], 'epsc');
    delete(fig);

    fig = figure(1);    dataShow = data4;
    set(gcf, 'Position', [0 70 1024 600]);
    for chr = 1:num_chr
        left   = chr_posX(chr);    bottom = chr_posY(chr);
        width  = chr_width(chr);   height = chr_height(chr);
        subplot('Position',[left bottom width height]);
        hold on;
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
        for i = 1:length(dataShow{chr});
            x_ = [i i i-1 i-1];
            if (dataShow{chr}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr}(i);    end;
            startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
	for edge = 1:length(locs1{chr})
	    plot([locs1{chr}(edge) locs1{chr}(edge)], [0 maxY],'color',[0 0 1]);
	end;
	for edge = 1:length(locs2{chr})
	    plot([locs2{chr}(edge) locs2{chr}(edge)], [0 maxY],'color',[1 0 0]);
	end;
        hold off;
        xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
        set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''});
    end;
    saveas(fig,[figureDir projectName '.ChARM_test.5.eps'], 'epsc');
    delete(fig)
end;


%% ###################################################################################################
%=====================================================================================================
% Implementation of ChARM algorithm #2, "an EM (Expectation Maximization) based edge=placement
%	algorithm that statistically optimizes these start and end locations".
%-----------------------------------------------------------------------------------------------------
if (Run_algorithms_2 == true)
    %% ===================================================================================================
    % Expectation-Maximization Edge-Placement Algorithm
    %-----------------------------------------------------------------------------------------------------
    for chr = 1:num_chr
	% Initial local peaks found in first section of algorithm.
        position = locs{chr};    % locations of edges for this chromosome.
	data     = data1{chr};   % data to be examined for this chromosome.

        for edge = 1:length(position)
	    roi_left(edge)                    = position(edge)-max_roi;
	    roi_right(edge)                   = position(edge)+max_roi;
	    left_post{edge}                   = zeros(1,ceil(roi_right(edge))-floor(roi_left(edge)));
	    left_post{edge}(1:max_roi)        = 1;
	    right_post{edge}                  = zeros(1,ceil(roi_right(edge))-floor(roi_left(edge)));
	    right_post{edge}((max_roi+1):end) = 1;
        end;

        for t = 1:num_permutations
	    % Calculate Conditional probabilities that a data point is in the left vs. right distributions around each edge.
	    for edge = 1
		% Left edge of chromosome.
	    end;
	    for edge = 2:(length(position)-1);
		% max_roi             = 20;
		% percent_window_size = 0.5;    % The radius of influence is determined as follows: min(percent_window_size*window_size,max_ROI).
		L_windowSize   = min(floor(percent_window_size*(position(edge  )-position(edge-1))),max_ROI);
		R_windowSize   = min(floor(percent_window_size*(position(edge+1)-position(edge  ))),max_ROI);
		L_dist         = data((position(edge)-L_windowSize):(position(edge)-1));   %
		R_dist         = data((position(edge)+1):(position(edge)+R_windowSize));   % anything outside these ranges has a zero chance of being on the wrong side of the edge.
		L_distAve      = ave(L_dist);
		R_distAve      = ave(R_dist);
		L_distStdev    = std(L_dist);
		R_distStdev    = std(R_dist);
		cP_L_dist_is_L = normcdf(L_dist,L_distAve,L_distStdev);
		cP_L_dist_is_R = normcdf(L_dist,R_distAve,R_distStdev);
		cP_R_dist_is_L = normcdf(R_dist,L_distAve,L_distStdev);
		cP_R_dist_is_R = normcdf(R_dist,R_distAve,R_distStdev);
	    end;
	    for edge = length(position{chr})
		% Right edge of chromosome.
	    end;

	    % Summarize all individual posterior (post) probabilities for left and right ends.
	    left_prev_post_tot  = zeros(1,length(CNVplot2{chr}));
	    right_prev_post_tot = zeros(1,length(CNVplot2{chr}));

	    left_prev_post  = left_post;
	    right_prev_post = right_post;
	    for edge = 1:length(position)
	        Lcount = 0;
	        for pos = max(1,floor(roi_left(edge))):min(length(CNVplot2{chr}),position(edge))
		    Lcount = Lcount+1;
		    left_prev_post_tot(pos) = left_prev_post_tot(pos)+left_post{edge}(Lcount);
	        end;

	        Rcount = 0;
	        for pos = position(edge):min(length(CNVplot2{chr}),ceil(roi_right(edge)))
		    Rcount = Rcount+1;
		    right_prev_post_tot(pos) = right_prev_post_tot(pos)+right_post{edge}(Rcount);
	        end;
	    end;

	    left_prev_post_ave  = sum(left_prev_post_tot) /length(position);
	    right_prev_post_ave = sum(right_prev_post_tot)/length(position);

	    %% -----------------------------------------------
	    % Update membership (E-step)                     |
	    %-------------------------------------------------
	    for edge = 1:length(position)
	        % Find the position of adjacent edges, either at the chr ends or internal edges.
	        if (edge == 1)
		    leftEdgePosition(edge) = 1;
	        else
		    leftEdgePosition(edge) = position(edge-1);
	        end;
	        if (edge == length(position))
		    rightEdgePosition(edge) = position(length(position));
	        else
		    rightEdgePosition(edge) = position(edge+1);
	        end;
	        position(edge) = position(edge);

	        % Sizes of windows at left and right.
	        left_window_size(edge)  = position(edge) - leftEdgePosition(edge);
	        right_window_size(edge) = rightEdgePosition(edge) - position(edge);

	        % Determine influence window, the sizes of regions of interest for left and right sides of current edge position.
	        % Initialize the ROI size as [max_roi], then set it as the shorter of the two legs to either an adjacent chr end or halfway to an adjacent edge.
	        roi_size(edge) = max_roi;
	        if ((position(edge) - left_window_size(edge)) <= 1)
		    half = left_window_size(edge);
	        else
		    half = left_window_size(edge)*percent_window_size;
	        end;
	        if (half < roi_size(edge))
		    roi_size(edge) = half;
	        end;
	        if ((position(edge) + right_window_size(edge)) >= length(CNVplot2{chr}))
		    half = right_window_size(edge)-1;
	        else
		    half = right_window_size(edge)*percent_window_size;
	        end;
	        if (half < roi_size(edge))
		    roi_size(edge) = half;
	        end;

	        % define left and right bounds of the ROI around the current edge.
	        roi_left(edge)     = position(edge) - roi_size(edge);
	        roi_right(edge)    = position(edge) + roi_size(edge) + 1;

%	        if (edge == 1)
%		    fprintf('left\tpos\tright\t[ROI]\n');
%	        end;
%	        fprintf([num2str(leftEdgePosition(edge)) '\t' num2str(position(edge)) '\t' num2str(rightEdgePosition(edge)) '\t[' num2str(roi_left(edge)) ':' num2str(roi_right(edge)) ']\n']);

	        % Determine the Means and variances for the data within the left and right ROIs.
	        leftROIdata{edge}  = CNVplot2{chr}( max(1,floor(roi_left(edge))) : min(length(CNVplot2{chr}),position(edge))        );
	        rightROIdata{edge} = CNVplot2{chr}( max(1,position(edge))        : min(length(CNVplot2{chr}),ceil(roi_right(edge))) );
	        left_mean(edge)    = mean(leftROIdata{edge});
	        right_mean(edge)   = mean(rightROIdata{edge});
	        left_var(edge)     = var(leftROIdata{edge});
	        right_var(edge)    = var(rightROIdata{edge});

	        % Calculate the conditional (cond) probability of observing the data points in the zone of influence in either distribution.
	        left_cond{edge}    = pdf('normal', CNVplot2{chr}( max(1,floor(roi_left(edge))) : min(length(CNVplot2{chr}),ceil(roi_right(edge))) ), left_mean(edge),  left_var(edge) );
	        right_cond{edge}   = pdf('normal', CNVplot2{chr}( max(1,floor(roi_left(edge))) : min(length(CNVplot2{chr}),ceil(roi_right(edge))) ), right_mean(edge), right_var(edge));

	        % Calculate the posterior (post) probability of observing the data poitn ins the zone of influence in either distribution.
	        left_post{edge}    = (left_cond{edge} .*left_prev_post_ave )./((left_cond{edge}*left_prev_post_ave)+(right_cond{edge}*right_prev_post_ave));
	        right_post{edge}   = (right_cond{edge}.*right_prev_post_ave)./((left_cond{edge}*left_prev_post_ave)+(right_cond{edge}*right_prev_post_ave));
	    end;

	    %% ----------------------------------------------
	    % Mean and variance computation (M-step 1)      |
	    %------------------------------------------------
	    for edge = 1:length(position)
	        left_post_tot  = zeros(1,length(CNVplot2{chr}));
	        right_post_tot = zeros(1,length(CNVplot2{chr}));
	    
	        Lcount = 0;
	        for pos = max(1,floor(roi_left(edge))):position(edge)
		    Lcount = Lcount+1;
		    left_post_tot(pos)  = left_post_tot(pos)+left_post{edge}(Lcount);
	        end;
	        Rcount = 0;
	        for pos = position(edge):min(length(CNVplot2{chr}),ceil(roi_right(edge)))
		    Rcount = Rcount+1;
		    right_post_tot(pos) = right_post_tot(pos)+right_post{edge}(Rcount);
	        end;

	        left_mean_numerator(edge)  = sum(left_post_tot .*CNVplot2{chr});
	        right_mean_numerator(edge) = sum(right_post_tot.*CNVplot2{chr});
	        left_denominator(edge)     = sum(left_post_tot);
	        right_denominator(edge)    = sum(right_post_tot);
	        ML_left_mean(edge)         = left_mean_numerator/left_denominator;
	        ML_right_mean(edge)        = right_mean_numerator/right_denominator;
	        left_var_numerator(edge)   = sum((CNVplot2{chr}-ML_left_mean(edge) ).^2.*left_post_tot );
	        right_var_numerator(edge)  = sum((CNVplot2{chr}-ML_right_mean(edge)).^2.*right_post_tot);
	        ML_left_var(edge)          = left_var_numerator /left_denominator;
	        ML_right_var(edge)         = right_var_numerator/right_denominator;
	    end;

	    %% ----------------------------------------------
	    % Edge adjustment (M-step 2)                    |
	    %------------------------------------------------
	    for edge = 1:length(position)
	        surprise = zeros(1,ceil(roi_right(edge))-floor(roi_left(edge)));
	        for i = 0:(ceil(roi_right(edge))-floor(roi_left(edge)))
		    left_surprise  = sum( log10(left_post_tot(floor(roi_left(edge)):(ceil(roi_left(edge))+i))) );
		    right_surprise = sum( log10(right_post_tot((floor(roi_left(edge))+i+1):ceil(roi_right(edge)))) );
		    surprise(i+1)  = -(left_surprise+right_surprise);
	        end;
	        [surprise_min,surprise_min_index] = min(surprise);
	        position_new(edge) = floor(roi_left(edge))+surprise_min_index-1;

	        % Find the position of adjacent edges, either at the chr ends or internal edges.
	        if (edge == 1)
		    leftEdgePosition(edge) = 0;
	        else
		    leftEdgePosition(edge) = position(edge-1);
	        end;
	        if (edge == length(position))
		    rightEdgePosition(edge) = position(length(position))
	        else
		    rightEdgePosition(edge) = position(edge+1);
	        end;

	        % Sizes of windows at left and right.
	        left_window_size(edge)  = position_new(edge)      - leftEdgePosition(edge);
	        right_window_size(edge) = rightEdgePosition(edge) - position_new(edge);

	        % Determine influence window, the sizes of regions of interest for left and right sides of current edge position.
	        % Initialize the ROI size as [max_roi], then set it as the shorter of the two legs to either an adjacent chr end or halfway to an adjacent edge.
	        roi_size(edge) = max_roi;
	        if ((position_new(edge) - left_window_size(edge)) <= 0)
		    half = left_window_size(edge)
	        else
		    half = left_window_size(edge)*percent_window_size;
	        end;
	        if (half < roi_size(edge))
		    roi_size(edge) = half;
	        end;
	        if ((position_new(edge) + right_window_size(edge)) >= length(CNVplot2{chr}))
		    half = right_window_size(edge)-1;
	        else
		    half = right_window_size(edge)*percent_window_size;
	        end;
	        if (half < roi_size(edge))
		    roi_size(edge) = half;
	        end;

	        % define left and right bounds of the ROI around the current edge. 
	        roi_left_new(edge)  = position_new(edge) - roi_size(edge);
	        roi_right_new(edge) = position_new(edge) + roi_size(edge) + 1;

	        % Determine the Means and variances for the data within the left and right ROIs.
	        leftROIdata_new{edge}  = CNVplot2{chr}( max(1,floor(roi_left_new(edge))) : min(length(CNVplot2{chr}),position_new(edge))        );
	        rightROIdata_new{edge} = CNVplot2{chr}( max(1,position_new(edge))        : min(length(CNVplot2{chr}),ceil(roi_right_new(edge))) );
	        left_mean_new(edge)    = mean(leftROIdata_new{edge});
	        right_mean_new(edge)   = mean(rightROIdata_new{edge});
	        left_var_new(edge)     = var(leftROIdata_new{edge});
	        right_var_new(edge)    = var(rightROIdata_new{edge});

	        % Edge adjustments.
	        edge_adjustment(edge) = position(edge) - position_new(edge);
	    end;

	    % Calculate average change in edge positions.
	    delta_edge_mean = mean(edge_adjustment);

	    %% ----------------------------------------------
            % Window similarity test                        |
            %------------------------------------------------
            for edge = 1:length(position)
	        % Find the position of adjacent edges, either at the chr ends or internal edges.
	        if (edge == 1)
		    leftEdgePosition(edge) = 0;
	        else
		    leftEdgePosition(edge) = position(edge-1);
	        end;
	        if (edge == length(position))
		    rightEdgePosition(edge) = position(length(position))
	        else
		    rightEdgePosition(edge) = position(edge+1);
	        end;

	        % Sizes of windows at left and right.
	        left_window_size(edge)  = position(edge)          - leftEdgePosition(edge);
	        right_window_size(edge) = rightEdgePosition(edge) - position(edge);

	        % Medians of windows at left and right.
	        leftWINDOWdata            = CNVplot2{chr}( max(1,floor(leftEdgePosition(edge))) : min(length(CNVplot2{chr}),position(edge))                );
	        rightWINDOWdata           = CNVplot2{chr}( max(1,position(edge))                : min(length(CNVplot2{chr}),ceil(rightEdgePosition(edge))) );
	        left_window_median(edge)  = median(leftWINDOWdata);
	        right_window_median(edge) = median(rightWINDOWdata);
	        left_window_mean(edge)    = mean(leftWINDOWdata);
	        right_window_mean(edge)   = mean(rightWINDOWdata);
	        left_window_var(edge)     = var(leftWINDOWdata);
	        right_window_var(edge)    = var(rightWINDOWdata);

	        % Calculate sums for denominator of window similarity test.
	        left_sum  = sum(abs( CNVplot2{chr}( max(1,floor(leftEdgePosition(edge))) : min(length(CNVplot2{chr}),position(edge))                )-left_window_median(edge)  ));
	        right_sum = sum(abs( CNVplot2{chr}( max(1,position(edge))                : min(length(CNVplot2{chr}),ceil(rightEdgePosition(edge))) )-right_window_median(edge) ));

	        % Calculate Signal-Noise-Ratio as defined in ChARM paper.
	        SNR_numerator   = abs(left_window_median(edge) - right_window_median(edge));
	        SNR_denominator = (left_sum + right_sum)/(left_window_size(edge) + right_window_size(edge));
	        SNR_1(edge)     = SNR_numerator/SNR_denominator;

	        % ChARM java source uses a different SNR function than is defined in the paper.
	        % Primarily, the source uses 'mean's instead of 'median's.
	        %
	        %    double snr = Math.abs(leftMean-rightMean)/Math.sqrt(leftVar/currLeft.getSize() + rightVar/currRight.getSize());
	        %
	        % The threshold in the source is defined as "MERGE_SNR_THRESH = 1.5"
	        %     SNR values less than this lead to the edge being removed from further consideration.
	        % The paper uses a threshold dependent on the current delta_edge_mean value.

	        % Calculate Signal-Noise-Ratio as defined in ChARM software.
	        SNR_2(edge) = abs(left_window_mean-right_window_mean)/sqrt( left_window_var/left_window_size(edge) + right_window_var/right_window_size(edge) );

	        % Remove edge if less if SNR is less than threshold.
	        position_new   = position;
	        roi_left_new   = roi_left;
	        roi_right_new  = roi_right;
	        left_post_new  = left_post;
	        right_post_new = right_post;
	        MERGE_SNR_THRESH = 1.5;
	        if (SNR_2(edge) < MERGE_SNR_THRESH)
		    position_new(edge)   = [];
		    roi_left_new(edge)   = [];
		    roi_right_new(edge)  = [];
		    left_post_new(edge)  = [];
		    right_post_new(edge) = [];
	        end;
	    end;

	    %% ----------------------------------------------
	    % Update edge and associated lists              |
	    %------------------------------------------------
	    position   = position_new;
	    roi_left   = roi_left_new;
	    roi_right  = roi_right_new;
	    left_post  = left_post_new;
	    right_post = right_post_new;
        end;

        locs{chr} = position;
    end;
end;


%% ###################################################################################################
%=====================================================================================================
% Implementation of ChARM algorithm #3, "a window significance test that determines whether predicted
%	amplifications and deletions are statistically significant or are artifacts of noise".
%-----------------------------------------------------------------------------------------------------
if (Run_algorithms_3 == true)
	%
	% This part of the ChARM algorithm examines the statistical significance of windows with CNV
	% changes.   For our purposes, this section is less useful... but may be implemented later.
	%

    %% ###################################################################################################
    %=====================================================================================================
    % Save chromosome breakpoints as 'links_dir/pileup_dir/'$projectName'_segments.txt'
    % File format is as follows :
    %	Any number of header lines starting with '#' characters.
    %	One or more lines of breakpoint data (two tab-delimited columns).
    %		1) chromosome id.
    %		2) breakpoint in percent along length of chromosome.
    %-----------------------------------------------------------------------------------------------------
    segmentFile = [workingDir '/pileup_dir/' projectName '_segments.txt'];
    segmentsFID = fopen(segmentFile,'w');
    fprintf(segmentsFID,'# chrID       : chromosome id defined in "chr" column of "figure_definitions.txt" for the in-use genome.\n');
    fprintf(segmentsFID,'# break-point : percentage of chromosome length where aneuploidy breakpoint is found.\n');
    fprintf(segmentsFID,'#-----------------------------------------------------------------------------------------------------------\n');
    fprintf(segmentsFID,'# chrID\t break-point\n');
    fprintf(['num edges = ' num2str(length(locs)) '\n']);
    for chr = 1:length(locs)
        for edge = 1:(length(locs{chr}))
	    fprintf(segmentsFID,[num2str(chr) '\t' num2str(locs{chr}(edge)) '\n']);
        end;
    end;
    fclose(segmentsFID);
end;

end
