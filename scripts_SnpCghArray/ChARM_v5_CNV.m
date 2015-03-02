function [] = ChARM_v5_CNV(projectName, workingDir)
%% =========================================================================================
% Analyze CNV information for copy number changes.
% Code based on algorithms described in:
% 	Accurate detection of aneuploidies in array CGH and gene expression microarray data 
%	Chad L. Myers, Maitreya J. Dunham, S.Y. Kung, Olga G. Troyanskaya (2004)
%===========================================================================================

fprintf(['\nGenerating ChARM figures from project "' projectName '" CNV data.\n']);

HistPlot          = true;
ChrNum            = true;
Chr_max_width     = 0.8;
show_annotations  = true;
   temp_figures   = true;

%%=========================================================================
% Load Common_CNV file for project : 'CNVplot2', 'genome_CNV'.
%--------------------------------------------------------------------------
dataFile = [workingDir 'Common_CNV.mat'];
fprintf(['\nLoading common_CNV file for "' projectName '" : ' dataFile '\n']);
load(dataFile);
vars          = who('-file',dataFile);
genome        = genome_CNV;
examined_data = CNVplot2;

%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines annotation locations in bp.
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1();

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
%bases_per_bin               = 5000;
bases_per_bin               = max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;


%% ###################################################################################################
%===================================================================================================== 
% Implementation of ChARM algorithms.
%-----------------------------------------------------------------------------------------------------
% default values from ChARMRunner.java.   Values in brackets represent values chosen for Seq data.
median_window_width   = 9;		% Median filter window size.        [9]
smooth_window_width   = 9;		% Smoothing filter window size.     [9]
smooth_gaussian_sigma = 4;		% Sigma factor for Gaussian-kernal smoothing.

adjacent_dist         = 2;		% The minimum index (in genes) between any two initial edge estimates from filtering stage.
init_snr              = 0.5;		% Edge removal SNR criteria BEFORE edge convergence.
max_snr               = 0.75;		% Edge removal SNR criteria AFTER edge convergence.
percent_edges         = 0.2;		% The number of initial edge estimates reported from filtering stage: num_estimates = percent_edges*chromosome_length.
num_permutations      = 20;		% Number of permutations required for mean permutation test (with Gaussian approximation). [200]
merge_snr_thresh      = 1.5;		% (?)

% default values from Edge.java
max_ROI               = 20;		% The maximum allowable radius of influence (i.e. the number of genes on either side of an edge that influence edge placement). [20]
percent_window_size   = 0.5;		% The radius of influence is determined as follows: min(percent_window_size*window_size,max_ROI).
start_likelihood      = -1000000;	% Starting likelihood value for EM iterations.
min_window            = 2;		% Minimum allowable window size (windows smaller than this are discarded during EM iterations).

edge_nearness_limit   = 2;		% Edge pairs closer than this are converted to a single edge.
SNR_threshold         = 4;		% Signal to Noise threshold for identifying side windows as the same.



%% ###################################################################################################
%===================================================================================================== 
% Implementation of ChARM algorithm #1, "an edge detection filter that identifies points on chromosomes
%	where potential aneuploidies start or end".
%-----------------------------------------------------------------------------------------------------
%% ===================================================================================================
% [sliding-box] median filter.
%-----------------------------------------------------------------------------------------------------
window_width = median_window_width;
window_halfwidth = (window_width-1)/2;

fprintf('\nMedian Filter');
for chr = 1:num_chr
	fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(examined_data{chr})) ]);
	for data = 1:length(examined_data{chr})
		window_start   = max(data-window_halfwidth, 1);
		window_end     = min(data+window_halfwidth, length(examined_data{chr}));
		window         = examined_data{chr}(window_start:window_end);
		if (window_start == 1)
			if (length(window) < window_width)
				for jj = 1:(window_width - length(window))
				window = [examined_data{chr}(1) window];
			end;
		end;
	elseif (window_end == length(examined_data{chr}))
		if (length(window) < window_width)
			for jj = 1:(window_width - length(window))
				window = [window examined_data{chr}(end)];
			end;
		end;
	end;
	CNV_median{chr}(data) = median(window);
	end;
end;

%%====================================================================================================
% [gaussian] smoothing filter
%     Original ChARM implementation used a sliding-box smoothing filter, but the use of a gaussian
%     smoothing filter reduces sharp-edge artifacts.
%-----------------------------------------------------------------------------------------------------
fprintf('\nSmoothing Filter after median');
window_width = smooth_window_width;
window_halfwidth = (window_width-1)/2;
for chr = 1:num_chr
	fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(CNV_median{chr})) ]);
	CNV_median_smoothed{chr} = smooth_gaussian(CNV_median{chr},smooth_gaussian_sigma,smooth_gaussian_sigma*16);
end;

%% =================================================================================================== 
% differentiation filter
%-----------------------------------------------------------------------------------------------------
fprintf('\nDifferentiation Filter');
for chr = 1:num_chr
	fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(CNV_median_smoothed{chr})) ]);
	for data = 1:length(CNV_median_smoothed{chr})
		window_start   = max(data-1, 1);
		window_end     = min(data+1, length(CNV_median_smoothed{chr}));
		if (window_start == 1)
			window(1) = CNV_median_smoothed{chr}(data  );
			window(2) = CNV_median_smoothed{chr}(data+1);
		elseif (window_end == length(CNV_median_smoothed{chr}))
			window(1) = CNV_median_smoothed{chr}(data-1);
			window(2) = CNV_median_smoothed{chr}(data  );
		else
			window(1) = CNV_median_smoothed{chr}(data-1);
			window(2) = CNV_median_smoothed{chr}(data+1);
		end;
		CNV_differentiated{chr}(data) = (window(2)-window(1))/2;
	end;
end;

%%====================================================================================================
% [gaussian] smoothing filter
%     Original ChARM implementation did not smooth the differentiated signal before peak finding, but
%     doing so reduces the number of spurious peaks found by local peak-finder algorithms.
%-----------------------------------------------------------------------------------------------------
fprintf('\nSmoothing Filter 2');
window_width = smooth_window_width;
window_halfwidth = (window_width-1)/2;
for chr = 1:num_chr
	fprintf(['\n\t' num2str(chr) ':' num2str(num_chr) ':' num2str(length(CNV_differentiated{chr})) ]);
	CNV_differentiated_smoothed{chr} = smooth_gaussian(CNV_differentiated{chr},smooth_gaussian_sigma,smooth_gaussian_sigma*16);
end;

%% ===================================================================================================
% Find local maxima/minima.
%.....................................................................................................
% If peaks are smaller than 0.01, delete them as not reliable.
%-----------------------------------------------------------------------------------------------------
fprintf('\nFinding local maxima & minima\n');
fprintf('\tPeak size minimum cutoff = 0.001\n');
for chr = 1:num_chr
	[pks1{chr},locs1{chr}]       = findpeaks( CNV_differentiated_smoothed{chr});
	[pks2{chr},locs2{chr}]       = findpeaks(-CNV_differentiated_smoothed{chr});
	locs1{chr}(pks1{chr} < 0.01) = [];
	locs2{chr}(pks2{chr} < 0.01) = [];
	pks1{chr}( pks1{chr} < 0.01) = [];
	pks2{chr}( pks2{chr} < 0.01) = [];
	[locs{chr},keys]             = sort([locs1{chr} locs2{chr}]);
	pks_temp                     = [pks1{chr} pks2{chr}];
	pks{chr}                     = pks_temp(keys);

	%% move first peak to start of chromosome.
	locs{chr}(1) = 1;

	%% move last peak to the end of the chromosome.
	locs{chr}(end) = length(CNV_differentiated_smoothed{chr});

	%% If there is only one edge, fix so end edges are present.
	if (length(locs{chr}) == 1)
		locs{chr}(1) = 1;
		locs{chr}(2) = length(CNV_differentiated_smoothed{chr});
	end;

	fprintf(['\nInternal peak positions on chr : "' num2str(chr) '"\n']);
	fprintf(['\t{' num2str(locs{chr}(1)) ', ']);
	num_edges = length(locs{chr});
	if (num_edges > 2)
		for edge = 2:(length(locs{chr})-1)
			fprintf(['[' num2str(locs{chr}(edge)) ':' num2str(pks{chr}(edge)) '], ']);
			if (mod(edge,30) == 0);
				fprintf('\n\t');
			end;
		end;
	end;
	fprintf([num2str(locs{chr}(end)) '}\n']);
end;

%=====================================================================================================
% Generate figures representing stages of algorithm.
%-----------------------------------------------------------------------------------------------------
if (temp_figures == true)
	%% Assigning raw and filtered data to convenient names for later use.
	data1 = examined_data;
	for chr = 1:num_chr
		data2{chr} = CNV_median{chr};
		data3{chr} = CNV_median_smoothed{chr};
		data4{chr} = CNV_differentiated_smoothed{chr};
	end;
	maxY = 2;

	fig = figure(1);    dataShow = data1;
	set(gcf, 'Position', [0 70 1024 600]);
	for chr = 1:num_chr
		left   = chr_posX(chr);    bottom = chr_posY(chr);
		width  = chr_width(chr);   height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		hold on;
		c_ = [0 0 0];
		fprintf(['chr' num2str(chr) ':' num2str(length(examined_data{chr})) '\n']);
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
	saveas(fig,[workingDir 'ChARM_test.1.eps'], 'epsc');
	saveas(fig,[workingDir 'ChARM_test.1.png'], 'png');
	delete(fig);

	fig = figure(2);    dataShow = data2;
	set(gcf, 'Position', [0 70 1024 600]);
	for chr = 1:num_chr
		left   = chr_posX(chr);    bottom = chr_posY(chr);
		width  = chr_width(chr);   height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		hold on;
		c_ = [0 0 0];
		fprintf(['chr' num2str(chr) ':' num2str(length(examined_data{chr})) '\n']);
		for i = 1:length(dataShow{chr});
			x_ = [i i i-1 i-1];
			if (dataShow{chr}(i) == 0);    CNVhistValue = 0;    else;    CNVhistValue = dataShow{chr}(i);    end;
			startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);   
			set(f,'linestyle','none');
		end;
		x2 = chr_size(chr)*chr_length_scale_multiplier;
		plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
		hold off;
		xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
		set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''}); 
	end;
	saveas(fig,[workingDir 'ChARM_test.2.eps'], 'epsc');
	saveas(fig,[workingDir 'ChARM_test.2.png'], 'png');
	delete(fig);

	fig = figure(3);    dataShow = data3;
	set(gcf, 'Position', [0 70 1024 600]);
	for chr = 1:num_chr
		left   = chr_posX(chr);    bottom = chr_posY(chr);
		width  = chr_width(chr);   height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		hold on;
		c_ = [0 0 0];
		fprintf(['chr' num2str(chr) ':' num2str(length(examined_data{chr})) '\n']);
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
	saveas(fig,[workingDir 'ChARM_test.3.eps'], 'epsc');
	saveas(fig,[workingDir 'ChARM_test.3.png'], 'png');
	delete(fig);

	fig = figure(4);
	dataShow = data4;
	set(gcf, 'Position', [0 70 1024 600]);
	for chr = 1:num_chr
		left   = chr_posX(chr);    bottom = chr_posY(chr);
		width  = chr_width(chr);   height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		hold on;
		c_ = [0 0 0];
		for i = 1:length(dataShow{chr});
			x_ = [i i i-1 i-1];
			CNVhistValue = dataShow{chr}(i);
			startY = 0;
			endY   = CNVhistValue*8;
			y_ = [startY endY endY startY];
			f = fill(x_,y_,c_);
			set(f,'linestyle','none');
		end;
		x2 = chr_size(chr)*chr_length_scale_multiplier;
		plot([0; x2], [0; 0],'color',[0 0 0]);  % 2n line.
		for edge = 1:length(locs1{chr})
			plot([locs1{chr}(edge) locs1{chr}(edge)], [-1 1],'color',[0 0 1]);
		end;
		for edge = 1:length(locs2{chr})
			plot([locs2{chr}(edge) locs2{chr}(edge)], [-1 1],'color',[1 0 0]);
		end;
		hold off;
		xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
		ylim([-1 1]);
		set(gca,'YTick',[-1 0 1]);
		set(gca,'YTickLabel',{'','',''});
	end;
	saveas(fig,[workingDir 'ChARM_test.4.eps'], 'epsc');
	saveas(fig,[workingDir 'ChARM_test.4.png'], 'png');
	delete(fig)
end;


%% ###################################################################################################
%=====================================================================================================
% Expectation-Maximization Edge-Placement Algorithm
%-----------------------------------------------------------------------------------------------------
%% Initialize initial (t=0) posterior probabilities that a data point is in the left vs. right distributions adjacent to each edge.
%  Posterior probabilities are likelihood of membership in left vs. right distributions for each edge.
for chr = 1:num_chr
	position  = locs{chr};       % locations of edges for this chromosome.
	num_edges = length(position);
	data      = examined_data{chr};
	for edge = 1
		pos                = position(edge);
		R_windowSize(edge) = min(ceil(percent_window_size*(position(edge+1)-position(edge  ))),max_ROI);
		pP_dist_is_L{edge} = zeros(1,length(data));
		pP_dist_is_R{edge} = zeros(1,length(data));
		pP_dist_is_R{edge}((pos):(pos+R_windowSize(edge))) = 1;
	end;
	if (num_edges > 2)
		for edge = 2:(num_edges-1)
			pos                = position(edge);
			L_windowSize(edge) = min(ceil(percent_window_size*(position(edge  )-position(edge-1))),max_ROI);
			R_windowSize(edge) = min(ceil(percent_window_size*(position(edge+1)-position(edge  ))),max_ROI);
			pP_dist_is_L{edge} = zeros(1,length(data));
			pP_dist_is_R{edge} = zeros(1,length(data));
			pP_dist_is_L{edge}((pos-L_windowSize(edge)):(pos)) = 1;
			pP_dist_is_R{edge}((pos):(pos+R_windowSize(edge))) = 1;
		end;
	end;
	for edge = num_edges
		pos                = position(edge);
		L_windowSize(edge) = min(ceil(percent_window_size*(position(edge  )-position(edge-1))),max_ROI);
		pP_dist_is_L{edge} = zeros(1,length(data));
		pP_dist_is_R{edge} = zeros(1,length(data));
		pP_dist_is_L{edge}((pos-L_windowSize(edge)):(pos)) = 1;
	end;
	old_Ppost_dist_is_L{chr} = pP_dist_is_L;
	old_Ppost_dist_is_R{chr} = pP_dist_is_R;
end;

%% This iterated section of the ChARM algoritm adjusts the local positions of the edges.
for t = 1:1; % num_permutations
	%% ###################################################################################################
	%=====================================================================================================
	% Update membership (E-step).
	%-----------------------------------------------------------------------------------------------------
	%% Calculate Conditional probabilities that a data point is in the left vs. right distributions adjacent to each edge.
	for chr = 1:num_chr
		position  = locs{chr};       % locations of edges for this chromosome.
		num_edges = length(position);
		if (num_edges > 2)
			data      = examined_data{chr};   % data to be examined for this chromosome.
			for edge = 1
				% Left edge of chromosome.
				position(edge) = 1;
			end;
			for edge = 2:(num_edges-1);
				pos                  = position(edge);
				% max_ROI             = 20;
				% percent_window_size = 0.5;
				% The radius of influence is determined as follows: min(percent_window_size*window_size,max_ROI).
				L_windowSize(edge)   = min(ceil(percent_window_size*(position(edge  )-position(edge-1))),max_ROI);
				R_windowSize(edge)   = min(ceil(percent_window_size*(position(edge+1)-position(edge  ))),max_ROI);
				L_dist{edge}         = data((pos-L_windowSize(edge)):(pos));
				R_dist{edge}         = data((pos):(pos+R_windowSize(edge)));
				% Anything outside these ranges has a zero chance of being on either side of the edge.
				% This means that though the conditional probability is calculated for all positions relative to each edge,
				%    only a small area around each edge has to be stored.
				L_distMean(edge)     = mean(L_dist{edge});
				R_distMean(edge)     = mean(R_dist{edge});
				L_distStdev(edge)    = std(L_dist{edge});
				R_distStdev(edge)    = std(R_dist{edge});

				cP_dist_is_L{edge}   = zeros(1,length(data));
				cP_dist_is_R{edge}   = zeros(1,length(data));
				for loc = (pos-L_windowSize(edge)):(pos)
					cP_dist_is_L{edge}(loc) = normpdf(L_dist{edge}(loc-(pos-L_windowSize(edge))+1),L_distMean(edge),L_distStdev(edge));
					cP_dist_is_R{edge}(loc) = normpdf(L_dist{edge}(loc-(pos-L_windowSize(edge))+1),R_distMean(edge),R_distStdev(edge));
				end;
				for loc = (pos):(pos+R_windowSize(edge))
					cP_dist_is_L{edge}(loc) = normpdf(R_dist{edge}(loc-(pos)+1)  ,L_distMean(edge),L_distStdev(edge));
					cP_dist_is_R{edge}(loc) = normpdf(R_dist{edge}(loc-(pos)+1)  ,R_distMean(edge),R_distStdev(edge));
				end;
			end;
			for edge = length(position)
				% Right edge of chromosome.
				position(edge) = length(data);
			end;
			Pcond_dist_is_L{chr} = cP_dist_is_L;
			Pcond_dist_is_R{chr} = cP_dist_is_R;
		end;
		locs{chr} = position;
	end;

	%% Calculate P(theta_(j,k)^(t-1)) terms used in calculating posterior probabilities.
	for chr = 1:num_chr
		position         = locs{chr};       % locations of edges for this chromosome.
		num_edges        = length(position);
		if (num_edges > 2)
			data             = examined_data{chr};
			chr_bins         = length(data);
			old_pP_dist_is_L = old_Ppost_dist_is_L{chr};
			old_pP_dist_is_R = old_Ppost_dist_is_R{chr};
			for edge = 1:num_edges
				L_term(edge) = sum(old_pP_dist_is_L{edge})/chr_bins;
				R_term(edge) = sum(old_pP_dist_is_R{edge})/chr_bins;
			end;
			L_terms{chr} = L_term;
			R_terms{chr} = R_term;
		end;
	end;

	%% Calculate posterior probabilities that a data point is in the left vs. right distributions adjacent to each edge.
	%  Posterior probabilities are likelihood of membership in left vs. right distributions for each edge.
	for chr = 1:num_chr
		position     = locs{chr};       % locations of edges for this chromosome.
		num_edges    = length(position);
		if (num_edges > 2)
			L_term       = L_terms{chr};
			R_term       = R_terms{chr};
			cP_dist_is_L = Pcond_dist_is_L{chr};
			cP_dist_is_R = Pcond_dist_is_R{chr};
			data         = examined_data{chr};
			for edge = 2:(num_edges-1);
				L_windowSize(edge)   = min(ceil(percent_window_size*(position(edge  )-position(edge-1))),max_ROI);
				R_windowSize(edge)   = min(ceil(percent_window_size*(position(edge+1)-position(edge  ))),max_ROI);
				pP_dist_is_L{edge} = cP_dist_is_L{edge}.*L_term(edge)/(sum(cP_dist_is_L{edge})*L_term(edge) + sum(cP_dist_is_R{edge})*R_term(edge));
				pP_dist_is_R{edge} = cP_dist_is_R{edge}.*R_term(edge)/(sum(cP_dist_is_L{edge})*L_term(edge) + sum(cP_dist_is_R{edge})*R_term(edge));
			end;

			Ppost_dist_is_L{chr} = pP_dist_is_L;
			Ppost_dist_is_R{chr} = pP_dist_is_R;
		end;
	end;

	%% ###################################################################################################
	%=====================================================================================================
	% Mean and variance computation (M-step 1).
	%-----------------------------------------------------------------------------------------------------
	for chr = 1:num_chr
		position     = locs{chr};       % locations of edges for this chromosome.
		num_edges    = length(position);
		if (num_edges > 2)
			data         = examined_data{chr};   % data to be examined for this chromosome.
			pP_dist_is_L = Ppost_dist_is_L{chr};
			pP_dist_is_R = Ppost_dist_is_R{chr};
			for edge = 2:(num_edges-1)
				L_mean(edge)     = sum(pP_dist_is_L{edge}.*data)/sum(pP_dist_is_L{edge});
				R_mean(edge)     = sum(pP_dist_is_R{edge}.*data)/sum(pP_dist_is_R{edge});
				L_sigma_sq(edge) = sum((data-L_mean(edge)).^2.*pP_dist_is_L{edge})/sum(pP_dist_is_L{edge});
				R_sigma_sq(edge) = sum((data-R_mean(edge)).^2.*pP_dist_is_R{edge})/sum(pP_dist_is_R{edge});
			end;

			L_means{chr}     = L_mean;
			R_means{chr}     = R_mean;
			L_sigmas_sq{chr} = L_sigma_sq;
			R_sigmas_sq{chr} = R_sigma_sq;
		end;
	end;

	%% ###################################################################################################
	%=====================================================================================================
	% Edge adjustment (M-step 2).
	%-----------------------------------------------------------------------------------------------------
	pos_change   = [];
	fprintf(['\nIteration : ' num2str(t) ]);
	new_locs = locs;
	for chr = 1:num_chr
		position     = locs{chr};       % locations of edges for this chromosome.
		num_edges    = length(position);
		if (num_edges > 2)
			chr_bins     = length(data);
			pP_dist_is_L = Ppost_dist_is_L{chr};
			pP_dist_is_R = Ppost_dist_is_R{chr};
			L_mean       = L_means{chr};
			R_mean       = R_means{chr};
			L_sigma_sq   = L_sigmas_sq{chr};
			R_sigma_sq   = R_sigmas_sq{chr};
			new_position = [];
			for edge = 2:(num_edges-1)
				old_pos   = position(edge);
				pos_start = position(edge) - min(ceil(percent_window_size*(position(edge  )-position(edge-1))),max_ROI);
				pos_end   = position(edge) + min(ceil(percent_window_size*(position(edge+1)-position(edge  ))),max_ROI);
				new_edge_pos = [];
				count = 0;
				for pos = pos_start:(pos_end-1)
					count = count+1;
					new_edge_pos(count) = -(sum(log10(pP_dist_is_L{edge}(pos_start:pos))) + sum(log10(pP_dist_is_R{edge}((pos+1):pos_end))));
				end;
				[minVal,minIndex]  = min(new_edge_pos);
				new_pos            = pos_start+minIndex-1;
				new_position(edge) = new_pos;
			end;
			new_position(1)         = 1;
			new_position(num_edges) = chr_bins;
			pos_change{chr}         = new_position - position;
			new_locs{chr}           = new_position;

			fprintf(['\n\tchr' num2str(chr) ' : ' num2str(pos_change{chr}) ]);

			old_Ppost_dist_is_L{chr} = Ppost_dist_is_L{chr};
			old_Ppost_dist_is_R{chr} = Ppost_dist_is_R{chr};
		end;
	end;
%	fprintf(['\n\tave : ' num2str(mean(abs([pos_change{1} pos_change{2} pos_change{3}]))) ]);
    
	%% Update positions of edges.
	locs = new_locs;
end;

%% ###################################################################################################
%=====================================================================================================
% Window similarity test.
%    Delete edges when side window distrubutions are too similar.
%        ChARM used median of data in each side window, but I've found performing a Gaussian fit to
%        the data in each side window and using the derived median is more useful.
%    Convert pairs of edges that are too close into single edges.
%    Force left and right edges to left and right of datasets.
%-----------------------------------------------------------------------------------------------------
fprintf('\nRemoving edges that are bound by too similar windows.\n');
for chr = 1:num_chr
	test_edge = 2;
	test = locs{chr}
	num_starting_edges = length(locs{chr});
	for t = 1:num_starting_edges
		position        = locs{chr};
		position_sorted = position;
		num_edges       = length(position);
		if (num_edges > 2)
			data         = examined_data{chr};
			%----------------------------------------------------------------------
			% Force leftmost and rightmost edges to left and right of datasets.
			%    This was not described in the ChARM paper, but is useful.
			%----------------------------------------------------------------------
			position(1)   = 1;
			position(end) = length(examined_data{chr});   

			%----------------------------------------------------------------------
			% Converts edge pairs that are too close into single edges.
			%    This was not described in the ChARM paper, but is needed.   If two
			%    edges are in the same position or adjacent, the second edge is
			%    removed. Otherwise, the SNR calculation produces excessively high
			%    values.
			%----------------------------------------------------------------------
			remove_edges = zeros(1,num_edges);
			for edge = 2:num_edges
				if (position(edge-1) == position(edge)) || (position(edge-1)+1 == position(edge))
					remove_edges(edge) = 1;
				end;
			end;
			position(remove_edges == 1) = [];
			num_edges                   = length(position);

			%----------------------------------------------------------------------
			% Delete invalid edges by examining relationship of median and stdev values for left and right windows.
			%    This was not described in the ChARM paper, but dramatically improves
			%    the results of the implemented ChARM algorithms.
			%----------------------------------------------------------------------
			signal       = [];
			fprintf(['&&&& chr = ' num2str(chr) '; t = ' num2str(t) '; num_edges = ' num2str(num_edges) '.\n']);
			for edge = 2:(num_edges-1)
				pos_L           = position(edge-1);
				pos             = position(edge);
				pos_R           = position(edge+1);
				data_j1         = data(pos_L:pos);
				data_j2         = data(pos:pos_R);
				stdev_j1        = std(data_j1);
				stdev_j2        = std(data_j2);
				med_j1          = median(data_j1);
				med_j2          = median(data_j2);
				delta_median    = abs(med_j2-med_j1);
				signal(edge)    = (stdev_j1+stdev_j2)/2/delta_median;
			end;
			% Define signal values for the end edges.
			signal(1)            = 0;
			signal(num_edges)    = 0;
			% Delete edges where the average of side standard deviations, divided by the difference in side medians, is greater than 1.
			position(signal > 1) = [];
			num_edges            = length(position);

			%----------------------------------------------------------------------
			% Calculate SNRs and medians of left and right data windows.
			%----------------------------------------------------------------------
			SNR          = [];
			signal       = [];
			delta_med    = [];
			for edge = 2:(num_edges-1)
				pos_L           = position(edge-1);
				pos             = position(edge);
				pos_R           = position(edge+1);
				data_j1         = data(pos_L:pos);
				data_j2         = data(pos:pos_R);
				length_j1       = length(pos_L:pos);
				length_j2       = length(pos:pos_R);
				stdev_j1        = std(data_j1);
				stdev_j2        = std(data_j2);
				noise_j1        = sum(abs(data_j1-med_j1));
				noise_j2        = sum(abs(data_j2-med_j2));
				med_j1          = median(data_j1);
				med_j2          = median(data_j2);
				delta_med(edge) = abs(med_j2-med_j1);
				SNR(edge)       = delta_med(edge)*(length_j1+length_j2)/(noise_j1+noise_j2);
			end;
			% define the SNR and delta_med for the end edges.
			SNR(1)               = 1000;
			SNR(num_edges)       = 1000;
			delta_med(1)         = 1000;
			delta_med(num_edges) = 1000;
			if (t == 1) || (t == num_starting_edges)
				fprintf(['%%%% SNR       = [' num2str(SNR      ) ']\n\n']);
				fprintf(['&&&& delta_med = [' num2str(delta_med) ']\n\n']);
			end;

			%----------------------------------------------------------------------
			% Sort edges by difference of side window medians.
			%----------------------------------------------------------------------
			[delta_med_sorted,delta_med_index] = sort(delta_med(2:(num_edges-1)));
			fixed_index                        = [1 (delta_med_index+1) num_edges];
			[position_sorted]                  = position(fixed_index);
			[SNR_sorted]                       = SNR(fixed_index);
		end;

		%----------------------------------------------------------------------
		% Examine the first edge and delete it if it is bound by windows which are too similar.
		%----------------------------------------------------------------------
		SNR_threshold   = 3;
		delta_threshold = 0.25;
		if (num_edges > test_edge)
			edge = test_edge;
		%	if (SNR_sorted(edge) < SNR_threshold)
		%		position_sorted(edge) = [];
		%	else
		%		test_edge = test_edge + 1;
		%	end;
			if (delta_med(edge) < delta_threshold)
				position_sorted(edge) = [];
			else
				test_edge = test_edge + 1;
			end;
		else
			% Chromosome only has two edges (lenft and right ends), so none need to be examined or removed.
		end;

		%----------------------------------------------------------------------
		% Resort and save edge positions.
		%----------------------------------------------------------------------
		locs{chr} = sort(position_sorted);

%		%=====================================================================================================
%		% Generate figures representing intermediate and final output of algorithm.
%		%-----------------------------------------------------------------------------------------------------
%		if (temp_figures == true)
%			fig = figure(1);    dataShow = data1;
%			set(gcf, 'Position', [0 70 1024 600]);
%			for chr_fig = 1:num_chr
%				left   = chr_posX(chr_fig);    bottom = chr_posY(chr_fig);
%				width  = chr_width(chr_fig);   height = chr_height(chr_fig);
%				subplot('Position',[left bottom width height]);
%				hold on;
%				c_ = [0 0 0];
%				for i = 1:length(dataShow{chr_fig});
%					x_ = [i i i-1 i-1];
%					if (dataShow{chr_fig}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr_fig}(i);    end;
%					startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);
%					set(f,'linestyle','none');
%				end;
%				x2 = chr_size(chr_fig)*chr_length_scale_multiplier;
%				plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
%				for edge = 1:length(locs{chr_fig})
%					plot([locs{chr_fig}(edge) locs{chr_fig}(edge)], [0 maxY],'color',[0 0 1]);
%				end;
%				hold off;
%				xlim([0,chr_size(chr_fig)*chr_length_scale_multiplier]);
%				ylim([0,maxY]);
%				set(gca,'YTick',[0 maxY/2 maxY]);
%				set(gca,'YTickLabel',{'','',''});
%			end;
%			saveas(fig,[workingDir 'ChARM_test.5.' num2str(chr_fig) '_' num2str(t,'%02d') '.eps'], 'epsc');
%			delete(fig)
%		end;
	end;
end;

%=====================================================================================================
% Generate figures representing intermediate and final output of algorithm.
%-----------------------------------------------------------------------------------------------------
if (temp_figures == true)
	fig = figure(1);    dataShow = data1;
	set(gcf, 'Position', [0 70 1024 600]);
	for chr = 1:num_chr
		left   = chr_posX(chr);    bottom = chr_posY(chr);
		width  = chr_width(chr);   height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		hold on;
		c_ = [0 0 0];
		fprintf(['chr' num2str(chr) ':' num2str(length(examined_data{chr})) '\n']);
		for i = 1:length(dataShow{chr});
			x_ = [i i i-1 i-1];
			if (dataShow{chr}(i) == 0);    CNVhistValue = 1;    else;    CNVhistValue = dataShow{chr}(i);    end;
			startY = maxY/2;    endY = CNVhistValue;    y_ = [startY endY endY startY];    f = fill(x_,y_,c_);
			set(f,'linestyle','none');   
		end;
		x2 = chr_size(chr)*chr_length_scale_multiplier;
		plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
		for edge = 1:length(locs{chr})
			plot([locs{chr}(edge) locs{chr}(edge)], [0 maxY],'color',[0 0 1]);
		end;
		hold off;
		xlim([0,chr_size(chr)*chr_length_scale_multiplier]);    ylim([0,maxY]);
		set(gca,'YTick',[0 maxY/2 maxY]);    set(gca,'YTickLabel',{'','',''});
	end;
	saveas(fig,[workingDir 'ChARM_test.5.eps'], 'epsc');
	saveas(fig,[workingDir 'ChARM_test.5.png'], 'png');
	delete(fig)
end;

%%=========================================================================
% Save common_ChARM file for project : 'segmental_aneuploidy'.
%--------------------------------------------------------------------------   
dataFile = [workingDir 'Common_ChARM.mat']
fprintf(['\nSaving common_ChARM file for "' projectName '" : ' dataFile '\n']);

i = 0;
segmental_aneuploidy = [];
for chr = 1:num_chr
	position  = locs{chr};
	num_edges = length(position);
	data      = examined_data{chr};
	chr_size  = length(data);
	for edge = 1:num_edges
		if (position(edge) == 1) || (position(edge) == chr_size)
			% nothing is added to file, as these edges are later assumed.
		else
			i = i+1;
			segmental_aneuploidy(i).chr     = chr;                       % chromosome being examined.
			segmental_aneuploidy(i).break   = position(edge)/chr_size;   % percent along chromosome of edge.
		end;
	end;
end;
save(dataFile, 'segmental_aneuploidy');

end
