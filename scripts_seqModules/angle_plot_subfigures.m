%% =========================================================================================
% Draw angleplots to left of main chromosome cartoons.
%-------------------------------------------------------------------------------------------

if (AnglePlot == true)
	width      = 0.075;
	height     = chr_height(chr);
	bottom     = chr_posY(chr);
	chr_length = chr_size(chr);
	for segment = 1:length(chrCopyNum{chr})
		fprintf(['^^^     segment#    = ' num2str(segment) ':' num2str(length(chrCopyNum{chr})) '\n']);

		if (segment == 1) % generate sublot for each segment.
			subplot('Position',[0.03 bottom width (height/length(chrCopyNum{chr}))]);
		else
			subplot('Position',[0.03 (bottom+height/length(chrCopyNum{chr})*(segment-1)) width (height/length(chrCopyNum{chr}))]);
		end;

		peaks                     = chrSegment_peaks{              chr}{segment};
		mostLikelyGaussians       = chrSegment_mostLikelyGaussians{chr}{segment};
		actual_cutoffs            = chrSegment_actual_cutoffs{     chr}{segment};
		segment_smoothedHistogram = chrSegment_smoothed{           chr}{segment};
		segment_copyNum           = round(chrCopyNum{              chr}(segment));
		segment_chrBreaks         = chr_breaks{                    chr}(segment);

		fprintf(['^^^     copyNum             = ' num2str(segment_copyNum)     '\n']);
		fprintf(['^^^     peaks               = ' num2str(peaks)               '\n']);
		fprintf(['^^^     mostLikelyGaussians = ' num2str(mostLikelyGaussians) '\n']);
		fprintf(['^^^     actual_cutoffs      = ' num2str(actual_cutoffs)      '\n']);

		copynum = round(chrCopyNum{chr}(segment));
		region_ = 0;
		hold on;
		for region = mostLikelyGaussians
			region_ = region_+1;
			% Define color of the histogram region.
			if (FillColors == true)
				fprintf(['region_ #= ' num2str(region_) '\n']);
				if (show_uncalibrated == true)
					color = colorAB;
				else
					fprintf(['    copyNum      = ' num2str(copynum) '\n']);
					if (copynum == 0) %deletion or error
					elseif (copynum == 1) %monosomy
						if (apply_phasing == true)
							if (region == 2);     color = colorB;
							else                  color = colorA;
							end;
						else
							if (region == 2);     color = unphased_color_1of1;
							else                  color = unphased_color_1of1;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',[0 200]);
							set(gca,'XTickLabel',{'a','b'});
						end;
					elseif (copynum == 2) %disomy
						if (apply_phasing == true)
							if (region == 3);     color = colorBB;
							elseif (region == 2); color = colorAB;
							else                  color = colorAA;
							end;
						else
							if (region == 3);     color = unphased_color_2of2;
							elseif (region == 2); color = unphased_color_1of2;
							else                  color = unphased_color_2of2;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:100:200);
							set(gca,'XTickLabel',{'a','ab','b'});
						end;
					elseif (copynum == 3) %trisomy
						if (apply_phasing == true)
							if (region == 4);     color = colorBBB;
							elseif (region == 3); color = colorABB;
							elseif (region == 2); color = colorAAB;
							else                  color = colorAAA;
							end;
						else
							if (region == 4);     color = unphased_color_3of3;
							elseif (region == 3); color = unphased_color_2of3;
							elseif (region == 2); color = unphased_color_2of3;
							else                  color = unphased_color_3of3;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',[0 66.667 133.333 200]);
							set(gca,'XTickLabel',{'a','aab','abb','b'});
						end;
					elseif (copynum == 4) %tetrasomy
						if (apply_phasing == true)
							if (region == 5);     color = colorBBBB;
							elseif (region == 4); color = colorABBB;
							elseif (region == 3); color = colorAABB;
							elseif (region == 2); color = colorAAAB;
							else                  color = colorAAAA;
							end;
						else
							if (region == 5);     color = unphased_color_4of4;
							elseif (region == 4); color = unphased_color_3of4;
							elseif (region == 3); color = unphased_color_2of4;
							elseif (region == 2); color = unphased_color_3of4;
							else                  color = unphased_color_4of4;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:50:200);
							set(gca,'XTickLabel',{'a', '3:1', '2:2', '1:3' 'b'});
						end;
					elseif (copynum == 5) %pentasomy
						if (apply_phasing == true)
							if (region == 6);     color = colorBBBBB;
							elseif (region == 5); color = colorABBBB;
							elseif (region == 4); color = colorAABBB;
							elseif (region == 3); color = colorAAABB;
							elseif (region == 2); color = colorAAAAB;
							else                  color = colorAAAAA;
							end;
						else
							if (region == 6);     color = unphased_color_5of5;
							elseif (region == 5); color = unphased_color_4of5;
							elseif (region == 4); color = unphased_color_3of5;
							elseif (region == 3); color = unphased_color_3of5;
							elseif (region == 2); color = unphased_color_4of5;
							else                  color = unphased_color_5of5;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:40:200);
							set(gca,'XTickLabel',{'a', '4:!', '3:2', '2:3', '1:4' 'b'});
						end;
					elseif (copynum == 6) %hexasomy
						if (apply_phasing == true)
							if (region == 7);     color = colorBBBBBB;
							elseif (region == 6); color = colorABBBBB;
							elseif (region == 5); color = colorAABBBB;
							elseif (region == 4); color = colorAAABBB;
							elseif (region == 3); color = colorAAAABB;
							elseif (region == 2); color = colorAAAAAB;
							else                  color = colorAAAAAA;
							end;
						else
							if (region == 7);     color = unphased_color_6of6;
							elseif (region == 6); color = unphased_color_5of6;
							elseif (region == 5); color = unphased_color_4of6;
							elseif (region == 4); color = unphased_color_3of6;
							elseif (region == 3); color = unphased_color_4of6;
							elseif (region == 2); color = unphased_color_5of6;
							else                  color = unphased_color_6of6;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:33.333:200);
							set(gca,'XTickLabel',{'a', '5:1', '4:2', '3:3', '2:4', '1:5' 'b'});
						end;
					elseif (copynum == 7) %heptasomy
						if (apply_phasing == true)
							if (region == 8);     color = colorBBBBBBB;
							elseif (region == 7); color = colorABBBBBB;
							elseif (region == 6); color = colorAABBBBB;
							elseif (region == 5); color = colorAAABBBB;
							elseif (region == 4); color = colorAAAABBB;
							elseif (region == 3); color = colorAAAAABB;
							elseif (region == 2); color = colorAAAAAAB;
							else                  color = colorAAAAAAA;
							end;
						else
							if (region == 8);     color = unphased_color_7of7;
							elseif (region == 7); color = unphased_color_6of7;
							elseif (region == 6); color = unphased_color_5of7;
							elseif (region == 5); color = unphased_color_4of7;
							elseif (region == 4); color = unphased_color_4of7;
							elseif (region == 3); color = unphased_color_5of7;
							elseif (region == 2); color = unphased_color_6of7;
							else                  color = unphased_color_7of7;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:28.571:200);
							set(gca,'XTickLabel',{'a', '', '5:2', '', '', '2:5', '' 'b'});
						end;
					elseif (copynum == 8)  %octasomy
						if (apply_phasing == true)
							if (region == 9);     color = colorBBBBBBBB;
							elseif (region == 8); color = colorABBBBBBB;
							elseif (region == 7); color = colorAABBBBBB;
							elseif (region == 6); color = colorAAABBBBB;
							elseif (region == 5); color = colorAAAABBBB;
							elseif (region == 4); color = colorAAAAABBB;
							elseif (region == 3); color = colorAAAAAABB;
							elseif (region == 2); color = colorAAAAAAAB;
							else                  color = colorAAAAAAAA;
							end;
						else
							if (region == 9);     color = unphased_color_8of8;
							elseif (region == 8); color = unphased_color_7of8;
							elseif (region == 7); color = unphased_color_6of8;
							elseif (region == 6); color = unphased_color_5of8;
							elseif (region == 5); color = unphased_color_4of8;
							elseif (region == 4); color = unphased_color_5of8;
							elseif (region == 3); color = unphased_color_6of8;
							elseif (region == 2); color = unphased_color_7of8;
							else                  color = unphased_color_8of8;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:22.222:200);
							set(gca,'XTickLabel',{'a', '', '6:2', '', '4:4', '', '2:6', '' 'b'});
						end;
					else % (copynum >= 9) %nonasomy
						if (apply_phasing == true)
							if (region == 10);    color = colorBBBBBBBBB;
							elseif (region == 9); color = colorABBBBBBBB;
							elseif (region == 8); color = colorAABBBBBBB;
							elseif (region == 7); color = colorAAABBBBBB;
							elseif (region == 6); color = colorAAAABBBBB;
							elseif (region == 5); color = colorAAAAABBBB;
							elseif (region == 4); color = colorAAAAAABBB;
							elseif (region == 3); color = colorAAAAAAABB;
							elseif (region == 2); color = colorAAAAAAAAB;
							else                  color = colorAAAAAAAAA;
							end;
						else
							if (region == 10);    color = unphased_color_9of9;
							elseif (region == 9); color = unphased_color_8of9;
							elseif (region == 8); color = unphased_color_7of9;
							elseif (region == 7); color = unphased_color_6of9;
							elseif (region == 6); color = unphased_color_5of9;
							elseif (region == 5); color = unphased_color_5of9;
							elseif (region == 4); color = unphased_color_6of9;
							elseif (region == 3); color = unphased_color_7of9;
							elseif (region == 2); color = unphased_color_8of9;
							else                  color = unphased_color_9of9;
							end;
						end;
						if (segment == 1)
							set(gca,'XTick',0:20:200);
							set(gca,'XTickLabel',{'a', '', '', '6:3', '', '', '3:6', '', '', 'b'});
						end;
					end;
				end;
			else
				color = colorAB;
			end;

			fprintf(['    mostLikelyGaussian   = ' num2str(region) '\n']);
			if (length(mostLikelyGaussians) <= 1)
				% draw entire smoothed histogram.
				area(1:200,segment_smoothedHistogram(1:200),'FaceColor',color,'EdgeColor',color);
			else
				% draw segment of smoothed histogram corresponding to region.
				if (region_ == 1) % first region in list.
					coord1 = max(min(round(actual_cutoffs(region_))+1, 200), 1);
					area(1:coord1, segment_smoothedHistogram(1:coord1), 'FaceColor',color,'EdgeColor',color);
					fprintf(['    angleplotCoordinates = 1:' num2str(coord1) '\n']);
				elseif (region_ == length(mostLikelyGaussians)) % last region in list.
					coord2 = max(min(round(actual_cutoffs(region_-1))+1, 200), 1);
					area(coord2:200, segment_smoothedHistogram(coord2:200), 'FaceColor',color,'EdgeColor',color);
					fprintf([' angleplotCoordinate = ' num2str(coord2) ':200\n']);
				else
					coord3 = max(min(round(actual_cutoffs(region_-1))+1, 200), 1);
					coord4 = max(min(round(actual_cutoffs(region_  ))+1, 200), 1);
					area(coord3:coord4, segment_smoothedHistogram(coord3:coord4), 'FaceColor',color,'EdgeColor',color);
					fprintf(['    angleplotCoordinates = ' num2str(coord3) ':' num2str(coord4) '\n']);
				end;
			end;
			fprintf(['    color = ' num2str(color) '   (colorA = [' num2str(colorA) ']; colorB = [' num2str(colorB) '])\n']);
		end;

		colorPeak   = [0.5 0.5 0.5]; % color of lines drawn at peak locations.
		colorCutoff = [1.0 0.0 0.0]; % color of lines drawn at cutoffs between Gaussian fits.
		for peak = 1:length(peaks)
			plot([peaks(peak); peaks(peak)],[0; 1],'color',colorPeak);
		end;
		for cutoff = 1:length(actual_cutoffs)
			plot([actual_cutoffs(cutoff); actual_cutoffs(cutoff)],[0; 1],'color',colorCutoff);
		end;
		set(gca,'FontSize',10);
		hold off;
		set(gca,'YTick',[]);
		if (segment ~= 1)
			set(gca,'XTick',[]);
		end;
		xlim([0,200]);
		ylim([0,1]);
	end;
end;
