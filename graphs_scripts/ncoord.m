function [ ] = RADSpaceGraphs()

	width = 4; 	% Width in inches
	height = 9;	% Height in inches
	alw = 1;	% AxesLineWidth
	afsz = 5.5;  	% AxesFontsize
	fsz = 5.5;  	% Fontsize
	lw = 0.7;  	% LineWidth
	msz = 0.8;   	% MarkerSize
	grey = 0.9; %grey lightness

	figure();

	t = 1;

        tempstr = {'acetic acid' 'acetone' 'acetonitrile' 'ammonia' 'aniline' 'benzene' 'benzyl alcohol' 'benzaldehyde' 'butane' 'butanol' '2-butoxyethanol' 'carbon dioxide' 'chloroform' 'cyclohexane' 'diazene' 'dichloromethane' 'diethanolamine' 'diethyl ether' 'DMFA' 'DMSO' '1,4-dioxane' 'ethane' 'ethanol' 'ethene' 'ethyl acetate' 'ethylamine' 'ethylene glycol' 'formamide' 'formic acid' 'furan' 'hexane' 'hexanol' 'hydrazine' 'hydrogen peroxide' 'hydrogen sulfide' 'methane' 'methanethiol' 'methanol' 'methylamine' 'NMA' 'octanol' 'pentane' 'pentanol' 'piperidine' 'propane' 'propanol' 'pyridine' 'styrene' 'TBA' 'tetrahydrofuran' 'TFE' 'toluene' 'triethylamine' 'm-xylene' 'o-xylene' 'p-xylene' ;};

	while t <= length(tempstr)

		ncFilename = [char(tempstr(t)) '/ncoord.txt'];
		ncArray = load(ncFilename);

                splot = subplot(7, 8, t);
		splotpos = get(splot,'position');
		splotpos(4) = splotpos(4)*0.7; % Take 20 percent off height (4);
		splotpos(3) = splotpos(3)*1.0; % Take 20 percent off height (4);
                if t > 45
                        splotpos(2) = splotpos(2)+0.0001;
                end
	        set(splot, 'position', splotpos);

                xlim([0 15]);
		ylim([0 0.39])
		
		hold on
		set(gca,'ticklength',2*get(gca,'ticklength'))
                plot(ncArray(1:15,1), ncArray(1:15,2), 'o','Color', 'k', 'linewidth', 1.3*lw, 'MarkerSize', 1*msz)
		hold off

		text(0.05,1.1,tempstr(t),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', fsz)
		set(gca, 'FontSize', afsz, 'LineWidth', alw);
                if (t == 1) || (t == 9) || (t == 17) || (t == 25) || (t == 33) || (t == 41) || (t == 49)
		  ylabel('{\itp}({\itN}{_c_o_o_r_d)}');
  		  set(gca, 'ytick', [0.0 0.1 0.2 0.3 0.4]);
	        else
		    set(gca,'YTickLabel',[]);
                end
		if t > 0
  	 	  xlabel('{{\it N}_c_o_o_r_d} ');
		  set(gca, 'xtick', [0 5 10 15]);
		  tix=get(gca,'xtick')';
		  set(gca,'xticklabel',num2str(tix,'%d'))
                else
                  set(gca,'XTickLabel',[]);
                end

		t = t + 1;
	end

	outputFilename = ['offset_ncoord.eps'];
%	set(gcf, 'PaperPositionMode', 'auto');
	print('-depsc', outputFilename);
	print('-djpeg', 'ncoord.jpg');
        print('-dtiff', '-r355', 'ncoord.tif');

	close all
