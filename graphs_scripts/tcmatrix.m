function [ ] = RADSpaceGraphs()

	width = 9; 	% Width in inches
	height = 3.23;	% Height in inches
	alw = 1;	% AxesLineWidth
	afsz = 5.5;  	% AxesFontsize
	fsz = 5.5;  	% Fontsize
	lw = 3;  	% LineWidth
	msz = 1;   	% MarkerSize

	x= [0.285 0.285];
        y= [0.285 0.285];
	grey = 0.8; %grey lightness

	figure();

	t = 1;

        tempstr = {'acetic acid' 'acetone' 'acetonitrile' 'ammonia' 'aniline' 'benzene' 'benzyl alcohol' 'benzaldehyde' 'butane' 'butanol' '2-butoxyethanol' 'carbon dioxide' 'chloroform' 'cyclohexane' 'diazene' 'dichloromethane' 'diethanolamine' 'diethyl ether' 'DMFA' 'DMSO' '1,4-dioxane' 'ethane' 'ethanol' 'ethene' 'ethyl acetate' 'ethylamine' 'ethylene glycol' 'formamide' 'formic acid' 'furan' 'hexane' 'hexanol' 'hydrazine' 'hydrogen peroxide' 'hydrogen sulfide' 'methane' 'methanethiol' 'methanol' 'methylamine' 'NMA' 'octanol' 'pentane' 'pentanol' 'piperidine' 'propane' 'propanol' 'pyridine' 'styrene' 'TBA' 'tetrahydrofuran' 'TFE' 'toluene' 'triethylamine' 'm-xylene' 'o-xylene' 'p-xylene' ;};


	while t <= length(tempstr)

		Filename = [char(tempstr(t)) '/tcmatrix.txt'];
		TCMatrix = load(Filename);
                [m,n] = size(TCMatrix);
logm = log(m);
if (( t == 0 ) || (t == 0) )
        clims = [-1 1];
 else
        clims = [-300 300];
   end

                splot = subplot(7, 8, t);
		splotpos = get(splot,'position');
		splotpos(4) = splotpos(4)*1.25; % Take 20 percent off height (4);
		splotpos(3) = splotpos(3)*1.15; % Take 20 percent off width (3);
                if t > 7
                        splotpos(2) = splotpos(2)+0.001;
                end
	        set(splot, 'position', splotpos);

                hold on;
                %box on;
                imagesc(TCMatrix, clims);
                colormap(gray(50));
%colormap(hot);
                yshift = m - 30; % Shift plot up to top
		xshift = 0 - 0.5*(30-m);
                x2shift = 30+xshift;
                axis([xshift x2shift yshift m]);
%                axis([0 34 yshift m]);
axis off;

                %xlim([-4 40]);
		%ylim([-4 40])
		
		set(gca,'ticklength',2*get(gca,'ticklength'))
		hold off

		text(0.52,1.16,tempstr(t),'Units', 'Normalized', 'VerticalAlignment', 'Top',  'HorizontalAlignment', 'center', 'FontSize', fsz)
		set(gca, 'ytick', []);
                set(gca, 'xtick', []);
		set(gca, 'FontSize', afsz, 'LineWidth', alw);

		t = t + 1;
	end

	outputFilename = ['offset_tcmatrix.eps'];
%	set(gcf, 'PaperPositionMode', 'auto');
	print('-depsc', outputFilename);
	print('-djpeg', 'tcmatrix.jpg');
        print('-dtiff', '-r355', 'tcmatrix.tif');

	close all
