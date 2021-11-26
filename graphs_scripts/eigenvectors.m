function [ ] = RADSpaceGraphs()

	width = 3; 	% Width in inches
	height = 3;	% Height in inches
	alw = 1;	% AxesLineWidth
	afsz = 8;  	% AxesFontsize
	fsz = 8;  	% Fontsize
	lw = 0.6;  	% LineWidth
	msz = 1;   	% MarkerSize

	x= [0.145 0.385];
        y= [0.145 0.385];
	grey = 0.8; %grey lightness

	figure();

	t = 1;

        tempstr = {'benzene';};

	while t <= length(tempstr)

		Filename = [char(tempstr(t)) '/teigenvectors.txt'];
		FTCMatrix = load(Filename);
                [m,n] = size(FTCMatrix);
        logm = log(m);
        clims = [0 0.5];

                splot = subplot(2, 7, t);
		splotpos = get(splot,'position');
		splotpos(4) = splotpos(4)*0.4; % Take 20 percent off height (4);
		splotpos(3) = splotpos(3)*1.2; % Take 20 percent off width (3);
                if t > 7
                        splotpos(2) = splotpos(2)+0.34;
                end
	        set(splot, 'position', splotpos);

                hold on;
                %box on;
                imagesc(FTCMatrix, clims);
                c = gray(50); % Change order of gray
                c = flipud(c);
                colormap(c);
%colormap(hot);
                yshift = m - 34; % Shift plot up to top
		xshift = 0 - 0.5*(34-m);
                x2shift = 34+xshift;
                axis([xshift x2shift yshift m]);
%                axis([0 34 yshift m]);
axis off;

                %xlim([-3 34]);
		%ylim([-3 34])
		
		set(gca,'ticklength',2*get(gca,'ticklength'))
		hold off

		text(0.52,1.2,tempstr(t),'Units', 'Normalized', 'VerticalAlignment', 'Top',  'HorizontalAlignment', 'center', 'FontSize', fsz)
		set(gca, 'ytick', []);
                set(gca, 'xtick', []);
		set(gca, 'FontSize', afsz, 'LineWidth', alw);

		t = t + 1;
	end

	outputFilename = ['eigenvectors.eps'];
%	set(gcf, 'PaperPositionMode', 'auto');
	print('-depsc', outputFilename);
	print('-djpeg', 'eigenvectors.jpg');

	close all
