function AX = subplotspace(varargin)
%SUBPLOTSPACE Spread subplots out.
%   SUBPLOTSPACE increases horizontal spacing between subplots      
%   by 1%.							    
%   SUBPLOTSPACE(DIRECTION) increases spacing between subplots      
%   by 1% where DIRECTION is either 'vertical' or 'horizontal'.     
%   SUBPLOTSPACE(DIRECTION,N) increases spacing between subplots    
%   by N%.  If N is negative then the space is reduced.					%   
%

	
%
%   Note, this function is useful for fixing the problem of 
%   axes labels and titles running into each other in figures 
%   using multiple subplots.
%
%   See also SUBPLOT
%

%Written by Daniel P. Dougherty 4/30/2004
%    

DIR = 'h';
N = 1;

if (nargin == 1)
	DIR = varargin{1};
elseif  (nargin == 2)
	DIR = varargin{1};
	N = varargin{2};	
elseif  (nargin > 2)
	error('Too many inputs.');
end

DIR = lower(DIR(1));

AX = findobj(gcf,'Type','axes');

if (isempty(AX) | (length(AX) <= 1))
	return;
end

UNITS = get(AX,'Units');

set(AX,'Units','pixels');

AXES = get(AX,'Position');
AXES = cat(1,AXES{:});

switch DIR
case 'h'
	UNODE = unique(AXES(:,1));
	FIGPOS = get(gcf,'Position');


	spacer = cumsum([0;diff(UNODE)]);
	UNODEnew = UNODE + (0.01*N)*spacer;

	for i = 1:length(UNODE)
		ind = find(AXES(:,1) == UNODE(i));
	%	ind
		for j = 1:length(ind)
			axes(AX(ind(j)));
			this_ax = get(gca,'Position');
			this_ax(1) = UNODEnew(i);
			set(gca,'Position',this_ax);
		end
	end

	FIGPOS(3) = FIGPOS(3)+ (0.01*N)*max(spacer);
	set(gcf,'Position',FIGPOS);

	for i = 1:length(AX)
		set(AX,'Units',UNITS{i});
	end
	
case 'v'
	UNODE = unique(AXES(:,2));
	FIGPOS = get(gcf,'Position');
	
	spacer = cumsum([0;diff(UNODE)]);
	UNODEnew = UNODE + (0.01*N)*spacer;

	for i = 1:length(UNODE)
		ind = find(AXES(:,2) == UNODE(i));
		for j = 1:length(ind)
			axes(AX(ind(j)));
			this_ax = get(gca,'Position');
			this_ax(2) = UNODEnew(i);
			set(gca,'Position',this_ax);
		end
	end
        
	FIGPOS(4) = FIGPOS(4)+ (0.01*N)*max(spacer);
	
	set(gcf,'Position',FIGPOS);

	for i = 1:length(AX)
		set(AX,'Units',UNITS{i});
	end
	
otherwise
	error('Unrecognized direction.');
end




