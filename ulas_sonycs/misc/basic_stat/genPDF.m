function [pdf,val] = genPDF(imSize,p0, params)
% [pdf,val] = ~(imSize,deg,p0 [params])
% property: sum(pdf(:))== numel(imSize)*p0
%
%	generates a pdf for a 1-3 dim. random sampling pattern
%	with polynomial variable density sampling
%
%	Input:
%		imSize - size of matrix or vector
%       p0 - partial sampling factor e.g. 0.5 for half
%		params.deg - power of polynomial (if 0 uniform distribution)
%		params.distType - 1 or 2 for L1 or L2 distance measure
%		params.radius - radius of fully sampled center
%		params.fig - display output
% 		params.box(=[-1,1]) - origin is center of box
%
%	Output:
%		pdf - the pdf
%		val - min sampling density
%
% 
%	Example:
%	[pdf,val] = genPDF([128,128],2,0.5,2,0,1);
%
%	(c) Michael Lustig 2007
%   G.Troll April 2014: -- extended to 3d, 
%                       -- automatic adaptation of degree
%                       -- change of signature (add. parameters moved to params)
%                       -- center of distribution made variable (via params.box)
%                       -- extendet to degree params.deg=0 (uniform distr)
%

if nargin < 3
    params=struct;
end
if ~isfield(params,'deg')
	params.deg = 5;
end
if ~isfield(params,'distType')
	params.distType = 2;
end
if ~isfield(params,'radius')
	params.radius = 0;
end
if ~isfield(params,'fig');
	params.fig= 0;
end
if ~isfield(params,'box');
    % origin is center of box
    params.box=[-1,1];
end
a=params.box(1);
b=params.box(2);
deg=params.deg;

imSize=reshape(imSize,1,[]);
if deg==0
    pdf = p0*ones(imSize);
    val=p0;
    return;
end

minval=0;
maxval=1;
val = 0.5;

L=length(imSize);
if L==1
	imSize = [imSize,1,1];
elseif L==2
    imSize = [imSize,1];
end

sx = imSize(1);
sy = imSize(2);
sz=  imSize(3);
PCTG = floor(p0*sx*sy*sz);


dims=3-sum(imSize==1);


if dims==3  % 3d
    [x,y,z] = meshgrid(linspace(a,b,sy),linspace(a,b,sx),linspace(a,b,sz));
	switch params.distType
		case 1
			r = max(abs(x),abs(y),abs(z));
		otherwise
			r = sqrt((x).^2+(y).^2+(z).^2);
			r = r/max(abs(r(:)));			
	end
elseif dims==2  % 2D
	[x,y] = meshgrid(linspace(a,b,sy),linspace(a,b,sx));
	switch params.distType
		case 1
			r = max(abs(x),abs(y));
		otherwise
			r = sqrt((x).^2+(y).^2);
			r = r/max(abs(r(:)));			
	end

else %1d
	r = abs(linspace(a,b,max(sx,sy)));
end

idx = find(r<params.radius);

ok=false;
it=0;
itmax=100;
while ~ok && it <itmax
    it=it+1;
    pdf = (1-r).^params.deg; pdf(idx) = 1;
    ok=floor(sum(pdf(:))) < PCTG;      
    if ~ok       
        params.deg=params.deg+1;
    end
end
deg=params.deg;
if ~ok
    error('infeasible without undersampling dc, increase deg');   
end

% begin bisection
N=0;
while N~=PCTG % optimal
	val = (minval+maxval)/2;
	pdf = (1-r).^deg + val; pdf(pdf>1) = 1; pdf(idx)=1;
	N = floor(sum(pdf(:)));
	if N > PCTG % infeasible
		maxval=val;
	end
	if N < PCTG % feasible, but not optimal
		minval=val;
    end	
end

if params.fig~=0
	prepfigure(params.fig);	
	if dims==2
        subplot(2,2,[1,2]); imshow(pdf); colorbar;
        title(['PDF with p_0= ',num2str(p0,'%3.2g')],'fontsize',12);
		subplot(223), plot(pdf(end/2+1,:)); 
        title('central cross section','fontsize',12);
        xlabel('x');
        subplot(224), plot(pdf(:,end/2+1)); 
        title('central cross section','fontsize',12);
        xlabel('y');
    elseif dims==1
		plot(pdf);  title(['PDF with p_0= ',num2str(p0,'%3.2g')],'fontsize',12);
    elseif dims==3
       s3=Signal3D(pdf);  s3.signalname=['PDF with p_0= ',num2str(p0,'%3.2g')];   
       
       subplot(221), s3.graph_signalScatter(false);
       subplot(222), plot(squeeze(pdf(:,floor(end/2)+1,floor(end/2)+1)));
       xlabel('y');
       title('central cross section','fontsize',12);
       
       subplot(223), plot(squeeze(pdf(floor(end/2)+1,:,floor(end/2)+1)));
       title('central cross section','fontsize',12);
       xlabel('x');
       
       subplot(224), plot(squeeze(pdf(floor(end/2)+1,floor(end/2)+1,:)));
       title('central cross section','fontsize',12);
       xlabel('z');
	end
end






