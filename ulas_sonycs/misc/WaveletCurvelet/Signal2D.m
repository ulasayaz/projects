classdef Signal2D<SignalClass
    % 2D signals wrapping image data
    % GT Berlin, 2013
    %
    % example:
    % omega=[10.1,6.3],; hN=101; sig=0.01; del0=[0.4,0.13];w=10;a=-1;
    % s2=Signal2D.make_sinus(omega,del0,hN,sig);
    % s2=Signal2D.make_fromImage('cameraman.bmp');
    % s2=Signal2D.make_zernike(48);
    % -- show signal:
    % s2.graph_signal;
    % -- change units:
    % s2.set_scale(PosScales.n); s2.graph_signal;
    % w=20; sig=0;
    % s2=Signal2D.make_rect(w,hN,sig);
    % s2.graph_signal;
    % hN=128;
    % s2=Signal2D.make_star(10,hN,0);
    % Delta signal:
    % s2=Signal2D(); s2=s2.make_delta([10,10],256);
    %
    
    properties
        xt      %@<function_handle> contin. 2D time signal (default is interpolation of xn)
        Rs      %@<real> signal sampling rectangle
        axis_type
    end
    
    properties (SetAccess=protected)
        
    end
    
    %% constructor and commands
    methods
        function obj=Signal2D(xdata,hN,hR)
            obj.requireOnly(nargin<1 || ismatrix(xdata) ... %&& ~isvector(xdata)) ...
                || ( isa(xdata,'function_handle')...
                && nargin>=2 && obj.isnaturalArray(hN) ),'local',...
                'input is empty or matrix or function plus sampling size & range');
            
            if nargin==0
                % needed to allow empty calls from subclasses
                % and prevent their premature testing of invariant
                xdata=[];
                hN=[0,0];
                hR=[0,0];
            end
            if nargin <3
                hR=1;
            end
            if isscalar(hR)
                hR=hR*[1,1];
            end
            if isnumeric(xdata) || islogical(xdata)
                obj.xn=xdata;
                if ~isempty(xdata)
                    obj.xt=@(x,y) interp2(obj.xn,x,y);
                else
                    obj.xt=@(x,y) 0;
                end
                obj.set_scale(PosScales.n);
            else % isa(xdata,'function_handle')
                obj.xt=xdata;
                [X,Y]=obj.meshgrid(hN,hR,false);
                % funktioniert nicht mit obj.xt=@(t) randn:
                % obj.xn=obj.xt(x,y);
                % statt dessen:
                obj.xn=arrayfun(obj.xt,X,Y);
                obj.set_scale(PosScales.x);
            end
            
            obj.Rs=hR;
            obj.scale=PosScales.n;
            obj.units={'',''};
            obj.colormap_active='default';  % e.g. 'default'
            obj.colormap_freeze=true;
            obj.winfun=[];
            obj.padd=0;
            obj.signalname=[];
            obj.origin=[0,0];
            obj.fig=1;
            obj.ison_colorbar =true;
            obj.axis_type='image';
            % turn warning for reading images off:
            warning('off', 'MATLAB:imagesci:readpnm:extraData');
            
            obj.ensure(true,'local','verify invariant');
            
        end % constructor
        
        function set_Rs(obj,hR,units)
            % set sampling positions and units of signal xn
            obj.requireOnly(length(hR)==2,'local',...
                'hR gives rectangle sides');
            obj.require(ischar(units)||iscellstr(units),'local',...
                'units is single string or cell string');
            
            obj.Rs=hR;
            obj.set_units(units);
        end
        
        function set_units(obj,units)
            % e.g. set_scale(PosScales.d)
            obj.require(ischar(units)||iscellstr(units),'local', ...
                'units is single string or cell string');
            if ~iscellstr(units)
                units={units,units};
            end
            obj.units=units;
        end
        
        function set_signal(obj,xn,hR)
            % ~(xn,[hR]) set signal to matrix xn
            obj.requireOnly(ismatrix(xn) ,'local','input is matrix');
            if nargin <3 || isempty(hR)
                hR=size(xn);
            end
            obj.reset_Signal;
            obj.signalname=[];
            
            obj.xn=xn;
            obj.set_scale(PosScales.n);
            obj.padd=[0,0];
            obj.Rs=hR;
            obj.xt=@(x,y) interp2(obj.xn,x,y);
            
            obj.ensureOnly(length(size(obj.xn))==2 ,...
                'local','xn is a true matrix');
            obj.ensureOnly(isa(obj.xt,'function_handle'),...
                'local','xt is function');
        end
        
        function set_origin(obj,z)
            % ~(z) set origin to z => domain D=[z-Ts/2, z+Ts/2]
            obj.requireOnly(length(z)==2,'local',' z is pair of numbers');
            obj.origin=z;
        end
        
        function set_signalfun(obj,xt,N,hR)
            % ~(xt,N,[hR]) set a signal function and sample it at N points
            obj.requireOnly( isa(xt,'function_handle')...
                && obj.isnaturualArray(N),'local',...
                'input is function plus sampling size');
            if nargin <4
                hR=2*[1,N(2)/N(1)];
            end
            
            obj.reset_Signal;
            obj.signalname=[];
            obj.Rs=hR;
            obj.xt=xt;
            obj.set_scale(PosScales.x);
            
            [X,Y]=obj.meshgrid(hN,hR,false);
            % funktioniert nicht mit obj.xt=@(t) randn:
            % obj.xn=obj.xt(x,y);
            % statt dessen:
            obj.xn=arrayfun(obj.xt,X,Y);
            
            obj.padd=[0,0];
            
            obj.ensureOnly(size(obj.xn,2)<=1,'local','xn is a row vector');
            obj.ensureOnly(isa(obj.xt,'function_handle'),...
                'local','xt is function');
        end
        
        function imwrite(obj,fn,fmt)
            % ~(fn,ext) write obj.xn to file in image fromat fmt
            imwrite(obj.xn,fn,fmt);
        end
        
        
    end
    
    %% transformation commands
    methods
        
        function set_0padding(obj,padd)
            % ~(padd) symm. 0-padding (left/right, bottom/top) by padd;
            % calling with padd=[0,0] removes old padding;
            
            obj.requireOnly(isequal(floor(padd),padd),'local',...
                'padd is pair of non-negative integers');
            if isscalar(padd)
                padd=padd*[1,1];
            end
            obj.reset_Signal;
            
            local_inv=obj.N()-obj.padd;
            if max(obj.padd)>0
                %remove old zero-padding:
                % proportion of unpadded part
                prop=(obj.N-obj.padd)./obj.N;
                % decrease sampling time interval to unpadded value:
                obj.Rs=obj.Rs.*prop;
                %remove old 0-padding:
                obj.remove_padding();
            end
            
            % increase sampling time interval to padded value
            obj.Rs=obj.Rs.*(obj.N+2*padd)./obj.N;
            % add new padding symmetrically
            obj.paddDir='both';
            pval=0;
            obj.xn=padarray(obj.xn,padd,pval, obj.paddDir);
            obj.padd=2*padd;
            
            obj.ensure(obj.N-obj.padd==local_inv,'local',...
                'padding ok');
        end
        
        function remove_padding(obj)
            % remove padding from signal
            if ~all(obj.padd==0)
                if strcmp(obj.paddDir,'post')
                    obj.xn = obj.xn(1:end-obj.padd(1),1:end-obj.padd(2) );
                elseif strcmp(obj.paddDir,'both')
                    obj.xn = obj.xn(1+obj.padd(1)/2:end-obj.padd(1)/2,...
                        1+obj.padd(2)/2:end-obj.padd(2)/2);
                elseif strcmp(obj.paddDir,'pre')
                    obj.xn = obj.xn(1+obj.padd(1):end,1+obj.padd(2):end);
                end
                obj.padd=0;
            end
        end
        
        function coarsen(obj,f,method)
            % ~(f,[method]) reduce to resolution f<1 (instead to a certain size)
            % method is in ['nearest','bilinear','bicubic'];
            obj.requireOnly(isnumeric(f) && isscalar(f) && f<=1,'local', ' f is rescaling factor');
            if nargin <3 || isempty(method)
                method='bilinear';
            end
            sn=obj.signalname;
            oldsize=obj.size;
            % reduce --> imresize introduces shift!
            obj.set_signal(imresize(obj.xn,f,method, 'Antialiasing',true));
            % magnify
            obj.set_signal(imresize(obj.xn,oldsize,method, 'Antialiasing',true));
            
            obj.signalname=[sn, ', coarsened by ', num2str(f,'%4.2g')];
            obj.ensure(isequal(oldsize,obj.size),'local','image size stays the same');
            
        end
        
        function tile(obj, m, n)
            % tile image
            obj.require(isequal(floor(m),m) && isequal(floor(n),n),'local',...
                'needs 2 integers');
            obj.xn=kron(ones(m,n),obj.xn);
        end
        
        function resize(obj,ss,method)
            % ~(ss) resize image by factor ss (if scalar) or to new number of rows and columns
            obj.require(nargin<3 || ischar(method),'local','method is string');
            if length(ss)<2 && ss >10
                % interpret as new number of columns and rows
                ss=ss*[1,1];
            end
            if nargin <3 || isempty(method)
                method='bilinear';
            end
            sn=obj.signalname;
            obj.set_signal(imresize(obj.xn,ss,method, 'Antialiasing',true));
            obj.signalname=sn;
        end
        
        
        function shift(obj,shiftsize)
            % ~(shiftsize) shift signal by shiftlen (subpixels)
            obj.requireOnly(nargin>1 && length(shiftsize)==2,'local',...
                'shift is pair of numbers');
            obj.xn=Signal2D.SubPixelShifter(obj.xn, shiftsize);
        end
        
        function obj=transpose(obj)
            % ~() transpose data matrix
            obj.xn=obj.xn';
        end
        
        
    end
    
    %% queries
    methods
        
        function [xaxis,yaxis]=get_gridnodes(obj)
            % sampling positions for signal xn
            if obj.scale==PosScales.n
                hN=obj.N;
                yaxis=1:hN(1);
                xaxis=1:hN(2);
            else
                % linspace(0,obj.Ts,obj.N); % is wrong!
                [xaxis, yaxis]=obj.gridnodes(obj.N,obj.Rs,obj.origin);
            end
        end
        
        function s2=make_like(obj)
            % s2=~() clone all but signal content
            s2=make_like@SignalClass(obj);
            s2.Rs=obj.Rs;
        end
        
        function mask=mask_partition(obj,X,Y,j,k)
            % create mask for partition (X,Y)
            kept=[Y(j,k),Y(j+1,k)-1;X(j,k),X(j,k+1)-1];
            mask=false(obj.N);
            mask(kept(1,1):kept(1,2),kept(2,1):kept(2,2))=true;
        end
        
        function [X,Y]= make_partition(obj,m,n)
            % [X,Y]=~(m,n) make a partition of image into rectangles
            [X,Y]=meshgrid(floor(linspace(1,obj.size(2),n+1)),floor(linspace(1,obj.size(1),m+1)));
        end
        
        function vals=apply2partition(obj,func,m,n)
            % vals=~(func,m,n) apply function func to partitions of image
            obj.require(isa(func,'function_handle'),'local','func is a function');
            [X,Y]=obj.make_partition(m,n);
            ss=size(X);
            vals=cell(ss-[1,1]);
            for j=1:ss(1)-1
                for k=1:ss(2)-1
                    vals{j,k}=func(obj.mask_partition(X,Y,j,k));
                end
            end
        end
        
        function weights= eval_winfun(obj)
            % calculate weights from windows function
            hN=obj.N;
            if ~isempty(obj.winfun)
                if any(obj.padd>0)
                    hN=hN-2*obj.padd;
                end
                weights=win2D(hN(1),hN(2),obj.winfun);
            else
                weights=ones(obj.N);
            end
        end
        
        function r=SSIM(obj,yn,K,window,L,do_rescaling)
            % Structural similarity index of comparing yn with the present image obj.xn
            % obj.xn is assumed to be the ideal signal
            % in windows x,y of default size of 8x8
            % SSIM(x,y)= (2*mux*muy+c1)/(mux^2+muy^2+c2) *
            % (2*sigxy+c2)/(sigx^2+sigy^2+c2);
            % with small deniminator stabilizers  ci=(ki*L)^2;
            % k1=0.01; k2=0.03; L ... dynamic range.
            % symmetries of SSI:
            %-------------------
            % SSIM(f*x, f*y)= SSIM(x,y)
            % BUT translation symmetry is missing:  SSIM(x-t, y-t) != SSIM(x,y)
            %
            obj.requireOnly(numel(yn)==numel(obj.xn),'local','compatible test signal');
            if ~exist('K','var') || isempty(K)
                K = [0.01 0.03];
            end
            if ~exist('window','var') || isempty(window)
                window = fspecial('gaussian', 11, 1.5);	%
            end
            if ~exist('L','var') || isempty(L)
                L=obj.pv; % L=1? dynamic range of image
            end
            if ~exist('do_rescaling','var')|| isempty(do_rescaling)
                do_rescaling=true;
            end
            
            if isa(yn,'Signal2D')
                yn=yn.xn;
            end
            
            if obj.borderSize>0
                ns=obj.size-2*obj.borderSize;
                yn=wkeep(yn,ns);
                xnc=wkeep(obj.xn,ns);
                Imin=min(xnc(:));
                if do_rescaling
                    % mean square error
                    % minimized by linear rescaling:
                    p = polyfit(xnc(:),yn(:),1);
                    yn=p(1)*yn+p(2);
                end
                % shift to positive values, because of missing
                % shift invariance of SSIM
                r=ssim(xnc-Imin,yn-Imin,K,window,L);
            else
                Imin=min(obj.xn(:));
                if do_rescaling
                    p = polyfit(yn(:),obj.xn(:),1);
                    yn=p(1)*yn+p(2);
                else  % remove overshooting values
                    Imax=max(obj.xn(:));
                    yn(yn>Imax)=Imax;
                    yn(yn<Imin)=Imin;
                end
                % shift to positive values, because of missing
                % shift invariance of SSIM
                r=ssim(obj.xn-Imin,yn-Imin,K,window,L);
            end
            
        end
        
        function [mv,maxcorr]= motionVector(obj,signal2,mask)
            % mv=~(signal2) global shift between signals (mv(1)=rows, mv(2)=cols);
            % --- should be applied to contour images rather than to the
            %     original images;
            % --- for supixel precision use independent function
            %     EstimateMotionVector;
            %
            % method uses phase correlation, because based on the correlation
            % theorem: F(Corr(f,g))=F^*(f).*F(g)
            if nargin < 3
                mask=[];
            end
            
            % cross-correlation measured as phase correlation:
            norm1=obj.norm;
            norm2=signal2.norm;
            normfac=min(norm1,norm2)^2;
            if normfac<1e-9
                normfac=max(norm1,norm2)^2;
            end
            if ~isempty(mask)   % mask part of the signal
                xn2=signal2.xn.*mask;
            else
                xn2=signal2.xn;
            end
            CC=abs(ifft2(fft2(obj.xn).*conj(fft2(xn2))))/normfac;
            % only values close to 1 are reliable
            CC(CC<0.6 | CC>1.001)=0;
            % figure; surf(fftshift(CC)); colorbar;
            % determine location of max. cross correlation:
            [maxcorr,mv]=max(CC(:));
            [mv(1),mv(2)]=ind2sub(size(CC),mv);
            % shift length more than half image size:
            if mv(1)>size(CC,1)/2
                mv(1)=mv(1)-size(CC,1);
            end
            if mv(2)>size(CC,2)/2
                mv(2)=mv(2)-size(CC,2);
            end
            mv=mv-[1,1];
            % subpixel precision cannot be reached by simple parabolic fit
            % either on the image or on some contour image.
            % Instead of using the following parabolic fit, use
            % EstimateMotionVector;
            % fpcontrol.npeaks=1;  % look for only 1 peak
            % [peaks,params] = findpeaks2D(CC, fpcontrol);
            % mv=[peaks.rows_fit,peaks.cols_fit]-[1,1];
        end
        
        function s3=make_VideofromSingleImage(obj,L,sig)
            % s3=~(L,sig) make 3d signal of L frames with noise std sig
            obj.requireOnly(nargin <2 || isscalar(L),'local','frame number is scalar');
            if ~exist('L','var') || isempty(L)
                L=3;
            end
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            s3=Signal3D();
            x3=zeros([obj.N,L]);
            for j=1:L
                x3(:,:,j)=obj.xn +sig*randn(obj.N);
            end
            s3.set_signal(x3);
            s3.signalname=[obj.signalname,' STATIC'];
            s3.colormap_active=obj.colormap_active;
            
        end
        
        
        
    end
    
    methods (Hidden)
        
    end
    
    %% statics and factories
    methods (Static)
        
        
        function ts=make()
            % contructor with universal name (callable by parent)
            ts=Signal2D();
        end
        
        function s2= make_standardExample()
            % standard example to be redefined in sublasses
            s2=Signal2D.make_star();
        end
        
        function [xaxis,yaxis]=gridnodes(hN,hR,horigin)
            % [yax,xax]=~() sampling positions for signal xn
            
            % yaxis=linspace(0,hR(1),hN(1));
            % xaxis=linspace(0,hR(2),hN(2));
            % that gives wrong result in fft;
            % therefore do not close periodic curves.
            yaxis=hR(1)*(0:hN(1)-1)/hN(1)-horigin(1);
            xaxis=hR(2)*(0:hN(2)-1)/hN(2)-horigin(2);
            
        end
        
        function [X,Y]=meshgrid(hN,hR,horigin)
            % [X,Y]=~(hN,hR,[horigin]) yields right meshgrid for this class
            % if iscentered=true, meshgrid is centered on (0,0).
            % The order of the arguments of meshgrid is important!
            % meshgrid: 1st argument is mapped to columns (x) (#cols=hN(2))
            %           2nd argument to rows (y) (#rows=hN(1)
            % The inverse order of these arguments would prevent
            % FourierTrafos2D.dominantFreq from reconstructing omega !
            if nargin <3
                horigin=0;
            end
            [xaxis,yaxis]=Signal2D.gridnodes(hN,hR,horigin);
            [X,Y]=meshgrid(xaxis,yaxis);
        end
        
        function [x , y] = SpectralCoordinates(M , N)
            % [x,y] =~(M,N) Create normalized spectral coordinates
            X = ifftshift(((0 : N - 1) - floor(N / 2)) / N);
            Y = ifftshift(((0 : M - 1) - floor(M / 2)) / M);
            [x , y] = meshgrid(X , Y);
            
        end
        
        function ShiftedImage = SubPixelShifter(img , shiftsize)
            % SI=~(img, shiftsize) shifts img by shiftsize(1) rows and
            % shiftsize(2) columns (same convention as circshift)
            % shiftsize values need not be integers (subpixel shift)
            % Method: Fourier interpolation
            assert(length(shiftsize)==2,'shiftsize is pair of numbers');
            % Define normalized spectral coordinates
            [M , N , P] = size(img);
            [x , y] = Signal2D.SpectralCoordinates(M , N);
            
            % Generate mask for eliminating negative Fourier components wo . positive counterpart
            Mask = ones(M , N);
            if M == 2 * floor(M / 2)
                Mask(1 + M / 2 , :) = 0;
            end
            if N == 2 * floor(N / 2)
                Mask(: , 1 + N / 2) = 0;
            end
            
            % Generate Fourier spectrum of image and shifted image
            ShiftedImage = zeros(M , N , P);
            for p = 1 : P
                Spectrum = Mask .* fft2(img(: , : , p));
                ShiftedImage(: , : , p) = ...
                    real(ifft2(exp(-2 * pi * 1i * (x * shiftsize(2) + ...
                    y * shiftsize(1))) .* Spectrum));
            end
        end
        
        
        function s2= make_zernike(n,hN,sig)
            % s2= ~(n,[hN,sig]) n-tes Zernike  setzen. Zernikes sind nur in
            % der Einheitskreisscheibe definiert!
            assert(floor(n)==n && n>0,'n is pos. integer');
            if nargin <2 || isempty(sig)
                sig=0;
            end
            if nargin <3 || isempty(hN)
                hN=256;
            end
            if isscalar(hN)
                hN=[hN,hN];
            end
            
            hR=2*[1,1]; horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            [X,Y]=Signal2D.meshgrid(hN,hR,horigin);
            r=X.^2+Y.^2;
            xf=@(x,y) zrnxy(x,y,n);
            
            xdata=xf(X,Y)+sig*randn(size(X));
            xdata(r>1)=0;
            
            s2=Signal2D(xdata,hN,hR);
            
            s2.signalname=['Zernike Z_{',num2str(n),'}(x,y)',s2.noise_str(sig)];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
        end
        
        function s2= make_sinus(omega,del0,hN,sig)
            % constructor s2=~([omega,del0,hN,sig])
            % omega ... analog. freq of signal
            % del0 ... phase shift
            % sig ... additive noise std
            % hN ...length of signal
            assert(nargin<4 || isempty(hN) || (all(hN>0) && all(floor(hN)==hN)),...
                'size N is integer');
            if ~exist('omega','var') || isempty(omega)
                omega=5;
            end
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('del0','var') || isempty(del0)
                del0=0;
            end
            
            if isscalar(hN)
                hN=[hN,hN];
            end
            if isscalar(del0)
                del0=[del0,del0];
            end
            if isscalar(omega)
                omega=[omega,omega];
            end
            hR=2*pi+del0; horigin=[0,0]; %  => D=[0,Rs]
            
            xf=@(x,y) sin(omega(1)*x+omega(2)*y);
            [X,Y]=Signal2D.meshgrid(hN,hR,horigin);
            
            xdata=xf(X,Y)+sig*randn(size(X));
            
            s2=Signal2D(xdata,hN,hR);
            
            s2.signalname=['sin(',num2str(omega(1),'%3.2f'),'\cdot{x}+',...
                num2str(omega(2),'%3.2f'),'\cdot{y})',s2.noise_str(sig)];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
            
        end
        
        function s2= make_rect(w,hN,c,sig)
            % constructor s2=~([w,N,c,sig]), rectangle size w, upper left c,
            % image size N, noise std sig.
            requireOnly(nargin <3 || all(c>=1),'local','c is upper left corner');
            
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('w','var') || isempty(w)
                w=max(floor(hN/5),1);
            end
            if ~exist('c','var') || isempty(c)
                c=hN/2-w/2;
            end
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            
            if isscalar(hN)
                hN=hN*[1,1];
            end
            if isscalar(w)
                w=w*[1,1];
            end
            if isscalar(c)
                c=c*[1,1];
            end
            
            hR=2*hN/hN(1); horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            
            xdata=zeros(hN);
            lrc1=min(hN(1),max(1,c(1)+w(1)-1));
            lrc2=min(hN(2),max(1,c(2)+w(2)-1));
            
            xdata(c(1):lrc1,c(2):lrc2)=1;
            xdata=xdata+sig*randn(size(xdata));
            
            s2=Signal2D(xdata,hN,hR);
            
            s2.signalname=['rect(x,y)',s2.noise_str(sig)];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
        end
        
        function s2= make_step(hN,sig)
            % constructor s2=~([hN,sig]), with noise sig, size hN
            
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if isscalar(hN)
                hN=hN*[1,1];
            end
            
            hR=2*[1,1]; horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            xdata=zeros(hN);
            [X,Y]=Signal2D.meshgrid(hN,hR,horigin);
            xdata(X>0 & Y>0)=1;
            xdata=xdata+sig*randn(size(xdata));
            
            s2=Signal2D(xdata,hN,hR);
            
            s2.signalname=['step(x,y)',s2.noise_str(sig)];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
        end
        
        function s2= make_exp_symm(a,hN,sig)
            % s2=~([a,hN,sig]) symmetric exp with discontin. 1st deriv at 0
            % with noise sig, size hN, exp a
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var')|| isempty(hN)
                hN=256;
            end
            if ~exist('a','var') || isempty(a)
                a=-1;
            end
            
            if isscalar(hN)
                hN=hN*[1,1];
            end
            if isscalar(a)
                a=a*[1,1];
            end
            
            hR=2*[1,1];horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            [X,Y]=Signal2D.meshgrid(hN,hR,horigin);
            
            xdata=exp(a(1)*abs(X)+a(2)*abs(Y));
            
            xdata=xdata+sig*randn(size(xdata));
            
            s2=Signal2D(xdata,hN,hR);
            
            sn=['exp(',num2str(a(1),'%3.1f'),'\cdot{|x|}+',...
                num2str(a(2),'%3.1f'),'\cdot{|y|})',s2.noise_str(sig)];
            sn=strrep(sn,'+-','-');
            s2.signalname=sn;
            s2.origin=horigin;
            s2.set_scale(PosScales.x);
            
        end
        
        function s2= make_exp(a, hN,sig)
            % constructor s2=~([a,hN,sig]), with noise sig, size hN, exp a
            
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var')|| isempty(hN)
                hN=256;
            end
            if ~exist('a','var') || isempty(a)
                a=-1;
            end
            
            if isscalar(hN)
                hN=hN*[1,1];
            end
            if isscalar(a)
                a=a*[1,1];
            end
            
            hR=2*[1,1];horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            [X,Y]=Signal2D.meshgrid(hN,hR,horigin);
            
            xdata=exp(a(1)*X+a(2)*Y);
            xdata(X<0 | Y<0)=0;
            xdata=xdata+sig*randn(size(xdata));
            
            s2=Signal2D(xdata,hN,hR);
            
            sn=['exp(',num2str(a(1),'%3.1f'),'\cdot{x}+',...
                num2str(a(2),'%3.1f'),'\cdot{y})',s2.noise_str(sig)];
            sn=strrep(sn,'+-','-');
            s2.signalname=sn;
            s2.origin=horigin;
            s2.set_scale(PosScales.x);
        end
        
        function s2=make_star(v,hN,ngaps,sig)
            % s2=~([v,hN,ngaps,sig]) v radial lines with noise sig, size hN
            % and number of gaps [e.g. s2=Signal2D.make_star(8,[],20); ]
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var')|| isempty(hN)
                hN=256;
            end
            if ~exist('v','var') || isempty(v)
                v=8;
            end
            if ~exist('ngaps','var') || isempty(v)
                ngaps=0;
            end
            if isscalar(hN)
                hN=hN*[1,1];
            end
            
            hR=2*[1,1]; horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            xdata=zeros(hN)+sig*randn(hN);
            
            center=hN/2;
            angles=2*pi*(0:v-1)/v;
            x1=1-round(cos(angles)*(min(hN/2)-1));
            y1=1-round(sin(angles)*(min(hN/2)-1));
            x2=1+round(cos(angles)*(min(hN/2)-1));
            y2=1+round(sin(angles)*(min(hN/2)-1));
            
            for j=1:v
                [x y]=bresenham(x1(j),y1(j),x2(j),y2(j));
                for k=1:length(x)
                    xdata(center(1)+x(k),center(2)+y(k))=1;
                end
            end
            
            % dashed lines
            if ngaps>0
                rseg=linspace(0,floor(min(hN/2)), 4*ngaps+2);
                del=rseg(2)-rseg(1);
                pixorigin=ceil(size(xdata)/2);
                [X,Y]=meshgrid(1:size(xdata,2),1:size(xdata,2));
                rads=sqrt((X-pixorigin(2)).^2+(Y-pixorigin(2)).^2);
                for j=1:ngaps
                    xdata(rads>=rseg(4*j)-1.8*del & rads<=rseg(4*j)+1.8*del)=0;
                end
            end
            
            s2=Signal2D(xdata,hN,hR);
            
            s2.signalname='star';
            s2.origin=horigin;
            s2.set_scale(PosScales.x);
            s2.colormap_active='gray';
        end
        
        function s2= make_triang(w,hN,sig)
            % constructor s2=~([w,sig,hN]), with noise sig, size hN, width w
            
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('w','var') || isempty(w)
                w=max(floor(hN/5),1);
            end
            if isscalar(hN)
                hN=hN*[1,1];
            end
            if isscalar(w)
                w=w*[1,1];
            end
            
            hR=2*[1,1]; horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            w=hR.*w./hN;
            [X,Y]=Signal2D.meshgrid(hN,hR,horigin);
            
            xdata=2-abs(X)/w(1)-abs(Y)/w(2);
            xdata(abs(X)>w(1) | abs(Y)>w(2))=0;
            xdata=xdata+sig*randn(size(xdata));
            
            s2=Signal2D(xdata,hN,hR);
            
            s2.signalname=['triang(x,y)',s2.noise_str(sig)];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
        end
        
        function s2=make_sparse(hN, ct)
            % ~(N,ct) make random signal of size NxN and sparsity N^2/ct
            assert(ct>1 && floor(hN)==hN && min(hN)>1,'size and sparsity compression >1');
            if isscalar(hN)
                hN=hN*[1,1];
            end
            hR=2*[1,1]; horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            density=1/ct;
            xdata=full(sprandn(hN(1),hN(2),density));
            s2=Signal2D(xdata,hN,hR);
            s2.signalname=['density ', num2str(density,'%3.2g')];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
            
        end
        
        function pix= shiftdist(L, unit)
            % pix=~(L,unit) multiples of shift unit
            assert(nargin<2 || isempty(unit)||isscalar(unit),...
                'unit is scalar');
            if ~exist('L','var') || isempty(L)
                L=32;
            end
            if ~exist('unit','var') || isempty(unit)
                unit=11;
            end
            
            j=1:L-1;
            pix=j*unit;
        end
        
        function s2=make_fromImage(fn,sig,scale)
            % s2=~(fn,[sig,scale]) normalized from filename fn with noise added
            % if sig>0. Rescale image if scale missing or ~=1
            assert(exist(fn,'file')~=0,'file exists');
            if ~exist('sig','var')
                sig=0;
            end
            if ~exist('scale','var')
                scale=1;
            end
            xdata=imread(fn);
            if length(size(xdata))>2
                xdata=rgb2gray(xdata);
            end
            % to double and normalise to 1
            xdata=double(xdata);
            xdata=xdata-min(xdata(:));
            xdata=xdata/max(xdata(:));
            
            if sig>0
                xdata=xdata+sig*randn(size(xdata));
            end
            
            if scale~=1
                xdata=imresize(xdata,scale);
            end
            s2=Signal2D(xdata);
            [~,s2.signalname,~] = fileparts(fn);
            s2.colormap_active='gray';
            
        end
        
        
    end
    
    %% graphics
    methods
        function graph_signal(obj, open_new)
            % ~([open_new, clims] show signal in the open window using caxis
            if nargin <2
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig);
            end
            sn=obj.get_signalname;
            if isempty(sn)
                sn='Signal';
            end
            if obj.padd>0
                sn=[sn, ' 0-padded'];
            end
            
            [x,y]=obj.get_gridnodes;
            colormap(obj.colormap_active);
            clims=obj.colorlimits;
            if isempty(clims)
                % adapt color limits to outliers
                if obj.hide_outliers
                    [thresh1, thresh2]=obj.outliers(1e-3, 3);
                    clims = [ thresh1, thresh2];
                elseif isempty(obj.repfun)
                    clims = [ min(obj.xn(:)), max(obj.xn(:))];
                end
            end
            if ~isempty(clims) && clims(1)==clims(2)
                clims(2)=(1+1e-9)*clims(1);
                if clims(2)==clims(1)
                    clims(2)=clims(1)+1e-9;
                end
            end
            
            if isempty(obj.repfun)
                imagesc(x,y,obj.xn,clims);
            else
                imagesc(x,y,obj.repfun(obj.xn));
            end
            axis(obj.axis_type);
            cblabel_str=[];
            if ~isempty(obj.repfun)
                cblabel_str=func2str(obj.repfun);
            end
            if obj.ison_colorbar
                colorbar;
                if ~open_new && obj.colormap_freeze
                    freezeColors;
                    cbfreeze;
                end
                if ~isempty(obj.repfun)
                    cblabel(cblabel_str);
                end
            elseif ~open_new && obj.colormap_freeze
                freezeColors;
            end
            
            xlabel([obj.scale.symbol{1},deblank([' ',strtrim(obj.units{1})])],...
                'units','normalized','Position',[0.95,-0.08]);
            ylabel([obj.scale.symbol{2},deblank([' ',strtrim(obj.units{2})])]);
            title(sn,'fontsize',12);
            
        end
        
        function graph_partition(obj,m,n,open_new)
            % ~(m,n) show signal together with the partition boundaries
            if ~exist('open_new','var')
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig);
            end
            s2=obj.clone;
            smax=obj.max;
            [X,Y]=obj.make_partition(m,n);
            for j=1:size(X,2)
                s2.xn(:,X(1,j))=smax-s2.xn(:,X(1,j));
            end
            for j=1:size(Y,1)
                s2.xn(Y(j,1),:)=smax-s2.xn(Y(j,1),:);
            end
            s2.graph_signal(false);
            
        end
        
        function graph_winfun(obj, open_new)
            % show window function winfun in the open window
            if nargin <2
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig);
            end
            fv=obj.eval_winfun;
            [x,y]=obj.get_gridnodes;
            imagesc(x, y,fv); colormap(obj.colormap_active);
            if obj.ison_colorbar
                colorbar;
            end
            if ~open_new && obj.colormap_freeze
                freezeColors;
                cbfreeze;
            end
            title(['Windowing: ',obj.winfun2str, ', [min,max]=',...
                vec2str([min(fv(:)), max(fv(:))],'%3.1e')],...
                'fontsize',12);
            xlabel([obj.scale.symbol{2},deblank([' ',strtrim(obj.units{1})])]);
            ylabel([obj.scale.symbol{1},deblank([' ',strtrim(obj.units{2})])]);
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 2d';
            ok= isempty(obj.xn) || ismatrix(obj.xn);
        end
    end
    
end

