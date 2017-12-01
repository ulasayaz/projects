classdef TimeSignal<SignalClass
    % time- (or 1d-) signals
    % example:
    % omega=10; hN=101; sig=0.01; del0=0.11;w=10;a=-1;
    % ts=TimeSignal.make_sinus(omega,del0,hN,sig);
    % -- show signal:
    % ts.graph_signal();
    % w=20; sig=0;
    % ts=TimeSignal.make_rect(w,hN,sig);
    % --- change scale form index to time:
    % ts.set_scale(PosScales.t);
    % ts.graph_signal();
    % -- Chebyshev
    % n=3; ts=TimeSignal.make_chebyshev(n);
    %
    
    properties
        xt      %@<function_handle> contin. 1D time signal (default is interpolation of xn)
        Ts      %@<real> sampling time interval [s]
        marker  %@<string> plot symbol for single points
        linestyle %@<string> plot line style
        markercolor %@<string> marker color
    end
    
    properties (SetAccess=protected)
        ok_license     %<logical> all necessary toolbox licenses exist
    end
    
    %% constructor and commands
    methods
        function obj=TimeSignal(xdata,N,hT)
            obj.requireOnly(nargin<1 ||  isvector(xdata) ...
                || ( isa(xdata,'function_handle')...
                && nargin>=3 && obj.isnatural(N) ),'local',...
                'input is empty or vector or function plus sampling length & time interval');
            
            if nargin==0
                % needed to allow empty calls from subclasses
                % and prevent their premature testing of invariant
                xdata=[];
                N=100;
                hT=0;
            end
            if isnumeric(xdata) || islogical(xdata)
                obj.xn=xdata(:);
                N=length(xdata);
                if nargin<3
                    hT=2;  % sampling time interval [-1,1] of length 2
                end
                if N>1
                    obj.xt=@(t) interp1(obj.Ts*(0:N-1)/N,obj.xn,t);
                else
                    obj.xt=@(t) obj.xn(1);
                end
            else % isa(xdata,'function_handle')
                obj.xt=xdata;
                % funktioniert mit micht obj.xt=@(t) randn:
                % obj.xn=reshape(obj.xt(Ts*(0:N-1)/N),[],1);
                obj.xn=reshape(arrayfun(obj.xt,(hT*(0:N-1)/N)),[],1);
            end
            
            obj.Ts=hT;
            obj.scale=PosScales.n;
            obj.units='';
            obj.marker='';  % 'x'
            obj.linestyle='-';
            obj.markercolor='k';
            obj.winfun=[];
            obj.padd=0;
            obj.signalname=[];
            obj.origin=0;
            obj.fig=1;
            obj.ok_license= license('test','image_toolbox');
            
            obj.ensureOnly(size(obj.xn,2)<=1,'local','xn is a row vector');
            obj.ensure(isa(obj.xt,'function_handle'),...
                'local','xt is function');
            
        end % constructor
        
        function reset_Signal(obj)
            % reset signal trafo (e.g. after signal has changed)
            % redefined in subclasses
            
        end
        
        function set_tn(obj,hT,units)
            % set sampling positions and units of signal xn
            obj.requireOnly(hT>0,'local',...
                'hT is sampling interval length');
            obj.require(ischar(units),'local', 'unit is string');
            
            obj.Ts=hT;
            obj.units=units;
        end
        
        function set_origin(obj,z)
            % ~(z) set origin to z => domain D=[z-Ts/2, z+Ts/2]
            obj.requireOnly(isnumeric(z) & isscalar(z),'local',' z is scalar number');
            obj.origin=z;
        end
        
        function set_signal(obj,xn,hT)
            % ~(xn,[Ts]) set signal to vector xn
            obj.requireOnly(isvector(xn) ,'local','input is vector');
            if nargin <3 || isempty(hT)
                hT=2;
            end
            obj.reset_Signal;
            obj.signalname=[];
            
            obj.xn=xn(:);
            obj.padd=0;
            N=length(xn);
            obj.Ts=hT;
            obj.xt=@(t) interp1((0:N-1),obj.xn,t);
            
            obj.ensureOnly(size(obj.xn,2)<=1,'local','xn is a row vector');
            obj.ensure(isa(obj.xt,'function_handle'),...
                'local','xt is function');
        end
        
        function set_automaticStyle(obj)
            % ~() set graphics style adapted to data;
            % depending on separation of consecutive data points.
            pointdist=median(diff(find(~isnan(obj.xn))));
            if pointdist>2
                obj.marker='.';
                obj.linestyle='';
            else
                obj.marker='';
                obj.linestyle='-';
            end
        end
        
        function set_signalfun(obj,xt,N,Ts)
            % ~(xt,N,[Ts]) set a signal function and sample it at N points
            obj.requireOnly( isa(xt,'function_handle')...
                && obj.isnatural(N),'local',...
                'input is function plus sampling length');
            if nargin <4
                Ts=2;
            end
            obj.reset_Signal;
            obj.signalname=[];
            
            obj.xt=xt;
            obj.Ts=Ts;
            
            obj.xn=reshape(arrayfun(xt,(obj.Ts*(0:N-1)/N)),[],1);
            
            obj.ensureOnly(size(obj.xn,2)<=1,'local','xn is a row vector');
            obj.ensure(isa(obj.xt,'function_handle'),...
                'local','xt is function');
        end
        
        function weights= eval_winfun(obj)
            % calculate weights from windows function
            if ~isempty(obj.winfun)
                if obj.padd>0
                    hN=length(obj.xn)-2*obj.padd;
                    weights=obj.winfun(hN);
                    if obj.ok_license
                        weights=padarray(weights,obj.padd,0);
                    else
                        warning('TimeSignal:eval_winfun','function failed because of missing license');
                    end
                else
                    weights=obj.winfun(obj.length);
                end
            else
                weights=ones(obj.length,1);
            end
        end
        
    end
    
    %% transformation commands
    methods
        
        function set_0padding(obj,padd)
            % ~(padd) symm. 0-padding (left and right) by padd;
            % calling with padd=0 removes old padding;
            
            obj.requireOnly(obj.isnatural0(padd),'local',...
                'padd is non-negative integer');
            obj.reset_Signal;
            
            local_inv=obj.length()-obj.padd;
            if obj.padd>0
                % proportion of unpadded part
                prop=(obj.length-obj.padd)/obj.length;
                % decrease sampling time interval to unpadded value:
                obj.Ts=obj.Ts*prop;
                %remove old 0-padding:
                obj.remove_padding();
            end
            
            % increase sampling time interval to padded value
            obj.Ts=obj.Ts*(obj.length+2*padd)/obj.length;
            % add padding symmetrically
            obj.paddDir='both';
            pval=0;
            obj.xn=padarray(obj.xn,padd,pval, obj.paddDir);
            obj.padd=2*padd;
            
            obj.ensure(obj.length-obj.padd==local_inv,'local',...
                'padding ok');
        end
        
        function remove_padding(obj)
            % remove padding from signal
            if ~all(obj.padd==0)
                if strcmp(obj.paddDir,'post')
                    obj.xn = obj.xn(1:end-obj.padd);
                elseif strcmp(obj.paddDir,'both')
                    obj.xn = obj.xn(1+obj.padd/2:end-obj.padd/2);
                elseif strcmp(obj.paddDir,'pre')
                    obj.xn = obj.xn(1+obj.padd:end);
                end
                obj.padd=0;
            end
        end
        
        function interp1(obj, upsamp, ipmethod)
            % ~(upsm, [ipm]) interpolate by upsampling by a factor of upsm
            % using method ipm (default is 'spline').
            % REMARK
            %-----------
            % interp1 expands the principal analog frequency range
            % but does not increase the resolution in analog frequency.
            % Only sampling of longer time intervals (more periods!)
            % increases analog freq. resolution (for 0-padding only
            % the appearance of increased resolution).
            % Observe that DFT and DTFT
            % approximate the integrand of CFT only with a step function.
            % A linear or cubic interpolation of the signal improves this
            % approximation provided the signal is smooth enough.
            %
            obj.requireOnly(obj.isnatural(upsamp),'local',...
                'upsamp is a natural number');
            obj.requireOnly(nargin<3 || ischar(ipmethod),'local',...
                'ipmethod is a string denoting a method');
            if nargin < 3 || isempty(ipmethod)
                ipmethod='spline';
            end
            osn=obj.signalname;
            htn=obj.get_timenodes;
            % linspace should be true upsampling, i.e. keep old htn
            tn=linspace(htn(1),htn(end),(obj.length-1)*upsamp+1);
            
            ipsig=interp1(htn,obj.xn,tn,ipmethod);
            obj.set_signal(ipsig,obj.Ts); % use old Ts
            obj.signalname=osn;
            
        end
        
    end
    %% queries
    methods
        
        function sL=length(obj)
            % L=~(): sample length with padding inmcluded
            sL=length(obj.xn);
        end
        
        function ss=size(obj,k)
            % s=~(): sample size (needs same signature as in general case)
            ss=length(obj.xn);  % need scalar for invariant
        end
        
        
        function ss=size_dyadic(obj)
            % s=~(): ceiling of log2 of sample size
            ss=ceil(log2(length(obj.xn)));
        end
        
        function ts=get_timenodes(obj)
            % sampling positions for signal xn
            if obj.scale==PosScales.n
                ts=1:obj.length;
            else
                % ts=linspace(0,obj.Ts,obj.length); % is wrong!
                ts=obj.timenodes(obj.length,obj.Ts,obj.origin);
            end
        end
        
        
        function ts2=make_like(obj)
            % clone all properties but signal content
            ts2=make_like@SignalClass(obj);
            ts2.Ts=obj.Ts;
            ts2.marker=obj.marker;
            ts2.ok_license=obj.ok_license;
        end
        
        
        function xdiff=diff(obj, signal2)
            % xd=~([signal2]) difference with signal2 (mean value corrected)
            % or between adjacent values (approx. derivative)
            if nargin <2
                signal2=[];
            end
            if ~isempty(signal2)  % difference between two signals
                xdiff =obj.diff@SignalClass(signal2);
            else  % approximate derivative of present signal
                xdiff=obj.clone;
                xdiff.xn=diff(obj.xn);
                xdiff.signalname=['derivative ',obj.signalname];
                xdiff.amp_quantity=['\Delta ',obj.amp_quantity];
            end
            
        end
        
        function xdiff=difflog(obj)
            % xd=~([signal2]) log-difference between adjacent values (approx. derivative)
            xdiff=obj.clone;
            xdiff.xn=log(abs(diff(obj.xn)));
            xdiff.signalname=['log-derivative ',obj.signalname];
            xdiff.amp_quantity=['log-\Delta ',obj.amp_quantity];
        end
        
    end
    
    methods (Hidden)
        
        function sL=N(obj)
            % ~() sample size (gives at least a pair
            % of numbers in contrast to obj.size)
            sL=size(obj.xn);
        end
        
    end
    
    %% factories
    
    methods
        
        function s2= make_delta(obj,pos,hN,sig)
            % constructor s2=~([pos,hN,sig]]) delta signal at pos
            % sig ... aditiive noise std, hN ... length.
            % used to compute representation matrices of linear maps.
            
            if nargin <4 || isempty(sig)
                sig=0;
            end
            if nargin <3 || isempty(hN)
                hN=256;
            end
            if nargin <2 || isempty(pos)
                pos= floor(hN(1)/2);
            end
            
            hT=2; horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            xdata=zeros(hN,1);
            xdata(pos)=1;
            if sig>0
                xdata=xdata+sig*randn(size(xdata));
            end
            
            s2=obj.make();
            s2.set_signal(xdata,hT);
            
            s2.signalname=['\delta(t)',s2.noise_str(sig)];
            s2.set_scale(PosScales.t);
            s2.origin=horigin;
        end
    end
    methods (Static)
        
        
        function ts=make()
            % contructor with universal name (callable by parent)
            ts=TimeSignal();
        end
        
        function s2= make_standardExample()
            % standard example to be redefined in sublasses
            s2=TimeSignal.make_sinus();
        end
        
        
        function X = SpectralCoordinates(N)
            % [x,y] =~(M,N) Create normalized spectral coordinates
            X= ifftshift(((0 : N - 1) - floor(N / 2)) / N);
            
        end
        
        function ShiftedSig = SubPixelShifter(data , shiftsize)
            % SI=~(data, shiftsize) shifts data by shiftsize(1) columns and shiftsize(2) rows
            % shiftsize values need not be integers (subpixel shift)
            % Method: Fourier interpolation
            assert(isscalar(shiftsize),'shiftsize is scalar');
            % Define normalized spectral coordinates
            [N , P] = size(data);
            x = TimeSignal.SpectralCoordinates(N);
            
            % Generate mask for eliminating negative Fourier components wo . positive counterpart
            Mask = ones(N,1);
            if N == 2 * floor(N / 2)
                Mask(1 + N / 2) = 0;
            end
            
            % Generate Fourier spectrum of image and shifted image
            ShiftedSig = zeros(N , P);
            for p = 1 : P
                Spectrum = Mask .* fft2(data(: , p));
                ShiftedSig(: , p) = ...
                    real(ifft2(exp( -2 * pi * 1i * x(:) * shiftsize(1)  ...
                    ) .* Spectrum));
            end
        end
        
        function ts=timenodes(hN,hT,horigin)
            % ts=~(hN,hT,[horigin]) sampling positions for signal xn
            if nargin <3
                horigin=0;
            end
            % ts=linspace(0,hT,hN); % gives wrong result in fft;
            % therefore do not close periodic curves.
            ts=(0:hN-1)*hT/hN-horigin;
        end
        
        function ts=make_chebyshev(n,type,hN,sig)
            % ts=~(n,[type,hN,sig]) Tschebyscheff Polynom Typ type, Ordn. n
            % sig ... Rausch std, hN ... #Punkte
            assert(nargin>0 && (nargin<2 ||ismember(type,[1,2])),'type is in [1,2]');
            if ~exist('type','var') || isempty(type)
                type=1;
            end
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            
            hT=2; horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            % D= [-1,1]:
            if type==1
                xf=@(t) cos(n*acos(t));
                symb='T';
            else
                xf=@(t) sin((n+1)*acos(t))./sin(acos(t));
                symb='U';
            end
            xdata=xf(nv)+sig*randn(size(nv));
            ts=TimeSignal(xdata,hN,hT);
            ts.signalname=['Cheby ',symb,'_{',num2str(n),'}(t)',...
                ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
            
        end
        
        function ts= make_superpos(omegas,amps,hN,sig)
            % ts=~([omegas,amps,hN,sig]) superpostion of oscillations
            % omegas ... analog. frequencies of signal
            % amps ... amplitiudes of freqs.
            % sig ... additive noise std
            % hN ...length of signal
            assert(nargin<2 || length(omegas)==length(amps),...
                'same number of freqs and amplitudes');
            if ~exist('omegas','var') || isempty(omegas)
                omegas=[2.3,5.2];
            end
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('amps','var') || isempty(amps)
                amps=[1,0.5];
            end
            
            xf=@(t) sum( amps(:).*sin(omegas(:)*t));
            hT=2*pi; horigin=0; %  => D=[0,Ts]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            % nv=hT*(0:hN-1)/hN;
            xdata=arrayfun(xf,nv)+sig*randn(size(nv));
            
            ts=TimeSignal(xdata,hN,hT);
            
            ts.signalname=['oscillations \omega=',vec2str(omegas,'%3.1f'),', ',...
                'a=',vec2str(amps,'%3.1f'), ...
                ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
            
        end
        
        function ts= make_sinus(omega,del0,hN,sig)
            % constructor ts=~([omega,del0,hN,sig])
            % omega ... analog. freq of signal
            % del0 ... phase shift
            % sig ... additive noise std
            % hN ...length of signal
            assert(nargin<4 || isempty(hN) || (hN>0 && floor(hN)==hN),...
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
            
            xf=@(t) sin(omega*t);
            hT=2*pi+del0; horigin=0; %  => D=[0,Ts]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            % nv=hT*(0:hN-1)/hN;
            xdata=xf(nv)+sig*randn(size(nv));
            
            ts=TimeSignal(xdata,hN,hT);
            
            del0_str=[];  % phase shift
            if del0~=0
                del0_str=['+',num2tex(del0,'%3.1e','none')];
            end
            ts.signalname=['sin(',num2str(omega,'%3.2f'),'\cdot{t}',...
                del0_str,') ', ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
        end
        
        function ts= make_sinusRamp(omega,del0,hN,sig)
            % constructor ts=~([omega,del0,hN,sig])
            % omega ... analog. freq of signal
            % del0 ... phase shift
            % sig ... additive noise std
            % hN ...length of signal
            assert(nargin<4 || isempty(hN) || (hN>0 && floor(hN)==hN),...
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
            
            xf=@(t) sin(omega*t)+t;  % ramp adds y=t
            hT=2*pi+del0; horigin=0; %  => D=[0,Ts]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            % nv=hT*(0:hN-1)/hN;
            xdata=xf(nv)+sig*randn(size(nv));
            
            ts=TimeSignal(xdata,hN,hT);
            
            del0_str=[];  % phase shift
            if del0~=0
                del0_str=['+',num2tex(del0,'%3.1e','none')];
            end
            ts.signalname=['sin(',num2str(omega,'%3.2f'),'\cdot{t}',...
                del0_str,')+t ', ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
        end
        
        function ts= make_rect(w,hN,sig)
            % constructor ts=~([w,hN,sig]), width w, noise sig, length hN
            requireOnly(nargin<1 || isempty(w) || (w>=1 && floor(w)==w),'local',...
                'width is natural number');
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('w','var') || isempty(w)
                w=max(floor(hN/5),1);
            end
            
            hT=2; horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            wt=w*hT/hN;
            
            xdata=zeros(hN,1);
            xdata(abs(nv)<=wt)=1;
            xdata=xdata+sig*randn(size(xdata));
            
            ts=TimeSignal(xdata,hN,hT);
            
            ts.signalname=['rect(t)',ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
        end
        
        function ts= make_interv(w,hN,LB,sig)
            % ts=~([w,hN,LB,sig]), width w, length hN, left LB, noise sig
            requireOnly(nargin<1 || isempty(w) || (w>=1 && floor(w)==w),'local',...
                'width is natural number');
            requireOnly(nargin<3 || isempty(LB) || floor(LB)==LB,'local',...
                'left border is integer (possibly neg.)');
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('w','var') || isempty(w)
                w=max(floor(hN/5),1);
            end
            if ~exist('LB','var') || isempty(LB)
                LB=min(hN,max(1,round(hN/2-w/2)));
            end
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            
            if LB<1
                w=w+LB-1;
                LB=1;
            end
            
            hT=2; horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            xdata=zeros(hN,1);
            if w>0
                lrc1=min(hN(1),max(1,LB(1)+w(1)-1));
                xdata(LB:lrc1)=1;
            end
            xdata=xdata+sig*randn(size(xdata));
            
            ts=TimeSignal(xdata,hN,hT);
            
            ts.signalname=['interval, left=',num2str(LB),',width=',num2str(w)];
            ts.set_scale(PosScales.n);
            ts.origin=horigin;
        end
        
        function ts= make_step(hN,sig)
            % constructor ts=~([hN,sig]), with noise sig, length hN
            
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            
            hT=2; horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            
            xdata=zeros(hN,1);
            xdata(nv>=0)=1;
            xdata=xdata+sig*randn(size(xdata));
            
            ts=TimeSignal(xdata,hN,hT);
            
            ts.signalname=['step(t)',ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
        end
        
        function ts= make_exp(a,hN,sig)
            % constructor ts=~([a,hN,sig]), with noise sig, length hN, exp a
            
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('a','var') || isempty(a)
                a=-1;
            end
            
            hT=2; horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            xdata=zeros(hN,1);
            ydata=exp(a*abs((1:hN)-floor(hN/2)));
            xdata(floor(hN/2):end)=ydata(floor(hN/2):end);
            xdata=xdata+sig*randn(size(xdata));
            
            ts=TimeSignal(xdata,hN,hT);
            
            ts.signalname=['exp(',num2tex(a,'%3.2f'),'\cdot{t})',ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
        end
        
        function ts= make_triang(w,hN,sig)
            % constructor ts=~([w,hN,sig]), with noise sig, length hN, width w
            
            if ~exist('sig','var') || isempty(sig)
                sig=0;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            if ~exist('w','var') || isempty(w)
                w=max(floor(hN/5),1);
            end
            
            hT=2;  horigin=hT/2; %  => D=[-Ts/2,Ts/2]
            % sampled nodes:
            nv=TimeSignal.timenodes(hN,hT,horigin);
            xdata=zeros(hN,1);
            wt=w*hT/hN;
            % args=linspace(-hN/2,hN/2,hN);
            interv= abs(nv)<=wt;
            ydata=1-abs(nv)/wt;
            xdata(interv)=ydata(interv);
            
            xdata=xdata+sig*randn(size(xdata));
            
            ts=TimeSignal(xdata,hN,hT);
            
            ts.signalname=['triang(t)',ts.noise_str(sig)];
            ts.set_scale(PosScales.t);
            ts.origin=horigin;
        end
        
        function ns=noise_str(sig)
            % text output for noisy signals
            if sig==0
                ns=[];
            else
                ns=['+\epsilon_\sigma, \sigma=',...
                    num2tex(sig,'%3.2e','none')];
            end
        end
        
    end
    
    %% graphics
    methods
        function graph_signal(obj,open_new)
            % ~([open_new]) show signal in the open window by default
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            sn=obj.get_signalname;
            if obj.padd>0
                sn=[sn, ' 0-padded'];
            end
            
            tnv=obj.get_timenodes;
            markersize=10;
            
            % do not set line color in plot, so that
            % one can still cylce through the colororder
            % by calling hold all.
            if isempty(obj.repfun)
                plot(tnv,obj.xn,[obj.marker,obj.linestyle],...
                    'MarkerEdgeColor',obj.markercolor,...
                    'MarkerFaceColor',obj.markercolor,...
                    'MarkerSize',markersize);
            else
                plot(tnv,obj.repfun(obj.xn),[obj.marker,obj.linestyle],...
                    'MarkerEdgeColor',obj.markercolor,...
                    'MarkerFaceColor',obj.markercolor,...
                    'MarkerSize',markersize);
                legend(func2str(obj.repfun),'location','best');
            end
            title(sn,'fontsize',12);
            xlabel([obj.scale.label,' ',obj.units]);
            ylabel([obj.amp_quantity,' ',obj.amp_unit]);
            axis tight;
            
            
        end
        
        
        function graph_winfun(obj, open_new)
            % ~([open_new]) show window function winfun in the open window by default
            if nargin < 2 || open_new
                prepfigure(obj.fig);
            end
            fv=obj.eval_winfun;
            plot(obj.get_timenodes,fv,[obj.marker,obj.linestyle],...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r');
            title(['Windowing: ',obj.winfun2str, ', [min,max]=',...
                vec2str([min(fv), max(fv)],'%3.1e')],...
                'fontsize',12);
            xlabel([obj.scale.label,' ',obj.units]);
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 1d';
            ok= isempty(obj.xn) || isvector(obj.xn);
        end
    end
    
end

