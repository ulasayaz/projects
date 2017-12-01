classdef Signal3D<SignalClass
    % 3D signals wrapping a stack of 2d images (video)
    % GT Stuttgart, 2014
    % example:
    %
    % --- read 32 frame starting from frame 33 of a video file:
    %  fn='1920x1080x50p_11_nn0497_tennismen.avi';
    %  signal3=Signal3D.make_fromVideo(fn,[33,32]);
    % --- reduce frames in size:
    %  signal3.signal2square(512);
    % --- show video as montage
    %  signal3.graph_signal(true,4);
    % -- play back video:
    %  signal3.play_signal;
    % -- show total variation per pixel in temporal direction:
    %  signal3.graph_totalVar;
    % -- save video:
    %  signal3.writeToVideoFile('tennis.avi');
    % --- cf. signal factory: Signal_Factory
    %
    
    properties
        Rs     %@<real> signal sampling rectangle
        video_sim  %@<integer> 1 is rotation (default) , 2 is shift
        zref
    end
    
    properties (SetAccess=protected)
        
    end
    
    %% constructor and commands
    methods
        function obj=Signal3D(xdata,hN,hR)
            obj.requireOnly(nargin<1 || ismember(length(size(xdata)),[2,3]) ...  ...
                || ( isa(xdata,'function_handle')...
                && nargin>=2 && obj.isnaturalArray(hN) ),'local',...
                'input is empty or 3d or function plus sampling size & range');
            
            if nargin==0
                % needed to allow empty calls from subclasses
                % and prevent their premature testing of invariant
                xdata=[];
                hN=[0,0,0];
                hR=[0,0,0];
            end
            if nargin <3
                hR=1;
            end
            if isscalar(hR)
                hR=hR*[1,1,1];
            end
            if isnumeric(xdata) || islogical(xdata)
                obj.xn=xdata;
                obj.set_scale(PosScales.n);
            else % isa(xdata,'function_handle')
                obj.xt=xdata;
                obj.set_scale(PosScales.x);
            end
            
            obj.Rs=hR;
            obj.scale=PosScales.n;
            obj.units={'','',''};
            obj.colormap_active='default';  % e.g. 'default'
            obj.colormap_freeze=true;
            obj.padd=0;
            obj.signalname=[];
            obj.origin=[0,0,0];
            obj.fig=1;
            obj.ison_colorbar =true;
            obj.video_sim=1;
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
                units={units,units,units};
            end
            obj.units=units;
        end
        
        function set_zref(obj,j)
            % ~(j) set reference frame
            obj.requireOnly(j>=1 && j<=obj.size(3),'local','valid index');
            obj.zref=j;
            obj.colorlimits=[];
        end
        
        function set_winfun(obj,wf)
            % ~(wf) set windowing function to wf
            obj.requireOnly( isempty(wf)|| isa(wf,'function_handle'),'local',...
                'input is empty or function handle');
            
            obj.winfun=wf;
            
        end
        
        function set_signal(obj,xn,hR)
            % ~(xn,[hR]) set signal to matrix xn
            obj.requireOnly(isnumeric(xn) || islogical(xn),'local','input is matrix');
            if nargin <3 || isempty(hR)
                hR=size(xn);
            end
            obj.reset_Signal;
            obj.signalname=[];
            
            obj.xn=xn;
            obj.set_scale(PosScales.n);
            obj.padd=[0,0];
            obj.Rs=hR;
            
            obj.ensureOnly(isempty(obj.xn) || length(size(obj.xn))==3 ,...
                'local','xn is a 3d-matrix');
        end
        
        
        function set_origin(obj,z)
            % ~(z) set origin to z => domain D=[z-Ts/2, z+Ts/2]
            obj.requireOnly(length(z)==3,'local',' z is triplet of numbers');
            obj.origin=z;
        end
        
        function writeToVideoFile(obj,fn)
            % ~(fn) write 3d-signal to video file
            % Prepare the new file.
            vidObj = VideoWriter(fn);
            open(vidObj);
            % transform to range [0,1] for play back!
            minval=min(obj.xn(:));
            maxval=max(obj.xn(:));
            for k = 1:obj.size(3)
                % Write each frame to the file.
                writeVideo(vidObj,(obj.xn(:,:,k)-minval)/(maxval-minval));
            end
            % Close the file.
            close(vidObj);
        end
        
        function repSignal(obj,L)
            % ~(L) replicate signal L times
            obj.requireOnly(obj.isnatural(L),'local','needs natrual number');
            obj.xn=repmat(obj.xn,[1,1,L]);
        end
        
        
    end
    
    %% transformation commands
    
    methods
        
        function crop(obj,ss)
            % ~([ss]) crop to size ss (default lower dyadic bound)
            if nargin <2
                ss=min(2.^floor((log2(obj.size))));
            end
            ssold=size(obj.xn);
            if isscalar(ss) && ~isvector(obj.xn)
                ss=ss*ones(1,length(ssold));
            end
            ss=min(ss,ssold);
            s2=obj.clone;
            obj.xn=zeros(ss,'single'); % saves memory
            wait_handle = waitbar(0,...
                ['cropping ',num2str(ss(3)),' video frames ...']);
            for j=1:ss(3)
                waitbar(j /ss(3));
                % central part
                obj.xn(:,:,j)=wkeep(s2.xn(:,:,j),[ss(1),ss(2)]);
            end
            clear('s2');
            obj.xn=double(obj.xn);
            close(wait_handle );
        end
        
        function resize(obj,ss,method)
            % ~(ss) resize image by by factor ss (if scalar) or to new number of rows and columns
            if length(ss)<2 && ss >10
                % interpret as new number of columns and rows
                ss=ss*[1,1];
            end
            if nargin <3 || isempty(method)
                method='bilinear';
            end
            nFrames = obj.size(3);
            if max(ss(1),ss(2))>10
                vidHeight=ss(1);
                vidWidth =ss(2);
            else
                vidHeight = ceil(obj.size(1)*ss(1));
                vidWidth = ceil(obj.size(2)*ss(2));
            end
            s2=obj.clone;
            obj.xn=zeros([vidHeight,vidWidth,nFrames],'single'); % saves memory
            wait_handle = waitbar(0,...
                ['resizing ',num2str(nFrames),' video frames ...']);
            for j=1:obj.size(3)
                waitbar(j / nFrames);
                obj.xn(:,:,j)=imresize(s2.xn(:,:,j),[vidHeight,vidWidth],method, 'Antialiasing',true);
            end
            clear('s2');
            obj.xn=double(obj.xn);
            close(wait_handle );
        end
        
        
        function signal2square(obj,ss)
            % ~([ss]) crop to next dyadic size and then resize to square
            % of length ss
            obj.requireOnly(ss(1)>0 && (isscalar(ss) || ss(1)==ss(2)),'local',...
                ' ss is length of new frame size');
            sdyadic=2.^floor(log2(min(obj.size(1),obj.size(2))));
            obj.crop([sdyadic,sdyadic,obj.size(3)]);
            obj.resize(ss);
        end
        
        function remove_padding(obj)
            % remove padding from signal
            if ~all(obj.padd==0)
                if strcmp(obj.paddDir,'post')
                    obj.xn = obj.xn(1:end-obj.padd(1),1:end-obj.padd(2),...
                        1:end-obj.padd(3));
                elseif strcmp(obj.paddDir,'both')
                    obj.xn = obj.xn(1+obj.padd(1)/2:end-obj.padd(1)/2,...
                        1+obj.padd(2)/2:end-obj.padd(2)/2,...
                        1+obj.padd(3)/2:end-obj.padd(3)/2);
                elseif strcmp(obj.paddDir,'pre')
                    obj.xn = obj.xn(1+obj.padd(1):end,1+obj.padd(2):end,...
                        1+obj.padd(3):end);
                end
                obj.padd=0;
            end
        end
        
        function invert_time(obj)
            % ~() video backward in time (3d dimension)
            obj.xn=flipdim(obj.xn,3);
        end
        
        function circshift_time(obj)
            % ~() circular shift in time
            obj.xn=circshift(obj.xn,[0,0,1]);
        end
        
    end
    
    %% queries
    methods
        
        function s3=make_like(obj)
            % s2=~() clone all but signal content
            s3=make_like@SignalClass(obj);
            s3.Rs=obj.Rs;
            s3.zref=obj.zref;
        end
        
        function s=frame(obj,j)
            % s=~(j) retrieves frame number j from signal
            obj.requireOnly(min(j)>=1 && max(j)<=obj.size(3),'local',...
                'j admissible');
            if isscalar(j)
                s=Signal2D(obj.xn(:,:,j));
            else
                s=Signal3D(obj.xn(:,:,j));
            end
            sn=obj.signalname;
            if isempty(sn)
                sn='signal';
            end
            s.signalname=[sn,', z=',vec2str(j)];
            s.colormap_active=obj.colormap_active;
            s.fig=obj.fig;
        end
        
        function s=frameDiff(obj,j,k)
            % s=~(j,k) difference between frame j and frame k
            obj.requireOnly(min(j,k)>=1 && max(j,k)<=obj.size(3),'local',...
                'j,k admissible');
            s=Signal2D(obj.xn(:,:,j)-obj.xn(:,:,k));
            s.signalname=['frame diff of ',obj.signalname,', z=',vec2str([j,k])];
            s.colormap_active='default';
        end
        
        function s=diffSignal(obj,shiftsize)
            % s=~(j,k) difference between between consecutive frames
            if nargin <2
                shiftsize=1;
            end
            s=obj.make_like();
            s.xn=diff(obj.xn,shiftsize,3);
            s.signalname=['diff. signal of ',obj.signalname,', z=',num2str(j)];
        end
        
        function s=quotientSignal(obj,shiftsize)
            % s=~(j,k) quotient between consecutive frames
            if nargin <2
                shiftsize=1;
            end
            s=obj.make_like();
            L=obj.size(3);
            yn=zeros(obj.N-[0,0,shiftsize]);
            for z1=1:L-shiftsize
                yn(:,:,z1)=obj.xn(:,:,z1+shiftsize)./obj.xn(:,:,z1);
            end
            s.xn=yn;
            s.signalname=['quot. signal of ',obj.signalname];
            s.colormap_active='default';
        end
        
        function s2=cross_section(obj,section)
            % s2=~(section) cross section defined by string section
            obj.requireOnly(ischar(section)&& length(strfind(section,':'))==2,...
                'local','needs section of form "n,:,:"');
            cmd=['squeeze(obj.xn(',section,'))'];
            try
                img=eval(cmd);
            catch
                error('value for section is invalid');
            end
            
            s2=Signal2D(img);
            s2.signalname=[obj.signalname,', section ',section];
            obj.ensureOnly(isa(s2,'Signal2D'),'local','returns 2d signal');
        end
        
        function r=PSNR(obj, sig2)
            % r=~(signal2) ... peak signal to noise ratio of reconstructed signal2 in dB
            % for each pair of frames of same z-value.
            % obj.xn is assumed to be the ideal signal
            % used to test quality of similarity of yn to signal xn
            obj.requireOnly(isa(sig2,'Signal3D') && sig2.numel==obj.numel,...
                'local','compatible test signal');
            
            L=obj.size(3);
            r=zeros(L,1);
            
            for j=1:L
                % compute separately for each plane
                x=Signal2D(obj.xn(:,:,j));
                y=sig2.xn(:,:,j);
                r(j)=x.PSNR(y);
            end
        end
        
        function r=PSNRcross(obj, sig2)
            % r=~(signal2) ... cross PSNR between all pairs of frames
            obj.requireOnly(nargin<2 || isempty(sig2) || ...
                (isa(sig2,'Signal3D') && sig2.numel==obj.numel),...
                'local','compatible test signal');
            
            if nargin<2 || isempty(sig2)
                sig2=obj;
            end
            L=obj.size(3);
            r=NaN+zeros(L,L);
            
            if sig2==obj
                offset=1; % PSNR of frame with itself would be outlier
            else
                offset=0;
            end
            for z1=1:L
                x=Signal2D(obj.xn(:,:,z1));
                for z2=z1+offset:L
                    % compute separately for each plane
                    y=sig2.xn(:,:,z2);
                    r(z1,z2)=x.PSNR(y);
                end
            end
        end
        
        
        function r=PSNRshift(obj, shiftsize)
            % r=~([shiftsize]) ...  PSNR between pairs of consecutive frames;
            % is a measure of the amount of global or local motion between
            % frames
            obj.requireOnly(nargin<2 || floor(shiftsize)==shiftsize,'local',...
                'shiftsize is integer');
            
            if nargin<2 || isempty(shiftsize)
                shiftsize=1;
            end
            shiftsize=abs(shiftsize);
            L=obj.size(3);
            r=zeros(L-shiftsize,1);
            for z1=1:L-shiftsize
                y=Signal2D(obj.xn(:,:,z1));
                x=Signal2D(obj.xn(:,:,z1+shiftsize));
                r(z1)=x.PSNR(y);
            end
        end
        
        function r=SSIM(obj,sig2,K,window)
            % r=~(signal2) Structural similarity index between pairs of
            % frames with the same z-value.
            % obj.xn is assumed to be the ideal signal
            obj.requireOnly(isa(sig2,'Signal3D') && sig2.numel==obj.numel,...
                'local','compatible test signal');
            if nargin <3 || isempty(K)
                K = [0.01 0.03];
            end
            if nargin <4 || isempty(window)
                window = fspecial('gaussian', 11, 1.5);	%
            end
            
            L=obj.size(3);
            r=zeros(L,1);
            
            for z=1:L
                % compute separately for each plane
                x=Signal2D(obj.xn(:,:,z));
                y=sig2.xn(:,:,z);
                r(z)=x.SSIM(y,K,window);
            end
            
        end
        
        function r=SSIMcross(obj,sig2,K,window)
            % r=~(signal2) cross SSIM between all pairs of frames.
            % obj.xn is assumed to be the ideal signal
            obj.requireOnly(nargin<2 || isempty(sig2) || ...
                (isa(sig2,'Signal3D') && sig2.numel==obj.numel),...
                'local','compatible test signal');
            if nargin <2 || isempty(sig2)
                sig2=obj;
            end
            if nargin <3 || isempty(K)
                K = [0.01 0.03];
            end
            if nargin <4 || isempty(window)
                window = fspecial('gaussian', 11, 1.5);	%
            end
            
            L=obj.size(3);
            r=NaN+zeros(L,L);
            if sig2==obj
                offset=1;  % SSIM of frame with itself would be outlier
            else
                offset=0;
            end
            for z1=1:L
                x=Signal2D(obj.xn(:,:,z1));
                for z2=z1+offset:L
                    y=sig2.xn(:,:,z2);
                    r(z1,z2)=x.SSIM(y,K,window);
                end
            end
            
        end
        
        function r=SSIMshift(obj, shiftsize,K,window)
            % r=~([shiftsize]) ...  PSNR between pairs of consecutive frames;
            % is a measure of the amount of global or local motion between
            % frames.
            obj.requireOnly(nargin<2 || floor(shiftsize)==shiftsize,'local',...
                'shiftsize is integer');
            
            if nargin<2 || isempty(shiftsize)
                shiftsize=1;
            end
            if nargin <3 || isempty(K)
                K = [0.01 0.03];
            end
            if nargin <4 || isempty(window)
                window = fspecial('gaussian', 11, 1.5);	%
            end
            shiftsize=abs(shiftsize);
            L=obj.size(3);
            r=zeros(L-shiftsize,1);
            
            for z1=1:L-shiftsize
                y=Signal2D(obj.xn(:,:,z1));
                x=Signal2D(obj.xn(:,:,z1+shiftsize));
                r(z1)=x.SSIM(y,K,window);
            end
        end
        
        function [v,l1]=totalVar(obj,dim)
            % [v,l1]=~(dim) pixelwise 1-d total variation along dimension dim
            % optional output: l1-norm
            if nargin<2
                dim=3;
            end
            v= sum(abs(diff(obj.xn,1,dim)),dim);
            if nargout>1
                l1=sum(v(:));
            end
        end
        
        
        function v=totalVarShift(obj,shiftsize)
            % v=~(dim) relative 1-d total variation between consecutive images in z
            % is a measure of the amount of global or local motion between
            % frames
            if nargin <2
                shiftsize=1;
            end
            L=obj.size(3);
            v=zeros(L-shiftsize,1);
            for z1=1:L-shiftsize
                v(z1)= sum(abs(reshape(obj.xn(:,:,z1+shiftsize)-obj.xn(:,:,z1),[],1)))/...
                    sum(abs(reshape(obj.xn(:,:,z1+shiftsize)+obj.xn(:,:,z1),[],1)));
            end
        end
        
        function [s2,params]= denoise(obj, mra,params)
            % s2=~([mra,params]) denoise signal; must be redefined!
            % mra is of class MultiResTransform, e.g. mra=Wavelet2D_mlab()
            % example for 2d signal: s2=obj.denoise(Wavelet2D_mlab());
            obj.requireOnly(isa(mra,'MultiResTransform2'),'local',...
                'denoising each 2d frame separately');
            if ~exist('params','var') || ~isstruct(params)
               params=struct;
            end
            L=obj.size(3);
            s2=obj.clone;
            for j=1:L
                [s2d,params]=denoise@SignalClass(obj.frame(j),mra,params);
                s2.xn(:,:,j)=s2d.xn;
            end
            s2.signalname=s2d.signalname;    
            
        end
        
        function z0=midplane(obj)
            % z0=~() index of central frame
            z0=max(1,ceil(obj.size(3)/2));
        end
        
        function clims=get_colorlimits(obj)
            % clims=~() get or compute colorlimits
            clims=obj.colorlimits;
            if isempty(clims)
                if isempty(obj.zref)
                    clims=[obj.min-eps,obj.max+eps];
                else
                    clims=[obj.frame(obj.zref).min-eps,obj.frame(obj.zref).max+eps];
                end
            end
        end
        
    end
    
    
    %% factories
    methods (Static)
        
        function [ssize,vid,err_videoreader]=VideoSize(fn)
            % ss=~(fn) size of video contained in file fn
            assert(ischar(fn) && exist(fn,'file')~=0,'file exists');
            try % video reader has problems with remote connections
                vid = VideoReader(fn);
                ssize=[vid.Height,vid.Width,vid.NumberOfFrames];
                err_videoreader=false;
            catch err
                err_videoreader=true;
                vid= mmread(fn);
                ssize=[vid.height,vid.width,vid.nrFramesTotal];
            end
            
        end
        
        function ts=make()
            % contructor with universal name (callable by parent)
            ts=Signal3D();
        end
        
        function s2=make_UsingVideoReader(vid, clip,scale,sig)
            % s2=~(fn,[clip,scale,sig]) create from Videoreader object vid
            assert(isa(vid,'VideoReader'),'vid is a VideoReader object');
            assert(nargin<2 || isnumeric(clip),'clip is #frames or start frame and #frames');
            if ~exist('sig','var')
                sig=0;
            end
            if ~exist('scale','var')
                scale=1;
            end
            if ~exist('clip','var')
                clip=[1,32];
            end
            if isscalar(clip)
                clip=[1,clip];
            end
            nFrames = min(clip(2),vid.NumberOfFrames);
            vidHeight = ceil(vid.Height*scale);
            vidWidth = ceil(vid.Width*scale);
            test=read(vid, 1);
            if isa(test,'integer')
                if isa(test,'uint8')
                    normfac=255;
                else
                    normfac=2^ceil(log2(double(max(vid.frames(1).cdata(:)))))-1;
                end
            else
                normfac=1;
            end
            
            
            % Preallocate movie structure.
            s2=Signal3D(zeros(vidHeight,vidWidth,nFrames));
            
            % Read one frame at a time.
            wait_handle = waitbar(0,...
                ['importing ',num2str(nFrames),' video frames ...']);
            for k = 1 : nFrames
                waitbar(k / nFrames);
                frameidx= clip(1)+k-1;
                if scale==1
                    s2.xn(:,:,k)= double(rgb2gray(read(vid, frameidx)))/normfac;
                else
                    s2.xn(:,:,k)=imresize(double(rgb2gray(read(vid, frameidx)))/normfac,...
                        [vidHeight,vidWidth]);
                end
                if sig>0
                    s2.xn(:,:,k)=s2.xn(:,:,k)+sig*randn([vidHeight,vidWidth]);
                end
            end
            close(wait_handle );
            drawnow;
            [~,s2.signalname,~] = fileparts(vid.name);
            s2.colormap_active='gray';
            
        end
        
        function s2=make_UsingMMReader(fn, clip,scale,sig)
            % s2=~(fn,[clip,scale,sig]) create from Videoreader object vid
            assert(ischar(fn) && exist(fn,'file')~=0,'file exists');
            assert(nargin<2 || isnumeric(clip),'clip is #frames or start frame and #frames');
            if ~exist('sig','var')
                sig=0;
            end
            if ~exist('scale','var')
                scale=1;
            end
            if ~exist('clip','var')
                clip=[1,32];
            end
            if isscalar(clip)
                clip=[1,clip];
            end
            
            vid= mmread(fn);
            if isa(vid.frames(1).cdata,'integer')
                if isa(vid.frames(1).cdata,'uint8')
                    normfac=255;
                else
                    normfac=2^ceil(log2(double(max(vid.frames(1).cdata(:)))))-1;
                end
            else
                normfac=1;
            end
            
            nFrames = min(clip(2),vid.nrFramesTotal);
            vidHeight = ceil(vid.height*scale);
            vidWidth = ceil(vid.width*scale);
            
            % Preallocate movie structure.
            s2=Signal3D(zeros(vidHeight,vidWidth,nFrames));
            
            % Read one frame at a time.
            wait_handle = waitbar(0,...
                ['importing ',num2str(nFrames),' video frames ...']);
            for k = 1 : nFrames
                waitbar(k / nFrames);
                frameidx= clip(1)+k-1;
                if scale==1
                    s2.xn(:,:,k)= double(rgb2gray(vid.frames(frameidx).cdata))/normfac;
                else
                    s2.xn(:,:,k)=imresize(double(rgb2gray(vid.frames(frameidx).cdata))/normfac,...
                        [vidHeight,vidWidth]);
                end
                
                if sig>0
                    s2.xn(:,:,k)=s2.xn(:,:,k)+sig*randn([vidHeight,vidWidth]);
                end
            end
            close(wait_handle );
            drawnow;
            [~,s2.signalname,~] = fileparts(fn);
            s2.colormap_active='gray';
            
        end
        
        function s2=make_fromVideo(fn,clip,scale,sig)
            % s2=~(fn,[clip,scale,sig]) create from file fn
            % sig ... std of additive Gaussian noise to be added if sig>0.
            % scale ... rescale factor for video frames
            % clip ... #frames or pair ([start frame, #frames]).
            
            assert(ischar(fn) && exist(fn,'file')~=0,'file exists');
            assert(nargin<2 || isnumeric(clip),'clip is #frames or start frame and #frames');
            if ~exist('sig','var')
                sig=0;
            end
            if ~exist('scale','var')
                scale=1;
            end
            if ~exist('clip','var')
                clip=[1,32];
            end
            if isscalar(clip)
                clip=[1,clip];
            end
            
            % read video file
            try
                vid = VideoReader(fn);
                err_videoreader=false;
            catch err
                err_videoreader=true;
            end
            
            if ~err_videoreader
                s2=Signal3D.make_UsingVideoReader(vid, clip,scale,sig);
            else
                s2=Signal3D.make_UsingMMReader(fn, clip,scale,sig);
            end
            
        end
        
    end
    
    %% graphics
    methods
        
        
        function play_signal(obj)
            % play 3d-signal as movie
            obj.requireOnly(~isempty(obj.xn),'local','non-empty signal');
            sz=obj.size;
            % transform to range [0,1] for play back!
            minval=min(obj.xn(:));
            maxval=max(obj.xn(:));
            implay((reshape(obj.xn,[sz(1),sz(2),1,sz(3)])-minval)/(maxval-minval));
        end
        
        
        function graph_signal(obj, open_new, zinterv)
            % ~([zinterv]) show 3d signal as a sequence of 2d sections in an new window
            sz=obj.size;
            if ~exist('zinterv','var') || isempty(zinterv)
                zinterv=20;
            end
            if length(zinterv)~=2
                % construct symmetric interval around center:
                Lmax=min(zinterv,sz(3));
                rmax=floor(Lmax/2);
                z1=max(1,ceil(sz(3)/2)-rmax);
                if mod(Lmax,2)==0
                    z2=z1+2*rmax-1;
                else
                    z2=z1+2*rmax;
                end
                zinterv=[z1,z2];
            end
            if ~exist('open_new','var')
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            z1=zinterv(1);
            z2=zinterv(2);
            L=max(1,z2-z1+1);
            sn=obj.get_signalname;
            if isempty(sn)
                sn='Signal';
            end
            
            linsize=sqrt(obj.numel/sz(3));
            limitsize=256;
            if linsize<=limitsize
                if ~isempty(obj.repfun)
                    % mdat requires 4 subscripts !
                    mdat=reshape(obj.repfun(obj.xn),[sz(1),sz(2),1,sz(3)]);
                else
                    mdat=reshape(obj.xn,[sz(1),sz(2),1,sz(3)]);
                end
            else
                f=limitsize/linsize;
                szneu=ceil(f.*[sz(1),sz(2)]);
                mdat=zeros(szneu(1),szneu(2),1,L);
                for j=1:L
                    % use method 'nearest' to avoid averaging over NaN
                    % values!
                    z=j+z1-1;
                    if ~isempty(obj.repfun)
                        mdat(:,:,1,j)=imresize(obj.repfun(obj.xn(:,:,z)),szneu,'nearest');
                    else
                        mdat(:,:,1,j)=imresize(obj.xn(:,:,z),szneu,'nearest');
                    end
                end
                z1=1;
                z2=L;
            end
            % transform to range [0,1] for function montage:
            % mdat=(mdat-min(mdat(:)))/(max(mdat(:))-min(mdat(:)));
            % warning('off','images:initSize:adjustingMag');
            clims=obj.get_colorlimits; 
            if any(isnan(clims))
                montage(mdat(:,:,:,z1:z2));
            else
                montage(mdat(:,:,:,z1:z2),'DisplayRange' ,clims);
            end
            colormap(obj.colormap_active);
            if obj.ison_colorbar
                colorbar;
            end
            
            if obj.colormap_freeze
                freezeColors;
                try
                    if obj.ison_colorbar
                        cbfreeze; % sometimes problems
                    end
                catch
                    clf;
                    montage(mdat(:,:,:,z1:z2));
                    colormap(obj.colormap_active);
                    colorbar;
                    freezeColors;
                end
            end
            if ~isempty(obj.repfun)
                cblabel_str=func2str(obj.repfun);
                cblabel(cblabel_str);
            end
            title(sn,'fontsize',12);
        end
        
        function show_signal(obj)
            % ~() show 3d signal as subplots of frames in new window
            prepfigure(obj.fig);
            L=obj.size(3);
            sd= factor_subplots(L);
            suptitle(obj.get_signalname,14);
            clims=obj.get_colorlimits;            
            for j=1:L
                subplot(sd(1),sd(2),j);
                if any(isnan(clims))
                    imagesc(obj.xn(:,:,j));
                else
                    imagesc(obj.xn(:,:,j),clims); % all subplots get same clims
                end
                if obj.ison_colorbar
                    colorbar;
                end
                title(['z= ',num2str(j)],'fontsize',12);
            end
            colormap(obj.colormap_active);
        end
        
        function show_partition(obj,m,n)
            % ~() show 3d signal and its partition as subplots of frames in new window
            obj.require(nargin>=3 && isnumeric(m) && isnumeric(n),'local',...
                'needs 2 natural numbers');
            prepfigure(obj.fig);
            L=obj.size(3);
            sd= factor_subplots(L);
            sn=[obj.get_signalname,' partitioned'];
            
            suptitle(sn,14);
            for j=1:L
                subplot(sd(1),sd(2),j);
                s=Signal2D(obj.xn(:,:,j));
                s.graph_partition(m,n,false);
                title(['z= ',num2str(j)],'fontsize',12);
            end
        end
        
        function graph_signalScatter(obj,open_new,filter_quantile)
            % show signal in the open window as scatter plot of non-zero
            % locations.
            if ~exist('filter_quantile','var') || isempty(filter_quantile)
                nmax=100;
                nvals=nnz(abs(obj.xn));
                vals=sort(nonzeros(abs(obj.xn)),'ascend');
                filter_quantile=max(0.99,(nvals-nmax)/nvals);
                minmodulus=vals(floor(filter_quantile*numel(vals)));
            else
                minmodulus=quantile(abs(obj.xn(:)),filter_quantile);
            end
            if ~exist('open_new','var')
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            
            coords=abs(obj.xn)>=minmodulus;
            [m n o]=size(obj.xn);
            [x,y,z] = meshgrid(1:m,1:n,1:o);
            
            scatter3(x(coords),y(coords),z(coords),[],'blue','filled');
            %scatter3(x(coords),y(coords),z(coords),[],obj.xn(coords),'filled');
            sn=obj.obj.get_signalname;
            if isempty(sn)
                sn='Signal';
            end
            title(sn,'fontsize',12);
            
        end
        
        
        function graph_totalVar(obj,open_new,dim)
            % ~([open_new, dim]) total variation along dimension dim
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('dim','var')
                dim=3;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            [img,l1]=obj.totalVar(dim);
            imagesc(img);
            colormap('default');
            if obj.ison_colorbar
                colorbar;
            end
            sn=obj.get_signalname;
            title({['TV_',num2str(dim),' of ',sn,' ',vec2str(obj.N)],...
                ['||.||_1=',num2str(l1,'%3.1e')]},'fontsize',12);
            
        end
        
        
        
        function graph_motionTV(obj,open_new, shiftsize)
            % ~([open_new]) temporal signal of total variation between
            % consecutive frames
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('shiftsize','var') || isempty(shiftsize)
                shiftsize=1;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            yn=obj.totalVarShift(shiftsize);
            
            s2=TimeSignal(yn);
            s2.amp_quantity='rel.TV';
            s2.marker='x';
            s2.signalname=['inter-frame rel. total variation, ',obj.signalname];
            s2.graph_signal(false);
            
        end
        
        function graph_diffshift(obj,open_new, shiftsize)
            % ~([open_new, dim]) difference between consecutive frames
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('shiftsize','var') || isempty(shiftsize)
                shiftsize=1;
            end
            shiftsize=abs(shiftsize);
            if open_new
                prepfigure(obj.fig);
            end
            
            yn=zeros(obj.N-[0,0,shiftsize]);
            L=obj.size(3);
            
            for z1=1:L-shiftsize
                yn(:,:,z1)=obj.xn(:,:,z1+shiftsize)-obj.xn(:,:,z1);
            end
            s2=Signal3D(yn);
            s2.signalname=['differences (shiftsize=',num2str(shiftsize),...
                ') ',obj.signalname];
            s2.graph_signal(false);
            
        end
        
        
        function graph_PSNRcross(obj,sig2, open_new)
            % ~([open_new]) cross PSNR of each frame pair
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('sig2','var') || isempty(sig2)
                sig2=obj;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            imagescnan(obj.PSNRcross(sig2));
            L=obj.size(3);
            if L<10
                set(gca,'XTick',1:L);
                set(gca,'XTickLabel',1:L);
                set(gca,'YTick',1:L);
                set(gca,'YTickLabel',1:L);
            end
            
            colormap('default');
            if obj.ison_colorbar
                colorbar;
            end
            sn=[obj.get_signalname, ' vs. ',sig2.get_signalname];
            title({['cross similarity PSNR of ',vec2str(obj.N)],sn},'fontsize',12);
        end
        
        function graph_SSIMcross(obj,sig2, open_new)
            % ~([open_new]) cross SSIM of each frame pair
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('sig2','var') || isempty(sig2)
                sig2=obj;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            imagescnan(obj.SSIMcross(sig2));
            L=obj.size(3);
            if L<10
                set(gca,'XTick',1:L);
                set(gca,'XTickLabel',1:L);
                set(gca,'YTick',1:L);
                set(gca,'YTickLabel',1:L);
            end
            colormap('default');
            if obj.ison_colorbar
                colorbar;
            end
            sn=[obj.get_signalname, ' vs. ',sig2.get_signalname];
            title({['cross similarity SSIM of ',vec2str(obj.N)],sn},'fontsize',12);
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 3d';
            ok= isempty(obj.xn) || ismember(length(size(obj.xn)),[2,3]);
        end
    end
    
end


