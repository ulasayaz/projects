classdef Signal_Factory <DC
    % make 3d or 2d signals and motion vector fields
    %
    % Example
    % ==========
    %{
    %------ ex1: create moving small image in larger one:
    s2=Signal2D.make_fromImage('cameraman.bmp');
    s2.resize(64);
    M=512; N=512; L=3; ulc=round([M,N]/4); shiftsize=25;
    sbig=Signal2D(zeros(M,N));
    [s3,motion]= Signal_Factory.make_VideofromSmallImage(sbig,s2,L,ulc,shiftsize);
    motion.show_field_TV;
    %
    %------ ex2: construct 3d signal from 2d signal by shifts:
    signal2=Signal2D.make_fromImage('cameraman.bmp');
    signal2.resize(128);
    signal3= Signal_Factory.make_CyclefromSignal2D(signal2,1,32);
    %
    % --- padd 3d signal to dyadic cube:
    % signal3.padd2dyadic; display(signal3.size);
    %  
    % ------ ex3:  camera panover big images: --------------  
    signalBig2=Signal2D.make_fromImage('parrots.tiff');
    L=8; pansize=256*[1,1]; origin=[700,300]; motionunit=[3,5];
    signal3= Signal_Factory.make_PanfromSignal2D(signalBig2,L,pansize,origin,motionunit);
    %       -- estimate motion
    vf=MotionVF.make_fromSignal3(signal3);
    vf.fig=5;
    %       -- all vector fields: 
                                  vf.show_field_TV();
    %       -- or only the first: vf.show_field_TV(1);
    %       -- show norm of all vector fields: 
                                   vf.show_norm;
    vf.show_motionCompensation();
    vf.test_motionCompensation(1);
    vf.test_motionCompensation();
    %---------------------------------------------------------
    %}
    
    properties
    end
    
    
    methods (Static)
        
        function [s3, motion]=make_VideofromSmallImage(sbig,s2,L,ulc,shiftsize)
            % [s3, motion]=~(sbig,s2,[L,ulc,shift]) move s2 over empty canvas
            requireOnly(isa(sbig,'Signal2D'),'local','needs 2D signal');
            requireOnly(isa(s2,'Signal2D'),'local','needs 2D signal');
            requireOnly(all(sbig.N>s2.N),'local','size requirements');
            requireOnly(nargin<5 ||isempty(ulc) || (isvector(ulc) && length(ulc)==2),...
                'local','upper left corner is pair of integers');
            if ~exist('L','var')
                L=3;
            end
            if ~exist('ulc','var')
                ulc=round(sbig.size/4);
            end
            if ~exist('shiftsize','var')
                shiftsize=round(s2.size/2);
            end
            if isscalar(shiftsize)
                shiftsize=shiftsize*[1,1];
            end
            
            s3=Signal3D();
            xn=zeros(sbig.size(1),sbig.size(2),L);
            N2=s2.N;
            N=sbig.N;
            motion2d=zeros(N(1),N(2),L-1,2);  % motion vector field
            cc=ulc;
            for k=1:L
                xn(:,:,k)=sbig.xn;
                xn(cc(1):cc(1)+N2(1)-1,cc(2):cc(2)+N2(2)-1,k)=s2.xn;
                if k<L
                    % define motion vector field
                    motion2d(cc(1):min(N(1),cc(1)+N2(1)-1),cc(2):min(N(2),cc(2)+N2(2)-1),...
                        k,1)=shiftsize(1);
                    motion2d(cc(1):min(N(1),cc(1)+N2(1)-1),cc(2):min(N(2),cc(2)+N2(2)-1),...
                        k,2)=shiftsize(2);
                end
                cc=cc+shiftsize;
            end
            s3.xn=xn;
            s3.signalname='moving small image';
            
            motion=MotionVF(motion2d);
            motion.set_signal(s3);
            motion.method='sim';
            
            assert(isa(motion,'MotionVF'),'returns a motion vector field');
            assert(isa(s3,'Signal3D'),'returns a 3d signal');
            assert(isequal(motion.size,[N(1),N(2),L-1,2]),'size matches input');
            
        end
        
        
        function [s3, motion]=make_MovingRectangles(countMovingRectangles,M,N,L,params)
            % [s2,vectorfield]= ~(countRectangles,M,N,L) signal size MxNxL
            % s2 is a 3d signal: s3.show_signal; s3.play_signal;
            % motion is a motion vector field: motion.show_field;
            requireOnly(countMovingRectangles>=1,'local', 'positive number');
            if ~exist('M','var') || ~exist('N','var')
                M=128;
                N=128;
            end
            if ~exist('L','var')
                L=3;
            end
            if ~exist('params','var')
                params=struct;
            end
            if ~isfield(params,'sig') || isempty(params.sig)
                params.sig=0; % noise
            end
            if ~isfield(params,'shiftsize')
                params.shiftsize=[];
            end
            if ~isfield(params,'w')
                params.w=[];  % interval width
            end
            if ~isfield(params,'ulc')
                params.ulc=[];  % interval width
            end
            
            if isscalar(params.shiftsize)
                params.shiftsize=params.shiftsize*[1,1];
            end
            
            R=uint32(countMovingRectangles);
            s3=Signal3D(zeros(M,N,L));
            s3.signalname=[num2str(R),' moving rectangles'];
            
            motion2d=zeros(M,N,L-1,2);
            
            maxshift=10;
            maxwidth=max(round(min(M,N)/8),10);
            minwidth=max(round(min(M,N)/32),5);
            
            
            for j=1:R
                
                % define movable rectangle
                ok=false;
                it=0;
                while ~ok
                    it=it+1;
                    if isempty(params.shiftsize)
                        shiftsizeA=randi(2*maxshift,[1,2])-maxshift*[1,1];
                    else
                        shiftsizeA=params.shiftsize;
                    end
                    if isempty(params.w)
                        w=randi(maxwidth-minwidth,[1,2])+minwidth; % moving rectangle size
                    else
                        w=params.w;
                        if isscalar(w)
                            w=w*ones(1,2);
                        end
                    end
                    ulc=params.ulc;
                    if isempty(params.ulc) || j>1
                        ulc=[randi(M),randi(N)];  % upper left corner of moving rectangle
                    end
                    % movability conditions
                    ulc2=ulc+(L-1)*shiftsizeA;
                    ok= it>1000 || (ulc2(1)>=1 && ulc2(1)<=M && ulc2(2)>=1 && ulc2(2)<=N);
                end
                
                for k=1:L   % simulate motion of rectangle
                    s2=Signal2D.make_rect(w,[M,N],ulc+(k-1)*shiftsizeA,params.sig);
                    s3.xn(:,:,k)=min(1,s3.xn(:,:,k)+s2.xn);
                    if k<L
                        % define motion vector field
                        ulcK=ulc+(k-1)*shiftsizeA;
                        motion2d(ulcK(1):min(M,ulcK(1)+w(1)-1),ulcK(2):min(N,ulcK(2)+w(2)-1),...
                            k,1)=shiftsizeA(1);
                        motion2d(ulcK(1):min(M,ulcK(1)+w(1)-1),ulcK(2):min(N,ulcK(2)+w(2)-1),...
                            k,2)=shiftsizeA(2);
                    end
                end
                
            end
            
            motion=MotionVF(motion2d);
            motion.set_signal(s3);
            motion.method_estim ='sim';
            
            assert(isa(motion,'MotionVF'),'returns a motion vector field');
            assert(isa(s3,'Signal3D'),'returns a 3d signal');
            assert(motion.sizes_match,'sizes of signal and motion field match');
            
        end
        
        function [s2, motion]=make_MovingIntervals(countMovingIntervals,hN,hL,params)
            % [s2,vectorfield]= ~(countIntervals,N,L,params) 2d signal size NxL
            % s2 is a 2d signal: s2.show_signal; s2.graph_signal;
            % motion is a 1d-motion vector field: motion.show_field;
            requireOnly(countMovingIntervals>=1,'local', 'positive number');
            if ~exist('hN','var') || ~exist('hN','var')
                hN=128;
            end
            if ~exist('hL','var')
                hL=32;
            end
            if ~exist('params','var')
                params=struct;
            end
            if ~isfield(params,'sig') || isempty(params.sig)
                params.sig=0; % noise
            end
            if ~isfield(params,'shiftsize')
                params.shiftsize=[];
            end
            if ~isfield(params,'w')
                params.w=[];  % interval width
            end
            
            R=uint32(countMovingIntervals);
            s2=Signal2D(zeros(hN,hL));
            s2.signalname=[num2str(R),' moving intervals'];
            
            motion1d=zeros(hL-1,hN);
            
            maxshift=10;
            maxwidth=max(round(hN/8),10);
            minwidth=max(round(hN/32),5);
            
            
            for j=1:R
                
                if isempty(params.shiftsize)
                    shiftsizeA=randi(2*maxshift)-maxshift;
                else
                    shiftsizeA=params.shiftsize;
                end
                if isempty(params.w)
                    w=randi(maxwidth-minwidth)+minwidth; % moving interval size
                else
                    w=params.w;
                end
                
                ulc=randi(hN);  % left corner of moving interval
                
                
                
                for k=1:hL   % simulate motion of interval
                    s1=TimeSignal.make_interv(w,hN,ulc+(k-1)*shiftsizeA,params.sig);
                    s2.xn(:,k)=min(1,s2.xn(:,k)+s1.xn(:));
                    if k<hL
                        % define motion vector field
                        ulcK=ulc+(k-1)*shiftsizeA;
                        if ulcK>=1 && ulcK<=hN
                            motion1d(k,ulcK:min(hN,ulcK+w(1)-1))=shiftsizeA(1);
                        end
                    end
                end
                
            end
            
            s2.transpose;
            motion=MotionVF1D(motion1d);
            motion.set_signal(s2);
            
            assert(isa(motion,'MotionVF1D'),'ensure violated: returns a motion vector field');
            assert(isa(s2,'Signal2D'),'ensure violated: returns a 3d signal');
            assert(isequal(motion.size,[hL-1,hN]),'ensure violated: size matches input');
            
        end
        
        function s2=make_CyclefromSignal1D(s1,L,motionunit)
            % s2=~(fn,[L,motionunit]) build 2D signal from 1D signal by moving image
            % periodically by multiples of motionunit 0 to L-1 times
            assert(isa(s1,'TimeSignal'),'needs 1D signal');
            assert(nargin<3 || isempty(motionunit)||isscalar(motionunit),...
                'unit is scalar');
            
            if ~exist('L','var')
                L=32;
            end
            if ~exist('motionunit','var') || isempty(motionunit)
                motionunit=11;
            end
            
            xdata=s1.xn;  % signal content
            
            data2d=zeros([length(xdata),L]);
            % build 2d image as sequence of shifted 1d signals
            
            data2d(:,1)=xdata;
            pix=Signal2D.shiftdist(L,motionunit);
            if isequal(pix, floor(pix))
                % all shifts by full pixels
                shifter=@circshift;
            else
                % subpixel shifts
                shifter=@TimeSignal.SubPixelShifter;
            end
            for j=2:L
                data2d(:,j)=shifter(xdata,pix(j-1));
            end
            
            
            s2=Signal2D(data2d);
            s2.signalname=[s1.signalname,' - shift:',num2str(motionunit,'%3.1f')];
            s2.colormap_active='gray';
            
        end
        
        function s3=make_PanfromSignal2D(s2,L,pansize,origin,motionunit)
            %s3=~(s2[,L,pansize,origin,motionunit]) camera pan over big image s2
            % difference to make_CyclefromLowerDimSig: instead of periodidicity
            % new image parts appear and others disappear during pan.
            requireOnly(isa(s2,'Signal2D'),'local','needs 2D signal');
            
            if ~exist('L','var') || isempty(L)
                L=8;
            end
            if ~exist('motionunit','var') || isempty(motionunit)
                motionunit=[11,17];
            end
            if ~exist('origin','var') || isempty(origin)
                origin=[1,1];
            end
            if ~exist('pansize','var') || isempty(pansize)
                pansize=[256,256];
            end
            
            xdata=s2.xn;  % signal content
            
            data3d=zeros([pansize,L]);
            pix=Signal_Factory.shiftdist(L,motionunit);
            
            % build 3d video as sequence of moved images
            % video sequence of shifts
            data3d(:,:,1)=xdata(origin(1):pansize(1)+origin(1)-1,...
                origin(2):pansize(2)+origin(2)-1);
            
            if isequal(pix, floor(pix))
                % all shifts by full pixels
                shifter=@circshift;
            else
                % subpixel shifts
                shifter=@Signal2D.SubPixelShifter;
            end
            for j=2:L
                % shift big image in the direction of the shift -> velocity
                % field gets same direction as shift
                h=shifter(xdata,pix(:,j-1));
                data3d(:,:,j)=h(origin(1):pansize(1)+origin(1)-1,...
                    origin(2):pansize(2)+origin(2)-1);
            end
            
            
            s3=Signal3D(data3d);
            s3.video_sim=2;
            motionType='shifted';
            s3.signalname=[s2.signalname,' - ',motionType,': ',vec2str(motionunit,'%3.1f')];
            s3.colormap_active='gray';
            
        end
        
        function s3=make_CyclefromSignal2D(s2,videoflag,L,motionunit,scale)
            % s3=~(s2,[videoflag,modesig,scale]) build 3d signal (pseudo-video sequence)
            % from 2d signal by by rotating or shifting periodically image by multiples
            % of motionunit;
            % rescale image if scale missing or ~=1
            requireOnly(isa(s2,'Signal2D'),'local','needs 2D signal');
            requireOnly(nargin<3 || isempty(videoflag) || ismember(videoflag,[1,2]),...
                'local','videoflag admitted values');
            
            if ~exist('scale','var') || isempty(scale)
                scale=1;
            end
            if ~exist('L','var') || isempty(L)
                L=8;
            end
            if ~exist('motionunit','var') || isempty(motionunit)
                motionunit=[11,17];
            end
            if ~exist('videoflag','var') || isempty(videoflag)
                videoflag=1;
            end
            
            xdata=s2.xn;  % signal content
            
            if scale~=1
                xdata=imresize(xdata,scale);
            end
            
            data3d=zeros([size(xdata),L]);
            % build 3d video as sequence of moved images
            if videoflag==1
                % video sequence of rotations
                % prepare image for rotation by cropping it to a disk:
                motionType='rotated';
                pixorigin=ceil(size(xdata)/2);
                [X,Y]=meshgrid(1:size(xdata,2),1:size(xdata,2));
                rads=sqrt((X-pixorigin(2)).^2+(Y-pixorigin(2)).^2);
                xdata(rads>min(size(xdata))/2)=0;
                
                data3d(:,:,1)=xdata;
                alphas=Signal_Factory.rotangles(L,motionunit);
                for j=2:L
                    data3d(:,:,j)=imrotate(xdata,alphas(j-1),'crop');
                end
            else
                % video sequence of shifts
                data3d(:,:,1)=xdata;
                pix=Signal_Factory.shiftdist(L,motionunit);
                motionType='shifted';
                if isequal(pix, floor(pix))
                    % all shifts by full pixels
                    shifter=@circshift;
                else
                    % subpixel shifts
                    shifter=@Signal2D.SubPixelShifter;
                end
                for j=2:L
                    data3d(:,:,j)=shifter(xdata,pix(:,j-1));
                end
            end
            
            s3=Signal3D(data3d);
            s3.video_sim=videoflag;
            s3.signalname=[s2.signalname,' - ',motionType,': ',vec2str(motionunit,'%3.1f')];
            s3.colormap_active='gray';
            
        end
        
        
    end
    
    
    methods (Static)
        
        function alphas= rotangles(L, rotunit)
            % rotation angles for pseudo-video sequence
            if ~exist('L','var') || isempty(L)
                L=8;
            end
            if ~exist('rotunit','var') || isempty(rotunit)
                rotunit=90/L;
            end
            j=2:L;
            alphas=rotunit(1)*(j-1);
        end
        
        function pix= shiftdist(L, unit)
            % pix=~(L,unit) multiples of shift unit
            if ~exist('L','var') || isempty(L)
                L=8;
            end
            if ~exist('unit','var') || isempty(unit)
                unit=[11;17];
            end
            if isscalar(unit)
                unit=unit*[1;1];
            end
            j=1:L-1;
            pix=kron(j,unit(:));
        end
        
    end
    
end

