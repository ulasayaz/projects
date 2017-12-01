classdef CSPMM <CompressiveSensing
    %% Class to implement PMM + NTVICT based recovery
%             signal=Signal2D.make_fromImage('cameraman.bmp');
%             Phi=Curvelet2_clab(signal);
%             csp=CSPMM();
%             csp.blocksize = 16;%16
%             N=256; csp.set_original(signal,N);
%             c=16; csp.set_compression(c);
%             ct=csp.s1w_Donoho;
%             csp.set_Phi(Phi);
%             f = double(imread('cameraman.bmp')); % use same i/p as 
%             csp.xorig.xn = f/max(f(:)); % NTVICT script for bit true
%             % Select one of the following measurement schemes (1 out of 4)
%             % csp.Mtype='Gaussian';           % Random Gaussian PMM
%             csp.Mtype='Hadamard';             % Hadamard Ensemble PMM
%             % csp.Mtype='sparse1';            % Random Binary PTM
%             % csp.Mtype='sparseplusminus1';   % Rademacher Matrix PTM
%             SNR=[]; csp.PMM(); % PMM based Encoder
%             csp.SolveNTVICT; % NTVICT based Solver
%      
    properties
        blocksize %@blocksize of measurement matrix
        Mtype     %@measurement matrix type
    end
    
    
    %% constructors and commands
    methods
        function obj=CSPMM()
            % constructor obj=~(blocksize) setting PMM block size
            obj = obj@CompressiveSensing();
                                            
        end
        function set_blocksize(obj)
            % ~(blocksize) setting PMM block size
            set_blocksize@CSSparseSignal(obj);
            obj.blocksize=[];%obj.A.nodes;
            %obj.blocksize=[];
        end
        
        function set_Projection(obj, block_size,subrate)
            %   This function generates the random projection matrix
            %   Phi for the given block size and subrate.
            %   Phi is returned as a M x N matrix, where N = block_size *
            %   block_size, and M = round(subrate * N)
            %% PMM Measurements Matrices
            if strcmp(obj.Mtype,'Gaussian'),
                N = block_size * block_size;
                M = round(subrate * N);
                Phi = orth(randn(N, N))';
                Phi = Phi(1:M, :);
                obj.A = Phi;
                obj.Atype = 'Gaussian';
            end
            if strcmp(obj.Mtype,'Hadamard'), 
                    %                 if subrate == 0.25 && block_size == 16,
                    %                     load('PhiSBHE.mat'); % 64x256 hadamard meas. matrix
                    %                     Phi = Phi/4;
                    %                     obj.A = Phi;
                    %                 elseif subrate == 0.125 && block_size == 16,
                    %                     load('PhiSBHE8.mat'); % 32x256 hadamard meas. matrix
                    %                     Phi = Phi/4; % check division factor by PhixPhi'???
                    %                     obj.A = Phi;
                    %                 elseif subrate == 0.0625 && block_size == 16,
                    %                     load('PhiSBHE16.mat'); % 16x256 hadamard meas. matrix
                    %                     Phi = Phi/4; % check division factor by PhixPhi'???
                    %                     obj.A = Phi;
                    %                 else
                    hm=sbhe(block_size,block_size,subrate*block_size*block_size); % problem with eval function
                    Phi = hm/sqrt(block_size);
                    obj.A = Phi;
                    obj.Atype = 'Hadamard';
                    %                 end
            end
            %% PTM Measurements Matrices
            if strcmp(obj.Mtype, 'sparse1') || strcmp(obj.Mtype, 'sparseplusminus1'),%            
                 % Random Binary or Rademacher Matrices for PTM type sampling
                 mbin1 = gen_matrix(subrate*block_size*block_size, block_size*block_size, obj.Mtype);
                 Phi = mbin1.A';
                 obj.A = Phi;
                 obj.Atype = 'RandomBinary';
            end
            
        end
        
        function sim_BCS_SPL_Encoder(obj)
            %   This function performs BCS projections of each block of
            %   current_image. The number of columns of the projection matrix,
            %   Phi, determines the size of the blocks into which current_image
            %   is partitioned.
            current_image = obj.xorig.xn;
            Phi = obj.A;
            N = size(Phi, 2);
            block_size = obj.blocksize;%sqrt(N);
            %[num_rows, num_cols] = size(current_image);
            x = im2col(current_image, [block_size block_size], 'distinct');
            obj.y = Phi * x;
                    
        end
        
        function PMM(obj)
            %% PMM Encoder
            % Simulate a measurement vector y by applying a Gaussian or Hadamard
            % matrix with measurement compression obj.c and encoding is
            % done in block by block fashion i.e. pixel muliplexing method
            % of CS sampling                
            
            %% Spatial CS Block based Encoding
            block_size = obj.blocksize; 
            subrate = 1.0/double(obj.c);%0.25;
            %% Measurement Matrix
            obj.set_Projection(block_size,subrate);
            
            %% CS PMM/PTM Encoder
            obj.sim_BCS_SPL_Encoder;%(f, Phi);%y = BCS_SPL_Encoder(f, Phi);
                                               
        end
        function [f,Ct,Midx]=curveletTh_m(obj,sigma,W, Mmu)
                %Curvelet thresholding
                img = obj.x;
                n = size(img,1);
                F = ones(n);
                X = fftshift(ifft2(F)) * sqrt(prod(size(F)));% try also w/ fftshift b/f ad also a/f fft2
                %C = fdct_wrapping(X,1,2);
                C = fdct_wrapping_m(X,0,1); % C = fdct_wrapping(u0,0,finest,nbscales,nbangles_coarse);
                if nargin < 4, % added by myself
                    Mmu = cell(size(C));
                    for s=1:length(C)
                      Mmu{s} = cell(size(C{s}));
                      for w=1:length(C{s})
                        %A = C{s}{w};
                        Mmu{s}{w} = 1;%sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
                      end
                    end
                end
                if nargin < 3, % added by myself
                    W = 1;
                end
                

                % Compute norm of curvelets (exact)
                E = cell(size(C));
                for s=1:length(C)
                  E{s} = cell(size(C{s}));
                  for w=1:length(C{s})
                    A = C{s}{w};
                    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
                  end
                end

                % Take curvelet transform
                %C = fdct_wrapping(img,1,2);%1 for curvelet, 2 for wavelet
                 C = fdct_wrapping_m(W.*img,0,1); 

                % Apply thresholding 
                Ct = C;
                %CtDiff = C;
                Midx = C;
                %%%% remove the form error %%%%% 
                
                for s =1:length(C)
                  thresh = 3*sigma + sigma*(s == length(C));%ww=length(C{s})
                  for w = 1:length(C{s})
               
                    Ct{s}{w} = Mmu{s}{w} .* C{s}{w}.* (abs(C{s}{w}) >thresh*E{s}{w});%0.5**0.85 % Mmu added by myself
                    
                    Midx{s}{w} = (Ct{s}{w} ~= 0);
                  end
                end
                % Inverse curvelet transform
                
                f= (1./W).*real(ifdct_wrapping_m(Ct,0,size(img,1),size(img,2)));
        end
%         function f = curveletTh_m(obj,sigma)
%             %Curvelet thresholding
%             img = obj.x;
%             n = size(img,1);
%             F = ones(n);
%             X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
%             %C = fdct_wrapping(X,1,2);
%             C = fdct_wrapping_m(X,0,1); % C = fdct_wrapping(u0,0,finest,nbscales,nbangles_coarse);
%             %   is_real     Type of the transform
%             %                   0: complex-valued curvelets
%             %                   1: real-valued curvelets
%             %               [default set to 0]
%             %   finest      Chooses one of two possibilities for the coefficients at the
%             %               finest level:
%             %                   1: curvelets
%             %                   2: wavelets
% 
%             % Compute norm of curvelets (exact)
%             E = cell(size(C));
%             for s=1:length(C)
%               E{s} = cell(size(C{s}));
%               for w=1:length(C{s})
%                 A = C{s}{w};
%                 E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
%               end
%             end
% 
%             % Take curvelet transform
%             %C = fdct_wrapping(img,1,2);%1 for curvelet, 2 for wavelet
%              C = fdct_wrapping_m(img,0,1); 
% 
%             % Apply thresholding 
%             Ct = C;
% 
%             %%%% remove the form error %%%%% 
%             %for w=1:length(C{1})
%             %    Ct{1}{w}=0;
%             %end
%             %%%%%
%             %length_c=length(C);
%             for s =1:length(C)
%               thresh = 3*sigma + sigma*(s == length(C));%ww=length(C{s})
%               for w = 1:length(C{s})
%                 Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) >thresh*E{s}{w});%0.5**0.85
%               end
%             end
%             % Inverse curvelet transform
%             % f= real(ifdct_wrapping(Ct,1));
%             f= real(ifdct_wrapping_m(Ct,0,size(img,1),size(img,2))); 
%         end
        
        function [f,Ct]=curveletUpdate_m(obj,Cu,Mmu,W) %curveletTh(img,sigma,Mmu)
            %Curvelet thresholding
            Uv = obj.x;
            img = Uv;
            n = size(img,1);
            F = ones(n);
            X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
            %C = fdct_wrapping(X,1,2);
            C = fdct_wrapping_m(X,0,1); % C = fdct_wrapping(u0,0,finest,nbscales,nbangles_coarse);

            % Compute norm of curvelets (exact)
            E = cell(size(C));
            for s=1:length(C)
              E{s} = cell(size(C{s}));
              for w=1:length(C{s})
                A = C{s}{w};
                E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
              end
            end

            % Take curvelet transform
            %C = fdct_wrapping(img,1,2);%1 for curvelet, 2 for wavelet
            Cv = fdct_wrapping_m(W.*Uv,0,1); 

            % Apply thresholding 
            Ct = C;

            %%%% remove the form error %%%%% 
            
            for s =1:length(C)
              %thresh = 3*sigma + sigma*(s == length(C));%ww=length(C{s})
              for w = 1:length(C{s})
            %       Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > sigma);% as mentioned in paper 2 of Ma et al
                Ct{s}{w} = Mmu{s}{w} .* Cu{s}{w} + (1-Mmu{s}{w}).* Cv{s}{w}; %.* (abs(C{s}{w}) >thresh*E{s}{w});%0.5**0.85 % Mmu added by myself
                
              end
            end
            % Inverse curvelet transform
            f = (1./W).*real(ifdct_wrapping_m(Ct,0,size(img,1),size(img,2))); 
        end
    end
    
    %% solvers   
    methods     
         
        function SolveNTVICT(obj) 
            %% NTVICT
            % NTVICT based Solver making use of Curvelet as sparsifying
            % tranform and explained in the papers of Prof. J. Ma et al.
            % Modify accordingly based on obj class
            % Maybe also consider writing curveletTH as a separate class
            % based function
            N1 = size(obj.xorig.xn,1);
            N2 = size(obj.xorig.xn,2);
            block_size = obj.blocksize;
            Phi = obj.A;
            y = obj.y;
            thre=0.3;%0.09;%0.4
            iteration=960;%60;%30;
            alf=1.2;bait=2*alf/(1+0.1);

            u0 = Phi' * y;
            u0im = col2im(u0,[block_size,block_size],[N1,N2],'distinct'); % tranforming coordinates to without block based
            u0im_ = u0im; % for TwICT

            applyWindowing = 0; % apply hamming type raised cosine window to reduce artifacts in Curvelet based Recovery
            if applyWindowing == 1,
                L=N1; w = hamming(L); W = (w*w'); %figure;imagesc((1./(w*w')).*((ifft2(fft2((w*w').*f)))))
            else
                W = 1;
            end

            for ii=1:iteration
                disp(ii);
                if ii == 1,
                    obj.x = u0im;
                    [u1im,Ctim, Midx]=obj.curveletTh_m(thre*(1-1/iteration)^(1/1), W); % Midx added chk if correct
                    u2im = u1im;
            %         u1im_ = u1im; % for TwICT
                end
                if ii > 1,
                    obj.x = u1im + col2im(Phi'*(y - Phi*im2col(u1im,[block_size,block_size],'distinct')),[block_size,block_size],[N1,N2],'distinct');
                    [u2im,Ctim, Midx] = obj.curveletTh_m(thre*(1-ii/iteration)^(1/1),W);% ,Midx at the end added chk if correct
                    %uim = (1-alf)*u0im + (alf-bait)*u1im + bait*u2im;
                    u0im = u1im;


                end

                %% introduce here now the projected NTV part acc. to paper 2 of prof. Ma
                % setting parameters acc. to paper 2 para III.A
                %lambda = 25;
                %beta = 2;
                %mu = 1;
                %d0 = 0;
                %b0 = 0;

                opts=[];
                if ~isfield(opts,'mu'),     opts.mu=10; end %% regularization term scale
                if ~isfield(opts,'delta'),  opts.delta=1; end %% delta<||A^*A||, since A here is a subsampled Fourier matrix, could be fixed to be 1.
                if ~isfield(opts,'nOuter'), opts.nOuter=200; end %% Outer Bregman iteration steps.
                if ~isfield(opts,'nDenoising'), opts.nDenoising=5; end%% 3st level: denoising/regularization level, in general could be fixed to be 10
                if ~isfield(opts,'type'),  opts.type=1; end %% 1: by bos,2:PBOS. 3:by Linearized Bregman If the without noise (or low), choose 1
                if ~isfield(opts,'bTol'), opts.bTol=10^-5; end %%stopping criterion on residual std2(Fmask.*fft2(u)/N-Fp), if the noise standard variation is known, can set btol to be sigma
                if ~isfield(opts,'xTol'), opts.xTol=10^-5; end %%stopping criterion, if the noise standard variation is known, can set btol to be sigma
                if ~isfield(opts,'verbose'), opts.verbose=0; end %% display message
                if ~isfield(opts,'h0'),  opts.h0=16; end %% weight filter parater, depends on noise and image standard variation. for example: for barbara [0, 255], h0=20; To be adapted for normalized image
                if ~isfield(opts,'nWin'), opts.nWin=2;  end    %% patch size [2*nwin+1, 2*nwin+1]
                if ~isfield(opts,'nBloc'), opts.nBloc=5; end   %% search window size [2*bloc+1, 2*bloc+1]
                if ~isfield(opts,'nNeigh'), opts.nNeigh=10; end   %% number of best neighbors (real neighbors size: 2*nNeigh+4)
                if ~isfield(opts,'nWeightupdate'), opts.nWeightupdate=20;  end %0: no weight update, otherwise update steps 
                if ~isfield(opts,'denoising_type'), opts.denoising_type=1; end %% NLTV denoising algorithm:1: Split bregman, 2: Projection in dual


                %% Projected NTV
                v1 = u2im; % init with TwICT o/p
                tau=1./(opts.delta*opts.mu);
                n=0;
                condition=1; % if condition '0' means TwICT elseif '1' means NTVICT
                %% OUTPUT M INDEX OF SIGNIFIACNT COMPONENTS  IN TWICT STAGE; THAT HAS TO BE EXPERIMENTED
                %% Main loop
                while (condition)  
                  if (opts.nWeightupdate) %% update weight
                      if(mod(n,opts.nWeightupdate)==0)
                        wopts=update_weight(real(v1),opts.h0,opts.nWin,opts.nBloc,opts.nNeigh);          
                      end
                  end

                %u2Pv = denoising_SBNLTV(v1,mu,beta./lambda,opts.nDenoising,wopts); % shud be acc. to paper 2 now
                u2Pv = denoising_SBNLTV(v1,tau,tau./4,opts.nDenoising,wopts); % is it acc. to step 2 of paper 2, pls. chk.
                 n = n+1; % increment for loop
                 condition=(n<1);%opts.nOuter&& energy(n)>opts.bTol); is this K???
                end
                obj.x = u2Pv;
                [u1im,Ct]=obj.curveletUpdate_m(Ctim,Midx,W); % updating significant coef with modified insignifiact coef. acc to paper 2, chk if ok
                obj.x = u1im; % storing intermediate/final recovered image
                %u1im = u2Pv; % step 6 of paper 2

                % writing every o/p of NTVICT algorithm
                %psnr1 = PSNR(double(uint8(f*256)), double(uint8(u1im*256)));
                %SSIMR = ssim(double(uint8(f*256)), double(uint8(u1im*256)));
            end
            %% Commented Out to Introduce Latest NTVICT Code
%             thre=0.6;%0.09;%0.4
%             iteration=960;%60;%30;
%             alf=1.2;bait=2*alf/(1+0.1);
%             Phi = obj.A;
%             y = obj.y;
%             % CHANGES DONE TILL HERE -> CONTINUE FROM HERE
%             u0 = Phi' * y;
%             u0im = col2im(u0,[block_size,block_size],[N1,N2],'distinct'); % tranforming coordinates to without block based
%             for ii=1:iteration
%                 disp(ii);
%                 if ii == 1,
%                     obj.x = u0im;
%                     u1im = obj.curveletTh_m(thre*(1-1/iteration)^(1/1)); % Midx added chk if correct
%                     u2im = u1im;
%             %         u1im_ = u1im; % for TwICT
%                 end
%                 if ii > 1,
%                     obj.x = u1im + col2im(Phi'*(y - Phi*im2col(u1im,[block_size,block_size],'distinct')),[block_size,block_size],[N1,N2],'distinct');
%                     u2im = obj.curveletTh_m(thre*(1-ii/iteration)^(1/1));% ,Midx at the end added chk if correct
%                     % uim = (1-alf)*u0im + (alf-bait)*u1im + bait*u2im;
%                     u0im = u1im;
%             
%                 end
%                 % introducing here now the projected NTV part acc. to paper 2 of prof. Ma
%                 % setting parameters acc. to paper 2 para III.A
%                 %lambda = 200;
%                 %Beta = 2;
%                 %d0 = 0;
%                 %b0 = 0;
%                 opts=[];
%                 if ~isfield(opts,'mu'),     opts.mu=10; end %% regularization term scale
%                 if ~isfield(opts,'delta'),  opts.delta=1; end %% delta<||A^*A||, since A here is a subsampled Fourier matrix, could be fixed to be 1.
%                 if ~isfield(opts,'nOuter'), opts.nOuter=200; end %% Outer Bregman iteration steps.
%                 if ~isfield(opts,'nDenoising'), opts.nDenoising=5; end%% 3st level: denoising/regularization level, in general could be fixed to be 10
%                 if ~isfield(opts,'type'),  opts.type=1; end %% 1: by bos,2:PBOS. 3:by Linearized Bregman If the without noise (or low), choose 1
%                 if ~isfield(opts,'bTol'), opts.bTol=10^-5; end %%stopping criterion on residual std2(Fmask.*fft2(u)/N-Fp), if the noise standard variation is known, can set btol to be sigma
%                 if ~isfield(opts,'xTol'), opts.xTol=10^-5; end %%stopping criterion, if the noise standard variation is known, can set btol to be sigma
%                 if ~isfield(opts,'verbose'), opts.verbose=0; end %% display message
%                 if ~isfield(opts,'h0'),  opts.h0=16; end %% weight filter parater, depends on noise and image standard variation. for example: for barbara [0, 255], h0=20; To be adapted for normalized image
%                 if ~isfield(opts,'nWin'), opts.nWin=2;  end    %% patch size [2*nwin+1, 2*nwin+1]
%                 if ~isfield(opts,'nBloc'), opts.nBloc=5; end   %% search window size [2*bloc+1, 2*bloc+1]
%                 if ~isfield(opts,'nNeigh'), opts.nNeigh=10; end   %% number of best neighbors (real neighbors size: 2*nNeigh+4)
%                 if ~isfield(opts,'nWeightupdate'), opts.nWeightupdate=20;  end %0: no weight update, otherwise update steps 
%                 if ~isfield(opts,'denoising_type'), opts.denoising_type=1; end %% NLTV denoising algorithm
%                 % Projected NTV
%                 v1 = u2im; % init with TwICT o/p
%                 tau=1./(opts.delta*opts.mu);
%                 n=0;
%                 condition=1; % if condition '0' means TwICT elseif '1' means NTVICT
%                 %% OUTPUT M INDEX OF SIGNIFIACNT COMPONENTS  IN TWICT STAGE; THAT HAS TO BE EXPERIMENTED
%                 %% Main loop
%                 while (condition)  
%                   if (opts.nWeightupdate) %% update weight
%                       if(mod(n,opts.nWeightupdate)==0)
%                         wopts=update_weight(real(v1),opts.h0,opts.nWin,opts.nBloc,opts.nNeigh);          
%                       end
%                   end
%                   u2Pv = denoising_SBNLTV(v1,tau,tau./4,opts.nDenoising,wopts); % is it acc. to step 2 of paper 2, pls. chk.
%                   n = n+1; % increment for loop
%                   condition=(n<1);%opts.nOuter&& energy(n)>opts.bTol); is this K???
%                 end
%                 u1im = u2Pv; % step 6 of paper 2
%                 obj.x = u2Pv; % saving final reconstructed image
%                 % writing every o/p of NTVICT algorithm
%                 %psnr1 = PSNR(double(uint8(f*256)), double(uint8(u2Pv*256)));
%                 %SSIMR = ssim(double(uint8(f*256)), double(uint8(u2Pv*256)));
%                 %imwrite(uint8(u2Pv*256),[outpath 'cameraman.tif' '_' num2str(b
%             end
%             %toc
%             %psnr1 = PSNR(double(uint8(f*256)), double(uint8(u2Pv*256)));
%             %SSIMR = ssim(double(uint8(f*256)), double(uint8(u2Pv*256)));
            
            
        end      
      
    end
    
    %% queries
    methods
        
        
    end
    
    
    %% tests
    methods
        
        
        
    end
    
    %% graphics
    methods
        
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='encoder and decoder are set ';
            ok= isa(obj.A,'FrameTrafo') && isa(obj.Phi,'FrameTrafo');
        end
    end
    
    
    
end