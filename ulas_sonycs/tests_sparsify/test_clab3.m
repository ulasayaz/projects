% test matlab's 3D wavelet transformation using Gunter's code. 

%fn='cam512.png';
%fn='cameraman.bmp';
%fn='cam256.png';
%fn='cam128.png';
fn='cam64.png';

cs=33; % 1/ratio

'start'
L=32; % number of frames

videoflag=2; % 1 for rotation, 2 for shifts

signal2=Signal2D.make_fromImage(fn);

pixshift=1;

%rot=90/L;
%rot=5; 
rot=pixshift; % to rotate as much as 'pixshift'

if videoflag==2  % shift
        units=pixshift*[1;1];
else % rotation
        units=rot;
end

signal3=Signal3D.make_CyclefromLowerDimSig(signal2,videoflag,L,units);

x=signal3.xn;

tic
xc=fdct3d_forward_mod(x);
%xc=fdct3d_forward(x);
toc

% Get the total size of the matrices
for s=1:length(xc)
    length(xc{s})
end

Le=0;
for s=1:length(xc)
    for w=1:length(xc{s})
        if mod(w,100)==0
            w
        end
        Le = Le+length(xc{s}{w}(:));
    end
end

% Get threshold value
%cfs =[];
% for s=1:length(xc)
%     fprintf('.')
%     for w=1:length(xc{s})
%         if mod(w,100)==0
%             w
%         end
%         cfs = [cfs; abs(xc{s}{w}(:))];
%     end
% end
cfs=zeros(Le,1);
ind=1;
for s=1:length(xc)
    for w=1:length(xc{s})
        if mod(w,100)==0
            w
        end
        ind2=ind+length(xc{s}{w}(:));
        cfs(ind:ind2-1) = abs(xc{s}{w}(:));
        ind=ind2;
    end
end
cfs = sort(cfs,'descend'); %cfs = cfs(end:-1:1);
nb = round(length(cfs)/cs);
cutoff = cfs(nb);

% set small coefficients to zero
for s=1:length(xc)
    for w=1:length(xc{s})
        xc{s}{w} = xc{s}{w} .* (abs(xc{s}{w})>cutoff);
    end
end
%%
inv_xc = real(fdct3d_inverse_mod(xc));
%inv_xc = real(fdct3d_inverse(xc));

'mean psnr and ssim'
PSNR3(x,inv_xc)
SSIM3(x,inv_xc)

y=normalize_frame(inv_xc);
playData(y);


