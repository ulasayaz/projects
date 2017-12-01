close all;

fn='cam512.png';
cs=33; % 1/ratio

'start'
L=512; % number of frames

videoflag=2; % 1 for rotation, 2 for shifts
compress=1; % 1 for mpeg4 video compression, 0 for raw video (no compression)

signal2=Signal2D.make_fromImage(fn); % this scales the image

pixshift=1;

%rot=90/L;
%rot=5; 
rot=pixshift; % to rotate as much as 'pixshift'

if videoflag==2  % shift
        units=pixshift*[1,1];
else % rotation
        units=rot;
end

signal3=Signal3D.make_fromLowerDimSig(signal2,videoflag,L,units);

data3d=signal3.xn;

%%

% D1=data3d(:,:,2)-data3d(:,:,1);
% 
% D1_sort = sort(abs(D1(:)),'descend');
% 
% thresh = D1_sort(floor(numel(D1)/cs));
% 
% s_D1 = D1 .* (abs(D1) > thresh);
% 
% inv_D1=wlab2D(D1,1/cs,1);

% subplot 3d video
% ver=ceil(sqrt(L));
% hor=ver-1;
% for k=1:L
%     k
%     subplot(ver,ver,k);
%     imgray(data3d(:,:,k));
% end

'finish'
        

%figure;

% x=imread('cam512.png');
% 
% %x = ReadImage('Daubechies');
% 
% x=double(x);
% 
% %subplot(2,2,1);
% %imgray(x); % displays pictures
% title('original');
% 
% ratio=0.05;
% 
% %% CREATE a VIDEO SEQ.
% 'start'
% L=32; % number of frames
% videoflag=1; % 1 for rotation, 2 for shifts
% 
% % Gunter's rescaling code from Signal3D.make_fromImage
% x=x-min(x(:));
% x=x/max(x(:));
% 
% data3d=zeros([size(x),L]);
% % build 3 d image as sequence of rotated images
% data3d(:,:,1)=x;
% 
% 
% pixshift=0.1;
% 
% rot=90/L;
% %rot=5; 
% %rot=pixshift; % to rotate as much as 'pixshift'
% 
% 
% j=2:L;
% alphas=rot*(j-1);
% 
% 
% unit=pixshift*[1,1];
% j=1:L-1;
% pix=kron(j,unit(:));
% 
% if videoflag==1
%     pixorigin=ceil(size(x)/2);
%     [X,Y]=meshgrid(1:size(x,2),1:size(x,2));
%     rads=sqrt((X-pixorigin(2)).^2+(Y-pixorigin(2)).^2);
%     x(rads>min(size(x))/2)=0;
%     for j=2:L
%         data3d(:,:,j)=imrotate(x,alphas(j-1),'crop');
%     end
% else
%     if isequal(pix, floor(pix))
%         % all shifts by full pixels
%         shifter=@circshift;
%     else
%         % subpixel shifts
%         shifter=@Signal2D.SubPixelShifter;
%     end
%     for j=2:L
%         data3d(:,:,j)=shifter(x,pix(:,j-1));
%     end
% end


