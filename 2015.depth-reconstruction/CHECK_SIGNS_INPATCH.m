if strcmp(settings.sampleMode, 'hvgrid')==1
    cc_mask = bwconncomp(1-mask,4);
    labels_mask = labelmatrix(cc_mask);
    RGB_label_mask = label2rgb(labels_mask, @spring, 'w', 'shuffle');
    figure
    imshow(RGB_label_mask)
    signFlag = 0;
    
    for i=1:cc_mask.NumObjects
        fprintf('checking area %d\n',i)
        Mi_indices = cc_mask.PixelIdxList{i};
        [Xi,Yi] = ind2sub([Nx,Ny],Mi_indices);
        Xmax = max(Xi);
        Xmin = min(Xi);
        Ymax = max(Yi);
        Ymin = min(Yi);
        % Yi = Yq(Mi_indices);
        Zi = ZGT(Xmin-2:Xmax+2, Ymin-2:Ymax+2);
        [Nxi,Nyi] = size(Zi);
        [Hii,Vii] = createFiniteDiff2(Nxi,Nyi);
        
        HZi = vec(Zi * Hii);
        J = find(abs(HZi) < tol);
        I = setdiff(1:length(HZi),J);
        sHZi = zeros(length(HZi),1);
        sHZi(I) = sign(HZi(I));
        
        VZi = vec(Vii * Zi);
        J = find(abs(VZi) < tol);
        I = setdiff(1:length(VZi),J);
        sVZi = zeros(length(VZi),1);
        sVZi(I) = sign(VZi(I));
        
        if ( max(sVZi) - min(sVZi) > 1 ) || ( max(sHZi) - min(sHZi)>1 )
            signFlag = 1;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure;
            subplot(1,2,1)
            surf(Zi*Hii); shading interp
            subplot(1,2,2)
            surf(Vii*Zi); shading interp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure(10)
            maskIm = mask;
            maskIm(Mi_indices) = ZGT(Mi_indices);
            subplot(222); imshow(maskIm);
            fprintf('sign change found in area %d\n',i)
        end
    end
    if signFlag==0
        disp('sign check passes successfully')
    end
    if signFlag == 0 && (norm( vec(ZGT * H), 1) + norm( vec(V * ZGT), 1 )) - (norm( vec(Z * H), 1) + norm( vec(V * Z), 1)) > 1e-5
        solMismatch = (norm( vec(ZGT * H), 1) + norm( vec(V * ZGT), 1 )) - (norm( vec(Z * H), 1) + norm( vec(V * Z), 1))
        error('sign change passed, but xo not in solution set, (solMismatch = %g)\n',solMismatch)
    end
end
%% compute mask and visualize
Fig = figure(11); clf
colormap hsv
cm = colormap;
hold on; set(gcf,'name','original VS reconstruction')
if settings.useInputImage
    mesh(ZGT); % black
    mesh(Z);
else
    mesh(1:Nx, 1:Ny, ZGT); % black
    mesh(1:Nx, 1:Ny, Z);
end
view(3)

%% Check the gradients
Fig = figure(4); set(Fig,'Position',[pos * horz, (height- vert), horz, vert]); pos = pos+1; hold on
subplot(331); imshow(zeroOneTh(abs(VZ_GT)));  title('V * img');
subplot(332); imshow(zeroOneTh(abs(ZH_GT)));  title('img * H')
subplot(333); imshow(zeroOneTh(normGrad_GT)); hold on; title('img * grad')
plot(X_sample,Y_sample,'xb','markersize',5)
subplot(334); imshow(zeroOneTh(abs(VZ)));  title('V * Z')
subplot(335); imshow(zeroOneTh(abs(ZH)));  title('Z * H')
subplot(336); imshow(zeroOneTh(normGrad)); hold on; title('est. grad')
plot(X_sample,Y_sample,'xb','markersize',5)
subplot(337); imshow(zeroOneTh(abs(VZ_naive)));  title('V * naive')
subplot(338); imshow(zeroOneTh(abs(ZH_naive)));  title('naive * H')
subplot(339); imshow(zeroOneTh(normGrad_naive)); hold on; title('naive grad')

%% surf figures
Fig = figure(8); set(Fig,'Position',[pos * horz, (height-vert), horz, vert]); pos = pos+1;
[Xq, Yq] = meshgrid(1:size(ZGT,2), 1:size(ZGT,1));
subplot(131); surf(Xq,Yq,ZGT); title('original');
subplot(132); surf(Xq,Yq,Z); title('Z');
set(gcf,'name','original with samples')
subplot(133); ERZ = zeroOne(abs(Z-ZGT)); imshow(ERZ);

%% surf figures
Fig = figure;
subplot(121); surf(Xq,Yq,ZGT); title('ZGT');
subplot(122); surf(Xq,Yq,Z); hold on; title('Z');
plot3(X_sample, Y_sample, Z(samples), 'or', 'markersize', 20)