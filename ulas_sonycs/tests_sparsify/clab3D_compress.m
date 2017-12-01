function [inv_xc] = clab3D_compress(x,cs,mode)
% performs curvelet transform with compression (rate=cs)
% returns the compressed image 'inv_xc'
% mode=0 uses default curvelet transform
% mode=1 uses modified curvelet transform functions

if mode==0
    forward=@fdct3d_forward;
    inverse=@fdct3d_inverse;
elseif mode==1
    forward=@fdct3d_forward_mod;
    inverse=@fdct3d_inverse_mod;
end

xc=forward(x);

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

inv_xc = real(inverse(xc));

end

