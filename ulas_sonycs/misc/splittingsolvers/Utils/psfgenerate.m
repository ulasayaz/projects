function h=psfgenerate(type,n)

if type == 1
% Guassian
      h = fftshift(normal_pdf(linspace(-2,2,n),0,0.0015))';
      h = real(ifft(fft(h)));
      h = h*h';
elseif type == 2
% Moving-Average 
      q = 9;
      h = shift([ ones(1,q) , zeros(1, n-q)]'/q,-floor(q/2));
      h = h*h';
elseif type == 3
% Exponential
      nu = 0.5;
      [Ii,Ji] = meshgrid([-n/2:n/2-1],[-n/2:n/2-1]);
      h = fftshift(exp(-(abs(Ii).^nu+abs(Ji).^nu)));
elseif type == 4
% Algebraic
      nu = 1;
      [Ii,Ji] = meshgrid([-n/2:n/2-1],[-n/2:n/2-1]);
      h = fftshift(1./(1 + abs(Ii) + abs(Ji)).^nu);
elseif type == 5
% Spherically-symmetric 1/f
      %q  = 7;
      %[Ii,Ji] = meshgrid([-q:q-1],[-q:q-1]);
      [Ii,Ji] = meshgrid([-n/2:n/2-1],[-n/2:n/2-1]);
      h = fftshift(1./(1 + Ii.^2 + Ji.^2));
      %h = zeros(n,n);
      %h(end/2-q+1:end/2+q,end/2-q+1:end/2+q) = 1./(1 + Ii.^2 + Ji.^2);
      %h = fftshift(h);
elseif type == 6
% The r used below is the same as that used in the Kalifa paper.
      r = 1;
      fft_row = [ones(1,n/4) (2^r) * (abs( 2*((n/4+1):(n/2))/(n)-1)).^r ] ; 
      sum1 = norm([ fft_row  fliplr(fft_row)]);
      fft_row = [ fft_row  fliplr(fft_row)];
      fft_row(~fft_row) = 1E-5;
      sum2 = norm(fft_row(:));
      fft_row = fft_row*sum1/sum2;  
      h_row   = real(ifft(fft_row));
      h_row = h_row(:);
      clear fft_col
      h_col = h_row;
      h = h_col * h_row.';
elseif type == 7
% Dirac
      h = zeros(n,n);h(1,1)=1;
elseif type == 8
% HST
      h=fitsread('~/Matlab/WaveLab802/MCA802/MCALab110/Two-D/Datasets/astro/simu_psf.fits');
      h = fftshift(h);
else
      error('You did not specify any valid PSF')
end
