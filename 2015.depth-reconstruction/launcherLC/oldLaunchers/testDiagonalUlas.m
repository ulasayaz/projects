clear all
close all
clc

Nx = 10;
Ny = 20;

for test = 1:100
    Z = rand(Nx,Ny);
    [H1,V1] = createFiniteDiff1(Nx,Ny);
    
    vecDD_expected = 2 * vec(V1 * Z * H1);
    TV2_xy = 2 * kron(H1',V1);
    vecDD_actual = TV2_xy * vec(Z);
    
    if norm(vecDD_expected - vecDD_actual) > 1e-7
        error('wrong vectorization')
    else
        disp('succcess!!')
    end
end



   
   
   