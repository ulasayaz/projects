function isCol = isCollinear3x3(indices, N)

if N ~= 3
    error('collinearity check only works for a 3x3 matrix')
end
[I,J] = ind2sub([N,N],indices);
[~,order] = sort(I);
I = I(order); J = J(order);
if ( I(1) == I(2) && I(2) == I(3) ) || ... % same row
        ( J(1) == J(2) && J(2) == J(3) ) || ... % same column
        ( I(1)==1 && I(2)==2 && I(3)==3 && J(1)==1 && J(2)==2 && J(3)==3 ) || ... % main diagonal
        ( I(1)==1 && I(2)==2 && I(3)==3 && J(1)==3 && J(2)==2 && J(3)==1 ) % opposite diagonal
    isCol = 1;
else
    isCol = 0;
end