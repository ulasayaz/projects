function w = ProjectOntoTranslatedL1Ball(v,y,b)
% ~project v onto l1 ball centered at y with radius b
% PROJECTONTOL1BALL Projects point onto a translated L1 ball of specified radius.
%
% w = ProjectOntoL1Ball(v,y,b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  ||w - y||_1 <= b.
%
%   That is, performs Euclidean projection of v to the 1-norm ball 
% with center y and radius b.
%
% v and y should be column vectors. 
%
% Author: John Duchi (jduchi@cs.berkeley.edu)

if (b < 0)
  error('Radius of L1 ball is negative: %2.3f\n', b);
end

w = ProjectOntoL1Ball(v-y,b) + y;