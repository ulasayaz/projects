function n = LengthList(list)
% LengthList -- length of a list data structure
% "list" is either a string where the elements (substrings) are
% separated by a blank or a vector where the subvectors are separated by inf
%  Usage
%	n = LengthList(list)
%  Inputs
%	list	the list data structure
%  Outputs
%	n	the length of list
%  See Also
%	MakeList, NthList
%

if isstr(list)
	n = sum(list == ' ');
else
	n = sum(list == inf);
end

if n == 0,
	n  = 1;
end