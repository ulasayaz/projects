function item = NthList(list, pos)
% Nthlist -- find the nth argument in a list data structure
% "list" is either a string where the elements (substrings) are
% separated by a blank or a vector where the subvectors are separated by inf
%  Usage
%	item = NthList(list, pos)
%  Inputs
%	list	the list data structure
%	pos	integer, the nth
%  Outputs
%	item	the data if the pos'th of the list
%  See Also
%	MakeList, LengthList
%

if LengthList(list) == 1,
	item = list;
else
	if isstr(list),
		ind = 1:length(list);
		ends = [0 ind(list == ' ')];
		item = list((ends(pos)+1):(ends(pos+1)-1));
	else
		ind = 1:length(list);
		ends = [0 ind(list == inf)];
		item = list((ends(pos)+1):(ends(pos+1)-1));
	end
end

