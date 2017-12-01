function str = arr2str(arr)
str = strjoin(cellstr(num2str(arr(:))),', ');
end