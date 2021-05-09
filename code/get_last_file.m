function last_file = get_last_file (path)
% This function gets the path of the last data file in the path 'path'
%  both 'last_file' and 'path' are strings

files = dir(fullfile(path, '*.mat'));
filenames = {files.name};
if isempty(filenames)
	error(['get_last_file: No .mat file in ', path]);
end
numbers = zeros(1, length(filenames));
for i = 1:length(filenames)
	name = filenames{i};
	if strcmp(name(1:3), 'dat') && strcmp(name(8:11), '.mat')
		numbers(i) = str2num(name(4:7));
	end
end
max_no = max(numbers);
if max_no == 0
	error(['get_last_file: No available data file in ', path]);
end
last_file = fullfile(path, ['dat', sprintf('%4.4d', max_no), '.mat']);
