% FUNCTION TO LOAD DATA FROM TEXT FILE - WORKS WITH BOTH MATLAB AND OCTAVE 
function C = textscan_mod_v1(file_name,string_read,delimiter,varargin)


if rem(length(varargin),2)>0
    error('Wrong number of inputs');
end
if strcmp(delimiter,'\t')
  delimiter=char(9);
end

headerlines_num=0;
for i=1:2:length(varargin)
    switch varargin{i}
        case 'headerlines'
            headerlines_num=varargin{i+1};
        otherwise
            error('Wrong inputs');
    end
end
	
fid = fopen(file_name,'rt');
for i=1:headerlines_num
	tLines = fgets(fid);
end

tLines = fgets(fid);

numCols = numel(strfind(tLines,delimiter)) + 1;
fclose(fid);

numColsNeeded=length(find(string_read=='%'));
if numCols>numColsNeeded
	fid=fopen(file_name);
	C= textscan(fid,[string_read '%*[^\n]'],'headerlines',headerlines_num,'delimiter',delimiter);
	fclose(fid);
elseif  numCols==numColsNeeded
	fid=fopen(file_name);
	C= textscan(fid,string_read,'headerlines',headerlines_num,'delimiter',delimiter);
	fclose(fid);
else
	error('Wrong number of columns');
end

end


