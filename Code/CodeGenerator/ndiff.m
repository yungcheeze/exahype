% preprocessing: replace in Fortran output all
% whitespaces with linebreaks
% sed -e 's/\s\+/\n/g' inputfile > outputfile
%
% now the two files have the same format
% compare entry by entry
id = fopen('result1','r');
A = cell2mat(textscan(id,'%f','headerlines',0));
fclose(id);
id = fopen('result2','r');
B = cell2mat(textscan(id,'%f','headerlines',0));
fclose(id);
display('Max error:')
err = A-B;
max(err)

% dangerous without tolerance? Worked for my examples
display('Line where error occurs:')
idx = find(max(err)==err)
