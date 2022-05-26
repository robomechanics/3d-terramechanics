function result = matrix2latex(A, dig)
%This function takes in a matrix and prints a string with the formatting
%to turn its values into a latex table, so that the result can be directly
%copy-pasted from the terminal. If "dig" is provided, it is the number
%of digits to include after the decimal

if ~exist('dig')
    dig = 2; %defaults to two digit precision if no option entered
end
formatspec = strcat('%', sprintf('0.%0.0ff', dig));

[n, m] = size(A);
result = '';
for i=1:n
    result = sprintf('%s %s ', result, num2str(A(i,1), formatspec));
    for j=2:m
        result = sprintf('%s & %s ', result, num2str(A(i,j), formatspec));
    end
    result = strcat(result, ' \\\\ \n');
end
fprintf(result)