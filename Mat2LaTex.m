function LaTexEx = Mat2LaTex(A)
% Copyright: ËóËøËø

Asize = size(A);
LaTexEx = ['\left[\begin{matrix}'];

for i=1:Asize(1)
    for j=1:Asize(2)
        if j < Asize(2)
            if mod(A(i,j),1)<1e-5 % if this number is integer
                LaTexEx = [LaTexEx num2str(A(i,j)) '&'];
            else % if it is a rational number, extract its numerator and denominator
                [N,D] = numden(sym(A(i,j)));
                LaTexEx = [LaTexEx, '\frac{', num2str(double(N)),'}{',num2str(double(D)) ,'}', '&'];
            end
        else
            if mod(A(i,j),1)<1e-5 % if this number is integer
                LaTexEx = [LaTexEx num2str(A(i,j))];
            else % if it is a rational number, extract its numerator and denominator
                [N,D] = numden(sym(A(i,j)));
                LaTexEx = [LaTexEx, '\frac{', num2str(double(N)),'}{',num2str(double(D)) ,'}'];
            end
        end
    end
    LaTexEx = [LaTexEx '\\'];
end
LaTexEx = [LaTexEx '\end{matrix}\right]'];
clipboard('copy', LaTexEx)
end

