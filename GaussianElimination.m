function Latexex = GaussianElimination(A,b, ifLU)
Latexex = ['\begin{eqnarray}&&'];
if ifLU==1 % it's a square matrix
    % factorize LU
    [P, L, U, Latexex] = FactLU(A, Latexex);
else
    G = [A b];
    Latexex = [Latexex '\left[\begin{matrix}\textbf{\textit{A}}&\textbf{\textit{b}}\\\end{matrix}\right]=' Mat2LaTex(G)];
    for i =1:size(G,1)-1
        % multiSet = zeros(size(G));
        [multi, tempG, Latexex] = eliminate(G, G(i:size(G,1), i:size(G,2)), i, Latexex); % i is column
        G(i:size(G,1), i:size(G,2)) = tempG;
        
        % multiSet(i:size(G,1), i:size(G,2)) = multi;
    end
    for i=1:size(G,1)
        if sum(G(i, :)) == 0
            continue
        end
        G(i,:) = G(i,:)/G(i,count_zeros(G(i,:)));
    end
    Latexex = arrow(Latexex, 'clean\ up', '');
    Latexex = [Latexex Mat2LaTex(G)];
    Latexex = [Latexex '\notag\\&&Then\ apply\ substitution:\notag\\&&'];
    for i = 1:min(size(G))
        [multi, tempG, Latexex] = Substitute(G, G(i:size(G,1), i:size(G,2)), i, Latexex);
        G(i:size(G,1), i:size(G,2)) = tempG;
    end
end
Latexex = [Latexex '\end{eqnarray}'];
end

function [P, L, U, Latexex] = FactLU(A, Latexex)
P = eye(size(A,1));
L = P;
    while 1
        [ifchange, P, A] = PA(P, A);
        if ifchange==0
            Latexex = [Latexex 'P=' Mat2LaTex(P) '\notag\\&&'];
            break
        end
    end
    for i = 1:size(A,1)-1
        [multi, tempA, Latexex] = eliminate(A, A(i:size(A,1), i:size(A,2)), i, Latexex);
        A(i:size(A,1), i:size(A,2)) = tempA;
        L(i:size(A,1), i:size(A,2)) = multi;
    end
    U=A;
end


function [multi, A, Latexex] = eliminate(G, A, column, Latexex) 
    %G(G(:,1)<-1e-5, :) = -1 * G(G(:,1)<-1e-5, :);
    %G = sortrows(G, -1);
    multi = eye(size(A));
    %if abs(A(1,1)) < 1e-5
    %    return
    %end
    for i = 2:size(A, 1)
        if A(i, 1) == 0
            continue
        else
            multi(i, 1) = A(i, 1)/A(1, 1);
            up = ['(',num2str(i), ',', num2str(column), ')'];
            down = ['l_{' num2str(i) num2str(column) '}=' num2str(multi(i,1))];
            Latexex = arrow(Latexex, up, down);
            A(i, :) = A(i, :) - A(1, :)*multi(i, 1);
            G(column:size(G,1), column:size(G,2)) = A;
            Latexex = [Latexex Mat2LaTex(G) '\notag\\&&'];
        end
    end
    % G(1, :) = G(1, :)/G(1,1);
end

function [ifchange, P, A] = PA(P, A)
    ifchange = 0;
    ifLUable = [];
    for i=1:min(size(A))
        if det(A(1:i,1:i))==0
            ifLUable = [ifLUable i];
        end
    end
    if sum(ifLUable)==0
        return
    else
        c = randi(size(P,1));
        P(ifLUable(1),ifLUable(1)) = 0;
        P(ifLUable(1),1+mod(ifLUable(1)+c,size(P,1))) = 1;
        P(1+mod(ifLUable(1)+c,size(P,1)),1+mod(ifLUable(1)+c,size(P,1))) = 0;
        P(1+mod(ifLUable(1)+c,size(P,1)),ifLUable(1)) = 1;
        A = P*A;
        ifchange = 1;
    end
    
end

function [multi, G, Latexex] = Substitute(OrgG, G, column, Latexex)
multi = eye(size(G));
    for i = 2:size(G, 1)
        multi(1,count_zeros(G(i, :))) = G(1, count_zeros(G(i, :)));
        up = ['(',num2str(column), ',', num2str(i), ')'];
        down = ['l_{' num2str(column) num2str(i) '}=' num2str(multi(1,count_zeros(G(i, :))))];
        Latexex = arrow(Latexex, up, down);
        G(1, :) = G(1, :) - G(i, :)*G(1, count_zeros(G(i, :)));
        OrgG(column:size(OrgG,1), column:size(OrgG,2)) = G;
        Latexex = [Latexex Mat2LaTex(OrgG) '\notag\\&&'];
    end
end

function Latexex = arrow(Latexex, up, down)
    Latexex = [Latexex '\xrightarrow[' down ']{' up '}'];
end

function ind = count_zeros(row)
    for i=1:size(row,2)
        if row(1,i) ~= 0
            break
        end
    end
    ind = i;
end