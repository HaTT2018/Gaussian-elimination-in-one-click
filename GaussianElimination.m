function Latexex = GaussianElimination(A,b, ifLU)
Latexex = ['\begin{eqnarray}&&'];
if size(A,1) == size(A,2) && ifLU==1 % it's a square matrix
    % factorize LU
    [P, L, U, Latexex] = FactLU(A, Latexex);
else
    G = [A b];
    Latexex = [Latexex Mat2LaTex(['A' 'b']) '=' Mat2LaTex(G)];
    for i =1:size(G,1)-1
        % multiSet = zeros(size(G));
        
        [multi, tempG, Latexex] = eliminate(G, G(i:size(G,1), i:size(G,2)), i, Latexex); % i is column
        G(i:size(G,1), i:size(G,2)) = tempG;
        
        % multiSet(i:size(G,1), i:size(G,2)) = multi;
    end
    for i=1:size(G,1)
        G(i,:) = G(i,:)/G(i,sum(G(i,:)==0)+1);
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
P = eye(size(A));
L = P;
    while 1
        [ifchange, P, A] = PA(P, A);
        if ifchange==0
            Latexex = [Latexex 'P=' Mat2LaTex(P) '\notag\\&&'];
            break
        end
    end
    for i = 1:size(A,1)
        [multi, tempA] = eliminate(A(i:size(A,1), i:size(A,2)));
        A(i:size(A,1), i:size(A,2)) = tempA;
        L(i:size(A,1), i:size(A,2)) = multi;
    end
    U=A;
end


function [multi, A, Latexex] = eliminate(G, A, column, Latexex) 
    %G(G(:,1)<-1e-5, :) = -1 * G(G(:,1)<-1e-5, :);
    %G = sortrows(G, -1);
    multi = eye(size(A));
    if abs(A(1,1)) < 1e-5
        return
    end
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
    for i=1:size(A,1)
        if det(A(1:i,1:i))==0
            ifLUable = [ifLUable i];
        end
    end
    if sum(ifLUable)==0
        return
    else
        P(ifLUable(1),ifLUable(1)) = 0;
        P(ifLUable(1),ifLUable(1)+1) = 1;
        P(ifLUable(1)+1,ifLUable(1)+1) = 0;
        P(ifLUable(1)+1,ifLUable(1)) = 1;
        A = P*A;
        ifchange = 1;
    end
    
end

function [multi, G, Latexex] = Substitute(OrgG, G, column, Latexex)
multi = eye(size(G));
    for i = 2:size(G, 1)
        multi(1,i) = G(1, i);
        up = ['(',num2str(column), ',', num2str(i), ')'];
        down = ['l_{' num2str(column) num2str(i) '}=' num2str(multi(1,i))];
        Latexex = arrow(Latexex, up, down);
        G(1, :) = G(1, :) - G(i, :)*G(1, i);
        OrgG(column:size(OrgG,1), column:size(OrgG,2)) = G;
        Latexex = [Latexex Mat2LaTex(OrgG) '\notag\\&&'];
    end
end

function Latexex = arrow(Latexex, up, down)
    Latexex = [Latexex '\xrightarrow[' down ']{' up '}'];
end