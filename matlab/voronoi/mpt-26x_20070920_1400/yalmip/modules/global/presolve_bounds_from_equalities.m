function p = presolve_bounds_from_equalities(p)
LU = [p.lb p.ub];
if p.K.f >0
    
    % Find bounds from sum(xi) = 1, xi>0
    for j = 1:p.K.f
        if p.F_struc(j,1)>0
            [row,col,val] = find(p.F_struc(j,:));
            if all(val(2:end) < 0)
                if all(p.lb(col(2:end)-1)>=0)
                    p.ub(col(2:end)-1) = min( p.ub(col(2:end)-1) , val(1)./abs(val(2:end)'));
                end
            end
        end
    end
        
    for j = 1:p.K.f
        % Simple x == y
        done = 0;
        if p.F_struc(j,1)==0
            [row,col,val] = find(p.F_struc(j,:));
            if length(row) == 2
                if val(1) == -val(2)
                    p.lb(col(1)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
                    p.lb(col(2)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
                    p.ub(col(1)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
                    p.ub(col(2)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
                    done = 1;
                elseif val(1) == val(2)
                    p.lb(col(1)-1) = max(p.lb(col(1)-1),-p.ub(col(2)-1));
                    p.lb(col(2)-1) = max(-p.ub(col(1)-1),p.lb(col(2)-1));
                    p.ub(col(1)-1) = min(p.ub(col(1)-1),-p.lb(col(2)-1));
                    p.ub(col(2)-1) = min(-p.lb(col(1)-1),p.ub(col(2)-1));
                    done = 1;
                end               
            end
        end
        if ~done
            b = p.F_struc(j,1);
            a = p.F_struc(j,2:end);
            ap = a.*(a>0);
            am = a.*(a<0);
            for k = find(a)
                L = p.lb;
                U = p.ub;
                L(k) = 0;
                U(k) = 0;
                if a(k) > 0 & (p.ub(k)-p.lb(k)) > 1e-8
                    newlower = (-b - ap*U - am*L)/a(k);
                    newupper = (-b - am*U - ap*L)/a(k);
                    p.ub(k) = min(p.ub(k),newupper);
                    p.lb(k) = max(p.lb(k),newlower);
                end
            end
            b = -p.F_struc(j,1);
            a = -p.F_struc(j,2:end);
            ap = a.*(a>0);
            am = a.*(a<0);
            for k = find(a)
                L = p.lb;
                U = p.ub;
                L(k) = 0;
                U(k) = 0;
                if a(k) > 0 & (p.ub(k)-p.lb(k)) > 1e-8
                    newlower = (-b - ap*U - am*L)/a(k);
                    newupper = (-b - am*U - ap*L)/a(k);
                    p.ub(k) = min(p.ub(k),newupper);
                    p.lb(k) = max(p.lb(k),newlower);
                end
            end    
        end      
    end
end
close = find(abs(p.lb - p.ub) < 1e-12);
p.lb(close) = (p.lb(close)+p.ub(close))/2;
p.ub(close) = p.lb(close);
if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end