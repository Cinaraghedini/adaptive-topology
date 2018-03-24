function varargout = iff_internal(varargin)
X = varargin{1};
Y = varargin{2};
switch class(varargin{1})

    case 'sdpvar'

        if length(X)>1
            error('IMPLIES not implemented for this case');
        end

        switch class(Y)
            case 'sdpvar'               % X <--> Y
                varargout{1} = set(Y == X);

            case {'lmi','constraint'}
                Y=set(Y,[],[],1);
                switch settype(Y)
                    case 'elementwise'  % X <--> Y(:)>=0
                        varargout{1} = binary_iff_lp(X,-sdpvar(Y));
                    case 'equality'     % X <--> Y(:)==0
                        varargout{1} = binary_iff_eq(X,sdpvar(Y));
                    otherwise
                        error('IFF not implemented for this case');
                end

            otherwise
                error('IFF not implemented for this case');
        end

    case {'lmi','constraint'}

        if isa(X,'constraint')
            X = set(X,[],[],1); % FIX: passes one to avoid pruning infeasible constraints
        end
        switch class(Y)
            case 'sdpvar'
                switch settype(X)
                    case 'elementwise'
                        varargout{1} = binary_iff_lp(Y,-sdpvar(X));
                    case 'equality'
                        binvar Z W

                        %                        varargout{1} = set(implies(Z&W,Y))+binary_iff_lp(Z,sdpvar(X)+eps)+binary_iff_lp(W,-sdpvar(X)+eps);
                        X = [sdpvar(X)-eps;eps-sdpvar(X)];
                        varargout{1} = binary_iff_lp(Y,X);%,sdpvar(X)+eps)+binary_iff_lp(W,-sdpvar(X)+eps);
                        %                        varargout{1} = binary_iff_eq(Y,sdpvar(X));
                    otherwise
                        error('IFF not implemented for this case');
                end

            case {'lmi','constraint'} % F(X) <--> F(Y)
                d = binvar(1,1);
                varargout{1} = iff_internal(X,d)+iff_internal(Y,d);

            otherwise
                error('IFF not implemented for this case');
        end

    otherwise
        error('IFF not implemented for this case');
end


function F = binary_iff_eq(X,f)
[M,m,infbound] = derivebounds(f);
if infbound
    warning('You have unbounded variables in IFF leading to a lousy big-M relaxation.');
end
eps = 1e-2;

[nf,mf]=size(f);
if mf>1
    f = reshape(f,nf*mf,1);
end

if 1%nf*mf>1
    di1 = binvar(nf*mf,1);
    di2 = binvar(nf*mf,1);
    F = set(M*(1-X) >= f >= m*(1-X));
    % F = F + set(-eps + (M+eps).*(1-di) >= f >= eps + (m-eps).*(1-di))+set(X>=sum(di)-length(di)+1);
    % F = F + set(-eps + (M+eps).*di >= f >= eps + (m-eps).*di)+set(X>=sum(di)-length(di)+1);

    F  = F + set(f>=eps+(m-eps).*di1)+set(-f>=-eps+(-M+eps).*di2)+set(X>=sum(di1)-length(di1)+1)+set(X>=sum(di2)-length(di2)+1);

else
    x1 = binvar(1,1);
    x2 = binvar(1,1);
    F = set(M*(1-x1) + eps >= f >= -eps + m*(1-x2));
    F = F + set(f  >=  eps + m.*x1);
    F = F + set(f  <= -eps + M.*x2);
    F = F + set(x1+x2-1 <=X);

    %   F = F + set(iff(~X,~W | ~Z));% >= 1-Z) + set(X >= 1-W);
    % F = F + set(f >= eps + (m-eps).*X)+set(-f >= eps + (-M-eps).*X);
    % F = F + set(f >= eps + (m-eps)*Z)+set(-f >= eps + (-M-eps)*W);
    % F = F + set(X == (Z | W));
end

function F = binary_iff_lp(X,f)
[M,m,infbound] = derivebounds(f);
if infbound
    warning('You have unbounded variables in IFF leading to a lousy big-M relaxation.');
end
eps = 1e-8;

[nf,mf]=size(f);
if nf*mf==1
    F  = set(f <= M*(1-X)) + set(f>=eps+(m-eps)*X);
else
    if mf>1
        f = reshape(f,nf*mf,1);
    end
    if nf*mf>1
        di = binvar(nf*mf,1);
        % di=0 means the ith hypeplane is violated
        % X=1 means we are in the polytope
        F  = set(f <= M*(1-X)) + set(f>=eps+(m-eps).*di)+set(X>=sum(di)-length(di)+1) + set(X <= di);
        
        % Add some cuts for c < a'x+b < d
        [bA] = getbase(f);
        b = bA(:,1);
        A = bA(:,2:end);
        S = zeros(0,length(di));
        for i = 1:length(b)
            j = findrows(abs(A),abs(A(i,:)));
            j = j(j > i);
            if length(j)==1               
                S(end+1,[i j]) = 1;             
            end
        end
        if size(S,1) > 0
            % Add cut cannot be outside both constraints
            F = F + set(S*di >= 1);           
        end              
    else
        F  = set(f <= M*(1-X)) + set(f>=eps+(m-eps)*X);
    end
end
