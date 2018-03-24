function dX = apply_recursive_differentiation(model,x,requested);
dX = [];
for variable = 1:length(model.linearindicies)
    dx = zeros(length(model.c),1);
    dx(model.linearindicies(variable)) = 1;
    for i = 1:length(model.evaluation_scheme)
        switch model.evaluation_scheme{i}.group
            case 'eval'
                for j = model.evaluation_scheme{i}.variables
                    k = model.evalMap{j}.variableIndex;
                    if any(dx(k))
                        if any(requested(model.evalMap{j}.computes))
                            if length(model.evalMap{j}.computes) == 1
                                dx(model.evalMap{j}.computes) = dx(k)'*model.evalMap{j}.properties.derivative(x(k));
                            else
                                dx(model.evalMap{j}.computes) = dx(k).*model.evalMap{j}.properties.derivative(x(k));
                            end
                        end
                    end
                end
            case 'monom'
                computed = model.monomials(model.evaluation_scheme{i}.variables);
                for j = computed
                    if requested(j)
                        dp = 0;
                        monomsj = model.monomtable(j,:);
                        for k = find(dx' & monomsj)
                            monoms = monomsj;
                            monoms(k) = 0;
                            r = model.monomtable(j,k);
                            s = find(monoms);
                            dp = dp + r*x(k)^(r-1)*dx(k)*prod((x(s)').^monoms(s));
                        end
                        dx(j) = real(dp);
                    end
                end
            otherwise
        end
    end
    dX = [dX dx];
end