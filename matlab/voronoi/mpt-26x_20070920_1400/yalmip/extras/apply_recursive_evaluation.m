function xevaled = apply_recursive_evaluation(p,xevaled)

xevaled = xevaled(:)';
for i = 1:length(p.evaluation_scheme)
    switch p.evaluation_scheme{i}.group
        case 'eval'
            xevaled = process_evals(p,xevaled,p.evaluation_scheme{i}.variables);
        case 'monom'
            xevaled = process_monomials(p,xevaled,p.evaluation_scheme{i}.variables);
            xevaled = real(xevaled);
        otherwise
    end
end
xevaled = xevaled(:);

function x = process_monomials(p,x,indicies);
indicies = p.monomials(indicies);
x(indicies) = prod(repmat(x,length(indicies),1).^p.monomtable(indicies,:),2);

function x = process_evals(p,x,indicies)
if isfield(p.evalMap{1},'prearg')
    for i = indicies
        arguments = p.evalMap{i}.prearg;
        arguments{2} = x(p.evalMap{i}.variableIndex);
        if isequal(arguments{1},'log') & (arguments{2}<=0)
            x(p.evalVariables(i)) = -1e4;
        else           
            x(p.evalMap{i}.computes(:)) = feval(arguments{:});
        end
    end
else
    for i = indicies
        arguments = {p.evalMap{i}.fcn,x(p.evalMap{i}.variableIndex)};
        arguments = {arguments{:},p.evalMap{i}.arg{2:end-1}};
        if isequal(arguments{1},'log') & (arguments{2}<=0)
            x(p.evalVariables(i)) = -1e4;
        else
            x(p.evalVariables(i)) = feval(arguments{:});
        end
    end
end