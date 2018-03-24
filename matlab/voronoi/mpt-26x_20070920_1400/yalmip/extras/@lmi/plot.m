function x_opt = plot(varargin)
%plot               plots feasible set
%
% p = plot(F,x,c,n,options)
%
% F:  set object
% x:  projected variables [At most three variables]
% c:  color [double] ([r g b] format)
% n:  #vertices [double ]
% options: options structure from sdpsettings


% Author Johan Löfberg
% $Id: plot.m,v 1.14 2007/04/02 14:26:00 joloef Exp $

% Get the onstraints
if nargin<1
    return
end

F = varargin{1};

if length(F)==0
    return;
end

if nargin < 5
    opts = sdpsettings('verbose',0);
else
    opts = varargin{5};
    if isempty(opts)
        opts = sdpsettings('verbose',0);
    end
end

if nargin < 3
    color = [1 0.1 0.1];
else
    color = varargin{3};
    if isempty(color)
        color = [1 0.1 0.1];
    end
end

% Plot onto this projection (at most in 3D)
if nargin < 2
    x = [];
else
    x = varargin{2};
    if ~isempty(x)
        x = x(:);
        x = x(1:min(3,length(x)));
    end
end

if isempty(F)
    return
end

% Create a model in YALMIPs low level format
% All we change later is the cost vector
%sol = solvesdp(F,sum(x),opts);
[model,recoverdata,diagnostic,internalmodel] = export(F,[],opts,[],[],0);
if isempty(internalmodel) | (~isempty(diagnostic) & diagnostic.problem)
    error('Could not create model. Can you actually solve problems with this model?')
end
internalmodel.options.saveduals = 0;
internalmodel.getsolvertime = 0;
internalmodel.options.dimacs = 0;

% Try to find a suitable set to plot
if isempty(x)
    if isempty(internalmodel.extended_variables)
        x = depends(F);
        x = x(1:min(3,length(x)));
        localindex = 1;
        localindex = find(ismember(recoverdata.used_variables,x));
    else
        % not extended variables
        candidates = setdiff(1:length(internalmodel.c),internalmodel.extended_variables);
        % Not nonlinear variables
        candidates = candidates(find(internalmodel.variabletype(candidates)==0));
        % Settle with this guess
        localindex = candidates(1:min(3,length(candidates)));
        x = localindex;
    end
else
    localindex = [];
    for i = 1:length(x)
        localindex = [localindex find(ismember(recoverdata.used_variables,getvariables(x(i))))];
    end
end

if nargin < 4
    if length(x)==3
        n = 100;
    else
        n = 25;
    end
else
    n = varargin{4};
    if isempty(n)
        if length(x)==3
            n = 100;
        else
            n = 25;
        end
    end
    if ~isa(n,'double')
        error('4th argument should be an integer>0');
    end
end


x_opt = [];
phi = [];
status = 0;
try % Try to ensure that we close h
    if length(x)==2
        mu =0.5;
    else
        mu=1;
    end
    n_ = n;
    n = ceil(mu*n);
    h = waitbar(0,['Please wait, solving ' num2str(n_) ' problems using ' internalmodel.solver.tag]);
    angles = (0:(n))*2*pi/n;
    if length(x)==2
        c = [cos(angles);sin(angles)];
    else
        c = randn(3,n);
    end
    i=1;
    while i<=n & status == 0
        xi = solvefordirection(c(:,i),internalmodel,localindex);
        x_opt = [x_opt xi];
        waitbar(i/n_,h)
        i=i+1;
    end

    if status==0 & length(x)==2
        % Close the set
        x_opt = [x_opt x_opt(:,1)];

        % Add points adaptively
        pick = 1;
        n = floor((1-mu)*n_);
        for i = 1:1:n
            for j= 1:(size(x_opt,2)-1)
                d = x_opt(:,j)-x_opt(:,j+1);
                distance(j,1) = d'*d;
            end
            [dist,pos]=sort(-distance);
            % Select insertion point
            phii=(angles(pos(pick))+angles(pos(pick)+1))/2;
            xi = solvefordirection([cos(phii);sin(phii)],internalmodel,localindex);
            d1=xi-x_opt(:,pos(pick));
            d2=xi-x_opt(:,pos(pick)+1);
            if d1'*d1<1e-3 | d2'*d2<1e-3
                pick = pick+1;
            else
                angles = [angles(1:pos(pick)) phii  angles((pos(pick))+1:end)];
                x_opt = [x_opt(:,1:pos(pick)) xi  x_opt(:,(pos(pick))+1:end)];
            end
            waitbar((ceil(n_*mu)+i)/n_,h);
        end
    end
    close(h);
catch
    %try
    close(h);
    %catch
    %end
end

if status
    if nargout==0
        plot(0);
    end
end

if nargout == 0 & status==0
    if length(x)==2
        patch(x_opt(1,:),x_opt(2,:),color);
    else
        K = convhulln(x_opt');
        p = patch('Vertices', x_opt', 'Faces', K, 'FaceVertexCData', color, 'FaceColor', color);%, 'FaceAlpha', 0.5);
        set(p,'LineStyle','none')
        lighting gouraud;
        view(3);
        camlight('headlight','local');
        camlight('headlight','local');
        camlight('right','local');
        camlight('left','local');
    end
end


function [xout,status] = solvefordirection(c,internalmodel,uv)
internalmodel.c = 0*internalmodel.c;
internalmodel.c(uv) = c;
sol  = feval(internalmodel.solver.call,internalmodel);
xout = sol.Primal;
xout = xout(uv(:));
status = sol.problem;




