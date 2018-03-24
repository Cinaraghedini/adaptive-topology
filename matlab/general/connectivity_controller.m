% filename: connectivity_controller.m
% Purpose: computes the weights for the algebraic connectivity maintenance
% - position - x and y coordinates for each network agent
% - param - parametrization struct
% - options - graph options (set the parametrization for graph properties
% computation
% Output: 
% data - weights for connectivity maintenance
% Reference:
% L. Sabattini, N. Chopra, and C. Secchi. Decentralized connectivity maintenance for cooperative control of mobile robotic sys-
% tems. The International Journal of Robotics Research (SAGE), 32(12):1411–1423, October 2013.

function  dotxy = connectivity_controller(position,param,options)


%% to compute global properties related to algebraic connectivity%%%%%%%%%%%%%%%%%%%%%

[data]=set_controlData(position,param,options); % compute global properties related to algebraic connectivity

% x position

nv=size(position,1);

pointer=1;
positionx=data(pointer:nv,1);

% y position

pointer=pointer+nv;
positiony=data(pointer:pointer+nv-1,1);

% algebraic connectivity estimated by node

pointer=pointer+nv;
ac=data(pointer:pointer+nv-1,1);

% eigenvector

pointer=pointer+nv;

eigV=data(pointer:pointer+nv-1,1);

[A] = initialize_matrixA(position,param,options);

dotx = zeros(nv,1);
doty = zeros(nv,1);

for i=1:nv
    for j=1:nv
        k=0; 
        if A(i,j)
            if ac(i)>param.epsilon
              k = (-(1/(param.sigma^2))*(csch(ac(i)-param.epsilon))^2 ) * (A(i,j)*((eigV(i)-eigV(j))^2)) ;
            else
              k = -(1/(param.sigma^2))*100 * (A(i,j)*((eigV(i)-eigV(j))^2));  % ORIGINAL
            end
        end
        
        dotx(i)=dotx(i)+ (k*(positionx(i)-positionx(j))); % X-axis component of the control effort
        doty(i)=doty(i)+ (k*(positiony(i)-positiony(j))); % Y-axis component of the control effort
    end
end

dotxy=[dotx doty] * param.gainConnectivity;