
% network model

param.networkSize = 20;  %number of nodes
param.numberNetworks =200; %number of networks
param.area = 50; %area
param.lowerB = 1; % area lower bound
param.upperB = param.area; % area upper bound
param.range = 16; % communication range

param.connected = 1; % if a connected graph is required
param.deals = 5000; % number of random networks will be generated at each iteration in order to achieve a connected graph 

param.network=0;
param.algebraic = 1; % if the algebraic threshold need to be fulfilled
param.normalized=0; %algebraic connectivity normalized 
param.fractionIteration=1/param.networkSize; %1 if is one node, < 1 if it is a fraction of the network n

%Algebraic Connectivity

param.threshold = .2;
param.epsilon = .2;
param.sigma = sqrt((-param.range^2)/(2*log(param.threshold)));
param.radius=.5*param.range;
param.gamma = -(1/(param.sigma^2))*100;

%persistence Information - you need to paramup the path where the experiment
%results will be saved

param.path = ['Z:\Project\database\' num2str(param.networkSize) '\' strcat(num2str(param.area),'_',strcat(num2str(param.range))) '\'];

if ~isequal(exist([param.path]),7)
    mkdir(param.path);
% else    
%    imageData = dir(fullfile(param.path,'*.*'));    
end
 
param.fileName = 'network'; % If you would like to define specific name prefix to result files 

param.obstacle = 'obstacle'; % If you would like to define specific name prefix to result files 

param.obstacleN = 100;

param.minDistObstacle=2;

param.lowerX= param.area + 5;
param.upperX= param.area * 10;

param.lowerY= param.lowerB;
param.upperY= param.upperB;

save([param.path 'setup'], 'param');
