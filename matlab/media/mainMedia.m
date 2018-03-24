% filename: mainMedia_Generation.m
% Purpose:  generating and persisting videos with the simulation evolution
% or images with network snapshots.
% Input:
% - parametrization file : setModelMedia_Video.m
% Output: videos or figure with network snapshots

% initialization setup

close all;
clear all;
clc;

setModelMedia; % parametrization setup

%network id initialization

if isempty(param.networkList) % if there is no specific networks set, sequencial numbers considering the number of network are generated
    param.networkList=transpose([1:param.numberNetworks].');
end

param.numberNetworks=length(param.networkList);  % set the number of networks

for k=1:size(param.labels,2)
    
    param.ref=k;
    
    param.color=param.colors(1,k);  %  set colors for identifying vulnerable nodes are
    
    % gain setting
    param.gainConnectivityController=param.gain(k,1);
    param.gainRobustnessControl=param.gain(k,2);
    param.gainCoverageController=param.gain(k,3);
    
    labelP=strcat(num2str(param.gainConnectivityController),'_',num2str(param.gainRobustnessControl),'_',num2str(param.gainCoverageController));
    
    mainworkingPath = [param.mainPath labelP '\']; % set main iteration data path
    
    % upload obstacles, failure time and nodes which failed
    
    aux=load([param.mainPath 'param.mat']);
    param.obstacles=aux.param.obstacles;
    param.failureT=aux.param.failureT;
    param.listNodeFailure=aux.param.listNodeFailure;
    
    for i=param.networkList % for each network
        
        workingPath = [mainworkingPath num2str(num2str(i)) '\']; % set the network data path
        
        resultsworkingPath = [workingPath param.dirName] ; % set path where video will be persisted according parametrization file
        
        % creating media folder
        
        if ~isequal(exist(workingPath, 'dir'),7)
            disp(sprintf('Network: %d, Connectivity Gain: %d, Robustness Gain: %d, Coverage Gain: %d - NOT FOUND', i, param.gain(k,1),param.gain(k,2),param.gain(k,3)));
            continue;
        else
            
            if ~isequal(exist(resultsworkingPath, 'dir'),7)
                mkdir(resultsworkingPath)
            end
            
            disp(sprintf('Network: %d, Connectivity Gain: %d, Robustness Gain: %d, Coverage Gain: %d - Gerenartion Video', i, param.gain(k,1),param.gain(k,2),param.gain(k,3)));
            
            % defining video settings
            
            if strcmp(param.typeMedia,'video')
                outputVideo = VideoWriter(fullfile(resultsworkingPath,strcat(param.videoName, '_', num2str(i), '_',labelP,'.mp4')), 'MPEG-4');
                outputVideo.FrameRate = 20;
                output.Quality=100;
                open(outputVideo);
            end
            
            param.idx=param.t0; % for controlling the iteration data to be loaded
            
            nodeF=[];
            
            % if prone failure scenario is on set the node label and its
            % position at the failure time
            
            if param.attack
                for ii=param.failureT % for each failure performed
                    if ii <= param.tf-2 % if is not the last simulation time
                        data1=load(strcat(workingPath,param.fileId,'_',num2str(ii),'.mat'));
                        data2=load(strcat(workingPath,param.fileId,'_',num2str(ii+1),'.mat'));
                        sizeIdx1=floor(size(data1.iterationData,2)/2);
                        sizeIdx2=floor(size(data2.iterationData,2)/2);
                        labels1=transpose(data1.iterationData(1,2:sizeIdx1+1));
                        labels2=transpose(data2.iterationData(1,2:sizeIdx2+1));
                        nodeF=[nodeF; [ii setdiff(labels1,labels2)]];
                    end
                end
            end
            
            node2obstacle=[];
            
            while param.idx < param.tf
                
                % loading file that contains all the data regarding the
                % simulation time generated during the experiments
                
                load(strcat(workingPath,param.fileId,'_',num2str(param.idx),'.mat'));
                nodeToFail=[];  % to persist nodes which failed
                    
                if param.attack % if fault-prone scenario is on
                    nodeTpos=(find(nodeF(:,1)==param.idx)); %verifies if some node failed at this iteration
                    if ~isempty(find(nodeF(:,1)==param.idx))
                        nodeToFail=nodeF(nodeTpos,2); % tracking node to fail
                    end
                end
                
                if ~isempty(iterationData) % if there is iteration data
                    
                    % first and second rows of iterationData contain
                    % information about the network - Notice that param.idx
                    % are those estimation points, where the network
                    % properties, vulnerability and coverage are evaluated.
                    % For instance, for param.idx = 10 there are 1 sec of
                    % simulation generated by the ODE exectuion
                    
                    sizeIdx=floor(size(iterationData,2)/2); % number of iteration performed
                    
                    labels=transpose(iterationData(1,2:sizeIdx+1)); % node labels
                    
                    vulnerability=transpose(iterationData(2,2:sizeIdx+1)); % node vulnerability status
                    
                    robustness=iterationData(2,1); % robustness level data
                    
                    % third to last rows contain each simulation time
                    % (continuos) generated by the ODE execution
                    
                    for l=3:size(iterationData,1)
                        
                        disp(sprintf('Network: %d, Gain: %s, ti: %d, t: %f', i, labelP,param.idx,iterationData(l,1)));
                        
                        timeI=repmat(iterationData(l,1),sizeIdx,1);
                        
                        data= [timeI labels  transpose(iterationData(l,2:sizeIdx+1)) transpose(iterationData(l,sizeIdx+2:(sizeIdx*2)+1)) vulnerability];
                        posFailure=[];
                        if ~isempty(nodeToFail)
                            pos=find(data(:,2)==nodeToFail);
                            posFailure=data(pos,3:4);
                        end
                        
                        % generate video for this time iteration
                        
                        currentFrame = plot_snapshot(data,robustness,posFailure,node2obstacle,param);
                        
                        param.pictureName=strcat(resultsworkingPath,param.fileName,'_Fig_',num2str(param.idx),'_',num2str(l-3),'_',labelP);
                                                
                        if strcmp(param.typeMedia,'video')
                            writeVideo(outputVideo, currentFrame);
                        else
                            export_fig(param.pictureName, '-eps','-png');
                            saveas(gcf,[param.pictureName,'.fig'],'fig');
                        end
                        
                        close all;
                    end
                    
                    if ~isempty(posFailure) % the node wich failure became an obstacle
                        node2obstacle=[node2obstacle; posFailure(size(posFailure,1),:)];
                    end
                end
                
                param.idx = param.idx + param.ti; % incrementing simulation for recovering iteration data
                
            end
            if strcmp(param.typeMedia,'video')
                close(outputVideo); % close media
            end
        end
    end
end
