%HIT_REATTRIBUTE_SUSPECTED Re-attribute inliers suspected to be mixed to the modes 
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% This is a batch file. The structure inliers is also updated if any
% reattribution took place.
% 
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% (none)
%
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% (none)

% Copyright is with the following author:
%
% (C) 2005 Giancarlo Ferrari Trecate,
%         giancarlo.ferrari@unipv.it
% -------------------------------------------------------------------------
% Legal note:
%     This program is free software; you can redistribute it and/or
%     modify it under the terms of the GNU General Public
%     License as published by the Free Software Foundation; either
%     version 2.1 of the License, or (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public
%     License along with this library; if not, write to the
%     Free Software Foundation, Inc.,
%     59 Temple Place, Suite 330,
%     Boston, MA  02111-1307  USA
%
%
% -------------------------------------------------------------------------

% Create the matrix exchange.
% Each row of exchange is a triplet (index of point , old set F , new set F)
%
% A point (X_j,y_j) is suspected to be mixed if it belongs to the i-th mode dataset
% but all other X-points in the j-th LD (i.e.the LD labeled with (X_j,y_i))
% are not in the same mode dataset
%
% Let I be the list of modes adjacent to i (including i). Then the point is
% assigned to the mode in I that minimizes the error

exchange=[];
for i=1:s, %cycle over the sets F
    np=length(F.pos{i});
    for k=1:np,
        index=F.pos{i}(k); %store the position of the point
        test_incluster=0;
        for j=1:idpar.c,
            if ~test_incluster
                % index of the point to be tested
                % note that the first element of the vector LDs.pos{i} is i
                % i.e. the point labeling the LD
                test=LDs.pos{index}(j); %index of the point to be tested
                % II is nonempty only if the point test belongs to F.pos{i}
                II=[];
                II=find(~(F.pos{i}-test));
                if isempty(II) % no point found -> possible outlier !
                    test_incluster=1; %avoid further testing in the same cluster
                    %
                    %TESTING COMMANDS:
                    % index
                    % LDs.X{index}
                    % LDs.pos{index}
                    % F.pos{i}
                    % idmodes1.par{:}
                    % END OF TESTING COMMANDS:
                    %
                    % find the most likely mode it belongs to among
                    % the modes adjacent to the one the points belonged to
                    distances=[];
                    [p,q]=find(adjacences==i);
                    % swap the columns where the i-th mode has been found
                    % i.e. if q(jjj)=1, set q(jjj)=2 and vice-cersa
                    q=q-1;
                    q=~q;
                    q=q+1;
                    % build the list of candidate modes the point can
                    % belong to
                    % this list contains the i-th mode and its adjacent
                    % modes
                    candidate_modes=[i];
                    for jjj=1:length(p),
                        candidate_modes=[candidate_modes;adjacences(p(jjj),q(jjj))];
                    end
                    for jjj=1:length(candidate_modes)
                        distances=[distances;(F.y{i}(k)-[F.X{i}(k,:) 1]*idmodes1.par{candidate_modes(jjj)}).^2];
                    end
                    [tempp,newi]=min(distances); %newi store the index of the most likely model
                    exchange=[exchange;[index,i,candidate_modes(newi)]]; % add the point to the exchange list
                end %ends for ii

            end %ends if isempty
        end
    end
end
% end of creation of exchange

% cut from exchange points that must be reattributed to the same set F
tempp=size(exchange,1);
temp=exchange;
exchange=[];
for i=1:tempp
    if temp(i,2)~=temp(i,3);
        exchange=[exchange;temp(i,:)];
    end
end
stat_reattr=size(exchange,1)/length(inliers.pos);
% Re-attribute the points in inliers.class
if ~isempty(exchange) %if there are points to be reattributed
    fprintf('\n Some datapoints have been reattributed: modes must be re-computed.');
    fprintf('\n Percentage of points reattributed: %2.0f %%.',stat_reattr*100);
    % re-attribute the points
    tempp=size(exchange,1);
    for i=1:tempp, %cycle over exchange
        npoint=exchange(i,1);
        oldF=exchange(i,2);
        newF=exchange(i,3);
        index_inl=find(~(inliers.pos-npoint));% find the index in inliers of the point to be removed
        % update inliers.class
        inliers.class(index_inl)=newF;
    end
else
  fprintf('\n No datapoint has been reattributed\n');
end

