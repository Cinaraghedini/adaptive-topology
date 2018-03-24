function [final_models,final_Variances,Weight_primal]=hit_compute_disc_modes(F,ndata,weig)
%HIT_COMPUTE_DISC_MODES Estimate the mode PVs in discontinuous PWA models 
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% [final_models,final_Variances,final_Variances_denormaliz,Weight_primal]=...
% hit_compute_disc_modes(F,ndata,weig)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% F: structure containing information about the mode datasets, see
% hit_regression.m
% ndata: total number of datapoints.
% weigh(i): weight of the i-th datapoint used in weighted LS for estimating
% the mode PVs.
% -------------------------------------------------------------------------
% OUTPUT                                                                                                   
% -------------------------------------------------------------------------
% final_models(i,:) estimted PV of the i-th mode. 
% final_Variances{i}: variance of the i-th mode PV.
% Weight_primal(i): if present returns the weight associated to the i-th
% datapoint that should be used in building LDs. This is an experimental
% feature currently not used for building LDs. 


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

m=size(F.X,1);
ndim=size(F.X{1}(1,:),2);


final_models=[]; 


final_Variances=cell(m,1);
final_Variances_denormaliz=cell(m,1);

% Weight_primal: the i-th row is a weight associated with the
%  				  i-th datapoint
%
% This vector is used only if a repeated primal-extended dual
% iteration is performed
% The use of such weights to estimate the local models in the data space
% should improve their identification against outliers in the clusters 
% in the data-space.
% 
% However the current rule for assigning such weights (the inverse of
% the distance from the
% center of the cluster of the correspnding vectors in the ext. dual)
% is not really effective .... some better strategy is needed !!
% 
% On the other hand this 'importance' measure works well for assigning the
% weights to be used in the identification of the final models !!



Weight_primal=zeros(ndata,1);

for i=1:m
   nC=size(F.X{i},1);
   if nC <= ndim+1
      fprintf('BE CAREFUL the number of the data attributed to mode %g are less than the parameters \n',i);
   end
   Phitemp=[F.X{i},ones(nC,1)];
   

   % if weig(i) is high a great importance to the sample should be assigned !
   % in Clust_primal.m Weight_primal is the variance of the noise

	 % extract from weig the elements corresponding to the i-th cluster
   wtemp=weig(F.pos{i}); 
   
   % Check if Weight_primal is required as output
   if nargout>4
  		 Weight_primal(F.pos{i})=wtemp;
   end

   
   %  uncomment the following line to NOT use the weights in WLS
   %  i.e. set them equal to 1
   %  W=eye(nC); 
   
   %Weighted LS formula
   [temp,stdx,mse,S]=hit_lscov(Phitemp,F.y{i},wtemp); 
   final_models=[final_models;temp'];
   final_Variances{i}=S;
   if ~all(all(isfinite(final_Variances{i}))) 
        fprintf('\n ======================================================== ');
        fprintf('\n ERROR IN ESTIMATING THE MODE PVs: the covariance'); 
        fprintf('\n of some mode do not have finite entries.');
        fprintf('\n A possible cause is that data are almost noiseless.');
        fprintf('\n Rough remedy: add a little noise to output data.');
        fprintf('\n hit_regression.m ends here. ');
        fprintf('\n ======================================================== \n');
        error('hit_regression:end','hit_regression ends here.')
    end
end

% % De-normalize the final model to obtain the coefficients 
% % in the data space 
% % -------------------------------------------------------------
% 
% % final_models_denormaliz is a matrix: the i-th row are the coefficients
% % of the i-th model in the data-space
% 
% final_models_denormaliz=[];
% 
%  % For the de-normalization formulas, see my handouts ... 
% for i=1:m
%    temp=final_models(i,:);
%    tempx=temp(1:ndim);
%    tempconst=temp(ndim+1:ndim+1);
%    stdy=inv(istdy);
%    tempconst=-istdy*tempx*istdX*mX'+istdy*tempconst+my;
%    tempx=stdy*tempx*istdX;
%    Trans_Variance=[stdy*istdX zeros(size(istdX,1),1);-mX*istdX*stdy stdy];
%    final_Variances_denormaliz{i}=Trans_Variance*final_Variances{i}*Trans_Variance';
% final_models_denormaliz(i,:)=[tempx tempconst];
%    end
% final_models_denormaliz
