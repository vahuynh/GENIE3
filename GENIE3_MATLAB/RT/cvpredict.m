
function [YPRED]=cvpredict(X,Y,cvparam,bl_learn,bl_param,verbose)
  
% Test by cross-validation

% parametres:
%  bl_learn: the learning function (e.g. 'rtenslearn_c')
%  bl_param: the parameters (for example the output of init_extra_trees())
%  cvparam: cross-validation parameters
%     cvparam.nbfolds: number of folds
%     cvparam.initseed: a random seed, to reproduce always the same splitting
%     cvparam.verbosebl: 1 to display messages

  N=size(X,1);
  nbfolds=cvparam.nbfolds;
  
  try
    rand('twister',cvparam.initseed);
  catch
  end

  try
    cvparam.verbosebl;
  catch
    cvparam.verbosebl=0;
  end

  ro=randperm(N);
  
  cvsize=round(N/nbfolds);
  
  YPRED=zeros(N,size(Y,2));
  
  for t=1:nbfolds
    
    if (verbose)
      fprintf('Fold %d...\n',t);
    end
    
    if (t==nbfolds) % on est au dernier, on prend tout ce qui reste
      cvts=int32(ro(((t-1)*cvsize+1):N));
      cvls=int32(ro(1:((t-1)*cvsize)));
    else
      cvts=int32(ro(((t-1)*cvsize+1):(t*cvsize)));
      cvls=int32(ro([[1:((t-1)*cvsize)],[(t*cvsize+1):N]]));
    end
    
    if (verbose)
      fprintf('Size LS=%d, TS=%d...\n',length(cvls),length(cvts));
    end   
    
    YPRED(cvts,:)=feval(bl_learn,X,Y,cvls,[],bl_param,X(cvts,:),cvparam.verbosebl);
      
  end
