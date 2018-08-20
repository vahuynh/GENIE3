
function [VI]=compute_rtens_variable_importance(treeensemble,XTS,YTS)

% Compute variable importances from an ensemble of (multiple output)
% regression trees.
  
% inputs:
%   treeensemble: a tree output by the function rtenslearn_c
%   XTS: the input data on which to estimate the importances
%   YTS: the output data on which to estimate the importance (XTS and YTS
%        have the same number of rows)
% Output:
%   VI: A matrix whose elements ij represents the importance of feature i for
%       output j
  
  Nts=size(XTS,1);
  Natt=size(XTS,2);
  Nout=size(YTS,2);
  VI=zeros(Natt,Nout);
  
  T=size(treeensemble.trees,1);
  
  for t=1:T
    
    VI=VI+compute_rt_variable_importance(treeensemble.trees(t),XTS,YTS);
    fprintf('.');
    
  end
  
  fprintf('\n');
  
  VI=VI/T;
  
function [VI]=compute_rt_variable_importance(tree,XTS,YTS)
  
  Nts=size(XTS,1);
  Natt=size(XTS,2);
  Nout=size(YTS,2);
  VI=zeros(Natt,Nout);
  
  totalvar=zeros(tree.nodenumber,Nout);
  totalweight=zeros(tree.nodenumber,1);
  assignednode=ones(Nts,1);
  
  opennodes = zeros(tree.nodenumber,1);
  opennodes(1)=1;
  indexopennodes=1;
  nextunusednode=2;
    
  totalvar(1,:)=Nts*var(YTS,1,1);
  
  while(indexopennodes > 0)
    
    currentnode=opennodes(indexopennodes);
    
    noderows = find(assignednode==currentnode);
    Nnode = length(noderows);
    testattribute=tree.testattribute(currentnode);
    
    if (testattribute==0 || Nnode==0) % une feuille, on ne fait rien
      indexopennodes=indexopennodes-1;
    else
      % split the TS
      x=XTS(noderows,tree.testattribute(currentnode));
      leftind=XTS(noderows,tree.testattribute(currentnode))<tree.testthreshold(currentnode);
      rightind=~leftind;
      leftnode=tree.children(currentnode,1);
      rightnode=tree.children(currentnode,2);
      assignednode(noderows(leftind)) = leftnode;
      assignednode(noderows(rightind)) = rightnode;
      
      % compute left/right statistics
      Nleft=sum(leftind);
      Nright=sum(rightind);
      
      totalvar(leftnode,:)=Nleft*var(YTS(noderows(leftind),:),1,1);
      totalvar(rightnode,:)=Nright*var(YTS(noderows(rightind),:),1,1);
      
      if (Nleft~=0 && Nright~=0)
        
        % add the score to the variable that is tested:
        TV=(totalvar(currentnode,:)-totalvar(leftnode,:)-totalvar(rightnode,:));
      
        VI(testattribute,1:Nout)=VI(testattribute,1:Nout)+TV;
      
        % put the nodes on the stack
      end
      indexopennodes=indexopennodes-1;
      if (Nleft~=0)
        indexopennodes=indexopennodes+1;
        opennodes(indexopennodes)=leftnode;
      end
      if (Nright~=0)
        indexopennodes=indexopennodes+1;
        opennodes(indexopennodes)=rightnode;
      end      
    end
    
  end
