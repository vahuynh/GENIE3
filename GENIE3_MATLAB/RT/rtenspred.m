
function [YTS]=rtenspred(treeensemble,XTSA)

% Make predictions with an ensemble of (multiple output) regression trees
  
% inputs:
%   tree: a tree output by the function rtenslearn_c
%   XTS: inputs for the test cases
%   YLS: outputs for the learning sample cases
% Output:
%   YTS: Predictions for the test cases

  verbose=0;
  
  global XTS;
  global tree;
  
  XTS=XTSA;

  Nts=size(XTS,1);
  T=length(treeensemble.trees);
  
  YTS=0;
  
  for t=1:T
    if (verbose)
      fprintf('t=%d\n',t);
    end
    tree=treeensemble.trees(t);
    YTS=YTS+tree.weight*rtpred();
  end
  
function [YTS]=rtpred()

% Test a multiple output regression tree
  
% inputs:
%   tree: a tree output by the function rtenslearn_c
%   YLS: outputs for the learning sample cases
%   XTS: inputs for the test cases
% Output:
%   YTS: output predictions for the test cases

  global assignednodets
  global XTS
  global tree
  
  Nts=size(XTS,1);
  assignednodets=zeros(Nts,1);
  
  verbose=0;
    
  YTS=zeros(Nts,size(tree.predictions,2));
  
  if (verbose)
    fprintf('computation of indexes\n');
  end
  
  getleafts(1,1:Nts);
  
  if (verbose)
    fprintf('computation of predictions\n');
  end  

  for i=1:Nts
    YTS(i,:)=tree.predictions(tree.indexprediction(assignednodets(i)),:);
  end    
  
function getleafts(currentnode,currentrows)

  global assignednodets
  global XTS
  global tree
  
  testattribute=tree.testattribute(currentnode);
  
  if testattribute==0 % a leaf
    assignednodets(currentrows)=currentnode;
  else
    testthreshold=tree.testthreshold(currentnode);
    leftind=(XTS(currentrows,testattribute)<testthreshold);
    rightind=~leftind;
    getleafts(tree.children(currentnode,1),currentrows(leftind));
    getleafts(tree.children(currentnode,2),currentrows(rightind));
  end
