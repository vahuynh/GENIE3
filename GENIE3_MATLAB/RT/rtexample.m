

% Set the parameters of the method
rtensparam=init_extra_trees() % extremely randomized trees
% rtensparam=init_rf() pour random forests
% rtensparam=init_mart() pour MART ou OKBOOST

% To change the number of trees in the ensemble
rtensparam.nbterms=100;

% To change the parameter K of the extra-trees method
rtensparam.rtparam.adjustdefaultk=0;
rtensparam.rtparam.extratreesk=10;
% or equivalently: rtensparam=init_extra_trees(10)

% Load a dataset (WARNING: the dataset must be converted in single format)
DATA=single(load('friedman1.csv'));

X=DATA(:,1:10); % inputs
Y=DATA(:,11); % output

ls=int32(1:300); % the learning sample is composed of the first 300 objects
                 % (WARNING: int32 format is mandatory)
w=[]; % We could attribute a different weight to each object

% Set the test sample (the last 2000 objects)
XTS=DATA(8001:10000,1:10);
YTS=DATA(8001:10000,11);

% Main function: rtenslearn_c:

% Build an ensemble of trees and test it on the test sample
% Input arguments:
%   - X: the inputs (a single float matrix of size nb_obj x nb_att)
%   - Y: the outputs (a single float matrix of size nb_obj x nb_outputs)
%   - ls: row indexes in X and Y of the learning sample objects (int32 vector)
%   - w: weight vector (when w=[], all objects have unit weights)
%   - rtensparam: parameters of the tree growing algorithm
%   - XTS: inputs for the test sample objects
%   - v: 1 to display informations when building the model
% Output argument:
%   - YPRED: Predictions for the test sample objects

[YPRED1]=rtenslearn_c(X,Y,ls,w,rtensparam,XTS,1);

% Computation of square error:

mean((YPRED1-YTS).^2)

% If no test sample is provided, the function outputs a tree structure and
% variable importances

[tree var_imp]=rtenslearn_c(X,Y,ls,w,rtensparam,1);

% Then, the tree ensemble can be tested on a test sample using the function
% rtenspred.m (but this is much slower than doing it directly from
% the function rtenslearn_c)

YPRED2=rtenspred(tree,XTS);
mean((YPRED2-YTS).^2)

