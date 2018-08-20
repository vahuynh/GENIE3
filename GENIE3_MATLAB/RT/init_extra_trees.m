
function rtensparam=init_extra_trees(k)
  
  rtensparam.nbterms=100;
  rtensparam.bootstrap=0;
  rtparam.nmin=1;
  rtparam.varmin=0;
  rtparam.savepred=1;
  rtparam.bestfirst=0;
  rtparam.extratrees=1;
  if nargin>0
    rtparam.adjustdefaultk=0;
    rtparam.extratreesk=k;
  else
    rtparam.adjustdefaultk=1;
  end 
  rtparam.savepred=1;
  rtensparam.rtparam=rtparam;