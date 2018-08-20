
function rtensparam=init_rf(k)
  
  rtensparam.nbterms=100;
  rtensparam.bootstrap=1;
  rtparam.nmin=2;
  rtparam.varmin=0;
  rtparam.savepred=1;
  rtparam.bestfirst=0;
  rtparam.rf=1;
  rtparam.extratrees=0;
  if nargin>0
    rtparam.adjustdefaultk=0;
    rtparam.extratreesk=k;
  else
    rtparam.adjustdefaultk=1;
  end 
  rtparam.savepred=1;
  rtensparam.rtparam=rtparam;