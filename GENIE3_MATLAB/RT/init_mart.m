
function rtensparam=init_mart(compl,mu)
  
  rtensparam.nbterms=500;
  rtensparam.mart=1;

  if nargin>1    
      rtensparam.martmu=mu;
  else
    rtensparam.martmu=0.2;
  end
  
  rtensparam.bootstrap=0;
  rtparam.nmin=1;
  rtparam.varmin=0;
  rtparam.savepred=1;
  rtparam.bestfirst=1;

  if nargin>0
    rtparam.maxnbsplits=compl;
  else
    rtparam.maxnbsplits=5;
  end  

  rtparam.extratrees=0;
  rtparam.savepred=1;
  rtensparam.rtparam=rtparam;
