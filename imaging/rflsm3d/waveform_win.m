function [win] = waveform_win(data,t,tstart,tend,td)
    % data:
    % t: time
    % tstart:
    % tend: 
    % td: begin to decrase (sec)
    
    win=zeros(size(data,1),1);
    
    dt = t(2) - t(1);
    Len = fix((tend - tstart)/dt)+1;
    
    % tend = tstart + wl;
    dlen = (td - tstart)/dt + 1;
    
    ratio = 2*dlen/Len;
    
    w = tukeywin(Len,ratio);
    
    n1 = fix((tstart-t(1))/dt)+1;
    
    w1 = n1;
    w2 = n1 + Len-1;
    
    win(w1:w2) = w;


end