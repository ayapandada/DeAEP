function [m2c,v2c]=RDC(Y,H,nTx,rou,pc,p,m2,v2)

 %gamma   

v = inv((1/rou)*H'*H + inv(v2));
m = v*((1/rou)*H'*Y + inv(v2)*m2);

  for i=1:nTx
   % Ratio Distribution Computation
    v2c(i) = ((v(i,i))^-1 - (v2(i,i))^-1)^-1;
    m2c(i) = v2c(i)*(v(i,i)^-1*m(i)- (v2(i,i))^-1*m2(i));
  end
% v = inv((1/rou)*H'*H + (1./v2));
% m = v*((1/rou)*H'*Y + (1./v2)*m2);
% 
%  for i=1:nTx
%    % Ratio Distribution Computation
%     v2c(i) = ((v(i,i))^-1 - (v2(i,i))^-1)^-1;
%     m2c(i) = v2c(i)*(v(i,i)^-1*m(i)- (v2(i,i))^-1*m2(i));
%  end


% v = inv((1/rou)*H'*H + v2);
% m = v*((1/rou)*H'*Y + v2*m2);
% 
%   for i=1:nTx
%    % Ratio Distribution Computation
%     v2c(i) = ((v(i,i))^-1 - v2(i,i))^-1;
%     m2c(i) = v2c(i)*(v(i,i)^-1*m(i)- (v2(i,i))^-1*m2(i));
%   end