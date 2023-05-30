function [intMatOut] = invIntMatrixAuto(intMatIn)
%INVMATRIX Summary of this function goes here
%   Detailed explanation goes here
intMatIn_int = intMatIn.int;
Acenter = (intMatIn_int.sup + intMatIn_int.inf)/2;
Adelta = (intMatIn_int.sup - intMatIn_int.inf)/2;


M = pinv(eye(size(pinv(Acenter)*Adelta)) - abs(pinv(Acenter))*Adelta);
mu = diag(M);
Tmu = diag(mu);
Tv =  pinv(2*Tmu - eye(size(Tmu)));
Blower = -M*abs(pinv(Acenter)) +Tmu*( pinv(Acenter) + abs(pinv(Acenter))) ;
Bupper = M*abs(pinv(Acenter)) +Tmu* (pinv(Acenter) - abs(pinv(Acenter))) ;

Blowerlower =min(Blower,Tv*Blower);
Bupperupper =max(Bupper,Tv*Bupper);

intMatOut = intervalMatrix((Bupperupper+Blowerlower)/2,(Bupperupper-Blowerlower)/2);

% M = pinv(eye(size(Acenter)) - abs(pinv(Acenter))*Adelta);
% mu = diag(M);
% Tmu = diag(mu);
% Tv =  pinv(2*Tmu - eye(size(Acenter)));
% Blower = -M*abs(pinv(Acenter)) +Tmu*( pinv(Acenter) + abs(pinv(Acenter))) ;
% Bupper = M*abs(pinv(Acenter)) +Tmu* (pinv(Acenter) - abs(pinv(Acenter))) ;
% 
% Blowerlower =min(Blower,Tv*Blower);
% Bupperupper =max(Bupper,Tv*Bupper);
end

