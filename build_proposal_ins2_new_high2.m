function [m,b,w,S] = build_proposal_ins2_new_high2(X,W,thetaW,thetaE,Mask,idw,S,t,type)

% BUILD_PROPOSAL builds a proposal from a set of support points 'S'.
%   f    target dist
% 	S	 set of support points
%   t    iterations on the chain
%   type type of proposal construction (0/1)

% Parameters for the initialization of the tails, to reduce the dependence of choice of the initial points.
beta= 0.95;	% 0<beta<1
tau= 0.01;	% tau>0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = length(S);			% number of points
b = zeros(1,N-1);       % intercepts
m = zeros(1,N-1);       % slopes
%f_S = eval_target(S);	% evaluate all points
% f_S = f(S);	% evaluate all points
[f_S] = f_estimate_ins2_high2(X,W,thetaW,thetaE,Mask,idw,S); % log function
xa=[-Inf S +Inf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the pieces			        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type 0 for approximations of order 0 (constant pieces)
if type == 0 | type~=1
    b(1)=0;
    %     aux(1,:)=log(f_S(1:end-1));
    %     aux(2,:)=log(f_S(2:end));
    aux(1,:)=f_S(1:end-1);
    aux(2,:)=f_S(2:end);
    b(2:N)=max(aux);
    m=zeros(1,N);
    
    % Type 1 for approximations of order 1 (linear pieces)
elseif type == 1
    m(1)=0;
    %     m(2:N) = zeros(1,N-1);
    m(2:N)=(exp(f_S(1:end-1))-exp(f_S(2:end)))./(S(1:end-1)-S(2:end));
    b(1)=0;
    %     b(2:N) = zeros(1,N-1);
    b(2:N)=(exp(f_S(1:end-1)))- m(2:N).*S(1:end-1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the tail pieces 			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The tails connect with the end pieces
% v(1) = log(f_S(1));
% v(2) = log(f_S(2));
% v(3) = log(f_S(end-1));
% v(4) = log(f_S(end));
v(1) = f_S(1);
v(2) = f_S(2);
v(3) = f_S(end-1);
v(4) = f_S(end);


% Control to avoid numerical problems
pos = find(v==-Inf);
if isempty(pos) == 0
    v(pos) = -50;
end

% first slope and intercept
m(1) = ( v(2) - v(1) ) / ( S(2) - S(1) );
m(1)= m(1)*(1-beta*exp(-tau*t));
b(1) = v(1) - m(1) * S(1);

% end slope and intercept
m(N+1) = ( v(end) - v(end-1) ) / ( S(end) - S(end-1));
m(N+1)= m(N+1)*(1-beta*exp(-tau*t));
b(N+1) = v(end) - m(N+1) * S(end);


% Control for 'suitable' slopes, to get a proper proposal pdf
if m(1) <= 0 | isnan(m(1))
    m(1) = 0.05;
    b(1) = v(1) - m(1) * S(1);
end
if m(N+1) >= 0 | isnan(m(N+1))
    m(N+1) = -0.05;
    b(N+1) = v(4) - m(N+1) * S(end);
end


%%%%%%%%%%%% COMPUTE AREA and weights! %%%%%%%%%%%%%
if type==0 | type~=1
    M = max([m(1)*xa(2)+b(1),m(N+1)*xa(N+1)+b(N+1),b(2:end-1)]);
    area(1)=(1/m(1)).*exp(m(1)*xa(2)+b(1)-M)-(1/m(1)).*exp(m(1)*xa(1)+b(1)-M);
    area(N+1)=(1/m(N+1)).*exp(m(N+1)*xa(N+2)+b(N+1)-M)-(1/m(N+1)).*exp(m(N+1)*xa(N+1)+b(N+1)-M);
    area(2:N)=(xa(3:end-1)-xa(2:end-2)).*exp(b(2:end-1)-M);
    %%%%%%%%%%%%%%%%%
elseif type==1
    %%%%%%%%%%%%%%%%%%
    M = max([m(1)*xa(2)+b(1),m(N+1)*xa(N+1)+b(N+1),f_S]);
    area(1)=(1/m(1)).*exp(m(1)*xa(2)+b(1)-M)-(1/m(1)).*exp(m(1)*xa(1)+b(1)-M);
    area(N+1)=(1/m(N+1)).*exp(m(N+1)*xa(N+2)+b(N+1)-M)-(1/m(N+1)).*exp(m(N+1)*xa(N+1)+b(N+1)-M);
    area(2:N)= exp(f_S(1:end-1)-M).*((xa(2:end-2)+xa(3:end-1)-2*S(1:end-1))./(S(1:end-1)-S(2:end))+2) ...
        - exp(f_S(2:end)-M).*((xa(2:end-2)+xa(3:end-1)-2*S(1:end-1))./(S(1:end-1)-S(2:end)));
    area(2:N) = area(2:N).*(xa(3:end-1)-xa(2:end-2))/2;
end
w=area./sum(area);
if(isnan(w))
    'warning'
end