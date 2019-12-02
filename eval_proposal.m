function  fp=eval_proposal(x,m,b,S,type)
% EVAL_PROPOSAL Evaluate point 'x' on the proposal defined by S,m,b
% 	x   the point to evaluate
% 	S	set of N points
%   m	set of slopes
%   b	set of intercepts
%   type type of proposal construction (0/1)

% 'pos' is the insertion order point for 'x' in 'S', i.e., the interval.
pos=find(sort([S x])==x);

% 'pos' could have > 1 dimensions (if 'x' falls between intervals, i.e., it is already a support point)
% we take the first ...  and take the exponential ... to avoid numerical problems.
if type == 0 | type~=1
    fp=m(pos(1))*x+b(pos(1));
elseif type == 1
    if pos(1)==1 | pos(1)==length(m)
        fp=m(pos(1))*x+b(pos(1));
    else
        fp=m(pos(1))*x+b(pos(1));
    end
end
