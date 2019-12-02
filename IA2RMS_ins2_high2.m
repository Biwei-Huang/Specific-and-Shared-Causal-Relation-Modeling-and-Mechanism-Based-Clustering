
function wij = IA2RMS_ins2_high2(X,W,thetaW,thetaE,Mask,iw,jw,N,type)
% N: number of samples we required
% type: 0/1


% if(iw==1 & jw==2)
%     S = [-0.02 0 0.02];
% else
%     S = [0,1,2];
% end
S = [-3,-2,-1,0,1,2,3,4];
S = sort(S); % Always sort support points

x = zeros(1,N);
if isempty(find(S==0))==0
    x(1)=0.02;
end

sign = 0;
alpha_before=1;
count=0;
k=2;
NewP_Yes=1;

while k<=N
    if NewP_Yes==1
        % BUILD PROPOSAL from a set of support points S
        %         [m,b,w] = build_proposal(f,S,k-1,type);
        [m,b,w,S] = build_proposal_ins2_new_high2(X,W,thetaW,thetaE,Mask,[iw,jw],S,k-1,type);
    end
    NewP_Yes=0;
    
    % SAMPLING AND EVALUATE POPOSAL: draw x'from and evaluate the proposal
    %         [x_prop,fp_now]=Sample_Eval_Proposal(f,m,b,S,type,w);
    [x_prop,fp_now]=Sample_Eval_Proposal_ins(X,W,thetaW,thetaE,[iw,jw],m,b,S,type,w);
    
    
    
    if count==0
        % EVALUATE x' first time
        fp_prev=eval_proposal(x(k-1),m,b,S,type);
        %                 ft_prev= f(x(k-1));
        [ft_prev] = f_estimate_ins2_high2(X,W,thetaW,thetaE,Mask,[iw,jw],x(k-1));
    end
    %fp_now=eval_proposal(x_prop,m,b,S,type);
    % EVALUATE x'
    %      ft_now= f(x_prop);
    [ft_now] = f_estimate_ins2_high2(X,W,thetaW,thetaE,Mask,[iw,jw],x_prop);
    
    
    u=rand(1,1);
    
    %     alpha1=ft_now/fp_now;
    alpha1=exp(ft_now-fp_now);
    
    if u > alpha1  %%%% RS test
        % reject step 1
        S=[S x_prop];
        S=sort(S);
        count=count+1;
        NewP_Yes=1;
    else
        % accept
        q_prev=min([ft_prev, fp_prev]);
        q_now=min([ft_now, fp_now]);
        %         rho=(ft_now*q_prev)/(ft_prev*q_now);
        rho = exp(ft_now - ft_prev + q_prev-q_now);
        alpha2=min([1,rho]);
        u2=rand(1,1);
        
        if u2<=alpha2 & abs(x_prop)<10 %%%% MH test  CHANGE to make it stable!!!!!!!!!!
            % accept
            x(k) = x_prop;
            y_aux=x(k-1);
            %%%%%%%%
            %%%%%%%%
            alpha3=1/alpha_before;
            %%%%%%%%
            %%%%%%%%
            fp_prev=fp_now;
            ft_prev=ft_now;
            
            sign = 1;
            
        else
            % reject MH test
            x(k) = x(k-1);
            y_aux=x_prop;
            %%%%%%%%
            %%%%%%%%
            alpha3=1/alpha1;
            %%%%%%%%
            %%%%%%%%
        end
        
        %%%% second control of IA2RMS
        alpha_before=alpha1;
        u3=rand(1,1);
        if u3> alpha3 & k~=N %%%% second condition just for plotting
            S=[S y_aux];
            S=sort(S);
            NewP_Yes=1;
        end
        
        k=k+1;
        count=count+1;
    end
%     if(~sign & k>N)
%         N=N+1;
%     end
    %     if(abs(x_prop)<0.1 & iw==2 & jw==1)
    % %         x(k) = -rand;
    %        display('warning')
    %     end
    if(count>50)
        count
        x(k) = rand;
        break;
    end
end
wij = x(end);
% figure,hist(x,50)
