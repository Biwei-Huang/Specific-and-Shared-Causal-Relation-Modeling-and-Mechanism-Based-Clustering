function [sample,fp]=Sample_Eval_Proposal_ins(X,W,thetaW,thetaE,idw,m,b,S,type,w)
% S		set of N points

N=length(S);
xa=[-Inf S +Inf];

piece=mnrnd(1,w)*[1:N+1]';

if(type==0)
    if piece==1
        sample = (m(1)*xa(2) + log(rand(1,1)))/m(1);
        fp = m(1)*sample+b(1);
    else
        if piece==N+1
            %             sample = (m(N+1)*xa(N+1) + log(1-rand(1,1)))/m(N+1);
            y=exp(m(N+1)*xa(N+1)+b(N+1))-exp(m(N+1)*xa(N+1)+b(N+1))*rand(1,1);
            if(y==0)
                y = 1e-300;
            end
            sample=(log(y)-b(N+1))/m(N+1);
            fp= m(N+1)*sample+b(N+1);
        else
            sample=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
            fp=b(piece);
        end
    end
end

if(type==1)
    if piece==1
        sample = (m(1)*xa(2) + log(rand(1,1)))/m(1);
        fp = m(1)*sample+b(1);
    else
        if piece==N+1
            sample = (m(N+1)*xa(N+1) + log(1-rand(1,1)))/m(N+1);
            fp= m(N+1)*sample+b(N+1);
        else
            u=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
            v=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
            fs1 = f_estimate_ins(X,W,mu,alpha,idw,xa(piece));
            fs2 = f_estimate_ins(X,W,mu,alpha,idw,xa(piece+1));

            %         auxvar=fs1/(fs1+fs2);
            M = max(fs1,fs2);
            auxvar=exp(fs1-M)/(exp(fs1-M)+exp(fs2-M));
            sample=mnrnd(1,[auxvar 1-auxvar])*[min([u,v]) max([u,v])]';
            fp=m(piece)*sample+b(piece);
        end
    end
end



% 
% 
% piece=mnrnd(1,w)*[1:N+1]';
% 
% if piece==1
%     y=exp(m(1)*xa(2)+b(1))*rand(1,1);
%     sample=(log(y)-b(1))/m(1);
%     %     fp=exp(m(1)*sample+b(1));
%     fp=m(1)*sample+b(1);
%     
% elseif piece==N+1
%     y=exp(m(N+1)*xa(N+1)+b(N+1))-exp(m(N+1)*xa(N+1)+b(N+1))*rand(1,1);
%     sample=(log(y)-b(N+1))/m(N+1);
%     %     fp=exp(m(N+1)*sample+b(N+1));
%     fp=m(N+1)*sample+b(N+1);
% else
%     if type == 0 | type~=1
%         sample=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
%         %         fp=exp(b(piece));
%         fp=b(piece);
%     elseif type == 1
%         %sample=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
%         u=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
%         v=xa(piece)+(xa(piece+1)-xa(piece))*rand(1,1);
%         fs1=f(xa(piece));
%         fs2=f(xa(piece+1));
%         auxvar=fs1/(fs1+fs2);
%         %sample=randsrc(1,1,[min([u,v]) max([u,v]);auxvar 1-auxvar]);
%         sample=mnrnd(1,[auxvar 1-auxvar])*[min([u,v]) max([u,v])]';
%         fp=m(piece)*sample+b(piece);
%         
%     end
% end
% 
% 
% 
