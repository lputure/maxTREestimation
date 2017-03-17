clc,close all,clear all

X = [-1 0; 1 0; 0 1]';
N = 1e4;
P = 2*rand(2,N)-1.5;
tic
[Y,V,F] = myLBWMEC(X,P);
toc
[c,r] = minboundcircle(X(1,:),X(2,:));

figure, hold on
for ii = 1:size(X,2)
    plot(X(1,ii),X(2,ii),'rx')
end
t = 0:0.01:2*pi;
xx = c(1) + r*cos(t);
yy = c(2) + r*sin(t);
plot(xx,yy)

id = 4;
plot(P(1,id),P(2,id),'ro')
plot(Y(1,id),Y(2,id),'bo')

% figure, hold on
% cnt = 0;
% for ii = 1:N
%     if(norm(Y(:,ii))>eps && F(ii)<0.5)
%         cnt = cnt + 1;
%         plot(Y(1,ii),Y(2,ii),'r.')
%     end
% end
% axis([-5 5 -5 5])
% %% some validattion
% figure, hold on
% for ii = 1:N
%     if(F(ii)>0)
%         plot(P(1,ii),P(2,ii),'r.')
%     else
%         if(norm(Y(:,ii))>eps)
%             plot(P(1,ii),P(2,ii),'b.')
%         else
%             plot(P(1,ii),P(2,ii),'y.')
%         end
%     end            
% end