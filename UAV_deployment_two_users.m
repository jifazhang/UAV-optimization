%% fixed wing UAV trajectory optimization
%% only optimize trajectory of UAV, two users case
 

clc,clear;
close all;

 %% paramters
qs = [0,0]'; % start position
qe = [100,100]';% end position
s1 = [-10,80]'; % BS station
s2 = [120,60]'; % BS station
dis  = norm(qe-qs); % distance between qs and qe
T =48; % times slots
delta_T = 0.5;% delta T
vmax =dis/(T*delta_T)*8; % maximial velocity
vmin = 3;% minimial velocity
amax = 5;% max accelerated velocity
noise_power = 1e-7;%50dbm
beta0 = 1e-3;% L0 pathloss
beta_hat=beta0/noise_power;

iter_num = 20; % number of iterations
q = zeros(T,2,iter_num+1);% positions variables
v = zeros(T,2,iter_num+1);% velocity variables
v(:,:,1) = sqrt(1/2*dis/(T*delta_T));% 
a = zeros(T,2,iter_num+1);% accelerated velocity
epsilon = 1e-3; % stop criterion
H =20;% hegith  
u = zeros(T,iter_num+1);% slack variables
w = zeros(T,iter_num+1);% slack variables
% straight line flying
q(:,1,1)=linspace(qs(1),qe(1),T);
q(:,2,1)=linspace(qs(2),qe(2),T);
val = zeros(1,iter_num+1);
for ii = 1:T
u(ii,1)= norm(q(ii,:,1)-s1)^2;
end
for ii = 1:T
w(ii,1)= norm(q(ii,:,1)-s2)^2;
end
%% optimization
for iter=1:iter_num
    
cvx_begin    quiet
variable q1(T,2)
variable u1(T,1)  nonnegative
variable w1(T,1)  nonnegative
variable v1(T,2,1)  
variable a1(T,2,1)   
rate = 0;
% calculate objective functions
for ii = 1:T
    rate = rate + log2(1+beta_hat/(u(ii,iter)+H^2))-1/log(2)*beta_hat/((u(ii,iter)+H^2)^2+beta_hat*(u(ii,iter)+H^2))*(u1(ii)-u(ii,iter))...
        +log2(1+beta_hat/(w(ii,iter)+H^2))-1/log(2)*beta_hat/((w(ii,iter)+H^2)^2+beta_hat*(w(ii,iter)+H^2))*(w1(ii)-w(ii,iter));
    
end
maximize (rate)
subject to
% start and end potisiton constraints
q1(1,1)==qs(1,1);
q1(1,2)==qs(2,1);

q1(T,1)==qe(1,1);
q1(T,2)==qe(2,1);
% start and end velocity constraints
v1(1,:) == [vmin*2,vmin*2]
v1(T,:) == [vmin*2,vmin*2]
% accelerated velocity constraint
norm(a1(1,:))<=amax
% slack variables constraints
norm([1-u1(1),2*(q1(1,:)-s1')])<=1+u1(1)
norm([1-w1(1),2*(q1(1,:)-s2')])<=1+w1(1)
for kk=2:T
    % maximal and minimal velocity constraints
    norm(v1(kk,:))<=vmax
%  
    v(kk,:,iter)*v(kk,:,iter)'+2*v(kk,:,iter)*(v1(kk,:)'-v(kk,:,iter)')>=vmin^2
%     % maximial accelerated constrains  
    norm(a1(kk,:))<=amax
    % position constrains
    q1(kk,:)==q1(kk-1,:)+v1(kk,:)*delta_T+1/2*a1(kk,:)*delta_T^2
    v1(kk,:)==v1(kk-1,:)+a1(kk,:)*delta_T
%     slack constrains
    norm([1-u1(kk),2*(q1(kk,:)-s1')])<=1+u1(kk)
    norm([1-w1(kk),2*(q1(kk,:)-s2')])<=1+w1(kk)
end
cvx_end
% save optimized results
q(:,:,iter+1)=q1;
v(:,:,iter+1)=v1;
a(:,:,iter+1)=a1;
u(:,iter+1)=u1;
w(:,iter+1)=w1;
val(iter+1)=cvx_optval;
disp(iter)
if abs(val(iter+1)-val(iter))<=epsilon
    break
end
end
%% plot trajectory
figure
hold on
%plot trajectory
plot(q(:,1,1),q(:,2,1),'r-.o','linewidth',1.5);
plot(q(:,1,iter+1),q(:,2,iter+1),'b-d','linewidth',1.5);
plot(s1(1),s1(2),'s','MarkerFaceColor', 'g','markersize',10)
plot(s2(1),s2(2),'s','MarkerFaceColor', 'g','markersize',10)
% plot positions
plot(qs(1),qs(2),'s','MarkerFaceColor', 'g','markersize',10)
plot(qe(1),qe(2),'s','MarkerFaceColor', 'g','markersize',10)
% add text
text(s1(1)+5,s1(2)+5,'User1')
text(s2(1)-17,s2(2)-5,'User2')
text(qs(1)-5,qs(2)-5,'Inital')
text(qe(1)-5,qe(2)+5,'Final')
grid on
% lim  & label & box
xlim([qs(1)-15,qe(1)+25])
ylim([qs(2)-10,qe(2)+10])
legend('straight flight','Optimized trajectory','Location','southeast')
box on
xlabel('x(m)')
ylabel('y(m)')
title(['T=',num2str(T*delta_T),'s'])
%% plot convergence

figure
plot(val(1:iter+1),'b-o','linewidth',1.5);
grid on

xlabel('iteration number')
ylabel('objective value')
box on
set(gca,'fontsize',10)
title('Convergence curve')
%% 
figure
hold on
for ii = 1:iter
plot(q(:,1,ii),q(:,2,ii),'linewidth',1.5);   
    
pause(1)    
end
grid on
