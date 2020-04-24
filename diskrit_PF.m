function [xvec,x_estvec]=diskrit_PF(A_bar,input1,input2,t,xo,disktrue,disknomi,p1,w,z,num_members,eps,stdnom)

xvec      = [];
x_estvec  = [];
C         = [0 1 0 0 0 0 0 0 0 0];

for i = 1:numel(t)-1
    disp(['Iterasi ke ',num2str(i)])
    time = [0 (t(i+1)-t(i))];
    if i == 1
        Pxxmin = zeros(p1,p1,num_members);
        xotr = xo;
        x_estmin = zeros(num_members,10)';
    end
    x_tr = lsim(disktrue,input1(i:i+1,:),time,xotr)+ w.*randn(2,p1);
    xotr = x_tr(end,:);
 
    for j=1:num_members
        W = w.*randn(2,p1);                          %create process noise
        Z = z.*randn(1,1);                          %create measurement noise
        x_est0 = lsim(disknomi,input2(i:i+1,:),time,x_estmin(:,j))+ W;
        x_est(j,:)= x_est0(end,:)+(A_bar*(normrnd(0,stdnom)*x_estmin(:,j)))';      %forecast state
        y(j,:)= C*xotr' + Z;                 %make measurement
        y_for(j,:)=C*x_est(j,:)';              %forecast measurement
    end
    x_estbar=mean(x_est,1);
    ybar=mean(y,1);
    y_forbar=mean(y_for,1);
    
    for j=1:p1
        Ex(:,j)= x_est(:,j)-x_estbar(j);
    end
    matm = eye(p1)+ [0 0 0 0 0 0 0 0 0 0;
                     0.9 0 0 0 0 0 0 0 0 0;
                     0.8 0.9 0 0 0 0 0 0 0 0;
                     0.7 0.8 0.9 0 0 0 0 0 0 0;
                     0.6 0.7 0.8 0.9 0 0 0 0 0 0;
                     0.5 0.6 0.7 0.8 0.9 0 0 0 0 0;
                     0.4 0.5 0.6 0.7 0.8 0.9 0 0 0 0;
                     0.3 0.4 0.5 0.6 0.7 0.8 0.9 0 0 0;
                     0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0 0;
                     0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0];
    for k=1:num_members
        Pxx(:,:,k) = (Ex(k,:)'*Ex(k,:).*matm)+(A_bar*((1+eps)*x_estmin(:,k)*x_estmin(:,k)'+(1+eps^-1)*Pxxmin(:,:,k))*A_bar');
        K(:,k)     = sum(Pxx(:,:,k)*C'*inv((((C*Pxx(:,:,k)*C'+(z.^2))'))),2);
    end

    perbaharui = y-y_for;
    for par = 1:num_members
        x_est(par,:) = x_est(par,:)+(K(:,par)*perbaharui(par))';
    end
    x_estbar=mean(x_est,1);
    % menghitung varian corrected state
    for j=1:p1
        Exmin(:,j)= x_est(:,j)-x_estbar(j);
    end
    for k=1:num_members
        Pxxmin(:,:,k) = Exmin(k,:)'*Exmin(k,:).*eye(p1);
    end
    x_estmin = x_est';
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    xvec = [xvec; xotr];
    x_estvec = [x_estvec; x_estbar];
    
end