function [gamma,L]= compLipConst(fun,U,R0,steps,initpoints,dim_x)
%compute Lipschitz constant

totalsamples = initpoints*steps;
%input random sample points
for i = 1:totalsamples
    u(:,i) = randPointExtreme(U);
end

%get state trajectories
index = 1;
for j = 1:dim_x:initpoints*dim_x
    x_free(j:j+dim_x-1,1) = randPoint(R0);
    for i=1:steps
        x_free(j:j+dim_x-1,i+1) = fun(x_free(j:j+dim_x-1,i),u(:,index));
        index=index+1;
    end
end


%combine trajectories
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1        
        x_free_vec_1(:,index_1) = x_free(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_free_vec_0(:,index_0) = x_free(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end



L(1) = [0];
gamma(1) = 0;
idx=1;
for i=1:totalsamples
    z1= [x_free_vec_0(:,i);u(:,i)];
    f1= x_free_vec_1(idx,i);
    for j=1:totalsamples
        z2= [x_free_vec_0(:,j);u(:,j)];    
        f2= x_free_vec_1(idx,j);
        if z1~=z2
            newnorm = norm(f1-f2)./norm(z1-z2);
            newgamma = norm(z1-z2);
            if newnorm > L(1)
                L(1) = newnorm(1);
            end
            if newgamma > gamma(1)
                gamma(1) = newgamma;
            end
           
        end
    end
end
L(2) = [0];
gamma(2) = 0;
idx=2;
for i=1:totalsamples
    z1= [x_free_vec_0(:,i);u(:,i)];
    f1= x_free_vec_1(idx,i);
    for j=1:totalsamples
        z2= [x_free_vec_0(:,j);u(:,j)];    
        f2= x_free_vec_1(idx,j);
        if z1~=z2
            newnorm = norm(f1-f2)./norm(z1-z2);
            newgamma = norm(z1-z2);
            if newnorm > L(2)
                L(2) = newnorm;
            end           
            if newgamma > gamma(2)
                gamma(2) = newgamma;
            end         
        end
    end
end

end

