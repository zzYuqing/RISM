function C= RISM(data,K)

    data=normalization(data);
    N=size(data,1);
    % Euclidean distance
    dist=pdist2(data,data);
    % Initial cluster number
    CP_num=ceil(sqrt(N));
    
    %% kNN
    k=10;%default nearest neighbor value
    % Distance sort
    [sort_dist,p]=sort(dist,2);
    % get kNN
    NK=p(:,2:k+1);
    % NK_V(:,:,i)Represents the kNN-Graph when k==i
    NK_V=zeros(N,N,k);
    for i=1:N
        for j=1:k
            NK_V(i,NK(i,1:j),j)=1;
        end
    end
    % kNN-Graph
    NK_VK=NK_V(:,:,k);
    % Shared kNN-Graph
    NK_VK_shared=NK_VK&NK_VK';
    
    %% Relative density and Relative distance

    % Absolute density: rho_abs
    rho_abs=1./mean(sort_dist(:,2:k+1),2);    
    % Relative density: rho_rela
    rho_abs_NK=rho_abs(NK);
    rho_rela=rho_abs./sum(rho_abs_NK,2)*k;
    rho_rela=normalization(rho_rela);
    % Relative distance: delta
    [~,ordrho]=sort(rho_rela,'descend');%Rho_rela sort
    delta=zeros(size(data,1),1);
    % The closest point of higher density: nec
    nec=zeros(size(data,1),1);
    % process the highest relative density point
	delta(ordrho(1))=max(dist(ordrho(1),:));
    nec(ordrho(1))=ordrho(1);
    for i=2:N
        b=ordrho(1:i-1);
        [v,p_rho]=min(dist(ordrho(i),b));
        delta(ordrho(i))=v;
        nec(ordrho(i))=b(p_rho);
    end
    % delta normalization
    delta=(delta-min(delta))./(max(delta)-min(delta));
    %% Split
    gamma=delta.*rho_rela;
    [~,ordgamma]=sort(gamma,'descend');
    % Initial clustering center points: centers
    centers=ordgamma(1:CP_num);
    centers=sort(centers);
    C=zeros(N,1);% Initialize clustering results
    [~,ordrho]=sort(rho_rela,'descend');
    C(centers)=1:size(centers,1);% Label the centers
    for i=2:N
        if(ismember(ordrho(i), centers))
            continue;
        end
        C(ordrho(i))=C(nec(ordrho(i)));%Label remaining points
    end
    C=clear_C(C);
    CP_num=length(unique(C));
    % merge single point clusters
    for i=1:CP_num
        if(CP_num<=K)
            C=clear_C(C);
            return
        end
        points=find(C==i);
        if(length(points)~=1)
            continue;
        end
        C(points)=C(p(points,2));
        CP_num=CP_num-1;
    end
    C=clear_C(C);
    CP_num=length(unique(C));
    
    %% Visualize small clusters
%     figure;
%     gscatter(data(:,1),data(:,2),C);
%     hold on;
%     gplot(NK_VK_shared,data);
    %% Distance factor
    P=zeros(CP_num,CP_num);
    for i=1:CP_num-1
        points_i=find(C==i);
        for j=i+1:CP_num
            points_j=find(C==j);
            V_ij=NK_VK_shared(points_i,points_j);
            % M(Ci,Cj) is empty
            if(length(find(V_ij==1))<1)
                P(i,j)=nan;
            else
                dist_ij=dist(points_i,points_j);
                dist_V_ij=dist_ij.*V_ij;
                dist_list=dist_V_ij(find(dist_V_ij~=0));
                P(i,j)=mean(dist_list);
            end
        end
    end
    P(find(P==0))=nan;
    %% Connectivity degree
    U=zeros(CP_num,CP_num);
    S=ones(CP_num,CP_num)*(-1);
    u_im=zeros(N,1);
    asP=zeros(N,1);
    MD=zeros(N,1);
    for i=1:CP_num
        % Ci to Cj
        points_i=find(C==i);
        Ni=length(points_i);
        if(Ni<k)% adapted k:Ki
            ki=Ni;
            NK_Vi=NK_V(:,:,ki);
        else
            ki=k;
            NK_Vi=NK_VK;
        end
        near_i=NK_Vi(points_i,:);
        MD_i=sum(near_i(:,setdiff(1:N,points_i)),2);
        % Ci's border points
        asP_i=find(MD_i~=0);
        asP(points_i(asP_i))=1;
        u_im(points_i)=(ki-MD_i);
        if(isempty(asP_i))
            U(i,:)=0;
        end
        MD_i=MD_i(find(MD_i~=0));
        MD(points_i(asP_i))=MD_i;
        for j=1:CP_num
            if(j==i)
                continue;
            end
            points_j=find(C==j);
            Nj=length(points_j);
            if(Nj<Ni)
                continue;
            end
            if(~isempty(asP_i))
                near_ij=NK_Vi(points_i(asP_i),points_j);%Nearest neighbor of Ci to Cj
                sum_near_ij=sum(near_ij,2);
                u_om=zeros(length(points_i),1);
                u_om(asP_i)=sum_near_ij./ki;% Unified connectivity degree of point
                C2C=u_om.*u_im(points_i);
                U(i,j)=mean(C2C(find(~isnan(C2C))));
            end
            if(i<j)
                S(i,j)=1-U(i,j);
                if(or(S(i,j)==1,isnan(S(i,j))))
                    S(i,j)=-1;
                end
            else
                S(j,i)=1-U(i,j);
                if(or(S(j,i)==1,isnan(S(j,i))))
                    S(j,i)=-1;
                end
            end

        end
    end
    S(find(S==-1))=nan;
    %% Inter-cluster distance:D
    D=P.*S;
    %% Merge stage 1
    min_disD_list=[];% save inter-cluster distance
    C_history=[];% save C
    if K<=1% Trigger the calculation of the number of clusters
        K=1;
    end
    while(CP_num>K)
        min_dis=min(D(:));
        min_disD_list=[min_disD_list,min_dis];% push inter-cluster distance
        [i,j]=find(D==min_dis(1));
        if(isempty(i))
            break;
        end
        if(length(i)>1)%void duplicate value
            i=i(1);
            j=j(1);
        end
        points=find(C==j);
        C(points)=i;
        C=clear_C(C);
        C_history=[C_history,C];% push C
        CP_num=length(unique(C));
        new_label=C(points(1));
        P(j,:)=[];
        P(:,j)=[];
        U(j,:)=[];
        U(:,j)=[];
        U(new_label,:)=0;
        U(:,new_label)=0;
        S(j,:)=[];
        S(:,j)=[];
        % update inter-cluster distance
        points_i=find(C==new_label);
        asP(points_i)=0;
        Ni=length(points_i);
        if(Ni<k)
            ki=Ni;
            NK_Vi=NK_V(:,:,ki);
        else
            ki=k;
            NK_Vi=NK_VK;
        end
        near_i=NK_Vi(points_i,:);
        MD_i=sum(near_i(:,setdiff(1:N,points_i)),2);
        inP=points_i(find(MD_i==0));
        asP_i=find(MD_i~=0);
        asP(points_i(asP_i))=1;
        u_im(points_i)=(ki-MD_i);
        if(isempty(asP_i))
            U(i,:)=0;
        end
        MD_i=MD_i(find(MD_i~=0));
        MD(points_i(asP_i))=MD_i;
        for j=setdiff(1:CP_num,new_label)
            points_j=find(C==j);
            Nj=length(points_j);
            V_ij=NK_VK_shared(points_i,points_j);
            % update distance factor
            if(length(find(V_ij==1))<1)
                P(i,j)=nan;
            else
                dist_ij=dist(points_i,points_j);
                dist_V_ij=dist_ij.*V_ij;
                dist_list=dist_V_ij(find(dist_V_ij~=0));
                if(new_label<j)
                    P(new_label,j)=mean(dist_list);
                else
                    P(j,new_label)=mean(dist_list);
                end
                
                % update connectivity degree
                if(Ni<Nj)
                    % Ci to Cj
                    if(~isempty(asP_i))
                        near_ij=NK_Vi(points_i(asP_i),points_j);
                        sum_near_ij=sum(near_ij,2);
                        u_om=zeros(length(points_i),1);
                        u_om(asP_i)=sum_near_ij./ki;
                        C2C=u_om.*u_im(points_i);
                        U(new_label,j)=mean(C2C(find(~isnan(C2C))));
                        U(j,new_label)=U(new_label,j);
                    end
                    if(new_label<j)
                        S(new_label,j)=1-U(new_label,j);
                        if(or(S(new_label,j)==1,isnan(S(new_label,j))))
                            S(new_label,j)=nan;
                        end
                    else
                        S(j,new_label)=1-U(j,new_label);
                        if(or(S(j,new_label)==1,isnan(S(j,new_label))))
                            S(j,new_label)=nan;
                        end
                    end
                else
                    % Cj to Ci
                    if(Nj<k)
                        kj=Nj;
                        NK_Vj=NK_V(:,:,kj);
                    else
                        kj=k;
                        NK_Vj=NK_VK;
                    end
                    asP_j=intersect(find(asP==1),points_j);
                    MD_j=MD(asP_j);
                    if(~isempty(asP_j))
                        near_ji=NK_Vj(asP_j,points_i);
                        sum_near_ji=sum(near_ji,2);
                        u_om=zeros(N,1);
                        u_om(asP_j)=sum_near_ji./kj;
                        C2C=u_om(points_j).*u_im(points_j);
                        U(new_label,j)=mean(C2C(find(~isnan(C2C))));
                        U(j,new_label)=U(new_label,j);
                    end
                    if(new_label<j)
                        S(new_label,j)=1-U(new_label,j);
                        if(or(S(new_label,j)==1,isnan(S(new_label,j))))
                            S(new_label,j)=nan;
                        end
                    else
                        S(j,new_label)=1-U(j,new_label);
                        if(or(S(j,new_label)==1,isnan(S(j,new_label))))
                            S(j,new_label)=nan;
                        end
                    end
                end
            end
        end
        % update inter-cluster distance
        D=P.*S;
    end
    C=clear_C(C);
    % update cluster number
    CP_num=length(unique(C));
    C_history=[C_history,C];% push C
    if(CP_num==K && K==1) % Termination state
        min_dis_diff=diff(min_disD_list);
        min_dis_diff(1:2)=0;
        [~,loc]=max(min_dis_diff);
        C=C_history(:,loc);% clustring result
        return;
    end
    %% Merge stage 2
    % if the merge continues
    if isnan(min_dis)
        % delete wrong results
        min_disD_list(end)=[];
        C_history(:,end)=[];
    end
    % calculate distance factor
    D=zeros(CP_num);
    for i=1:CP_num
        points_i=find(C==i);
        for j=i+1:CP_num
            points_j=find(C==j);
            dist_ij=dist(points_i,points_j);
            dist_list=dist_ij(find(dist_ij~=0));
            D(i,j)=mean(dist_list(~isnan(dist_list)));
        end
    end
    D(find(D==0))=nan;
    while(CP_num>K)
        min_dis=min(D(:));
        min_disD_list=[min_disD_list,min_dis];% push inter-cluster distance
        [i,j]=find(D==min_dis(1));
        if(length(i)>1)
            i=i(1);
            j=j(1);
        end
        points=find(C==j);
        C(points)=i;
        C=clear_C(C);
        C_history=[C_history,C];% push C
        CP_num=length(unique(C));
        new_label=C(points(1));
        points_i=find(C==new_label);       
        D(j,:)=[];
        D(:,j)=[];
        for j=setdiff(1:CP_num,new_label)
            % update distance factor
            points_j=find(C==j);
            dist_ij=dist(points_i,points_j);
            dist_list=dist_ij(find(dist_ij~=0));
            if(new_label<j)
                D(new_label,j)=mean(dist_list(~isnan(dist_list)));
            else
                D(j,new_label)=mean(dist_list(~isnan(dist_list)));
            end
        end
        D(find(D==0))=nan;
    end

    if(K==1) % Termination state
        min_dis_diff=diff(min_disD_list);
        min_dis_diff(1)=0;
        [~,loc]=max(min_dis_diff);
        C=C_history(:,loc);
        return;
    end
end

