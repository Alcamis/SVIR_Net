function [laplacian_L,K]=Laplacian_matrix(L)

    size_L=size(L,1);
    sum_L=sum(L);
    K=zeros(size_L,size_L);
    for i=1:size_L
        K(i,i)=sum_L(i);
    end
    
    laplacian_L=K-L;
end