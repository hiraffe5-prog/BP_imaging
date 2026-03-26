function output=taylor_windows(n,sll,N)
% n为泰勒阶数
% sll为第一旁瓣与主瓣的比值分贝数
%  N为输出点数
% 例如 taylor(5,-35,4096)
b=10.^(-sll/20);
a=(log(b+sqrt(b.^2-1)))/pi;
der=sqrt((n.^2)/(a.^2+(n-0.5).^2));
F=zeros(1,n-1);
sum1=1;
for i=1:n-1
    sum1=1;
    sum2=1;
    for j=1:n-1
        sum1=sum1*(1-((i/der).^2)/(a.^2+(j-0.5).^2));
        if j~=i
        sum2=sum2*(1-(i.^2)/(j.^2));
        end
    end
    F(i)=((-1).^(i+1))*sum1/(2*sum2);
end
output=zeros(1,N);
i=0:N-1;
for j=1:n-1
    output=output+F(j)*cos(2*pi*j*(i-N/2+1/2)/N);
end
output=ones(1,N)+2*output;
