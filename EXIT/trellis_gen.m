clear
clc
f=13;
g=15;

str=num2str(f);
len=length(str);
backward=[];
for i=1:len
    backward=[backward fliplr(dec2binvec(str2num(str(i)),3))];
end

str=num2str(g);
len=length(str);
forward=[];
for i=1:len
    forward=[forward fliplr(dec2binvec(str2num(str(i)),3))];
end

t=backward+forward;
ind=find(t>0);
if ind(1)>1
    backward(1:ind(1)-1)=[];
    forward(1:ind(1)-1)=[];
end
assert( length(backward) == length(forward) );

mu=length(backward)-1;
nstate=2^mu;
allstate=zeros(nstate,mu);
for i=1:nstate
    allstate(i,:)=dec2binvec(i-1,mu);
end
allstate=fliplr(allstate);

nxstate=zeros(2,nstate);
output=zeros(2,nstate);

for x=0:1
    temp=zeros(size(allstate));
    for i=1:nstate
        s=allstate(i,:);
        s1new=rem(dot(backward(2:end),s)+x,2);
        output(x+1,i)=rem( dot(forward,[s1new,s]) ,2);
        temp(i,1)=s1new;
        temp(i,2:end)=s(1:end-1);
    end
    
    col=zeros(nstate,1);
    for i=1:mu
        col=col+2^(mu-i)*temp(:,i);
    end
    nxstate(x+1,:)=col';
end

tailbits=zeros(nstate,mu);
for i=1:nstate
    s=allstate(i,:);
    for j=1:mu
        tailbits(i,j)=rem(dot(backward(2:end),s),2);
        s=[0 s(1:end-1)];
    end
end
col=zeros(nstate,1);
for i=1:mu
    col=col+2^(i-1)*tailbits(:,i);
end
tail=col';


fid=fopen(['(',num2str(f),',',num2str(g),')trellis_data.txt'],'wt');
fprintf(fid,'%d\n',mu);
for i=1:2
    for j=1:nstate
        fprintf(fid,'%d ',nxstate(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:2
    for j=1:nstate
        fprintf(fid,'%d ',output(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:length(tail)
    fprintf(fid,'%d ',tail(i));
end
fprintf(fid,'\nformat:\nmu\nnext state \noutput \ntail bits');
fclose(fid)

