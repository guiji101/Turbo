clear
clc

allfiles=dir;
allfiles=struct2cell(allfiles);
allfilenames=allfiles(1,:);
cnt=0;
filenames={};
len=length(allfilenames);
expr1='result.*\.txt$';
xlabel='EbN0';
ylabel='BER';
h=[];
for i=1:len
    if(~isempty(regexp(allfilenames{i},expr1)))
        cnt=cnt+1;
        filenames{1,cnt}=allfilenames{i};
    end
end
filenames
Nfiles=length(filenames);
colors='brmcgk';
expr2='EbN0=';
expr3='(\w+).*?=(\S+).*?';
figure
for i=1:Nfiles
    fid=fopen(filenames{1,i});
    foundPattern = 0;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline)   
            break;   
        end
        if (~isempty(regexp(tline,expr2)))
            [tok mat] = regexp(tline, expr3, 'tokens', 'match');
            if ~foundPattern
                foundPattern = 1;
                for j=1:length(tok)
                    eval([mat{1,j} ';'])
                end
            else
                for j=1:length(tok)
                    name=tok{j}{1,1};
                    val=tok{j}{1,2};
                    eval([name '=[' name ',' val '];']);
                end
            end
        end
    end
    fclose(fid);
    eval(['[' xlabel ',' 'ix]=sort(' xlabel ');']);
    eval([ylabel '=' ylabel '(ix);']);
    eval(['htmp=semilogy(' xlabel ',' ylabel ',''-o'',''color'',colors(i),''linewidth'',2),hold on;']);
    h=[h htmp];
    grid on
end
line([-0.4956 -0.4956],[1e-3 1e-5],'color','k','linewidth',2)

legend(h,filenames)