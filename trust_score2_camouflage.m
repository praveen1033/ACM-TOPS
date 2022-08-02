houseData = csvread('dataport-export(2).csv', 1,1);
C=unique(houseData(:,1));
temp=[];
for i=1:length(C)
    for j=1:length(houseData)
        if houseData(j,1) == C(i)
            temp=[temp;houseData(j,:)];
        end
    end
    f=num2str(i);
    data4{i}=temp;
    temp=[];
end




allPt1=[];
count=1;
for p=1:length(data4)
   tempo=data4{p};
   if length(tempo)>8640
      allPt1(:,count)=tempo(1:8640,2);
   end
   count=count+1;
end
allPt1=allPt1*1000;

allPt1(allPt1>6000)=1;
allPt1(allPt1<50)=1;

allPt2=reshape(allPt1',1,numel(allPt1));
allPt3 = boxcox(0.25,allPt2');
for i=1:length(allPt3)
allPt(ceil(i/175.00001),ceil(mod(i,175.00001))) = allPt3(i);
end

F=30;




for i=1:8640
frame=ceil(i/(F*24));
am1(frame,mod(i,F*24)+1)=mean(allPt(i,:));
tempi=allPt(i,:);
tempi(tempi<=0)=[];
hm1(frame,mod(i,F*24)+1)=hmn(tempi);
% gm1(frame,mod(i,F*24)+1)=geomean(tempi);
std1(frame,mod(i,F*24)+1)=std(tempi);
mad1(frame,mod(i,F*24)+1)=mad(tempi);
end


% for j=1:175
% for i=1:12
% k1(j,i)=kurtosis(mean(allPt((i-1)*720+1:720*i,j),2));
% sk1(j,i)=skewness(mean(allPt((i-1)*720+1:720*i,j),2));
% std1(j,i)=std(mean(allPt((i-1)*720+1:720*i,j),2));
% um1(j,i)=mean(mean(allPt((i-1)*720+1:720*i,j),2));
% mad1(j,i)=mad(mean(allPt((i-1)*720+1:720*i,j),2),0);
% end
% end
% 
% for j=1:175
% for i=1:12
% std1(j,i)=std(allPt((i-1)*720+1:720*i,j));
% m1(j,i)=mean(allPt((i-1)*720+1:720*i,j));
% mad1(j,i)=mad(allPt((i-1)*720+1:720*i,j));
% end
% end



for i=1:8640
for j=1:175
if allPt(i,j) <= (mad1(ceil(i/(F*24)),mod(i,F*24)+1) + am1(ceil(i/(F*24)),mod(i,F*24)+1)) && allPt(i,j) >= (-mad1(ceil(i/(F*24)),mod(i,F*24)+1) + am1(ceil(i/(F*24)),mod(i,F*24)+1))
X(j,ceil(i/(F*24)),mod(i,F*24)+1) = 1;
else
X(j,ceil(i/(F*24)),mod(i,F*24)+1) = 0;
end
end
end




M=16;                                                           %Number of attacked meters is assigned to M
N=175;


dmin=750;                                                        %Loading Delta min
dmax=850;                                                       %Loading Delta max



for i=1:M
            d=dmin+(dmax-dmin).*rand(8640,1);          %Creating an array of numbers b/w dmin and dmax to add them to affected meters
            if mean(allPt1(:,1)) < 200
                for j=1:8640
                    allPtadd1(j,i)=allPt1(j,i)+d(j,1);
                end
            elseif mean(allPt1(:,1)) > 800
                for j=1:8640
                    allPtadd1(j,i)=allPt1(j,i)+d(j,1);
                end
            else
                for j=1:8640
                    if mod(j,3)==1
                        allPtadd1(j,i)=allPt1(j,i)+d(j,1);
                    else
                        allPtadd1(j,i)=allPt1(j,i)-d(j,1);
                    end
                end
            end
end


for i=M+1:175
            for j=1:8640
                allPtadd1(j,i)=allPt1(j,i);
            end
end

allPtadd2=reshape(allPtadd1',1,numel(allPtadd1));
%[allPtadd3,l1] = boxcox(allPtadd2');
allPtadd3 = boxcox(0.25,allPtadd2');
for i=1:length(allPtadd3)
allPtadd(ceil(i/175.00001),ceil(mod(i,175.00001))) = allPtadd3(i);
end

% for j=1:175
% for i=1:12
% k2(j,i)=kurtosis(mean(allPtadd((i-1)*720+1:720*i,j),2));
% sk2(j,i)=skewness(mean(allPtadd((i-1)*720+1:720*i,j),2));
% std2(j,i)=std(mean(allPtadd((i-1)*720+1:720*i,j),2));
% um2(j,i)=mean(mean(allPtadd((i-1)*720+1:720*i,j),2));
% mad2(j,i)=mad(mean(allPtadd((i-1)*720+1:720*i,j),2),0);
% end
% end
% for j=1:175
% for i=1:12
% std2(j,i)=std(allPtadd((i-1)*720+1:720*i,j));
% um2(j,i)=mean(allPtadd((i-1)*720+1:720*i,j));
% mad2(j,i)=mad(allPtadd((i-1)*720+1:720*i,j));
% end
% end


for i=1:8640
frame=ceil(i/(F*24));
am2(frame,mod(i,F*24)+1)=mean(allPtadd(i,:));
tempi=allPtadd(i,:);
tempi(tempi<=0)=[];
hm2(frame,mod(i,F*24)+1)=hmn(tempi);
% gm2(frame,mod(i,F*24)+1)=geomean(tempi);
std2(frame,mod(i,F*24)+1)=std(tempi);
mad2(frame,mod(i,F*24)+1)=mad(tempi);
end

for i=1:12
for j=1:720
rm(i,j) = hm2(i,j) - (am2(i,j) - hm2(i,j));
end
end

for j=1:175
for i=1:8640
if allPtadd(i,j) <=  (mad1(ceil(i/(F*24)),mod(i,F*24)+1) + am1(ceil(i/(F*24)),mod(i,F*24)+1)) && allPtadd(i,j) >= (-mad1(ceil(i/(F*24)),mod(i,F*24)+1) + am1(ceil(i/(F*24)),mod(i,F*24)+1))
Y(j,ceil(i/(F*24)),mod(i,F*24)+1) = 1;
else
Y(j,ceil(i/(F*24)),mod(i,F*24)+1) = 0;
end
end
end





 for i=1:175
for j=1:12
r(i,j)=(sum(X(i,j,:))+1)/(720+2);
q(i,j)=(sum(Y(i,j,:))+1)/(720+2);
end
 end
 
 
    for i=1:175
    D(i) = (1-r(i,1)) * log((1-r(i,1))/(1-q(i,2))) + r(i,1) *  log((r(i,1))/(q(i,2)));
    end


    for i=1:175
    Q(i) = 1/(1+D(i)^0.5);
    end
    
    
idx=kmeans(Q',2);

count1=0;
count2=0;
for i=1:M
if idx(i)==1
count1 = count1+1;
end
end
for i=M+1:N
if idx(i)==2
count2 = count2+1;
end
end
disp(count1);
disp(count2);

hold on;
for i=1:N
if idx(i)==1
LH(1)=plot(i,Q(i),'r*');
else
LH(2)=plot(i,Q(i),'bo');
end
end

% [label,score] = predict(SVM,Q');
% hold on;
% for i=1:N
% if label(i)==1
% LH(1)=plot(i,Q(i),'r*');
% else
% LH(2)=plot(i,Q(i),'bo');
% end
% end

% for i=1:40
% if i<=16
% idx(i)=1;
% else
% idx(i)=2;
% end
% end
% SVM = fitcsvm(res',idx');