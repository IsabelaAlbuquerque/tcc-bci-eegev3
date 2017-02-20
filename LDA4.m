function [W1, W2, W3, acc1, acc2, acc3, acc11, acc22, acc33] = LDA4(dataf1,baseline,dataf2)


% dataf1 = datah671;
% baseline = hbaseline1;
% dataf2 = datah861;
Fs = 500;
n = length(dataf1);
L = 2;
nw = n/(500*L);

%% f1 x Baseline

% Stimulus: freq=6.7Hz
%dataf1=data67;
a=round(6.4*500*L/Fs);  %k=round(f*N/Fs), N=length/(60*ncalib) and f=6.7Hz
b=round(7*500*L/Fs);
c=round(13.1*500*L/Fs);
d=round(13.7*500*L/Fs);

e=round(8.3*500*L/Fs);
f=round(8.9*500*L/Fs);
g=round(16.9*500*L/Fs);
h=round(17.5*500*L/Fs);
% f1
% Channel 1   
featf1_1=zeros(nw,4);
for i=1:nw
    x=dataf1((i-1)*500*L+1:i*(500*L),1);
    X=fft(x);
    absX = abs(X).^2;
    featf1_1(i,1)=mean(absX(a:b));
    featf1_1(i,2)=mean(absX(c:d));
    featf1_1(i,3)=mean(absX(e:f));
    featf1_1(i,4)=mean(absX(g:h));
end;

% % Channel 2
% featf1_2=zeros(nw,4);
% for i=1:nw
%     x=dataf1((i-1)*500*L+1:i*(500*L),2);
%     X=fft(x);
%     absX = abs(X).^2;
%     featf1_2(i,1)=mean(absX(a:b));
%     featf1_2(i,2)=mean(absX(c:d));
%     featf1_2(i,3)=mean(absX(e:f));
%     featf1_2(i,4)=mean(absX(g:h));
% end;

% Channel 3
featf1_3=zeros(nw,4);
for i=1:nw
    x=dataf1((i-1)*500*L+1:i*(500*L),3);
    X=fft(x);
    absX = abs(X).^2;
    featf1_3(i,1)=mean(absX(a:b));
    featf1_3(i,2)=mean(absX(c:d));
    featf1_3(i,3)=mean(absX(e:f));
    featf1_3(i,4)=mean(absX(g:h));
end;

% Stimulus: freq=0Hz

% Channel 1
feat0f1_1=zeros(nw,4);
for i=1:nw
    x=baseline((i-1)*500*L+1:i*(500*L),1);
    X=fft(x);
    absX = abs(X).^2;
    feat0f1_1(i,1)=mean(absX(a:b));
    feat0f1_1(i,2)=mean(absX(c:d));
    feat0f1_1(i,3)=mean(absX(e:f));
    feat0f1_1(i,4)=mean(absX(g:h));
end;

% Channel 2
% feat0f1_2=zeros(nw,4);
% for i=1:nw
%     x=baseline((i-1)*500*L+1:i*(500*L),2);
%     X=fft(x);
%     absX = abs(X).^2;
%     feat0f1_2(i,1)=mean(absX(a:b));
%     feat0f1_2(i,2)=mean(absX(c:d));
%     feat0f1_2(i,3)=mean(absX(e:f));
%     feat0f1_2(i,4)=mean(absX(g:h));
% end;

% Channel 3
feat0f1_3=zeros(nw,4);
for i=1:nw
    x=baseline((i-1)*500*L+1:i*(500*L),3);
    X=fft(x);
    absX = abs(X).^2;
    feat0f1_3(i,1)=mean(absX(a:b));
    feat0f1_3(i,2)=mean(absX(c:d));
    feat0f1_3(i,3)=mean(absX(e:f));
    feat0f1_3(i,4)=mean(absX(g:h));
end;

featf1=[featf1_1 featf1_3];
feat0f1=[feat0f1_1 feat0f1_3];
feat1=vertcat(featf1(1:nw-10,:),feat0f1(1:nw-10,:));
target1=[ones(length(feat1)/2,1);zeros(length(feat1)/2,1)];       % Class 1 = f1. Class 0 = baseline.

W1=LDA(feat1,target1);

L1 = [ones(length(target1),1) feat1] * W1';
P1 = exp(L1)./repmat(sum(exp(L1),2),[1 2]);
[~,d1] = max(P1,[],2);

trainingLabels = target1;
for i=1:length(trainingLabels)
    if trainingLabels(i)==0
        trainingLabels(i)=1;
    else
        trainingLabels(i)=2;
    end;
end;    
TP = sum(d1==trainingLabels);
acc1 = TP/length(d1);
disp(['Training Accuracy Baseline1 vs. f1 : ' num2str(acc1)]);
if acc1 < 0.7
     warndlg(['f1 calibration failed. Accuracy = ' num2str(acc1)],'!! Warning !!')
end

 
%% f2 x Baseline
% Stimulus: freq=8.6Hz


% Channel 1
featf2_1=zeros(nw,4);
for i=1:nw
    x=dataf2((i-1)*500*L+1:i*(500*L),1);
    X=fft(x);
    absX = abs(X).^2;
    featf2_1(i,1)=mean(absX(a:b));
    featf2_1(i,2)=mean(absX(c:d));
    featf2_1(i,3)=mean(absX(e:f));
    featf2_1(i,4)=mean(absX(g:h));
end;

% Channel 2
% featf2_2=zeros(nw,4);
% for i=1:nw
%     x=dataf2((i-1)*500*L+1:i*(500*L),2);
%     X=fft(x);
%     absX = abs(X).^2;
%     featf2_2(i,1)=mean(absX(a:b));
%     featf2_2(i,2)=mean(absX(c:d));
%     featf2_2(i,3)=mean(absX(e:f));
%     featf2_2(i,4)=mean(absX(g:h));
% end;

% Channel 3   
featf2_3=zeros(nw,4);
for i=1:nw
    x=dataf2((i-1)*500*L+1:i*(500*L),3);
    X=fft(x);
    absX = abs(X).^2;
    featf2_3(i,1)=mean(absX(a:b));
    featf2_3(i,2)=mean(absX(c:d));
    featf2_3(i,3)=mean(absX(e:f));
    featf2_3(i,4)=mean(absX(g:h));
end;

% Stimulus: freq=0Hz
n=length(baseline);
% Channel 1
feat0f2_1=zeros(nw,4);
for i=1:nw
    x=baseline((i-1)*500*L+1:i*(500*L),1);
    X=fft(x);
    absX = abs(X).^2;
    feat0f2_1(i,1)=mean(absX(a:b));
    feat0f2_1(i,2)=mean(absX(c:d));
    feat0f2_1(i,3)=mean(absX(e:f));
    feat0f2_1(i,4)=mean(absX(g:h));
end;

% Channel 2
% feat0f2_2=zeros(nw,4);
% for i=1:nw
%     x=baseline((i-1)*500*L+1:i*(500*L),2);
%     X=fft(x);
%     absX = abs(X).^2;
%     feat0f2_2(i,1)=mean(absX(a:b));
%     feat0f2_2(i,2)=mean(absX(c:d));
%     feat0f2_2(i,3)=mean(absX(e:f));
%     feat0f2_2(i,4)=mean(absX(g:h));
% end;

% Channel 3
feat0f2_3=zeros(nw,4);
for i=1:nw
    x=baseline((i-1)*500*L+1:i*(500*L),3);
    X=fft(x);
    absX = abs(X).^2;
    feat0f2_3(i,1)=mean(absX(a:b));
    feat0f2_3(i,2)=mean(absX(c:d));
    feat0f2_3(i,3)=mean(absX(e:f));
    feat0f2_3(i,4)=mean(absX(g:h));
end;

featf2=[featf2_1 featf2_3];
feat0f2=[feat0f2_1 feat0f2_3];
feat2=vertcat(featf2(1:nw-10,:),feat0f2(1:nw-10,:)); 
target2=[ones(length(feat2)/2,1);zeros(length(feat2)/2,1)];  % Class 1 = f2. Class 0 = baseline.

W2=LDA(feat2,target2);

L2 = [ones(length(target2),1) feat2] * W2';
P2 = exp(L2)./repmat(sum(exp(L2),2),[1 2]);
[~,d1] = max(P2,[],2);

trainingLabels = target2;
for i=1:length(trainingLabels)
    if trainingLabels(i)==0
        trainingLabels(i)=1;
    else
        trainingLabels(i)=2;
    end;
end;    
TP = sum(d1==trainingLabels);
acc2 = TP/length(d1);
disp(['Training Accuracy Baseline2 vs. f2 : ' num2str(acc2)]);
if acc2 < 0.7
     warndlg(['f2 calibration failed. Accuracy = ' num2str(acc2)],'!! Warning !!')
end

%% f1 x f2


% Stimulus: freq=6.7Hz
%dataf1=data67;

% Channel 1   
featf12_1=zeros(nw,4);
for i=1:nw
    x=dataf1((i-1)*500*L+1:i*(500*L),1);
    X=fft(x);
    absX = abs(X).^2;
    featf12_1(i,1)=mean(absX(a:b));
    featf12_1(i,2)=mean(absX(c:d));
    featf12_1(i,3)=mean(absX(e:f));
    featf12_1(i,4)=mean(absX(g:h));
end;

% Channel 2
% featf12_2=zeros(nw,4);
% for i=1:nw
%     x=dataf1((i-1)*500*L+1:i*(500*L),2);
%     X=fft(x);
%     absX = abs(X).^2;
%     featf12_2(i,1)=mean(absX(a:b));
%     featf12_2(i,2)=mean(absX(c:d));
%     featf12_2(i,3)=mean(absX(e:f));
%     featf12_2(i,4)=mean(absX(g:h));
% end;

% Channel 3
featf12_3=zeros(nw,4);
for i=1:nw
    x=dataf1((i-1)*500*L+1:i*(500*L),3);
    X=fft(x);
    absX = abs(X).^2;
    featf12_3(i,1)=mean(absX(a:b));
    featf12_3(i,2)=mean(absX(c:d));
    featf12_3(i,3)=mean(absX(e:f));
    featf12_3(i,4)=mean(absX(g:h));
end;


% Stimulus: freq=8.6Hz

% Channel 1
featf21_1=zeros(nw,4);
for i=1:nw
    x=dataf2((i-1)*500*L+1:i*(500*L),1);
    X=fft(x);
    absX = abs(X).^2;
    featf21_1(i,1)=mean(absX(a:b));
    featf21_1(i,2)=mean(absX(c:d));
    featf21_1(i,3)=mean(absX(e:f));
    featf21_1(i,4)=mean(absX(g:h));
end;

% Channel 2
% featf21_2=zeros(nw,4);
% for i=1:nw
%     x=dataf2((i-1)*500*L+1:i*(500*L),2);
%     X=fft(x);
%     absX = abs(X).^2;
%     featf21_2(i,1)=mean(absX(a:b));
%     featf21_2(i,2)=mean(absX(c:d));
%     featf21_2(i,3)=mean(absX(e:f));
%     featf21_2(i,4)=mean(absX(g:h));
% end;

% Channel 3   
featf21_3=zeros(nw,4);
for i=1:nw
    x=dataf2((i-1)*500*L+1:i*(500*L),3);
    X=fft(x);
    absX = abs(X).^2;
    featf21_3(i,1)=mean(absX(a:b));
    featf21_3(i,2)=mean(absX(c:d));
    featf21_3(i,3)=mean(absX(e:f));
    featf21_3(i,4)=mean(absX(g:h));
    
end;

featf12=[featf12_1 featf12_3];
featf21=[featf21_1 featf21_3];
feat3=vertcat(featf21(1:nw-10,:), featf12(1:nw-10,:));
target3=[ones(length(feat3)/2,1);zeros(length(feat3)/2,1)];  % Class 1 = f2. Class 0 = baseline.

W3=LDA(feat3,target3);

L3 = [ones(length(target3),1) feat3] * W3';
P3 = exp(L3)./repmat(sum(exp(L3),2),[1 2]);
[~,d1] = max(P3,[],2);

trainingLabels = target3;
for i=1:length(trainingLabels)
    if trainingLabels(i)==0
        trainingLabels(i)=1;
    else
        trainingLabels(i)=2;
    end;
end;    
TP = sum(d1==trainingLabels);
acc3 = TP/length(d1);
disp(['Training Accuracy f1 vs. f2 : ' num2str(acc3)]);
if acc3 < 0.7
     warndlg(['f1xf2 calibration failed. Accuracy = ' num2str(acc3)],'!! Warning !!')
end

%% Testing

dec1 = vertcat(featf1((nw-9):nw,:),feat0f1((nw-9):nw,:));
target11 = [ones(length(dec1)/2,1);zeros(length(dec1)/2,1)];

dec2 = vertcat(featf2((nw-9):nw,:),feat0f2((nw-9):nw,:));
target22 = [ones(length(dec2)/2,1);zeros(length(dec2)/2,1)];

dec3 = vertcat(featf21((nw-9):nw,:),featf12((nw-9):nw,:));
target33 = [ones(length(dec3)/2,1);zeros(length(dec3)/2,1)];


L11 = [zeros(length(dec1),1) dec1] *W1';
P11= exp(L11)./repmat(sum(exp(L11),2),[1 2]);
[~,d1] = max(P11,[],2);

trainingLabels = target11;
for i=1:length(trainingLabels)
    if trainingLabels(i)==0
        trainingLabels(i)=1;
    else
        trainingLabels(i)=2;
    end;
end;    
TP = sum(d1==trainingLabels);
acc11 = TP/length(d1);
disp(['Testing Accuracy Baseline1 vs. f1 : ' num2str(acc11)]);
if acc11 < 0.7
     warndlg(['f1 training failed. Accuracy = ' num2str(acc11)],'!! Warning !!')
end


L22 = [zeros(length(dec2),1) dec2] *W2';
P22 = exp(L22)./repmat(sum(exp(L22),2),[1 2]);
[~,d1] = max(P22,[],2);

trainingLabels = target22;
for i=1:length(trainingLabels)
    if trainingLabels(i)==0
        trainingLabels(i)=1;
    else
        trainingLabels(i)=2;
    end;
end;    
TP = sum(d1==trainingLabels);
acc22 = TP/length(d1);
disp(['Testing Accuracy Baseline1 vs. f2 : ' num2str(acc22)]);
if acc22 < 0.7
     warndlg(['f2 training failed. Accuracy = ' num2str(acc22)],'!! Warning !!')
end


L33 = [zeros(length(dec3),1) dec3] *W3';
P33 = exp(L33)./repmat(sum(exp(L33),2),[1 2]);
[~,d1] = max(P33,[],2);

trainingLabels = target33;
for i=1:length(trainingLabels)
    if trainingLabels(i)==0
        trainingLabels(i)=1;
    else
        trainingLabels(i)=2;
    end;
end;    
TP = sum(d1==trainingLabels);
acc33 = TP/length(d1);
disp(['Testing Accuracy f2 vs. f1 : ' num2str(acc33)]);
if acc33 < 0.7
     warndlg(['f2 x f1 training failed. Accuracy = ' num2str(acc33)],'!! Warning !!')
end

end
