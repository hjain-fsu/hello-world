
function N=control(x,p) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control growth: Exponential
%
% p: growth rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = p*x(1);
end

function F1 = control_soln(pars,time,x0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control growth: Solve
%
% pars: growth rate
% time: time points
% x0: initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t,F1]=ode45(@control,time,x0,[],pars); 

end

function value = cutoff(months,x0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the size an untreated
% tumor reaches after a certain number of
% months
%
% months: time tumor grows
% x0: initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Growth rates considered
alphas = 0:0.001:1.55;

%%% Transforms months to days
time = 0:0.1:months*30;

%%% For each growth rate find final tumor size
for i = 1:length(alphas)
    y1 = control_soln(alphas(i)/(70),time,x0);
    F1(i) = y1(end);
end

%%% Calculate average value
value = sum(F1(:))/length(alphas);

end



function [Time_death,censorFlag] = transform_to_Kaplan(SOL,t0,Npatients,cutoff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reshapes data to perform
% survival analysis
%
% SOL: structure with two components
%      - SOL.x: time vector
%      - SOL.y: tumor growth vector
% Npatients: number of patients
% cutoff: tumor size considered as deadly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Npatients

    time = SOL(i).x;
    thisPatientTumor = SOL(i).y;
 
    % Find when patients dies
    indx_death = find(thisPatientTumor >= cutoff);
    
    if isempty(indx_death)  % patient never reached cutoff
        Time_death(i) = max(time);
        censorFlag(i) = 1;
    else
        indx_death = indx_death(1); % take first occurrence
        Time_death(i) = time(indx_death) - t0;
        censorFlag(i) = 0;
    end
end
end



function varargout=kmplot6(x1,x2,x3,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification of  kmplot.m available at
% MathWorks
%
% x1: patients under conventional treatment
% x2: patients under TMZ monotherapy optimized
% x3: patients under best combined therapy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Confidence Intervals
P_alpha = 0.05;
P_color = parula(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       CONVENTIONAL TREATMENT        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x1;

%Input Error handling
p=inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'alpha',P_alpha, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'cflag',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'flag',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x,varargin{:});
alpha=p.Results.alpha; cflag=p.Results.cflag; flag=p.Results.flag;
clear p
assert(all(x(:,2)==0 | x(:,2)==1),'Warning: all x(:,2) values must be 0 or 1')


%sort data by survival time
x=sortrows(x,1);

%table of patients observed for each survival time
%the TABULATE function sets up this matrix:
%table1=[time count percent(on total)]

table1=[0 size(x,1) 1; tabulate(x(:,1))];

%if all observed time are integers remove not observed time added by
%TABULATE function
table1(table1(:,3)==0,:)=[];

%Table of censored data
table12=tabulate(x(x(:,2)==1));
if ~isempty(table12)
    % remove not observed time added by TABULATE function
    table12(table12(:,3)==0,:)=[];
    % setup the vector of the censored data
    [cens,loc]=ismember(table1(:,1),table12(:,1)); %find censored data
end    

%the percents stored in the the third column are unuseful;
%so, place in the third column how many subjects are still alive at the
%beginning of the i-th interval.
a1=[table1(1,2); -1.*table1(2:end,2)];
table1(:,3)=cumsum(a1); table1(2:end,3)=table1(1:end-1,3);

%number of deaths in the intervals (don't take in account the censored
%data)
if ~isempty(table12)
    table1(cens,2)=table1(cens,2)-table12(loc(cens),2);
end
%finally, delete the first row that is now useless
table1(1,:)=[];


t1=[0;table1(:,1)]; %this is the x variable (time);
T1=[1;cumprod(1-(table1(:,2)./table1(:,3)))]; %this is the y variable (survival function)


if flag %if this function was not called by LOGRANK function
    
    %compute the standard error of the survival function
    SE=[0;T1(2:end).*sqrt(cumsum(table1(:,2)./(table1(:,3).* ...
        (table1(:,3)-table1(:,2)))))];
    
end

%censored data plotting
if ~isempty(table12) 
    %if there are censored data after max(t1), add a new cell into the t1,
    %T1 and SE arrays
    if table12(end,1)>=t1(end,1)
        t1(end+1,1)=table12(end,1)+1;
        T1(end+1,1)=T1(end,1);
        if flag %if this function was not called by LOGRANK function
            SE(end+1,1)=SE(end,1);
        end
    end
    if ~cflag
        %vectors preallocation
        xcg=zeros(1,sum(table12(:,2))); ycg=xcg; J=1;
        %for each censored data into the i-th time interval...
        for I=1:size(table12,1)
            %compute how many position into the array they must occupy
            JJ=J+table12(I,2)-1;
            %find the correct time interval in which censored data must be
            %placed
            A=find(t1<=table12(I,1),1,'last');
            B=find(t1>table12(I,1),1,'first');
            %equally divide this interval
            int=linspace(table12(I,1),t1(B,1),table12(I,2)+2);
            %put all in the vectors of the plotting variables
            xcg(J:JJ)=int(2:end-1);
            ycg(J:JJ)=T1(A);
            %update the counter
            J=JJ+1;
        end
    else
        xcg=table1(table1(:,2)==0,1);
        ycg=T1(table1(:,2)==0);
    end
else
    if ~flag %if this function was called by LOGRANK function
        xcg=[]; ycg=[];
    end
end
%compute the hazard rate
c1=T1.*numel(x);
c2=-(diff(log(c1(1:end-1)))./diff(t1(1:end-1)));
lambda1=mean(c2(c2~=0));

if flag %if this function was not called by LOGRANK function
    %compute the (1-alpha)*100% confidence interval curves
    cv=realsqrt(2)*erfcinv(alpha); %critical value
    %lower curve (remember that: the lower curve values can't be negative)
    lowc=max(0,T1-SE.*cv);
    %if the lower curve reaches the 0 earlier than survival function, trim the
    %data.
    if isequal(lowc(end-1:end),[0; 0])
        lowcend=find(lowc==0,1,'first');
    else
        lowcend=length(lowc);
    end
    %upper curve (remember that the upper curve values can't be >1)
    upc=min(1,T1+SE.*cv);
    %eventually, correct the data.
    if isequal(upc(end),1) 
        cupend=find(upc<1,1,'last');
        upc(cupend:end)=upc(cupend);
    end

    %compute the median survival time (if exist...)
    if isempty(T1(T1==0.5)) %if there is not a point where T=0.5...
        I=find(T1>0.5,1,'last'); %find the first point where T>0.5
        J=find(T1<0.5,1,'first'); %find the first point where T<0.5
        if isempty(J) %if all points are >0.5...
            mt=0; %...there is no median time
        else 
            %compute the median time by linear interpolation.
            p=polyfit([t1(I) t1(J)],[T1(I) T1(J)],1);
            mt=(0.5-p(2))/p(1);
            str2=['Median time ' num2str(mt)]; %string for LEGEND function
        end
    else
        mt=t1(T1==0.5);
        str2=['Median time ' num2str(mt)]; %string for LEGEND function
    end
    
    
 %%%%%%%%%% PLOTS %%%%%%%%%%%
    hFig=gcf;
    hold on
    S2 = stairs(t1(1:lowcend),lowc(1:lowcend),'k:','LineWidth',1); %lower confidence interval curve
    S2u = stairs(t1,upc,'k:','LineWidth',1); %upper confidence interval curve

    
    S1=stairs(t1,T1,'LineWidth',3); %Kaplan-Meier survival function
    S1.Color = P_color(1,:);
    [xb,yb] = stairs(t1,lowc);
    [xa,ya] = stairs(t1,upc);
    xx = [xb; flipud(xa)];
    inBetween = [yb; flipud(ya)];
    S2a = fill(xx, inBetween, P_color(1,:),'facealpha',.1,'EdgeColor','none');  

    if ~isempty(table12) %if there are censored data...
        S4=[];%plot(xcg,ycg,'rd','MarkerFaceColor','r','Markersize',7);
    else
        S4=[];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        OPTIMAL TMZ TREATMENT      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
 x = x2;

%Input Error handling
p=inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'alpha',P_alpha, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'cflag',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'flag',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x,varargin{:});
alpha=p.Results.alpha; cflag=p.Results.cflag; flag=p.Results.flag;
clear p
assert(all(x(:,2)==0 | x(:,2)==1),'Warning: all x(:,2) values must be 0 or 1')


%sort data by survival time
x=sortrows(x,1);

%table of patients observed for each survival time
%the TABULATE function sets up this matrix:
%table1=[time count percent(on total)]

table1=[0 size(x,1) 1; tabulate(x(:,1))];

%if all observed time are integers remove not observed time added by
%TABULATE function
table1(table1(:,3)==0,:)=[];

%Table of censored data
table12=tabulate(x(x(:,2)==1));
if ~isempty(table12)
    % remove not observed time added by TABULATE function
    table12(table12(:,3)==0,:)=[];
    % setup the vector of the censored data
    [cens,loc]=ismember(table1(:,1),table12(:,1)); %find censored data
end    

%the percents stored in the the third column are unuseful;
%so, place in the third column how many subjects are still alive at the
%beginning of the i-th interval.
a1=[table1(1,2); -1.*table1(2:end,2)];
table1(:,3)=cumsum(a1); table1(2:end,3)=table1(1:end-1,3);

%number of deaths in the intervals (don't take in account the censored
%data)
if ~isempty(table12)
    table1(cens,2)=table1(cens,2)-table12(loc(cens),2);
end
%finally, delete the first row that is now useless
table1(1,:)=[];


t1=[0;table1(:,1)]; %this is the x variable (time);
T1=[1;cumprod(1-(table1(:,2)./table1(:,3)))]; %this is the y variable (survival function)


if flag %if this function was not called by LOGRANK function
    
    %compute the standard error of the survival function
    SE=[0;T1(2:end).*sqrt(cumsum(table1(:,2)./(table1(:,3).* ...
        (table1(:,3)-table1(:,2)))))];
    
end

%censored data plotting
if ~isempty(table12) 
    %if there are censored data after max(t1), add a new cell into the t1,
    %T1 and SE arrays
    if table12(end,1)>=t1(end,1)
        t1(end+1,1)=table12(end,1)+1;
        T1(end+1,1)=T1(end,1);
        if flag %if this function was not called by LOGRANK function
            SE(end+1,1)=SE(end,1);
        end
    end
    if ~cflag
        %vectors preallocation
        xcg=zeros(1,sum(table12(:,2))); ycg=xcg; J=1;
        %for each censored data into the i-th time interval...
        for I=1:size(table12,1)
            %compute how many position into the array they must occupy
            JJ=J+table12(I,2)-1;
            %find the correct time interval in which censored data must be
            %placed
            A=find(t1<=table12(I,1),1,'last');
            B=find(t1>table12(I,1),1,'first');
            %equally divide this interval
            int=linspace(table12(I,1),t1(B,1),table12(I,2)+2);
            %put all in the vectors of the plotting variables
            xcg(J:JJ)=int(2:end-1);
            ycg(J:JJ)=T1(A);
            %update the counter
            J=JJ+1;
        end
    else
        xcg=table1(table1(:,2)==0,1);
        ycg=T1(table1(:,2)==0);
    end
else
    if ~flag %if this function was called by LOGRANK function
        xcg=[]; ycg=[];
    end
end
%compute the hazard rate
c1=T1.*numel(x);
c2=-(diff(log(c1(1:end-1)))./diff(t1(1:end-1)));
lambda2=mean(c2(c2~=0));

if flag %if this function was not called by LOGRANK function
    %compute the (1-alpha)*100% confidence interval curves
    cv=realsqrt(2)*erfcinv(alpha); %critical value
    %lower curve (remember that: the lower curve values can't be negative)
    lowc=max(0,T1-SE.*cv);
    %if the lower curve reaches the 0 earlier than survival function, trim the
    %data.
    if isequal(lowc(end-1:end),[0; 0])
        lowcend=find(lowc==0,1,'first');
    else
        lowcend=length(lowc);
    end
    %upper curve (remember that the upper curve values can't be >1)
    upc=min(1,T1+SE.*cv);
    %eventually, correct the data.
    if isequal(upc(end),1) 
        cupend=find(upc<1,1,'last');
        upc(cupend:end)=upc(cupend);
    end

    %compute the median survival time (if exist...)
    if isempty(T1(T1==0.5)) %if there is not a point where T=0.5...
        I=find(T1>0.5,1,'last'); %find the first point where T>0.5
        J=find(T1<0.5,1,'first'); %find the first point where T<0.5
        if isempty(J) %if all points are >0.5...
            mt2=0; %...there is no median time
        else 
            %compute the median time by linear interpolation.
            p=polyfit([t1(I) t1(J)],[T1(I) T1(J)],1);
            mt2=(0.5-p(2))/p(1);
            str2=['Median time ' num2str(mt)]; %string for LEGEND function
        end
    else
        mt2=t1(T1==0.5);
        str2=['Median time ' num2str(mt)]; %string for LEGEND function
    end
    
    
 %%%%%%%%%% PLOTS %%%%%%%%%%%
    hFig=gcf;
    hold on
    S2 = stairs(t1(1:lowcend),lowc(1:lowcend),'k:','LineWidth',1); %lower confidence interval curve
    S2u = stairs(t1,upc,'k:','LineWidth',1); %upper confidence interval curve


    S1b=stairs(t1,T1,'LineWidth',3); %Kaplan-Meier survival function
    S1b.Color = P_color(2,:);
    [xb,yb] = stairs(t1,lowc);
    [xa,ya] = stairs(t1,upc);
    xx = [xb; flipud(xa)];
    inBetween = [yb; flipud(ya)];
    S2a = fill(xx, inBetween, P_color(2,:),'facealpha',.1,'EdgeColor','none');  
    

    if ~isempty(table12) %if there are censored data...
        S4=[];%plot(xcg,ycg,'rd','MarkerFaceColor','r','Markersize',7);
    else
        S4=[];
    end
   
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           BEST TREATMENT        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
x = x3;

%Input Error handling
p=inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'alpha',P_alpha, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'cflag',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'flag',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x,varargin{:});
alpha=p.Results.alpha; cflag=p.Results.cflag; flag=p.Results.flag;
clear p
assert(all(x(:,2)==0 | x(:,2)==1),'Warning: all x(:,2) values must be 0 or 1')

%string for LEGEND function

% Sort data by survival
x=sortrows(x,1);

% table of patients observed for each survival time
table1=[0 size(x,1) 1; tabulate(x(:,1))];
% remove not observed time added 
table1(table1(:,3)==0,:)=[];
%Table of censored data
table12=tabulate(x(x(:,2)==1));

if ~isempty(table12)
    table12(table12(:,3)==0,:)=[];                  % remove not observed time
    [cens,loc]=ismember(table1(:,1),table12(:,1));  % find censored data
end    

% 3 column = # subjects are still alive at the beginning of the i-th interval.
a1=[table1(1,2); -1.*table1(2:end,2)];

% number of deaths in the intervals (don't take in account the censored data)
table1(:,3)=cumsum(a1); table1(2:end,3)=table1(1:end-1,3);

if ~isempty(table12)
    table1(cens,2)=table1(cens,2)-table12(loc(cens),2);
end

table1(1,:)=[];                               % delete the first row
t1=[0;table1(:,1)];                           % x variable (time);
T1=[1;cumprod(1-(table1(:,2)./table1(:,3)))]; % y variable (survival function)

if flag 
    %compute the standard error of the survival function
    SE=[0;T1(2:end).*sqrt(cumsum(table1(:,2)./(table1(:,3).* ...
        (table1(:,3)-table1(:,2)))))];
end

%censored data plotting
if ~isempty(table12) 
    %if there are censored data after max(t1), add a new cell into the t1,
    %T1 and SE arrays
    if table12(end,1)>=t1(end,1)
        t1(end+1,1)=table12(end,1)+1;
        T1(end+1,1)=T1(end,1);
        if flag 
            SE(end+1,1)=SE(end,1);
        end
    end
    if ~cflag
        %vectors preallocation
        xcg=zeros(1,sum(table12(:,2))); ycg=xcg; J=1;
        
        
        for I=1:size(table12,1) %for each censored data into the i-th time interval
            
            %compute how many position into the array they must occupy
            JJ=J+table12(I,2)-1;
            
            % find the correct time interval in which censored data must be placed
            A=find(t1<=table12(I,1),1,'last');
            B=find(t1>table12(I,1),1,'first');
            
            %equally divide this interval
            int=linspace(table12(I,1),t1(B,1),table12(I,2)+2);
            
            %put all in the vectors of the plotting variables
            xcg(J:JJ)=int(2:end-1);
            ycg(J:JJ)=T1(A);
            
            %update the counter
            J=JJ+1;
        end
    else
        xcg=table1(table1(:,2)==0,1);
        ycg=T1(table1(:,2)==0);
    end
else
    if ~flag 
        xcg=[]; ycg=[];
    end
end

%  Hazard rate
c1=T1.*numel(x);
c2=-(diff(log(c1(1:end-1)))./diff(t1(1:end-1)));
lambda=mean(c2(c2~=0));

if flag 
    % compute confidence interval
    cv=realsqrt(2)*erfcinv(alpha); % critical value
    
    % lower curve
    lowc=max(0,T1-SE.*cv);
    % if the lower curve reaches the 0 earlier than survival function, trim the data
    if isequal(lowc(end-1:end),[0; 0])
        lowcend=find(lowc==0,1,'first');
    else
        lowcend=length(lowc);
    end
    
    % upper curve 
    upc=min(1,T1+SE.*cv);
    % eventually, correct the data.
    if isequal(upc(end),1) 
        cupend=find(upc<1,1,'last');
        upc(cupend:end)=upc(cupend);
    end

    % compute the median survival time (if exist...)
    if isempty(T1(T1==0.5)) %if there is not a point where T=0.5...
        I=find(T1>0.5,1,'last'); %find the first point where T>0.5
        J=find(T1<0.5,1,'first'); %find the first point where T<0.5
        if isempty(J) %if all points are >0.5...
            mt3=0; %...there is no median time
        else 
            %compute the median time by linear interpolation.
            p=polyfit([t1(I) t1(J)],[T1(I) T1(J)],1);
            mt3=(0.5-p(2))/p(1);
            str2=['Median time ' num2str(mt)]; %string for LEGEND function
        end
    else
        mt3=t1(T1==0.5);
        str2=['Median time ' num2str(mt)]; %string for LEGEND function
    end
 
    %%%%%% PLOTS %%%%%%%%%%
    hFig=gcf;

    hold on
    S2b=stairs(t1(1:lowcend),lowc(1:lowcend),'k:','LineWidth',3); %lower confidence interval curve
    stairs(t1,upc,'k:','LineWidth',3); %upper confidence interval curve
    S1c=stairs(t1,T1,'LineWidth',3); %Kaplan-Meier survival function
    
    S1c.Color = P_color(3,:);
    [xb,yb] = stairs(t1,lowc);
    [xa,ya] = stairs(t1,upc);
    xx = [xb; flipud(xa)];
    inBetween = [yb; flipud(ya)];
    S2ba = fill(xx, inBetween, P_color(3,:),'facealpha',.1,'EdgeColor','none');  
    

    if ~isempty(table12) %if there are censored data...
        S4=[];%plot(xcg,ycg,'rd','MarkerFaceColor','r','Markersize',7);
    else
        S4=[];
    end
    hold off

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET PLOTS DETAILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %set the axis properly
    xmax=max(t1)+1;
    axis([0 xmax 0 1.1]);
    axis square
    %add labels and legend
    ylabel('Estimated Survival Probability','FontName','Arial','FontSize',30,'FontWeight','Bold'); 
    xlabel('Time','FontName','Arial','FontSize',30,'FontWeight','Bold'); 
    
    str1=[num2str((1-alpha)*100) '% confidence interval'];
    

    str_hazard1 = ['Conventional (Mean survival: ' num2str(mt,'%4.2f') ' months)'];
    str_hazard2 = ['Optimal TMZ (Mean survival: ' num2str(mt2,'%4.2f') ' months)'];
    str_hazard3 = ['Strategy D     (Mean survival: ' num2str(mt3,'%4.2f') ' months)'];

    if mt
        if isempty(S4)
           ll = legend([S1 S1b S1c S2],str_hazard1, str_hazard2,str_hazard3,str1)
        else
           ll = legend([S1 S1b S1c S2],str_hazard1, str_hazard2,str_hazard3,str1)
        end
    else
        if isempty(S4)
           ll = legend([S1 S1b S1c S2],str_hazard1, str_hazard2,str_hazard3,str1)
        else
            ll =legend([S1 S1b S1c S2],str_hazard1, str_hazard2, str_hazard3,str1)
        end
    end
end
set(ll,'Location','EastOutside')
if nargout
    varargout(1)={table1};
    varargout(2)={table12};
    varargout(3)={t1};
    varargout(4)={T1};
    varargout(5)={xcg};
    varargout(6)={ycg};
    varargout(7)={lambda};
end

hold off
end


end

function Kaplan_Plots()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-clinical in silico trial
% Kaplan-Meier plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Import Data using Import Tool: select numeric matrix

load('Human14.mat')
[n, m] = size(humanbest14);

%%% Length of study
months = 7+17;

initial_Tumor = humanbest14(1,1);
time = linspace(0,months,n);

%%% Calculate cutoff based on a tumor growing for 4 months
cutoff_value = cutoff(4,initial_Tumor);

%%% Set solution for each patient
for pat = 1:m  
SOL_con(pat).x = time;
SOL_con(pat).y = humancon14(:,pat);

SOL_opt(pat).x = time;
SOL_opt(pat).y = humanopt14(:,pat);

SOL_best(pat).x = time;
SOL_best(pat).y = humanbest14(:,pat);
end


%%% Get data in Kaplan Form

[Time_death_con,censor_con] = transform_to_Kaplan(SOL_con,time(1),m,cutoff_value);
[Time_death_opt,censor_opt] = transform_to_Kaplan(SOL_opt,time(1),m,cutoff_value);
[Time_death_best,censor_best] = transform_to_Kaplan(SOL_best,time(1),m,cutoff_value);


%%% Plot Results using modified kmplot6
figure(1) 
set(gca,'LineWidth',3,'FontSize',30);
kmplot6([Time_death_con;censor_con]',[Time_death_opt;censor_opt]',[Time_death_best;censor_best]');
ax1 = gca; 

end


function Phenotypic_evolution()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-clinical in silico trial
% Phenotypic evolution plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load cell composition data
load('CellInfo7.mat')


%%% Get data dimentions
[n, m] = size(humanapngcon7);

time = [0 1 2 3 4 5 6 7];

%%% Color scheme
P_color = parula(6);

%%% Rename data
APNGconN = humanapngcon7;
MGMTconN = humanmgmtcon7;
alphaconN = humanalphacon7;

APNGoptN = humanapngopt7;
MGMToptN = humanmgmtopt7;
alphaoptN = humanalphaopt7;

APNGbestN = humanapngbest7;
MGMTbestN = humanmgmtbest7;
alphabestN = humanalphabest7;

for k = 1:m %re-scale by maximum value
    APNGconN(:,k) = APNGconN(:,k)/humaninfocon7(4,k);
    MGMTconN(:,k) = MGMTconN(:,k)/humaninfocon7(5,k);
    alphaconN(:,k) = alphaconN(:,k)/humaninfocon7(1,k);
    
    APNGoptN(:,k) = APNGoptN(:,k)/humaninfoopt7(4,k);
    MGMToptN(:,k) = MGMToptN(:,k)/humaninfoopt7(5,k);
    alphaoptN(:,k) = alphaoptN(:,k)/humaninfoopt7(1,k);

    APNGbestN(:,k) = APNGbestN(:,k)/humaninfobest7(4,k);
    MGMTbestN(:,k) = MGMTbestN(:,k)/humaninfobest7(5,k);
    alphabestN(:,k) = alphabestN(:,k)/humaninfobest7(1,k);
end

%%% Quantile normalization of data
APNGconN = quantilenorm(APNGconN);
MGMTconN = quantilenorm(MGMTconN);
alphaconN = quantilenorm(alphaconN);

APNGoptN = quantilenorm(APNGoptN);
MGMToptN = quantilenorm(MGMToptN);
alphaoptN = quantilenorm(alphaoptN);

APNGbestN = quantilenorm(APNGbestN);
MGMTbestN = quantilenorm(MGMTbestN);
alphabestN = quantilenorm(alphabestN);


for i= 1:n 
    
%%% Calculate mean values
APNGcon(i) = mean(APNGconN(i,:));
MGMTcon(i) = mean(MGMTconN(i,:));
alphacon(i) = mean(alphaconN(i,:));

APNGopt(i) = mean(APNGoptN(i,:));
MGMTopt(i) = mean(MGMToptN(i,:));
alphaopt(i) = mean(alphaoptN(i,:));

APNGbest(i) = mean(APNGbestN(i,:));
MGMTbest(i) = mean(MGMTbestN(i,:));
alphabest(i) = mean(alphabestN(i,:));

%%% Calculate standard deviation
sdAPNGcon(i) = std2(APNGconN(i,:));
sdMGMTcon(i) = std2(MGMTconN(i,:));
sdalphacon(i) = std2(alphaconN(i,:));

sdAPNGopt(i) = std2(APNGoptN(i,:));
sdMGMTopt(i) = std2(MGMToptN(i,:));
sdalphaopt(i) = std2(alphaoptN(i,:));

sdAPNGbest(i) = std2(APNGbestN(i,:));
sdMGMTbest(i) = std2(MGMTbestN(i,:));
sdalphabest(i) = std2(alphabestN(i,:));

end


%%% Plot results

figure(1)
hold on
set(gca,'LineWidth',7,'FontSize',50,'FontWeight','normal','FontName','Helvetica')
ax = gca;

errorbar(time,APNGcon,sdAPNGcon,sdAPNGcon,'k>','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(1,:))
errorbar(time,APNGopt,sdAPNGopt,sdAPNGopt,'kd','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(2,:))
errorbar(time,APNGbest,sdAPNGbest,sdAPNGbest,'ks','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(3,:))

plot(time,APNGcon,'-','color', P_color(1,:),'LineWidth',5)
plot(time,APNGopt,':','color', P_color(2,:),'LineWidth',5)
plot(time,APNGbest,'--','color', P_color(3,:),'LineWidth',5)

errorbar(time,APNGcon,sdAPNGcon,sdAPNGcon,'k>','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(1,:))
errorbar(time,APNGopt,sdAPNGopt,sdAPNGopt,'kd','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(2,:))
errorbar(time,APNGbest,sdAPNGbest,sdAPNGbest,'ks','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(3,:))


hold off
xlabel('Treatment Cycle') % x-axis label

axis([0  7 0.3 0.65])
ylabel('Scaled Average Value')
title('APNG')
legend({'Conventional', 'Optimal TMZ', 'Strategy D'},'FontSize',50,'Location','northeast')
legend
legend boxoff;


figure(2)
hold on
set(gca,'LineWidth',7,'FontSize',50,'FontWeight','normal','FontName','Helvetica')
ax = gca;

errorbar(time,MGMTcon,sdMGMTcon,sdMGMTcon,'k>','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(1,:))
errorbar(time,MGMTopt,sdMGMTopt,sdMGMTopt,'kd','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(2,:))
errorbar(time,MGMTbest,sdMGMTbest,sdMGMTbest,'ks','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(3,:))

plot(time,MGMTcon,'-','color', P_color(1,:),'LineWidth',5)
plot(time,MGMTopt,':','color', P_color(2,:),'LineWidth',5)
plot(time,MGMTbest,'--','color', P_color(3,:),'LineWidth',5)

errorbar(time,MGMTcon,sdMGMTcon,sdMGMTcon,'k>','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(1,:))
errorbar(time,MGMTopt,sdMGMTopt,sdMGMTopt,'kd','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(2,:))
errorbar(time,MGMTbest,sdMGMTbest,sdMGMTbest,'ks','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(3,:))



hold off
xlabel('Treatment Cycle') % x-axis label
axis([0  7 0.3 0.65])
ylabel('Scaled Average Value')
title('MGMT')
legend({'Conventional', 'Optimal TMZ', 'Strategy D'},'FontSize',50,'Location','southeast')
legend
legend boxoff;


figure(3)
hold on
set(gca,'LineWidth',7,'FontSize',50,'FontWeight','normal','FontName','Helvetica')
ax = gca;

errorbar(time,alphacon,sdalphacon,sdalphacon,'k>','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(1,:))
errorbar(time,alphaopt,sdalphaopt,sdalphaopt,'kd','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(2,:))
errorbar(time,alphabest,sdalphabest,sdalphabest,'ks','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(3,:))

plot(time,alphacon,'-','color', P_color(1,:),'LineWidth',5)
plot(time,alphaopt,':','color', P_color(2,:),'LineWidth',5)
plot(time,alphabest,'--','color', P_color(3,:),'LineWidth',5)

errorbar(time,alphacon,sdalphacon,sdalphacon,'k>','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(1,:))
errorbar(time,alphaopt,sdalphaopt,sdalphaopt,'kd','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(2,:))
errorbar(time,alphabest,sdalphabest,sdalphabest,'ks','MarkerSize',25,'LineWidth',3,'MarkerFaceColor', P_color(3,:))



hold off
xlabel('Treatment Cycle') % x-axis label
axis([0  7 0.5 1])
ylabel('Scaled Average Value')
title('Proliferation')
legend({'Conventional', 'Optimal TMZ', 'Strategy D'},'FontSize',50,'Location','southeast')
legend
legend boxoff;

end
