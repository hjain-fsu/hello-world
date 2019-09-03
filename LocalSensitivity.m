function LocalSens_alpha()

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local values of alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = linspace(0.0001,1.55,4);

%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%

K = 5.8*10^6;
TMZ2 = 500;

Sa = 1000;
Sm = Sa;

Na = 10;
Mt = round(1/ht);

lambdaA = 0.0057;
lambdaM = 0.0043;

t = linspace(0,2*60,Mt);

clear y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when both APNG and MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = Population_Model(Na,Mt,t,TMZ2,alpha(1),K,Sa,Sm,lambdaA,lambdaM);
    y1(1,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(2),K,Sa,Sm,lambdaA,lambdaM);
    y1(1,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(3),K,Sa,Sm,lambdaA,lambdaM);
    y1(1,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(4),K,Sa,Sm,lambdaA,lambdaM);
    y1(1,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only APNG is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha(1),K,Sa,0,lambdaA,lambdaM); 
    y1(2,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(2),K,Sa,0,lambdaA,lambdaM);
    y1(2,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(3),K,Sa,0,lambdaA,lambdaM);
    y1(2,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(4),K,Sa,0,lambdaA,lambdaM);
    y1(2,4) = F(end,2)/10^5;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only MGMT is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    F = Population_Model(Na,Mt,t,TMZ2,alpha(1),K,0,Sm,lambdaA,lambdaM); 
    y1(3,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(2),K,0,Sm,lambdaA,lambdaM);
    y1(3,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(3),K,0,Sm,lambdaA,lambdaM);
    y1(3,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(4),K,0,Sm,lambdaA,lambdaM); 
    y1(3,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when neither APNG nor MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha(1),K,0,0,lambdaA,lambdaM);
    y1(4,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(2),K,0,0,lambdaA,lambdaM);
    y1(4,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(3),K,0,0,lambdaA,lambdaM);
    y1(4,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha(4),K,0,0,lambdaA,lambdaM);
    y1(4,4) = F(end,2)/10^5;

%%% Transform data to percentage
data = y1*100; 

%%% x-axis for histogram
xval = 1:4;

%%% Create color scheme
cmap = hot(6);


%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%%% Plot data
hBar = bar(xval,data); 

%%% Assign color scheme
set(gca,'Color',[0.91 0.91 0.91])
set(hBar(1), 'FaceColor', cmap(2,:))       
set(hBar(2), 'FaceColor', cmap(3,:))      
set(hBar(3), 'FaceColor', cmap(4,:))  
set(hBar(4), 'FaceColor', cmap(1,:))    

%%% Set fonts
set(gca,'LineWidth',1.25,'FontSize',30,'FontWeight','normal','FontName','Helvetica')

%%% Set Y-Axis Limits
set(gca, 'YLim', [0 100]) 


%%% Set Labels
ylabel({'Survival';'Percentage'})
legend('\alpha=0.0001','\alpha=0.5167','\alpha=1.0334','\alpha=1.55', 'Location', 'NE')
xtklbl{1} = sprintf('+\n+');
xtklbl{2} = sprintf('+\n-');
xtklbl{3} = sprintf('-\n+');
xtklbl{4} = sprintf('-\n-');

set(gca,'XTickLabel', {'','','',''});
[hx,hy] = format_ticks(gca,xtklbl);
text(0.3, -5, 'APNG', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')
text(0.3, -11, 'MGMT', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')

hold off
end



function LocalSens_S_MGMT()

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local values of S_MGMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sm = [10 100 1000 10000];

%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.453;

K = 5.8*10^6;
TMZ2 = 500;

Sa = 1000;

Na = 10;
Mt = round(1/ht);

lambdaA = 0.0057;
lambdaM = 0.0043;


t = linspace(0,2*60,Mt);

clear y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when both APNG and MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm(1),lambdaA,lambdaM);
    y1(1,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm(2),lambdaA,lambdaM);
    y1(1,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm(3),lambdaA,lambdaM);
    y1(1,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm(4),lambdaA,lambdaM);
    y1(1,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only APNG is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM); 
    y1(2,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM);
    y1(2,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM);
    y1(2,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM);
    y1(2,4) = F(end,2)/10^5;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only MGMT is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm(1),lambdaA,lambdaM); 
    y1(3,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm(2),lambdaA,lambdaM);
    y1(3,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm(3),lambdaA,lambdaM);
    y1(3,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm(4),lambdaA,lambdaM); 
    y1(3,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when neither APNG nor MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,4) = F(end,2)/10^5;

%%% Transform data to percentage
data = y1*100; 

%%% x-axis for histogram
xval = 1:4;

%%% Create color scheme
cmap = hot(6);


%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%%% Plot data
hBar = bar(xval,data); 

%%% Assign color scheme
set(gca,'Color',[0.91 0.91 0.91])
set(hBar(1), 'FaceColor', cmap(2,:))       
set(hBar(2), 'FaceColor', cmap(3,:))      
set(hBar(3), 'FaceColor', cmap(4,:))  
set(hBar(4), 'FaceColor', cmap(1,:))    

%%% Set fonts
set(gca,'LineWidth',1.25,'FontSize',30,'FontWeight','normal','FontName','Helvetica')

%%% Set Y-Axis Limits
set(gca, 'YLim', [0 100]) 


%%% Set Labels
ylabel({'Survival';'Percentage'})
legend('S_{MGMT}=10','S_{MGMT}=100','S_{MGMT}=1000','S_{MGMT}=10000', 'Location', 'NE')
xtklbl{1} = sprintf('+\n+');
xtklbl{2} = sprintf('+\n-');
xtklbl{3} = sprintf('-\n+');
xtklbl{4} = sprintf('-\n-');

set(gca,'XTickLabel', {'','','',''});
[hx,hy] = format_ticks(gca,xtklbl);
text(0.3, -5, 'APNG', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')
text(0.3, -11, 'MGMT', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')

hold off
end


function LocalSens_S_APNG()

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local values of S_APNG
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sa = [10 100 1000 10000];

%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.453;

K = 5.8*10^6;
TMZ2 = 500;

Sm = 1000;

Na = 10;
Mt = round(1/ht);

lambdaA = 0.0057;
lambdaM = 0.0043;


t = linspace(0,2*60,Mt);

clear y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when both APNG and MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(1),Sm,lambdaA,lambdaM);
    y1(1,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(2),Sm,lambdaA,lambdaM);
    y1(1,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(3),Sm,lambdaA,lambdaM);
    y1(1,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(4),Sm,lambdaA,lambdaM);
    y1(1,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only APNG is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(1),0,lambdaA,lambdaM); 
    y1(2,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(2),0,lambdaA,lambdaM);
    y1(2,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(3),0,lambdaA,lambdaM);
    y1(2,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa(4),0,lambdaA,lambdaM);
    y1(2,4) = F(end,2)/10^5;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only MGMT is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM); 
    y1(3,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM);
    y1(3,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM);
    y1(3,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM); 
    y1(3,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when neither APNG nor MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM);
    y1(4,4) = F(end,2)/10^5;

%%% Transform data to percentage
data = y1*100; 

%%% x-axis for histogram
xval = 1:4;

%%% Create color scheme
cmap = hot(6);


%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%%% Plot data
hBar = bar(xval,data); 

%%% Assign color scheme
set(gca,'Color',[0.91 0.91 0.91])
set(hBar(1), 'FaceColor', cmap(2,:))       
set(hBar(2), 'FaceColor', cmap(3,:))      
set(hBar(3), 'FaceColor', cmap(4,:))  
set(hBar(4), 'FaceColor', cmap(1,:))    

%%% Set fonts
set(gca,'LineWidth',1.25,'FontSize',30,'FontWeight','normal','FontName','Helvetica')

%%% Set Y-Axis Limits
set(gca, 'YLim', [0 100]) 


%%% Set Labels
ylabel({'Survival';'Percentage'})
legend('S_{APNG}=10','S_{APNG}=100','S_{APNG}=1000','S_{APNG}=10000', 'Location', 'NE')
xtklbl{1} = sprintf('+\n+');
xtklbl{2} = sprintf('+\n-');
xtklbl{3} = sprintf('-\n+');
xtklbl{4} = sprintf('-\n-');

set(gca,'XTickLabel', {'','','',''});
[hx,hy] = format_ticks(gca,xtklbl);
text(0.3, -5, 'APNG', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')
text(0.3, -11, 'MGMT', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')

hold off
end



function LocalSens_lambda_MGMT()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local values of lambda_MGMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambdaM = [1e-5 0.005 0.05 0.5];

%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.453;

K = 5.8*10^6;
TMZ2 = 500;

Sm = 1000;
Sa = 1000;

Na = 10;
Mt = round(1/ht);

lambdaA = 0.0057;


t = linspace(0,2*60,Mt);

clear y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when both APNG and MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA,lambdaM(1));
    y1(1,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA,lambdaM(2));
    y1(1,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA,lambdaM(3));
    y1(1,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA,lambdaM(4));
    y1(1,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only APNG is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM(1));
    y1(2,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM(2));
    y1(2,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM(3));
    y1(2,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA,lambdaM(4));
    y1(2,4) = F(end,2)/10^5;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only MGMT is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM(1));
    y1(3,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM(2));
    y1(3,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM(3));
    y1(3,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA,lambdaM(4));
    y1(3,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when neither APNG nor MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM(1));
    y1(4,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM(2));
    y1(4,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM(3));
    y1(4,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA,lambdaM(4));
    y1(4,4) = F(end,2)/10^5;

%%% Transform data to percentage
data = y1*100; 

%%% x-axis for histogram
xval = 1:4;

%%% Create color scheme
cmap = hot(6);


%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%%% Plot data
hBar = bar(xval,data); 

%%% Assign color scheme
set(gca,'Color',[0.91 0.91 0.91])
set(hBar(1), 'FaceColor', cmap(2,:))       
set(hBar(2), 'FaceColor', cmap(3,:))      
set(hBar(3), 'FaceColor', cmap(4,:))  
set(hBar(4), 'FaceColor', cmap(1,:))    

%%% Set fonts
set(gca,'LineWidth',1.25,'FontSize',30,'FontWeight','normal','FontName','Helvetica')

%%% Set Y-Axis Limits
set(gca, 'YLim', [0 100]) 


%%% Set Labels
ylabel({'Survival';'Percentage'})
legend('\lambda_{MGMT}=1e-5','\lambda_{MGMT}=0.005','\lambda_{MGMT}=0.05','\lambda_{MGMT}=0.5', 'Location', 'NE')
xtklbl{1} = sprintf('+\n+');
xtklbl{2} = sprintf('+\n-');
xtklbl{3} = sprintf('-\n+');
xtklbl{4} = sprintf('-\n-');

set(gca,'XTickLabel', {'','','',''});
[hx,hy] = format_ticks(gca,xtklbl);
text(0.3, -5, 'APNG', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')
text(0.3, -11, 'MGMT', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')

hold off
end


function LocalSens_lambda_APNG()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local values of lambda_APNG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambdaA = [1e-5 0.005 0.05 0.5];

%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.453;

K = 5.8*10^6;
TMZ2 = 500;

Sm = 1000;
Sa = 1000;

Na = 10;
Mt = round(1/ht);

lambdaM = 0.0043;


t = linspace(0,2*60,Mt);

clear y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when both APNG and MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA(1),lambdaM);
    y1(1,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA(2),lambdaM);
    y1(1,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA(3),lambdaM);
    y1(1,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,Sm,lambdaA(4),lambdaM);
    y1(1,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only APNG is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA(1),lambdaM);
    y1(2,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA(2),lambdaM);
    y1(2,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA(3),lambdaM);
    y1(2,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,Sa,0,lambdaA(4),lambdaM);
    y1(2,4) = F(end,2)/10^5;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when only MGMT is expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA(1),lambdaM);
    y1(3,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA(2),lambdaM);
    y1(3,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA(3),lambdaM);
    y1(3,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,Sm,lambdaA(4),lambdaM);
    y1(3,4) = F(end,2)/10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values when neither APNG nor MGMT are expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA(1),lambdaM);
    y1(4,1) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA(2),lambdaM);
    y1(4,2) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA(3),lambdaM);
    y1(4,3) = F(end,2)/10^5;
    F = Population_Model(Na,Mt,t,TMZ2,alpha,K,0,0,lambdaA(4),lambdaM);
    y1(4,4) = F(end,2)/10^5;

%%% Transform data to percentage
data = y1*100; 

%%% x-axis for histogram
xval = 1:4;

%%% Create color scheme
cmap = hot(6);


%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%%% Plot data
hBar = bar(xval,data); 

%%% Assign color scheme
set(gca,'Color',[0.91 0.91 0.91])
set(hBar(1), 'FaceColor', cmap(2,:))       
set(hBar(2), 'FaceColor', cmap(3,:))      
set(hBar(3), 'FaceColor', cmap(4,:))  
set(hBar(4), 'FaceColor', cmap(1,:))    

%%% Set fonts
set(gca,'LineWidth',1.25,'FontSize',30,'FontWeight','normal','FontName','Helvetica')

%%% Set Y-Axis Limits
set(gca, 'YLim', [0 100]) 


%%% Set Labels
ylabel({'Survival';'Percentage'})
legend('\lambda_{APNG}=1e-5','\lambda_{APNG}=0.005','\lambda_{APNG}=0.05','\lambda_{APNG}=0.5', 'Location', 'NE')
xtklbl{1} = sprintf('+\n+');
xtklbl{2} = sprintf('+\n-');
xtklbl{3} = sprintf('-\n+');
xtklbl{4} = sprintf('-\n-');

set(gca,'XTickLabel', {'','','',''});
[hx,hy] = format_ticks(gca,xtklbl);
text(0.3, -5, 'APNG', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')
text(0.3, -11, 'MGMT', 'HorizontalAlignment','center','FontSize',20,'FontWeight','normal','FontName','Helvetica')

hold off
end