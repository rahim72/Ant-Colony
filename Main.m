%%
%  meghdar iteration be sorat pesh farz 300 dar nazar grefte shode ast
%  tedad run baraber 25 run
% Data Set :[  kroA100 - kroB100- kroB150- kroA150- kroB200 -kroA200- XQF131 ]
%   Data Set >> XQF131 = 131 city
%
%%
clc;
clear all;
warning off;
close all;
%% Choice Algorithm
%%
disp('Start Ant Colony Algorithm :');
%% Qusetion For Data set Selection
choice2 = questdlg('Select  a DataSet!', ...
    'Select  Dataset', ...
    'Select Data Set','XQF131','No thank you');
% Handle response
switch choice2
%     case 'Tsp Point'
%         disp([choice2 ' , Selected.'])
%         load Mydata
%         dataset=Mydata;
%         sys.numcity=20;
    case 'Select Data Set'
        disp([choice2 ' , Selected.'])
        choice = selectdataset;
        %% Call Function Selectdataset And Choice Type Data Set
        if strcmp(choice,'kroA100')
            load kroA100
            disp(['Select Data Set :', num2str(choice)])
            dataset=kroA100;
            sys.numcity=100;
        elseif strcmp(choice,'kroB100')
            load kroB100
            disp(['Select Data Set :', num2str(choice)])
            dataset=kroB100;
            sys.numcity=100;
        elseif strcmp(choice,'kroB150')
            load kroB150
            disp(['Select Data Set :', num2str(choice)])
            dataset=kroB150;
            sys.numcity=150;
        elseif strcmp(choice,'kroA150')
            load kroA150
            disp(['Select Data Set :', num2str(choice)])
            dataset=kroA150;
            sys.numcity=150;
        elseif strcmp(choice,'kroB200')
            load kroB200
            disp(['Select Data Set :', num2str(choice)])
            dataset=kroB200;
            sys.numcity=200;
        elseif strcmp(choice,'kroA200')
            load kroA200
            disp(['Select Data Set :', num2str(choice)])
            dataset=kroA200;
            sys.numcity=200;
        end
        
    case 'XQF131'
        disp([choice2 ' , Selected.'])
        load XQF131
        dataset=XQF131;
        sys.numcity=131;
        fig=imread('city_point.JPG');
        
end

%% ACO Parameters
%%
prompt={ ' Maximum Number of Iterations', ' Number of Ants (Population Size)', 'hromone Exponential Weight', ' Heuristic Exponential Weight','Evaporation Rate'...
    'Maximum Number of Run Algorithm','Value Q'};
title='Enter a Value';
nline=([1 40;1 40;1 40;1 40;1 40;1 40;1 40]);
default={ '300', '90', '1', '1', '0.6','25','1'};
ev=inputdlg(prompt,title,nline,default);
aa=ev(1,:);aa=str2num(aa{:});                        Iteration=aa;                                                                              %  Maximum Number of Iterations
bb=ev(2,:);bb=str2num(bb{:});                   popant=bb;                                                                                %  Number of Ants (Population Size)
cc=ev(3,:);cc=str2num(cc{:});                       alpha=cc;                                                                                      % Phromone Exponential Weight
dd=ev(4,:);dd=str2num(dd{:});                     beta=dd;                                                                                        % Heuristic Exponential Weight
ff=ev(5,:);ff=str2num(ff{:});                             rho=ff;                                                                                              % Evaporation Rate
ee=ev(6,:);ee=str2num(ee{:});                     max_num_run=ee;
rr=ev(7,:);rr=str2num(rr{:});                          Q=rr;

%% Create Cost Functon And Set Parameter
%%
Matcost=zeros(size(dataset,1),size(dataset,1));
numbercity=sys.numcity;
numrow=size(dataset,1);
Matcost2=Matcost;
xi=dataset(:,1)';
yi=dataset(:,2)';

for ii=1:numrow-1
    for jj=ii+1:numrow
        Matcost2(ii,jj)=sqrt((xi(ii)-xi(jj))^2+(yi(ii)-yi(jj))^2);
        Matcost2(jj,ii)=Matcost2(ii,jj);
    end
end
sys.x=xi;                                                                                                                                                                                      % X point For City
sys.y=yi;                                                                                                                                                                                      % Y Point For City
sys.D=Matcost2;
eta=1./sys.D;                                                                                                                                                                           % Heuristic Information Matrix
tau0=10*Q/(numbercity*mean(sys.D(:)));
%% 25 Loop >>>>  Averge Response
%.%
tic
for num_run=1:max_num_run
    
    CostFunction=@(tour) Cost_Function(tour,sys);                                                                                              % Cost Function
    % Initial Phromone
    %% Initialization
    tau=tau0*ones(numbercity,numbercity);                                                                                        % Phromone Matrix
    BestCost=zeros(Iteration,1);                                                                                                            % Array to Hold Best Cost Values
    em_cell.Tour=[];
    em_cell.Cost=[];
    ant=repmat(em_cell,popant,1);                                                                                                       % Repeat copies of array
    BestSoFar.Cost=inf;                                                                                                                           % Best Ant
    
    for iter=1:Iteration
        for pop=1:popant
            ant(pop).Tour=randi([1 numbercity]);
            for l=2:numbercity
                
                ii=ant(pop).Tour(end);
                
                P=tau(ii,:).^alpha.*eta(ii,:).^beta;
                P(ant(pop).Tour)=0;
                P=P/sum(P);
                %%
                r=rand;
                sum_p=cumsum(P);                                                                                                                   % Cumulative sum
                jj=find(r<=sum_p,1,'first');
                %
                ant(pop).Tour=[ant(pop).Tour jj];
                
                %%
            end
            
            ant(pop).Cost=CostFunction(ant(pop).Tour);
            
            if ant(pop).Cost<BestSoFar.Cost
                BestSoFar=ant(pop);
            end
            
        end
        
        % Update Phromones
        for pop=1:popant
            
            tour=ant(pop).Tour;
            
            tour=[tour tour(1)]; %#ok
            
            for l=1:numbercity
                
                ii=tour(l);
                jj=tour(l+1);
                
                tau(ii,jj)=tau(ii,jj)+Q/ant(pop).Cost;
                
            end
            
        end
        tau=(1-rho)*tau;
        BestCost(iter)=BestSoFar.Cost;                                                                                                           % Store Best Cost
        
        
        disp(['Iteration ' num2str(iter) ': Best Cost = ' num2str(BestCost(iter))]);                                        % Show Iteration Information
        if strcmp(choice2,'XQF131')
            figure(11)
            ifig=imshow(fig);
            % Plot  TSP Map Proble
        end
        figure(num_run);
        tour1=[BestSoFar.Tour BestSoFar.Tour(1)];
        plot(sys.x(tour1),sys.y(tour1),'r-o',...
            'MarkerSize',4,...
            'MarkerFaceColor','k',...
            'LineWidth',1.5);
        axis equal;                                                                                                                                                                   % Set axis limits and aspect ratios
        grid on;                                                                                                                                                                          % Display or hide axes grid lines
        %% X
        %%
        xmin = min(sys.x);
        xmax = max(sys.x);
        dx = xmax - xmin;
        xmin = floor((xmin - 0.01*dx)/2)*2;                                                                                                           % Round toward negative infinity
        xmax = ceil((xmax + 0.01*dx)/2)*2;                                                                                                           % Round toward positive infinity
        xlim([xmin xmax]);                                                                                                                                                 % Set or query x-axis limits
        %% Y
        %%
        ymin = min(sys.y);
        ymax = max(sys.y);
        dy = ymax - ymin;
        ymin = floor((ymin - 0.1*dy)/5)*5;
        ymax = ceil((ymax + 0.1*dy)/5)*5;
        ylim([ymin ymax]);
        %
        pause(0.001);
        
    end
    ListBestCost(num_run,1)=BestCost(iter);
    %% Results
end


Average=mean(ListBestCost);
disp (' ' )
disp (' ' )
disp ('--------------------------------------------------------------------------------------------' )
disp ('List Best Cost : ')
disp (ListBestCost)
disp ('--------------------------------------------------------------------------------------------' )
disp (' ' )
disp ('--------------------------------------------------------------------------------------------' )
disp (['Average : ', num2str(Average)])
disp (['Time : ', num2str(toc)])
disp ('--------------------------------------------------------------------------------------------' )

bar(ListBestCost,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)

for ii=1:max_num_run
    
    figure('color','b','menuBar','none','position',[350 200 450 400],'Resize','off')
    plot(BestCost,'LineWidth',2);
    xlabel('Iteration');
    ylabel('Best Cost');
    grid on;
    
end


clear j k l Matcost nn title xi yi rr i it ev nline Iteration max_num_run num_run numbercity numrow popant prompt bb aa cc ee ff choice2 dd default




