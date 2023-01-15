%-------------------------------------------------------------------------%
% Individual Parameters
%-------------------------------------------------------------------------%
N=1000; %Population size
p0=.01; %Probabilities of spread upon contact
p1=.01;
m1=1;
ms=5;

mr=20;
pd=.9;
sim_days=150; %Simulation length (days)

%-------------------------------------------------------------------------%
% Population Parameters
%-------------------------------------------------------------------------%
ages=round(90*rand(N,1)); %Age distribution
has_home=100*rand(N,1)<=100; %Homeowners
stayshome=100*rand(N,1)<=80; %Individuals who largely isolate
goestowork=100*rand(N,1)<=60; % 60% working
commutes=100*rand(N,1)<30; % 30% commuting
goestoschool=100*rand(N,1)<=25; % 25% in school
isinsmallgroup=100*rand(N,1)<=40; % 40% in some small, in-person group

%-------------------------------------------------------------------------%
% Group Parameters
%-------------------------------------------------------------------------%
homemean=6; %average home size
homesdv=2; %deviation in home size
homegroups=generate_groups(homemean,homesdv,has_home,N);
%no social distancing condition for home groups bc home groups unaffected

workmean=25;
worksdv=10;
workgroups=generate_groups(workmean,worksdv,goestowork,N);
mask_workgroups=(workgroups>max(workgroups)/2); %mask used to eliminate portion of all workgroups. Proportion eliminated varies with group
social_distanced_workgroups=workgroups; 
social_distanced_workgroups(mask_workgroups)=0;

schoolmean=50;
schoolsdv=40;
schoolgroups=generate_groups(schoolmean,schoolsdv,goestoschool,N);
mask_schoolgroups=(schoolgroups>max(schoolgroups)/500);
social_distanced_schoolgroups=schoolgroups;
social_distanced_schoolgroups(mask_schoolgroups)=0;

smallrecmean=10;
smallrecsdv=5;
smallrecgroups=generate_groups(smallrecmean,smallrecsdv,isinsmallgroup,N);
mask_smallrecgroups=(smallrecgroups>max(smallrecgroups)/50);
social_distanced_smallrecgroups=smallrecgroups;
social_distanced_smallrecgroups(mask_smallrecgroups)=0;

commutemean=50;
commutesdv=10;
commutegroups=generate_groups(commutemean,commutesdv,commutes,N);
mask_commutegroups=(commutegroups>max(commutegroups)/15);
social_distanced_commutegroups=commutegroups;
social_distanced_commutegroups(mask_commutegroups)=0;

%-------------------------------------------------------------------------%
% Simulation Initial Conditions
%-------------------------------------------------------------------------%
isinfected=100*rand(N,1)<=100*p0; %establishes case zeros 
infectionday=zeros(N,1);
issick=zeros(N,1);
isimmune=zeros(N,1);
isdead=zeros(N,1);
adj=zeros(N,N);

ninfected=zeros(sim_days,1);
%establish disease metric vectors for later assignment
nsick=zeros(sim_days,1);
nimmune=zeros(sim_days,1);
ndead=zeros(sim_days,1);
nvulnerable=zeros(sim_days,1);
health_capacity=zeros(sim_days,1);

%-------------------------------------------------------------------------%
% Simulation Loop
%-------------------------------------------------------------------------%
for currentday=1:sim_days

%condition for authorities to advise social distancing/stay-home orders
%set groups to socially distanced
if sum(issick)>.005*N 
    social_distancing='[Active]';
workgroups=social_distanced_workgroups;
schoolgroups=social_distanced_schoolgroups;
smallrecgroups=social_distanced_smallrecgroups;
commutegroups=social_distanced_commutegroups;
largerecgroups=zeros(N,1);
else
    social_distancing='[Inactive]';
isinbiggroup=100*rand(N,1)<=50;
largerecmean=40;
largerecsdv=15;
largerecgroups=generate_groups(largerecmean,largerecsdv,isinbiggroup,N);
end

%calculate infection status and time of infection based on group exposures
[isinfected,infectionday]=infect(commutegroups,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups);
[isinfected,infectionday]=infect(schoolgroups,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups);
doeserrands=100*rand(N,1)<=14.3; %errans unaffected no matter what bc errands are essential
errandmean=10;
errandsdv=5;
errandgroups=generate_groups(errandmean,errandsdv,doeserrands,N);
[isinfected,infectionday]=infect(errandgroups,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups);
[isinfected,infectionday]=infect(smallrecgroups,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups);
[isinfected,infectionday]=infect(largerecgroups,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups);
[isinfected,infectionday]=infect(homegroups,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups);

[isinfected,issick,isimmune,isdead]=disease(isinfected,infectionday,issick,isimmune,isdead,currentday,N,p1,m1,ms,mr,pd,ages);

%update individual exposure adjacency
adj=update_adjacency(homegroups,adj);
adj=update_adjacency(workgroups,adj);
adj=update_adjacency(schoolgroups,adj);
adj=update_adjacency(smallrecgroups,adj);
adj=update_adjacency(largerecgroups,adj);
adj=update_adjacency(errandgroups,adj);

nvulnerable(currentday)=N-sum(isinfected)-sum(isimmune); %core metrics for graph
ninfected(currentday)=sum(isinfected);
nsick(currentday)=sum(issick);
nimmune(currentday)=sum(isimmune);
ndead(currentday)=sum(isdead);

network=graph(adj); %updating infection web for animation
drawnow limitrate nocallbacks
h=plot(network,'Linewidth',.000000000001);
title(['Day #' num2str(currentday) '    Total Cases: ' num2str(ninfected(currentday)) '    Healthcare System Capacity: ' num2str(round(100*.125*nsick(currentday)/(.003*N))) '%' '  Stay Home Orders: ' num2str(social_distancing)]);
nodelist_infected=find(isinfected);
highlight(h,nodelist_infected,'NodeColor','r')
nodelist_immune=find(isimmune);
highlight(h,nodelist_immune,'NodeColor','y')
end

figure
plot(nvulnerable)
title('Population Vulnerability')
legend('Vulnerable')
ylabel('Vulnerable Individuals')
xlabel('Days')

figure
plot(ninfected) 
hold on
plot(nsick)
plot(nimmune)
plot(ndead)
hold off
legend('Infected','Symptomatic','Immune','Dead');
ylabel('Individuals')
xlabel('Days')

%-------------------------------------------------------------------------%
% Function for Calculating Infection Status
%-------------------------------------------------------------------------%
function [isinfected_new,infectionday_new]=infect(group,isinfected,infectionday,issick,isimmune,stayshome,currentday,N,p1,m1,ms,mr,pd,homegroups)

subject_group=group;
isinfected_new=isinfected;
infectionday_new=infectionday;

infect_vulnerable=(100*rand(N,1)<=100*p1 & infectionday==0); %infection vulnerable: unlucky dice roll + nonreal infection day
infect_capable=(isinfected==1 & currentday-infectionday_new>=m1) & (issick==0 | stayshome==0 | subject_group==homegroups); %infection capable: infected + past latency period + (either asymptomatic, doesnt stay home, or infecting housemate)

for i=1:max(subject_group)
    group_members=subject_group==i; 
    
    if max(infect_capable+group_members)==2 %identify members in same group
    infect_mask=(group_members+infect_vulnerable==2);
    isinfected_new(infect_mask)=1;
    infectionday_new(infect_mask)=currentday;
    end
end
end

function [isinfected_new,issick_new,isimmune_new,isdead_new]=disease(isinfected,infectionday,issick,isimmune,isdead,currentday,N,p1,m1,ms,mr,pd,ages)
isinfected_new=isinfected;
issick_new=issick;
isimmune_new=isimmune;
isdead_new=isdead;

sick_mask=(currentday-infectionday==ms & isinfected_new==1); %identify individuals becoming symptomatic
issick_new(sick_mask)=1;

recover_mask=(currentday-infectionday==mr & isinfected_new==1); %identify individuals at end of infection 
isimmune_new(recover_mask)=1;
issick_new(recover_mask)=0;
isinfected_new(recover_mask)=0;
    if .125*sum(issick)<.003*N %account for limited medical resources and 12.5% avg hospitalization rate - based on 53,000 hospital beds and 19.45 million citizens in new york state
    isdead_new=isdead_new+((rand(N,1).*(recover_mask))>=(pd./(ages./60)));
    else
    isdead_new=isdead_new+((rand(N,1).*(recover_mask))>=((2*pd/3)./(ages./60)));
    end
end

%-------------------------------------------------------------------------%
% Function for Sorting Population into Distinct Groups
%-------------------------------------------------------------------------%
function groupnumbers=generate_groups(groupmeansize,groupsdv,participation,population)

% The size of the groups is approximately Gaussian, with prescribed mean and standard deviation
%       groupmeansize = average size of the groups
%       groupsdv = standard deviation of group size
%       participation(i=1:population) Logical flag indicating whether ith member of the population should be placed in a group
%       Population = total # people in the population
%
%       groupnumbers(i)   The group occupied by the ithe member of the population (zero if member is not in a group)
%
nmembers = sum(participation); % Number of people to be placed in groups
ngroups = round(nmembers/groupmeansize); % Number of groups
groupsize = round(normrnd(groupmeansize,groupsdv,ngroups,1)); % group sizes (approx Gaussian distribution)
mask = groupsize>0;
groupsize = groupsize(mask); % Remove groups with zero or fewer members (distribution no longer perfectly Gaussian)
groupnumbers = zeros(population,1); % Initialize group number for whole population
nassigned = 0;
assignmentorder = 3*population*ones(population,1); % Initialize vector used to place participants in queue
assignmentorder(participation) = randi(2*population,nmembers,1); % Each member of the population is put in the queue in random order
[~,index] = sort(assignmentorder);   % Generate index table pointing to group members in order of priority
% Fill the groups
for i = max(groupsize):-1:1
   availablegrouplist = find(groupsize>=i);  % Add a member to all groups with size >= i
   navailablegroups = length(availablegrouplist); % No. groups to receive a new member
   if (navailablegroups<1) break; end % Break out of loop if no groups left
   groupnumbers(index(nassigned+1:min(nassigned+navailablegroups,nmembers))) = availablegrouplist(1:min(navailablegroups,nmembers-nassigned)); % Assign people to groups
   nassigned = nassigned + min(navailablegroups,nmembers-nassigned);
   if (nassigned >= nmembers) break; end % Break out of loop if everyone has been assigned to a group
end
end

function [adj] = update_adjacency(groupnumbers,adj)

n_groups = max(groupnumbers);
for i=1:n_groups
    indices = find(groupnumbers==i);
    adj(indices,indices) = true;
end
% Remove self interactions
[n,~] = size(adj);
adj(1:n+1:end) = false;

end