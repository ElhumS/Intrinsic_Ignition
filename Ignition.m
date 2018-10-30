% WRITTEN BY: GUSTAVO DECO AND EDITED BY: ELHUM A SHAMSHIRI (Elhum.Shamshiri@unige.ch) AND MANEL VILA-VIDAL
% CODE SHOULD NOT BE DISTRIBUTED AND ANALYSIS SHOULD NOT BE CONDUCTED WITHOUT PRIOR CONSENT FROM  THE AUTHORS
% If you would like to use this software for publication please contact

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

cd I:\Data\Patients\subject1\fMRI_Analysis\IIM
load parcellation_matrix_reduced.mat
xs= y; % number of seeds x timecourse


flp = 0.04;                        % lowpass frequency of filter
fhi = 0.07;                        % highpass
delt = 1.90;                      % sampling interval
k=2;                               % 2nd order butterworth filter
fnq=1/(2*delt);                    % Nyquist frequency: half of the sampling rate
Wn=[flp/fnq fhi/fnq];              % butterworth bandpass non-dimensional frequency (normalization)
[bfilt,afilt]=butter(k,Wn);      % construct the filter


N=size(xs,1);
T=size(xs,2);



%% Obtain events for each region

for seed = 1:N % for each region
  x = demean(detrend(xs(seed,:)));
  timeseriedata = filtfilt(bfilt,afilt,x);    % zero phase filter the data
  Xanalytic = hilbert(demean(timeseriedata));
  Phases(seed,:) = angle(Xanalytic); % obtain the phases for each time point
  tise = detrend(demean(timeseriedata)); % clean the filtered signal
  ev1 = tise>(std(tise))+mean(tise); %there is an event if the signal is higher than 2*SD
  ev2 = [0 ev1(1:end-1)];
  events(seed,:) = (ev1-ev2)>0; % extract the timeseries of events (1=event, 0=no event)
end

%% INTEGRATION

 for t=1:T
     for i=1:N 
         for k=1:N
             phasematrix(i,k)=exp(-3*adif(Phases(i,t),Phases(k,t)));
         end
     end
     cc=phasematrix;
     cc=cc-eye(N); %returns an n-by-n identity matrix with ones on the main diagonal and zeros elsewhere
     pp=1;
     PR=0:0.01:0.99;
     for p=PR
         A=abs(cc)>p;
         [comps csize]=get_components(A);
         cs(pp)=max(csize);
         pp=pp+1;
     end
     integ(t)=sum(cs)*0.01/N;
 end

 % size of the largest connected subcomponent averaged across all possible
 % binarisations (whether it passes threshold or not --> you do it across subgroups) 
 % of the functional connectivity matrix
 
 %% Save the variables integration and events (optional)
 
 save ('SS_Ignition_subject1_ptprocess_1sd.mat', 'integ','events');
 
 %% Plot the integration and mark the events
 
plot(integ)
hold on
for seed=1:N
    t_ev=find(events(seed,:));
    scatter(t_ev,integ(t_ev),'o')
end
title('There is an event if the signal is higher than: 1SD')
xlabel('time (scans)')
ylabel('Integration')
ylim([0 1])
xlim([0 size(integ,2)])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'ptprocess_events_1sd.jpeg')  

%% Obtain IDMI for each brain region

figure(2)
TRwindow=3; % number of TRs after the event, adjust if needed
for seed=1:N
    t_ev=find(events(seed,:));
    t_ev=t_ev(t_ev<T-TRwindow);
    aux=[];
    for t=t_ev
	    aux=[aux mean(integ(t:t+TRwindow))];
    end
    IDMI(seed)=mean(aux);
end

save('IDMI_1sd','IDMI')
plot(IDMI)

figure(3)
imagesc(phasematrix)
saveas(gcf,'phasematrix_1sd.jpeg')
