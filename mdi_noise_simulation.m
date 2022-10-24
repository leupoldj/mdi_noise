% Parameters
GM_T1=1.6;  %T1 =1.6s for grey matter
TR=0.012;   %TR=12ms
angle1=4;   %angles in degree
angle2=16;
M0=1000;    %to go for high SNR regime 

GM_amp_a1=M0*ernst(GM_T1,TR,angle1); %computes the signal according to the Ernst equation (Eq. 1 in Ye et a. with E2*=exp(-i*phi)=1)
GM_amp_a2=M0*ernst(GM_T1,TR,angle2);

stan=1;             %standard dev. for the noise from one channel
no_bins=100;        %no of bins in the resulting histogram (i.e. the probability distribution function)
no_channels=32;
no_stat = 5000000;   %for the statistics (to repeat all no_stat times with independent noise instances). This value can be lowered to, e.g., 100000 to be faster but will result in uglier curves
TR_vec=[0.01:0.0001:0.019 0.02:0.01:0.2];   %vector for TR variation
a2_vec=[1:90];                              %vector for alpha2 variaiton

% initialize result vectors
eq4_vs_TR=zeros(1,length(TR_vec));
eq4_vs_a2=zeros(1,length(a2_vec));

bin_mat_TR=zeros(length(TR_vec),no_bins);
pdf_mat_TR=zeros(length(TR_vec),no_bins);   %pdf means "probability distribution function
pdf_mean_vs_TR=zeros(1,length(TR_vec));
pdf_std_vs_TR=zeros(1,length(TR_vec));

bin_mat_a2=zeros(length(a2_vec),no_bins);
pdf_mat_a2=zeros(length(a2_vec),no_bins);
pdf_mean_vs_a2=zeros(1,length(a2_vec));
pdf_std_vs_a2=zeros(1,length(a2_vec));

for k=1:length(TR_vec)      %Eq.4 in dependency of TR
 
    GM_amp_a1=M0*ernst(GM_T1,TR_vec(k),angle1); %caculate GM amplitude for angle1
    GM_amp_a2=M0*ernst(GM_T1,TR_vec(k),angle2); %caculate GM amplitude for angle2
    
    [pdf,bin_vec,pdf_std,pdf_mean]=distrib(GM_amp_a1,GM_amp_a2,stan,no_stat,no_bins,no_channels);
    bin_mat_TR(k,:)=bin_vec;     %bins for the histogram that displays the pdf
    pdf_mat_TR(k,:)=pdf;         %the pdf (probability distribution function)
    pdf_std_vs_TR(k)=pdf_std;    %standard deviation of the pdf   
    pdf_mean_vs_TR(k)=pdf_mean;  %mean of the pdf
    
    SD=eq4(GM_amp_a1,GM_amp_a2,stan,no_channels);   %Standard deviation according to Eq. 4
    eq4_vs_TR(k)=SD;

end

%% plot result of std vs TR
figure;
subplot(1,2,2);
plot(TR_vec,pdf_std_vs_TR,'LineWidth',2)
hold on; plot(TR_vec,eq4_vs_TR,'LineWidth',2)
lgd=legend('simulation',' Eq.4','Location','northwest');
lgd.FontSize = 8;
strx = ['TR [s]'];
xlabel(strx)
ylabel('MDI std')
strt = ['T1=',num2str(GM_T1),'s, \alpha_1=',num2str(angle1),'°, \alpha_2=',num2str(angle2),'°, \sigma=',num2str(stan)];
title(strt);
ax=gca;
ax.TitleFontSizeMultiplier = 0.8;
set(gca,'FontSize',18)
axis square;
grid on;

subplot(1,2,1);
plot(bin_mat_TR(21,:)-pdf_mean_vs_TR(21),pdf_mat_TR(21,:)/M0,'LineWidth',2)
hold on; plot(bin_mat_TR(100,:)-pdf_mean_vs_TR(100),pdf_mat_TR(100,:)/M0,'LineWidth',2)
lgd=legend('simulation TR=12ms','simulation TR=100ms','Location','northwest');
lgd.FontSize = 8;
strx = ['S_{aT1W}-<S_{aT1W}>'];
xlabel(strx)
ylabel('probability distribution')
strt = ['T1=',num2str(GM_T1),'s, \alpha_1=',num2str(angle1),'°, \alpha_2=',num2str(angle2),'°, \sigma=',num2str(stan)];
title(strt);
ax=gca;
ax.TitleFontSizeMultiplier = 0.8;
set(gca,'FontSize',18)
axis square;
set(gcf,'Position',[400 400 1200 500]);
grid on;


%%
for k=1:length(a2_vec)      %Eq.4 in dependency of alpha_2
    
    GM_amp_a1=M0*ernst(GM_T1,TR,angle1);
    GM_amp_a2=M0*ernst(GM_T1,TR,a2_vec(k));
    
    [pdf,bin_vec,pdf_std,pdf_mean]=distrib(GM_amp_a1,GM_amp_a2,stan,no_stat,no_bins,no_channels);
    bin_mat_a2(k,:)=bin_vec;       %bins for the histogram that displays the pdf
    pdf_mat_a2(k,:)=pdf;           %the pdf (probability distribution function)
    pdf_mean_vs_a2(k)=pdf_mean;    %mean of the pdf
    pdf_std_vs_a2(k)=pdf_std;      %standard deviation of the pdf  
    
    SD=eq4(GM_amp_a1,GM_amp_a2,stan,no_channels);  %Standard deviation according to Eq. 4
    eq4_vs_a2(k)=SD;
    
end

%% plot result of std vs alpha2
figure;
subplot(1,2,2);
plot(a2_vec,pdf_std_vs_a2,'LineWidth',2)
hold on;plot(a2_vec,eq4_vs_a2,'LineWidth',2)
lgd=legend('simulation',' Eq.4','Location','northwest');
lgd.FontSize = 8;
strx = ['\alpha_2 [deg]'];
xlabel(strx)
ylabel('MDI std')
strt = ['T1=',num2str(GM_T1),'s, \alpha_1=',num2str(angle1),'°, TR=',num2str(TR),'s, \sigma=',num2str(stan)];
title(strt);
ax=gca;
ax.TitleFontSizeMultiplier = 0.8;
set(gca,'FontSize',18)
axis square;
grid on;

subplot(1,2,1);
plot(bin_mat_a2(16,:)-pdf_mean_vs_a2(16),pdf_mat_a2(16,:)/M0,'LineWidth',2)
hold on; plot(bin_mat_a2(40,:)-pdf_mean_vs_a2(40),pdf_mat_a2(40,:)/M0,'LineWidth',2)
lgd=legend('simulation a2=16°','simulation a2=40°','Location','northwest');
lgd.FontSize = 8;
strx = ['S_{aT1W}-<S_{aT1W}>'];
xlabel(strx)
ylabel('probability distribution')
strt = ['T1=',num2str(GM_T1),'s, \alpha_1=',num2str(angle1),'°, TR=',num2str(TR*1000),'ms, \sigma=',num2str(stan)];
title(strt);
ax=gca;
ax.TitleFontSizeMultiplier = 0.8;
set(gca,'FontSize',18)
axis square;
set(gcf,'Position',[400 400 1200 500]);
grid on;


%% functions

%% compute the probability distribution function (pdf)
function  [f_pdf,f_bin_vec,f_pdf_std,f_pdf_mean] = distrib(f_x0, f_y0, f_std, f_no_stat, f_no_bins, f_channels);

    % Input parameters
    sigma_x = f_std; % std of noise in x   (x and y correspond to the signals S(a1) and S(a2)
    sigma_y = f_std; % std of noise in y

    % Calculations
    x_re = f_x0 + randn(f_channels,f_no_stat)*sigma_x;  % noisy data for x (real part)
    y_re = f_y0 + randn(f_channels,f_no_stat)*sigma_y;  % noisy data for y (real part)
    x_im = randn(f_channels,f_no_stat)*sigma_x;         % noisy data for x (imaginary part)
    y_im = randn(f_channels,f_no_stat)*sigma_y;         % noisy data for y (imaginary part)
    x=x_re+1i*x_im;
    y=y_re+1i*y_im;

    R = abs(sum(conj(x).*y,1))./sum(abs(conj(x).*x),1); %Eq.(3) in Ye et al. 

    %compute the histogram
    [n,b] = hist(R,f_no_bins); 
    dR = mean(diff(b)); 

    %output
    f_pdf_std=std(R);
    f_pdf_mean=mean(R);
    f_bin_vec=b;
    f_pdf=n/f_no_stat/dR; %for normalization

end

%% Ernst amplitude (with M0=Cn=E2*=exp(iphi)=1 ), Eq. 1 in Ye et al.
function f_ampl=ernst(f_T1,f_TR,f_alpha);

    f_ampl=sind(f_alpha)*(1-exp(-f_TR/f_T1))/(1-exp(-f_TR/f_T1)*cosd(f_alpha));

end

%% Eq.4 in Ye et al. 
function f_SD=eq4(f_Sa1,f_Sa2,f_sigma,f_channels);

    f_SD=sqrt(2)*f_sigma*((abs(1-f_Sa2/f_Sa1))*abs(f_Sa1)*f_channels) / (( f_channels*abs(f_Sa1*f_Sa1)) + 2*f_channels*f_sigma*abs(f_Sa1));

end
