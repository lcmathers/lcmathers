clear all;

%% Test of channel estimation


%% Defintion of the channel

fs=20e6;
Ts=1/fs;
max_delay=10e-6; % 最大时延
edge_time=1e-6;  % 边沿时间
N_delay=20;  % 时延抽头数目 

[delay_actual,delay_actual_dB]=delay_pdp1(N_delay,Ts,max_delay,edge_time);

delay_tap=floor(delay_actual./(Ts));

PowerdB=[0 -8 -17 ];
% PowerdB=delay_actual_dB;
% Delay=delay_tap;
Delay=[0 3 5];

Power=10.^(PowerdB/10);



channel_length=Delay(end)+1;

Ntap=length(Delay);

channel=sqrt(Power/2).*(randn(1,Ntap)+1j*randn(1,Ntap));
h=zeros(1,Ntap+1);
h(Delay+1)=channel;

H=fft(h,2048);
H_power_dB=10*log10(abs(H.*conj(H)));



%% OFDM

Nfft=2048;
Ng=Nfft/4;
Nvc=400;
Ndata=Nfft-Nvc;
Nframe=3;

M_mod=16;
M_bit=log2(M_mod);

sym_perframe=Ndata*Nframe/2;
bits_perframe=sym_perframe*M_bit;

Nps=2; % 导频间隔
Np=Ndata/2;

max_iter=1;
SNR_dB=30;

for i=1:length(SNR_dB)
    for iter=1:max_iter

        %% Tx

        Xp=2*(randn(1,Np)>0)-1;
        msgint=randi([0 1],(Ndata-Np)*M_bit,1);
        Data=qammod(msgint,M_mod,"gray",'InputType','bit',"UnitAveragePower",true);
        
        ip=0;
        pilot_loc=[1:2:2047];
        
        pilot_loc_outvc=pilot_loc;
        pilot_loc_outvc(pilot_loc_outvc==1)=[];
        pilot_loc_outvc((pilot_loc_outvc>=(Ndata/2+2)) & (pilot_loc_outvc<=(Ndata/2+Nvc)))=[];
        pilot_loc_vc=[1 pilot_loc(find(pilot_loc>=(Ndata/2+2) & (pilot_loc<=(Ndata/2+Nvc))))];

        loc_Data=[2:2:2048];
        loc_Data_outvc=loc_Data;
        loc_Data_outvc(loc_Data_outvc==1)=[];
        loc_Data_outvc(loc_Data_outvc>=(Ndata/2+2) & loc_Data_outvc<=(Ndata/2+Nvc))=[];


% 
%         for k=1:Ndata
%             if mod(k,Nps)==1
%                 X(k)=Xp(floor(k/Nps)+1);
%                 ip=ip+1;
%             else
%                 X(k)=Data(k-ip);
%             end
%         end
% 
% 
%         if Nvc==0
%             X_tx= [X(Ndata/2+1:end) X(1:Ndata/2)];
%         else
%             X_tx=[0 X(Ndata/2+1:end) zeros(1,Nvc-1) X(1:Ndata/2)];
%         end

        
        X_tx=zeros(1,Nfft);

        X_tx(pilot_loc_outvc)=[Xp(Ndata/4+1:end) Xp(1:Ndata/4)];
        X_tx(loc_Data_outvc)=[Data(Ndata/4+1:end) Data(1:Ndata/4)];

        x_tx=ifft(X_tx,Nfft);
        xt=[x_tx(Nfft-Ng+1:end) x_tx];


        %% Channel
        y=conv(xt,h);
        yt=awgn(y,SNR_dB(i),'measured');

        %% Rx

        y=yt(Ng+1:Nfft+Ng);
        Y=fft(y);


       if Nvc~=0
        H_LS_est_outvc= Y(pilot_loc_outvc)./[Xp(Ndata/4+1:end) Xp(1:Ndata/4)];
       end
       

       % 
       H_LS_est=zeros(1,Nfft);
       H_LS_est(pilot_loc_outvc)=H_LS_est_outvc;
       H_LS_est(pilot_loc_vc)=1e3;
       H_LS_est(H_LS_est==0)=[];
       
       index_vc=find(H_LS_est==1e3);
       H_LS_est(H_LS_est==1e3)=0;

       
       
       H_initial=H_LS_est;

       %% 时域处理
       for idk=1:20

       h_time=ifft(H_initial);

       h_delete=[h_time(1:channel_length) zeros(1,length(h_time)-channel_length)];


       H_LS_revice=fft(h_delete);

       %% 

       H_LS_est_modifed=H_initial;

       H_LS_est_modifed(index_vc)=H_LS_revice(index_vc);

       H_initial=H_LS_est_modifed;

       end

       
       % H_LS_est_modifed=H_initial;
       
       H_final=channel_interp(H_LS_est_modifed,pilot_loc,Nfft);

       % DFT 处理

       h_DFT=ifft(H_final);
       h_DFT_modified=[h_DFT(1:channel_length) zeros(1,length(h_DFT)-channel_length)];
       H_DFT_modified=fft(h_DFT_modified);


       H_final_dB= 10*log10(abs(H_DFT_modified.*conj(H_DFT_modified)));

       % H_final_dB= 10*log10(abs(H_final.*conj(H_final)));




    end
end

plot(H_final_dB);
hold on;
plot(H_power_dB);
legend('Est','Real');