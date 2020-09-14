clear all
clc


Modulasi = 2;
indeks_modulasi=log2(Modulasi);

frame= 1;

SNR=-6;

%LDPC
iterasi=50;

H=[1,0,1,0,1,0,0,0,0;1,0,0,1,1,1,0,0,0;1,0,1,0,0,1,1,0,0;0,1,0,1,0,0,1,1,0;0,1,1,0,0,0,0,1,1;0,1,0,1,0,0,0,0,1];

rows = size(H, 1);
cols = size(H, 2);

k = cols-rows;

G_awal = H;

r = 1;
for c = cols - rows + 1:cols
    if G_awal(r,c) == 0
        %         Swap needed
        for r2 = r + 1:rows
            if G_awal(r2,c) ~= 0
                tmp = G_awal(r, :);
                G_awal(r, :) = G_awal(r2, :);
                G_awal(r2, :) = tmp;
            end
        end
        
        %         Ups...
        if G_awal(r,c) == 0
            error('H is singular');
        end
    end
    
    %     Forward substitute
    for r2 = r + 1:rows
        if G_awal(r2, c) == 1
            G_awal(r2, :) = xor(G_awal(r2, :), G_awal(r, :));
        end
    end
    
    %     Back Substitution
    for r2 = 1:r - 1
        if G_awal(r2, c) == 1
            G_awal(r2, :) = xor(G_awal(r2, :), G_awal(r, :));
        end
    end
    
    %     Next row
    r = r + 1;
end

G=[eye(3) G_awal(:,1:3)'];







% for ii=1:length(SNR)
snr=10.^(SNR/10);
%     snr=10.^(SNR(ii)/10);
%
sigma = sqrt(1./(2*snr));
%
jum_error = 0;

for i=1:frame
    
    bit_informasi = [0 1 1];
    
    %encode LDPC
    %     c=xor(bit_informasi,G);
    
    c = mod(bit_informasi*G,2);
    
    x=1-c*2;
    
    
    
    h=1;
    
    noise=sigma*(randn(1,length(x)) + sqrt(-1)*randn(1,length(x)));
    
    y=h*x+noise;
    
    u=real(y);
    sebelum = u
    
    
    %soal
%     u=[-2.22184625490862,-0.906173357250190,0.165037758136963,2.11471620605673,-3.13753683350914,0.335291598543385,0.610903498910586,-3.66177320541728,0.351329372210964]
%     u=[5.6 -10.2 0.7 0.5 -7.5 12.2 -8.5 6.9 -7.7];
%       u=[1.19476685021654,-2.00275018457453,0.5962459254459502,1.87811287842705,-0.0866339309187617,-1.60050820969972,2.47940433564292,-0.0678324789579936,2.53953554523820]
      %soal
      u=[1.2, -2, 0.6, 1.9, -0.1, -1.6, 2.5, -0.1, -2.5];
    %       Soft Decoding
    
    Le_VND= zeros(rows,cols);
    La_VND= zeros(rows,cols);
    Le_CND= zeros(rows,cols);
    La_CND= zeros(rows,cols);
    
    %masuk LLR chanel
    %iterasi ke-0
    for k=1:cols
        el1=find(H(:,k));
        La_VND(el1,k)=u(k);
    end
    
    Le_VND=La_VND;
    
    for iter = 1 : iterasi
        
        %check node
        for m=1:rows
            el2=find(H(m,:));
            for ec=1:length(el2)
                el2_baru=setdiff(el2,el2(ec));                
                for ee=1:length(el2_baru)
                    if ee==1
                        operasi_CND=Le_VND(m,el2_baru(ee));
                    else
                        operasi_CND=sign(operasi_CND)*sign(Le_VND(m,el2_baru(ee)))*min(abs([operasi_CND Le_VND(m,el2_baru(ee))]));
                    end 
                end
                La_CND(m,el2(ec))=operasi_CND;       
            end
        end
        
        Le_CND=La_CND;
        
        
        %variable node
        for k=1:cols
            eel1=find(H(:,k));
            for ev=1:length(eel1)
                La_VND(eel1(ev),k)=u(k)+sum(Le_CND(eel1,k))-Le_CND(eel1(ev),k);
            end
        end
        
        Le_VND=La_VND;
        
        %mencari nilai per baris(representasi per check nodes)
        %         for ui = 1:rows
        %             mencari nilai 1 pada baris matriks tertentu
        %             EdgeI = find(H(ui,:));
        %             sumLrji =  Lrji;
        %             menghilangkan baris feedback sendiri
        %             sumLrji(ui,:) = [];
        %             for iterI = 1:numel(EdgeI)
        %                 ei = EdgeI(iterI);
        %                 perhitungan per nodes pada baris tertentu
        %                 EdgeJ = EdgeI;
        %                 menghilangkan nodes sendiri agar tidak masuk dalam perhitungan boxplus
        %                 EdgeJ(iterI) = [];
        %                 update nilai feedback
        %                 Lqij = u(ui) + sum(sumLrji);
        %                 state posisi dan nilai eh pada perhitungan boxplus
        %                 Lrji_h(ei) = Lqij(EdgeJ(1));
        %                 for iterJ = 1:numel(EdgeJ)-1
        %                     Lrji_h(ei);
        %                     Lqij(EdgeJ(iterJ+1));
        %                     Lrji_h(ei) = (-1)*sign(Lrji_h(ei))*sign(Lqij(EdgeJ(iterJ+1)))*min(abs([Lrji_h(ei) Lqij(EdgeJ(iterJ+1))]));
        %                 end
        %             end
        %             update boxplus untuk dilanjutkan pada baris selanjutnya
        %             Lrji(ui,:) = Lrji_h;
        %         end
        %update nilai LLR berdasarkan iterasi
    end
    
    for in=1:length(bit_informasi)
        ub(in)=sum(Le_CND(:,in))+u(in);
    end
    
    
    u=ub;
    
    for j=1:length(u)
        if u(j)>0
            b(j)=0;
        else
            b(j)=1;
        end
    end
    
    error = sum(b~=bit_informasi); %soft Demapper
    
    jum_error=jum_error+error;
    
end

ber=jum_error/(length(bit_informasi)*frame);


% ber=jum_error/(info*frame);

%     berplot(ii)=ber;
%
% end

% semilogy(SNR,berplot,'-b*','linewidth',1)

