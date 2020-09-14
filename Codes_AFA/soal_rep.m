clear all
clc


Modulasi = 2;
indeks_modulasi=log2(Modulasi);

frame= 1;

SNR=-5;

rate = 1/3;


% for ii=1:length(SNR)
snr=10.^(SNR/10);
%     snr=10.^(SNR(ii)/10);
%
sigma = sqrt(1./(2*snr));
%
jum_error = 0;

for i=1:frame
    
    bit_informasi = [0 1 1];
    
    %rep
    c=[];
    for r=1:length(bit_informasi)
        for re=1:1/(rate)
            c=[c bit_informasi(r)];
        end
    end
    
    x=1-c*2;
    
    
    
    h=1;
    
    noise=sigma*(randn(1,length(x)) + sqrt(-1)*randn(1,length(x)));
    
    y=h*x+noise;
    
    u=real(y);
    
    
    %soal
    u=[-0.911 3.153 1.646 -1.594 1.571 -1.394 -1.571 -1.271 -1.726]
    
    
    %       Soft Decoding
    ub=[];
    for xx=1:1/(rate):length(u)
        x=0;
        for j=1:1/(rate)
            x= x + u(xx+j-1);
        end
        ub=[ub x];
    end
    
    u=ub;
    
    for j=1:length(u)
        if u(j)>0
            b(j)=0;
        else
            b(j)=1;
        end
    end
    
    
%     %       Hard Decoding
    eb=[];
    m=1;
    for d=1:1/(rate):length(b)
        eb=[eb mode(b(d:1/(rate)*m))];
        m=m+1;
    end
    b=eb;
%soft Demapper

    
    error = sum(b~=bit_informasi);     
    jum_error=jum_error+error;
    
end

ber=jum_error/(info*frame);

%     berplot(ii)=ber;
%
% end

% semilogy(SNR,berplot,'-b*','linewidth',1)

