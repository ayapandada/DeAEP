function [de_EP]=de_EP_function_C(N,M,nTx,C,M_mod,H_c,Y_c,rou,Tmax,lenC)
%% Expectation Propagation algo
% Initialization
pc = ones(M*N*nTx,1);
p=1;
m2l = zeros(M*N*nTx,1);
v2l = diag(p*pc);
Th=1;
H=[];Y=[];
for c=1:1:length(H_c)
    H=[H;H_c{1,c}];
    Y=[Y;Y_c{1,c}];
end
mq=zeros(M*N*nTx,1);
Eq_record=[];
B_record=[];
B=norm(Y-H*mq,'fro');
Eq_record=[Eq_record;mq];

B_record=[B_record,B];


% Compute ratio distribution
for niter=1:1:Tmax

        % Ratio Distribution Computation
        va_ext_temp=0;
        xa_ext_temp=0;
        va_ext_l=[];
        xa_ext_l=[];
        for c=1:1:C
            if niter>1
             m2l=m2_l(:,c);
             v2l=v2_l{1,c};
            end
            [m2c,v2c]=RDC(Y_c{c},H_c{c},M*N*nTx,rou,pc,p,m2l,v2l); %pc-gamma m2-lambda  
            va_extl=v2c;
            xa_extl=m2c;

            va_ext_l=[va_ext_l,va_extl.'];
            xa_ext_l=[xa_ext_l,xa_extl.'];
             
            factorC=lenC(c)/sum(lenC);
            va_ext_temp=va_ext_temp+factorC*(1./va_extl)*4;
            xa_ext_temp= xa_ext_temp+factorC*(xa_extl./va_extl)*4;
           
        end
            va_ext=1./va_ext_temp;
            xa_ext=va_ext.*xa_ext_temp;
        % Moment Computation
    
            m2c=xa_ext;
            v2c=va_ext;


          for i=1:1:M*N*nTx
            m2ci=m2c(i);
            v2ci=v2c(i);

            G0 = (1-p)*cmplx_pdf(0,m2c(i),v2c(i)) + p*cmplx_pdf(0,m2c(i),pc(i) + v2c(i)); %first moment
            G1 = p*cmplx_pdf(0,m2c(i),pc(i) + v2c(i))*(m2c(i)*pc(i)/(pc(i) + v2c(i))); %second moment 
            G2 = p*cmplx_pdf(0,m2c(i),pc(i) + v2c(i))*((abs(m2c(i)*pc(i)/(pc(i) + v2c(i))))^2 + pc(i)*v2c(i)/(pc(i) + v2c(i))); 
            if G0==0
              G0=1;
            end
            Eq = G1/G0;
            Eq_all(i)=Eq;
            Vq = G2/G0 - abs(Eq)^2;
            Vq_all(i)=Vq;
          end
            % VB_mean=Vq_all;
            
            if niter==1
                Eq_all1=Eq_all;
            end


            % Threshold
             B=norm(Y-H*(Eq_all).','fro');
             B_temp=norm(Y-H*(Eq_all).','fro');
             if norm(B_temp-B,'fro')>Th
                 B=B_temp
                 B_record=[B_record,B];
                 Eq_record=[Eq_record,Eq_all.'];
             else
                 B=B_temp;
                 B_record=[B_record,B];
                 Eq_record=[Eq_record,Eq_all.'];
                 break;
             end

            % if B_temp>=B
            %     T_count=niter;
            %     break;
            % else
            %     T_count=niter;
            %     B=B_temp;
            % end

        % Moment Matching: Update parameters for next iteration
          v2_l=cell(1,C);
          m2_l=[];

            for c=1:1:C
             va_ecxl=va_ext_l(:,c);
             xa_extl=xa_ext_l(:,c);

             v2c=va_ecxl;
             m2c=xa_extl;
             for i=1:1:M*N*nTx
                % v2l(i,i) = (1/(mean(Vq_all)) - 1/v2c(i))^-1;
                v2l(i,i) = (1/(Vq_all(i)) - 1/v2c(i))^-1;
               %  m2l(i) = v2l(i,i)*(Eq_all(i)/(mean(Vq_all)) - m2c(i)/v2c(i));
                m2l(i) = v2l(i,i)*(Eq_all(i)/(Vq_all(i)) - m2c(i)/v2c(i));
             end
             
             v2_l{1,c}=v2l;
             m2_l=[m2_l,m2l]; 
            end

            % if G0==1
            %    de_EP=Eq_all.';
            %    return;
            % end

end

Eq_nan=isnan(Eq_all); %NAN Debug
findnan=find(Eq_nan==1);

if length(findnan)>0
  Eq_all=Eq_all1;  
else
    [value,index]=min(B_record);
    Eq_all=Eq_record(:,index).';
end

 de_EP=Eq_all.';
