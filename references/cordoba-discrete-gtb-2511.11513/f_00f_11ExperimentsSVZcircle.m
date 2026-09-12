
clear;
%%%%%%%%%%%%%%%%%%%
% numerical setup %
%%%%%%%%%%%%%%%%%%%
SVZ=0;
MeV=1/140;
GeV=1000*MeV;
%s0=(1.2*GeV)^2;
s0=(2*GeV)^2;
power=2; %2
normaval=2;
Mline=200;
Mcircle=Mline;

numpoints=Mcircle-1 ;
numpointsfin=0;
P = load('DATA/pointsnear.mat').Expression1;
disp(P);

for iM=3:3

    M=(2*iM+1)^2;
    for inmax=46:50
        nmax=inmax;
        %nmax=3*inmax-1;
        for lmax=11:11
            multiplier=20;
            for multnmaxi=0:0

                multnmax=M-multnmaxi*multiplier;

                %lmax='uniform<>15';
                vnu = importdata(join(['Experiments/vnuvalsM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
                %vnu = importdata(join(['Experiments/vnuvalsM=' num2str(M) 'nmax=' num2str(1) 'lmax=uniform<>15.mat'],""));
                %amplitude part
                H = importdata(join(['Experiments/hmatyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
                %rescaled amplitudes
                reH = importdata(join(['Experiments/rescaleReHM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
                imH = importdata(join(['Experiments/rescaleImHM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
                imH2 = importdata(join(['Experiments/rescaleImH2M=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));

                %multnmax="uniform81"

                %import veczz and vecsin
                veczz=importdata(join(['Experiments/veczzM=' num2str(M) 'nmaxtilde=' num2str(multnmax) '.mat'],""));
                vecsin=importdata(join(['Experiments/vecsinM=' num2str(M) 'nmaxtilde=' num2str(multnmax) '.mat'],""));








                %functional
                F=importdata( join(['Experiments/FvecyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));

                %functional unphysical
                Fun=importdata(join(['Experiments/FunvecyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));
                %Funcomplement=importdata(join(['Experiments/FuncomplementvecyifeiM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=15.mat'],""));

                %define factor for cursive F_1^1
                FF3factor=importdata(join(["Experiments/GTB/FF3factorM=" num2str(M) '.mat'],""));

                QCDs=importdata(join(["Experiments/GTB/SVZSRat2M=" num2str(M) '.mat'],""));
                %int coefficients at 2 GeV

                intrho01=importdata( join(["Experiments/GTB/intS01at2M=" num2str(M) '.mat'],""));
                intrho00=importdata( join(["Experiments/GTB/intS00at2M=" num2str(M) '.mat'],""));
                intrho02=importdata( join(["Experiments/GTB/intS02at2M=" num2str(M) '.mat'],""));

                intrho1m1=importdata( join(["Experiments/GTB/intP1m1at2M=" num2str(M) '.mat'],""));
                intrho10=importdata(  join(["Experiments/GTB/intP10at2M=" num2str(M) '.mat'],""));
                intrho11=importdata(  join(["Experiments/GTB/intP11at2M=" num2str(M) '.mat'],""));

                intrho2m2=importdata( join(["Experiments/GTB/intD0m2at2M=" num2str(M) '.mat'],""));
                intrho2m1=importdata(  join(["Experiments/GTB/intD0m1at2M=" num2str(M) '.mat'],""));
                intrho20=importdata(  join(["Experiments/GTB/intD00at2M=" num2str(M) '.mat'],""));








                dims=size(vnu);

                % multnmax Has to be less than M


                Mvnu=dims(1);
                Mreg= 10;


               MtolCSB=2*10^(-3);

                %1*10^-7
                %new 0.08*10^(-7);
                MSVZ0=1*10^(-7);
                %3*10^-8
                %better if we go higher
                MFF0=3*10^(-8);


                %5*10^-6
                %4*10^(-6);
                MSVZ=4*10^(-6);
                %2*10^-6
                %0.7*10^(-6);
                MFF=2*10^(-6);


                %6*10^-6
                MSVZ2=6*10^(-6);
                %4*10^-2
                MFF2=4*10^(-2);


                disp(["M" "nmax" "lmax" "MtolCSB" "MSVZ1" "MFF1" "MSVZ0" "MFF0" "MSVZ2" "MFF2" "nmax FF" "Mreg" ; M nmax lmax MtolCSB MSVZ MFF MSVZ0 MFF0 MSVZ2 MFF2  multnmax Mreg*10^power]);
                indexes=find(vnu>s0);
                is0=indexes(1)-1;
                disp(is0);
                %is0=indexes(1);


                fstore=MtolCSB;






                dims=size(H);
                i = dims(1); l = dims(2);m=dims(3); v=dims(4);


                F00=1:v;
                F00(:)=F(1,1,:);
                F11=1:v;
                F11(:)=F(2,1,:);
                Funcomplement=Fun;

                dims=size(Fun);
                MunphySB=dims(1);
                vsunphys=(1:MunphySB)/2;




                %define ratios CSB
                rS0P1chi=3*(2*vsunphys-1)./(vsunphys-4);
                rS2P1chi= 3*(2-vsunphys)./(vsunphys-4);


                %define functional unphysical region coefficients
                Fun00=zeros(MunphySB,v);
                Fun00(:,:)=Fun(:,1,:);

                Fun11=zeros(MunphySB,v);
                Fun11(:,:)=Fun(:,3,:);

                Fun02=zeros(MunphySB,v);
                Fun02(:,:)=Fun(:,2,:);

                Fun20=zeros(MunphySB,v);
                Fun20(:,:)=Funcomplement(:,1,:);

                Fun22=zeros(MunphySB,v);
                Fun22(:,:)=Funcomplement(:,2,:);

                Fun13=zeros(MunphySB,v);
                Fun13(:,:)=Funcomplement(:,3,:);


                %sqrt(4*pi/3)*(1/(8*pi^3))*


                cons1=FF3factor(2,:).';
                cons0=FF3factor(1,:).';
                cons2=FF3factor(3,:).';







                %define QCD sum rules values
                %equation 2.17g
                Nf=2;
                %QCD10=(1/(2*pi)^(4))*(1+0.4/pi)/(4*pi*(0+2));
                %QCD1m1=(1/(2*pi)^(4))*(1+0.4/pi)/(4*pi*(-1+2));
                %QCD11=(1/(2*pi)^(4))*(1+0.4/pi)/(4*pi*(1+2));
                %QCDs=importdata("Experiments/GTB/SVZSR.mat");


                QCD00=QCDs(1);
                QCD01=QCDs(2);
                QCD02=QCDs(3);
                QCD1m1=QCDs(4);
                QCD10=QCDs(5);
                QCD11=QCDs(6);
                QCD2m2=QCDs(7);
                QCD2m1=QCDs(8);
                QCD20=QCDs(9);



                %intrho1m1=importdata('Experiments/GTB/intP1m1.mat');
                %intrho10=importdata('Experiments/GTB/intP10.mat');
                %intrho11=importdata('Experiments/GTB/intP11.mat');
                %intrho01=importdata('Experiments/GTB/intS01.mat');
                %intrho00=importdata('Experiments/GTB/intS00.mat');
                %intrho02=importdata('Experiments/GTB/intS02.mat');




                % define vectors for storing maximization results
                sAA = zeros(Mline,5);
                % define vectors for storing f0 rho sigma
                sfrs=zeros(v,Mline);
                %define vector to compare with minimum distance to P

                normaup=10^10;
                normadown=10^10;
                sdup=zeros(1,2);
                sddown=zeros(1,2);
                vecprovup=zeros(1,2);
                vecprovdown=zeros(1,2);
                sdsvfup=zeros(v,1);
                sdsvfdown=zeros(v,1);

                sFormfactorup0=zeros(m,1);
                sFormfactorup1=zeros(m,1);
                sFormfactorup2=zeros(m,1);

                sspectralup0=zeros(m,1);
                sspectralup1=zeros(m,1);
                sspectralup2=zeros(m,1);

                sFormfactordown0=zeros(m,1);
                sFormfactordown1=zeros(m,1);
                sFormfactordown2=zeros(m,1);

                sFF0=zeros(2,multnmax);
                sFF1=zeros(2,multnmax);
                sFF2=zeros(2,multnmax);
                sspec0=zeros(2,multnmax);
                sspec1=zeros(2,multnmax);
                sspec2=zeros(2,multnmax);

                sspectraldown0=zeros(m,1);
                sspectraldown1=zeros(m,1);
                sspectraldown2=zeros(m,1);

                smatrixdown=zeros(6,m);
                smatrixup=zeros(6,m);

                phaseup1=zeros(m,1);
                phasedown1=zeros(m,1);

                phaseup2=zeros(m,1);
                phasedown2=zeros(m,1);

                phaseup3=zeros(m,1);
                phasedown3=zeros(m,1);

                phaseup4=zeros(m,1);
                phasedown4=zeros(m,1);

                phaseup5=zeros(m,1);
                phasedown5=zeros(m,1);

                phaseup6=zeros(m,1);
                phasedown6=zeros(m,1);


                elasticityup1=zeros(m,1);
                elasticitydown1=zeros(m,1);

                elasticityup2=zeros(m,1);
                elasticitydown2=zeros(m,1);

                elasticityup3=zeros(m,1);
                elasticitydown3=zeros(m,1);

                elasticityup4=zeros(m,1);
                elasticitydown4=zeros(m,1);

                elasticityup5=zeros(m,1);
                elasticitydown5=zeros(m,1);

                elasticityup6=zeros(m,1);
                elasticitydown6=zeros(m,1);

                tic
                cvx_clear

                cvx_begin quiet

                cvx_solver mosek

                variable f0
                variable rho1(nmax,nmax)
                variable rho2(nmax,nmax) symmetric
                variable sigma(2*nmax)
                variable ImFl1(multnmax)
                variable rhol1(multnmax)
                variable ImFl0(multnmax)
                variable rhol0(multnmax)
                variable ImFl2(multnmax)
                variable rhol2(multnmax)



                x1=reshape(rho1,nmax*nmax,1);
                x2=reshape(rho2,nmax*nmax,1);
                x=transpose(cat(2,f0 ,transpose(x1),transpose(x2),transpose(sigma)));

                % functional f2
                f2S0=F00*x;
                f2P1=F11*x;
                fS0unphys=Fun00*x;
                fP1unphys=Fun11*x;
                fS2unphys=Fun02*x;

                fD0unphys=Fun20*x;
                fF1unphys=Fun13*x;
                fD2unphys=Fun22*x;



                %SVZ coefficients l=1 I=1
                intrho10s=intrho10.'*(vecsin(1:is0,:)*rhol1);
                intrho1m1s=intrho1m1.'*(vecsin(1:is0,:)*rhol1);
                intrho11s=intrho11.'*(vecsin(1:is0,:)*rhol1);


                %SVZ coefficients l=0 I=0
                intrho00s=intrho00.'*(vecsin(1:is0,:)*rhol0);
                intrho01s=intrho01.'*(vecsin(1:is0,:)*rhol0);
                intrho02s=intrho02.'*(vecsin(1:is0,:)*rhol0);

                %SVZ coefficients l=0 I=2
                intrho2m2s=intrho2m2.'*(vecsin(1:is0,:)*rhol2);
                intrho2m1s=intrho2m1.'*(vecsin(1:is0,:)*rhol2);
                intrho20s=intrho20.'*(vecsin(1:is0,:)*rhol2);




                maximize(f2S0);
                subject to



                % regularization
                norm(x1,4) <= Mreg*10^(power);
                norm(x2,4) <= Mreg*10^(power);



                % quadratic cone
                tensor_mult(imH(: ,1:lmax,:,:),x,[4],[1]).^2+tensor_mult(reH(: ,1:lmax,:,:),x,[4],[1]).^2 <= 2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1])

                %chiSB

                norm([fS0unphys-rS0P1chi'.*fP1unphys; fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

                %optional
                %norm([fP1unphys])<=MtolCSB;
                %norm([fP1unphys,fD0unphys,fD2unphys])<=MtolCSB;

                %Positive definite matrix
                nh2=zeros(m,v);
                nh3=zeros(m,v);
                nh4=zeros(m,v);
                nh2(:,:)=1i*H(2,1,:,:);
                nh3(:,:)=1i*H(1,1,:,:);
                nh4(:,:)=1i*H(1,2,:,:);
                for iterpo=1:is0
                    FF1=1+veczz(iterpo,:)*ImFl1;
                    spec1=vecsin(iterpo,:)*rhol1;

                    FF0=1+veczz(iterpo,:)*ImFl0;
                    spec0=vecsin(iterpo,:)*rhol0;

                    FF2=1+veczz(iterpo,:)*ImFl2;
                    spec2=vecsin(iterpo,:)*rhol2;


                    %positivesemidefinite l=1, I=1
                    [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1;1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1);cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1] == hermitian_semidefinite( 3 );
                    %positivesemidefinite l=0, I=0
                    [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0;1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0);cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0]  == hermitian_semidefinite( 3 );
                    %positivesemidefinite l=0, I=2
                    [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2;1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2);cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2]   == hermitian_semidefinite( 3 );


                end



                norm([intrho1m1s-QCD1m1; intrho10s-QCD10; intrho11s-QCD11].' )<= MSVZ;
                norm([intrho01s-QCD01; intrho00s-QCD00; intrho02s-QCD02].')<= MSVZ0;
                norm([intrho2m2s-QCD2m2; intrho2m1s-QCD2m1; intrho20s-QCD20].')<= MSVZ2;

                %High energy limit
                real((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2+imag((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2<= MFF0;
                real((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2+imag((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2<= MFF;
                real((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2+imag((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2<= MFF2;


                cvx_end
                toc


                disp(["intrho1m1" abs(intrho1m1s) abs(QCD1m1);"intrho10"  abs(intrho10s) abs(QCD10);"intrho11"  abs(intrho11s) abs(QCD11)]);
                disp(["intrho01" abs(intrho01s) abs(QCD01);"intrho00"  abs(intrho00s) abs(QCD00);"intrho02"  abs(intrho02s) abs(QCD02)]);
                disp(["intrho2m2" abs(intrho2m2s) abs(QCD2m2);"intrho2m1"  abs(intrho2m1s) abs(QCD2m1);"intrho20"  abs(intrho20s) abs(QCD20)]);

                res1=rS0P1chi'.*fP1unphys;
                res2=rS2P1chi'.*fP1unphys;

                disp(["CSB S0" abs(fS0unphys(1)) abs(res1(1));"CSB S2"  abs(fS2unphys(1)) abs(res2(1))]);

                tic

                f2S0max=f2S0;
                disp(f2S0max);
                beep on; beep;
                % for f2S0bound = flip(((Mline-numpoints):(Mline-numpointsfin))*f2S0max/(Mline))
                for f2S0bound = (1:1)*P(1)
                    %ni=int32(f2S0bound*Mline/(f2S0max));
                    f2S0bound
                    ni=1;
                    sAA(ni,1) = f2S0bound;
                    lower=0.0725;
                    %lower=0;
                    upper=0.074;
                    %upper=100;

                    if f2S0max >= lower
                        disp("More");

                        cvx_clear
                        cvx_begin  quiet

                        cvx_solver mosek

                        variable f0
                        variable rho1(nmax,nmax)
                        variable rho2(nmax,nmax) symmetric
                        variable sigma(2*nmax)
                        variable ImFl1(multnmax)
                        variable rhol1(multnmax)
                        variable ImFl0(multnmax)
                        variable rhol0(multnmax)
                        variable ImFl2(multnmax)
                        variable rhol2(multnmax)


                        x1=reshape(rho1,nmax*nmax,1);
                        x2=reshape(rho2,nmax*nmax,1);
                        x=transpose(cat(2,f0 ,transpose(x1),transpose(x2),transpose(sigma)));

                        % functional f2
                        f2S0=F00*x;
                        f2P1=F11*x;
                        fS0unphys=Fun00*x;
                        fP1unphys=Fun11*x;
                        fS2unphys=Fun02*x;

                        fD0unphys=Fun20*x;
                        fF1unphys=Fun13*x;
                        fD2unphys=Fun22*x;





                        %SVZ coefficients l=1 I=1
                        intrho10s=intrho10.'*(vecsin(1:is0,:)*rhol1);
                        intrho1m1s=intrho1m1.'*(vecsin(1:is0,:)*rhol1);
                        intrho11s=intrho11.'*(vecsin(1:is0,:)*rhol1);


                        %SVZ coefficients l=0 I=0
                        intrho00s=intrho00.'*(vecsin(1:is0,:)*rhol0);
                        intrho01s=intrho01.'*(vecsin(1:is0,:)*rhol0);
                        intrho02s=intrho02.'*(vecsin(1:is0,:)*rhol0);

                        %SVZ coefficients l=0 I=2
                        intrho2m2s=intrho2m2.'*(vecsin(1:is0,:)*rhol2);
                        intrho2m1s=intrho2m1.'*(vecsin(1:is0,:)*rhol2);
                        intrho20s=intrho20.'*(vecsin(1:is0,:)*rhol2);


                        minimize(f2P1);
                        subject to

                        f2S0==f2S0bound;


                        % regularization
                        norm(x1,4) <= Mreg*10^(power);
                        norm(x2,4) <= Mreg*10^(power);


                        % quadratic cone
                        tensor_mult(imH(: ,1:lmax,:,:),x,[4],[1]).^2+tensor_mult(reH(: ,1:lmax,:,:),x,[4],[1]).^2 <= 2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1])

                        %chiSB

                        norm([fS0unphys-rS0P1chi'.*fP1unphys; fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

                        %optional
                        %norm([fP1unphys])<=MtolCSB;
                        %norm([fP1unphys,fD0unphys,fD2unphys])<=MtolCSB;

                        %Positive definite matrix
                        nh2=zeros(m,v);
                        nh3=zeros(m,v);
                        nh4=zeros(m,v);
                        nh2(:,:)=1i*H(2,1,:,:);
                        nh3(:,:)=1i*H(1,1,:,:);
                        nh4(:,:)=1i*H(1,2,:,:);
                        for iterpo=1:is0
                            FF1=1+veczz(iterpo,:)*ImFl1;
                            spec1=vecsin(iterpo,:)*rhol1;

                            FF0=1+veczz(iterpo,:)*ImFl0;
                            spec0=vecsin(iterpo,:)*rhol0;

                            FF2=1+veczz(iterpo,:)*ImFl2;
                            spec2=vecsin(iterpo,:)*rhol2;


                            %positivesemidefinite l=1, I=1
                            [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1;1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1);cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1]  == hermitian_semidefinite( 3 );
                            %positivesemidefinite l=0, I=0
                            [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0;1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0);cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0]  == hermitian_semidefinite( 3 );
                            %positivesemidefinite l=0, I=2
                            [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2;1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2);cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2]   == hermitian_semidefinite( 3 );


                        end



                        norm([intrho1m1s-QCD1m1; intrho10s-QCD10; intrho11s-QCD11].' )<= MSVZ;
                        norm([intrho01s-QCD01; intrho00s-QCD00; intrho02s-QCD02].')<= MSVZ0;
                        norm([intrho2m2s-QCD2m2; intrho2m1s-QCD2m1; intrho20s-QCD20].')<= MSVZ2;

                        %High energy limit
                        real((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2+imag((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2<= MFF0;
                        real((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2+imag((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2<= MFF;
                        real((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2+imag((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2<= MFF2;


                        cvx_end





                        vecprovdown(2)=f2P1;
                        vecprovdown(1)=f2S0;
                        if norm(vecprovdown-P)<normadown
                            sddown(:)=vecprovdown(:);
                            sdsvfdown(:)=x;
                            normadown=norm(sddown-P);
                            for iterform=1:m
                                FF1=1+veczz(iterform,:)*ImFl1;
                                spec1=vecsin(iterform,:)*rhol1;

                                sFormfactordown1(iterform)=cons1(iterform)*FF1;
                                sspectraldown1(iterform)=spec1;

                                FF0=1+veczz(iterform,:)*ImFl0;
                                spec0=vecsin(iterform,:)*rhol0;

                                sFormfactordown0(iterform)=cons0(iterform)*FF0;
                                sspectraldown0(iterform)=spec0;


                                FF2=1+veczz(iterform,:)*ImFl2;
                                spec2=vecsin(iterform,:)*rhol2;

                                sFormfactordown2(iterform)=cons2(iterform)*FF2;
                                sspectraldown2(iterform)=spec2;
                            end
                            sFF0(1,:)=ImFl0;
                            sFF1(1,:)=ImFl1;
                            sFF2(1,:)=ImFl2;
                            sspec0(1,:)=rhol0;
                            sspec1(1,:)=rhol1;
                            sspec2(1,:)=rhol2;

                            nh=zeros(m,v);
                            nh(:,:)=H(1,1,:,:);

                            smatrixdown(1,:)=1+1i*nh*x;


                            nh=zeros(m,v);
                            nh(:,:)=H(2,1,:,:);

                            smatrixdown(2,:)=1+1i*nh*x;


                            nh=zeros(m,v);
                            nh(:,:)=H(1,2,:,:);

                            smatrixdown(3,:)=1+1i*nh*x;


                            nh=zeros(m,v);
                            nh(:,:)=H(3,1,:,:);

                            smatrixdown(4,:)=1+1i*nh*x;


                            nh=zeros(m,v);
                            nh(:,:)=H(2,2,:,:);

                            smatrixdown(5,:)=1+1i*nh*x;

                            nh=zeros(m,v);
                            nh(:,:)=H(3,2,:,:);

                            smatrixdown(6,:)=1+1i*nh*x;


                        end

                         sdsvfdown(:)=x;

                    for iterform=1:m
                        FF1=1+veczz(iterform,:)*ImFl1;
                        spec1=vecsin(iterform,:)*rhol1;

                        sFormfactordown1(iterform)=cons1(iterform)*FF1;
                        sspectraldown1(iterform)=spec1;

                        FF0=1+veczz(iterform,:)*ImFl0;
                        spec0=vecsin(iterform,:)*rhol0;

                        sFormfactordown0(iterform)=cons0(iterform)*FF0;
                        sspectraldown0(iterform)=spec0;


                        FF2=1+veczz(iterform,:)*ImFl2;
                        spec2=vecsin(iterform,:)*rhol2;

                        sFormfactordown2(iterform)=cons2(iterform)*FF2;
                        sspectraldown2(iterform)=spec2;
                    end
                    sFF0(1,:)=ImFl0;
                    sFF1(1,:)=ImFl1;
                    sFF2(1,:)=ImFl2;
                    sspec0(1,:)=rhol0;
                    sspec1(1,:)=rhol1;
                    sspec2(1,:)=rhol2;

                    nh=zeros(m,v);
                    nh(:,:)=H(1,1,:,:);

                    smatrixdown(1,:)=1+1i*nh*x;
                    phasedown1=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown1=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,1,:,:);

                    smatrixdown(2,:)=1+1i*nh*x;
                    phasedown2=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown2=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(1,2,:,:);

                    smatrixdown(3,:)=1+1i*nh*x;
                    phasedown3=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown3=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,1,:,:);

                    smatrixdown(4,:)=1+1i*nh*x;
                    phasedown4=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown4=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,2,:,:);

                    smatrixdown(5,:)=1+1i*nh*x;
                    phasedown5=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown5=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,2,:,:);

                    smatrixdown(6,:)=1+1i*nh*x;
                    phasedown6=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown6=abs(1+1i*nh*x);



                        sAA(ni,2) = f2P1;

                        cvx_clear
                        cvx_begin  quiet

                        cvx_solver mosek

                        variable f0
                        variable rho1(nmax,nmax)
                        variable rho2(nmax,nmax) symmetric
                        variable sigma(2*nmax)
                        variable ImFl1(multnmax)
                        variable rhol1(multnmax)
                        variable ImFl0(multnmax)
                        variable rhol0(multnmax)
                        variable ImFl2(multnmax)
                        variable rhol2(multnmax)



                        x1=reshape(rho1,nmax*nmax,1);
                        x2=reshape(rho2,nmax*nmax,1);
                        x=transpose(cat(2,f0 ,transpose(x1),transpose(x2),transpose(sigma)));

                        % functional f2
                        f2S0=F00*x;
                        f2P1=F11*x;
                        fS0unphys=Fun00*x;
                        fP1unphys=Fun11*x;
                        fS2unphys=Fun02*x;

                        fD0unphys=Fun20*x;
                        fF1unphys=Fun13*x;
                        fD2unphys=Fun22*x;



                        %SVZ coefficients l=1 I=1
                        intrho10s=intrho10.'*(vecsin(1:is0,:)*rhol1);
                        intrho1m1s=intrho1m1.'*(vecsin(1:is0,:)*rhol1);
                        intrho11s=intrho11.'*(vecsin(1:is0,:)*rhol1);


                        %SVZ coefficients l=0 I=0
                        intrho00s=intrho00.'*(vecsin(1:is0,:)*rhol0);
                        intrho01s=intrho01.'*(vecsin(1:is0,:)*rhol0);
                        intrho02s=intrho02.'*(vecsin(1:is0,:)*rhol0);

                        %SVZ coefficients l=0 I=2
                        intrho2m2s=intrho2m2.'*(vecsin(1:is0,:)*rhol2);
                        intrho2m1s=intrho2m1.'*(vecsin(1:is0,:)*rhol2);
                        intrho20s=intrho20.'*(vecsin(1:is0,:)*rhol2);



                        maximize(f2P1);
                        subject to

                        f2S0==f2S0bound;



                        % regularization
                        norm(x1,4) <= Mreg*10^(power);
                        norm(x2,4) <= Mreg*10^(power);


                        % quadratic cone
                        tensor_mult(imH(: ,1:lmax,:,:),x,[4],[1]).^2+tensor_mult(reH(: ,1:lmax,:,:),x,[4],[1]).^2 <= 2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1])

                        %chiSB

                        norm([fS0unphys-rS0P1chi'.*fP1unphys; fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

                        %optional
                        %norm([fP1unphys])<=MtolCSB;
                        %norm([fP1unphys,fD0unphys,fD2unphys])<=MtolCSB;

                        %Positive definite matrix
                        nh2=zeros(m,v);
                        nh3=zeros(m,v);
                        nh4=zeros(m,v);
                        nh2(:,:)=1i*H(2,1,:,:);
                        nh3(:,:)=1i*H(1,1,:,:);
                        nh4(:,:)=1i*H(1,2,:,:);
                        for iterpo=1:is0
                            FF1=1+veczz(iterpo,:)*ImFl1;
                            spec1=vecsin(iterpo,:)*rhol1;

                            FF0=1+veczz(iterpo,:)*ImFl0;
                            spec0=vecsin(iterpo,:)*rhol0;

                            FF2=1+veczz(iterpo,:)*ImFl2;
                            spec2=vecsin(iterpo,:)*rhol2;


                            %positivesemidefinite l=1, I=1
                            [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1;1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1);cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1]  == hermitian_semidefinite( 3 );
                            %positivesemidefinite l=0, I=0
                            [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0;1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0);cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0]  == hermitian_semidefinite( 3 );
                            %positivesemidefinite l=0, I=2
                            [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2;1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2);cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2]   == hermitian_semidefinite( 3 );


                        end
                        norm([intrho1m1s-QCD1m1; intrho10s-QCD10; intrho11s-QCD11].' )<= MSVZ;
                        norm([intrho01s-QCD01; intrho00s-QCD00; intrho02s-QCD02].')<= MSVZ0;
                        norm([intrho2m2s-QCD2m2; intrho2m1s-QCD2m1; intrho20s-QCD20].')<= MSVZ2;


                        %High energy limit
                        real((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2+imag((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2<= MFF0;
                        real((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2+imag((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2<= MFF;
                        real((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2+imag((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2<= MFF2;

                        cvx_end




                        vecprovup(2)=f2P1;
                        vecprovup(1)=f2S0;
                        if norm(vecprovup-P)<normaup
                            sdup(:)=vecprovup(:);
                            sdsvfup(:)=x;
                            normaup=norm(sdup-P);
                            for iterform=1:m
                                FF1=1+veczz(iterform,:)*ImFl1;
                                spec1=vecsin(iterform,:)*rhol1;

                                sFormfactorup1(iterform)=cons1(iterform)*FF1;
                                sspectralup1(iterform)=spec1;

                                FF0=1+veczz(iterform,:)*ImFl0;
                                spec0=vecsin(iterform,:)*rhol0;

                                sFormfactorup0(iterform)=cons0(iterform)*FF0;
                                sspectralup0(iterform)=spec0;


                                FF2=1+veczz(iterform,:)*ImFl2;
                                spec2=vecsin(iterform,:)*rhol2;

                                sFormfactorup2(iterform)=cons2(iterform)*FF2;
                                sspectralup2(iterform)=spec2;
                            end
                            sFF0(2,:)=ImFl0;
                            sFF1(2,:)=ImFl1;
                            sFF2(2,:)=ImFl2;
                            sspec0(2,:)=rhol0;
                            sspec1(2,:)=rhol1;
                            sspec2(2,:)=rhol2;

                            nh=zeros(m,v);
                            nh(:,:)=H(1,1,:,:);

                            smatrixup(1,:)=1+1i*nh*x;


                            nh=zeros(m,v);
                            nh(:,:)=H(2,1,:,:);

                            smatrixup(2,:)=1+1i*nh*x;


                            nh=zeros(m,v);
                            nh(:,:)=H(1,2,:,:);

                            smatrixup(3,:)=1+1i*nh*x;

                            nh=zeros(m,v);
                            nh(:,:)=H(3,1,:,:);

                            smatrixup(4,:)=1+1i*nh*x;

                            nh=zeros(m,v);
                            nh(:,:)=H(2,2,:,:);

                            smatrixup(5,:)=1+1i*nh*x;

                            nh=zeros(m,v);
                            nh(:,:)=H(3,2,:,:);

                            smatrixup(6,:)=1+1i*nh*x;





                        end



                    sdsvfup(:)=x;

                    for iterform=1:m
                        FF1=1+veczz(iterform,:)*ImFl1;
                        spec1=vecsin(iterform,:)*rhol1;

                        sFormfactorup1(iterform)=cons1(iterform)*FF1;
                        sspectralup1(iterform)=spec1;

                        FF0=1+veczz(iterform,:)*ImFl0;
                        spec0=vecsin(iterform,:)*rhol0;

                        sFormfactorup0(iterform)=cons0(iterform)*FF0;
                        sspectralup0(iterform)=spec0;


                        FF2=1+veczz(iterform,:)*ImFl2;
                        spec2=vecsin(iterform,:)*rhol2;

                        sFormfactorup2(iterform)=cons2(iterform)*FF2;
                        sspectralup2(iterform)=spec2;
                    end

                    sFF0(2,:)=ImFl0;
                    sFF1(2,:)=ImFl1;
                    sFF2(2,:)=ImFl2;
                    sspec0(2,:)=rhol0;
                    sspec1(2,:)=rhol1;
                    sspec2(2,:)=rhol2;

                    nh=zeros(m,v);
                    nh(:,:)=H(1,1,:,:);

                    smatrixup(1,:)=1+1i*nh*x;
                    phaseup1=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup1=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,1,:,:);

                    smatrixup(2,:)=1+1i*nh*x;
                    phaseup2=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup2=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(1,2,:,:);

                    smatrixup(3,:)=1+1i*nh*x;
                    phaseup3=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup3=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,1,:,:);

                    smatrixup(4,:)=1+1i*nh*x;
                    phaseup4=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup4=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,2,:,:);

                    smatrixup(5,:)=1+1i*nh*x;
                    phaseup5=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup5=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,2,:,:);

                    smatrixup(6,:)=1+1i*nh*x;
                    phaseup6=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup6=abs(1+1i*nh*x);








                        sAA(ni,3) = f2P1;
                        fprintf('%f %f %f %f \n', ni, sAA(ni,1), sAA(ni,2),sAA(ni,3));
                        sfrs(:,ni)=x;











                    end

                end


                f2S0bound=sddown(1);
                f2P1bound=sddown(2);
            


                if f2S0max<=lower
                    disp("Less");
                    ni=1;
                    cvx_clear
                    cvx_begin  quiet

                    cvx_solver mosek

                    variable f0
                    variable rho1(nmax,nmax)
                    variable rho2(nmax,nmax) symmetric
                    variable sigma(2*nmax)
                    variable ImFl1(multnmax)
                    variable rhol1(multnmax)
                    variable ImFl0(multnmax)
                    variable rhol0(multnmax)
                    variable ImFl2(multnmax)
                    variable rhol2(multnmax)


                    x1=reshape(rho1,nmax*nmax,1);
                    x2=reshape(rho2,nmax*nmax,1);
                    x=transpose(cat(2,f0 ,transpose(x1),transpose(x2),transpose(sigma)));

                    % functional f2
                    f2S0=F00*x;
                    f2P1=F11*x;
                    fS0unphys=Fun00*x;
                    fP1unphys=Fun11*x;
                    fS2unphys=Fun02*x;

                    fD0unphys=Fun20*x;
                    fF1unphys=Fun13*x;
                    fD2unphys=Fun22*x;




                    %SVZ coefficients l=1 I=1
                    intrho10s=intrho10.'*(vecsin(1:is0,:)*rhol1);
                    intrho1m1s=intrho1m1.'*(vecsin(1:is0,:)*rhol1);
                    intrho11s=intrho11.'*(vecsin(1:is0,:)*rhol1);


                    %SVZ coefficients l=0 I=0
                    intrho00s=intrho00.'*(vecsin(1:is0,:)*rhol0);
                    intrho01s=intrho01.'*(vecsin(1:is0,:)*rhol0);
                    intrho02s=intrho02.'*(vecsin(1:is0,:)*rhol0);

                    %SVZ coefficients l=0 I=2
                    intrho2m2s=intrho2m2.'*(vecsin(1:is0,:)*rhol2);
                    intrho2m1s=intrho2m1.'*(vecsin(1:is0,:)*rhol2);
                    intrho20s=intrho20.'*(vecsin(1:is0,:)*rhol2);



                    minimize(f2P1);
                    subject to

                    f2S0==f2S0max;


                    % regularization
                    norm(x1,4) <= Mreg*10^(power);
                    norm(x2,4) <= Mreg*10^(power);


                    % quadratic cone
                    tensor_mult(imH(: ,1:lmax,:,:),x,[4],[1]).^2+tensor_mult(reH(: ,1:lmax,:,:),x,[4],[1]).^2 <= 2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1])

                    %chiSB

                    norm([fS0unphys-rS0P1chi'.*fP1unphys; fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

                    %optional
                    %norm([fP1unphys])<=MtolCSB;
                    %norm([fP1unphys,fD0unphys,fD2unphys])<=MtolCSB;

                    %Positive definite matrix
                    nh2=zeros(m,v);
                    nh3=zeros(m,v);
                    nh4=zeros(m,v);
                    nh2(:,:)=1i*H(2,1,:,:);
                    nh3(:,:)=1i*H(1,1,:,:);
                    nh4(:,:)=1i*H(1,2,:,:);
                    for iterpo=1:is0
                        FF1=1+veczz(iterpo,:)*ImFl1;
                        spec1=vecsin(iterpo,:)*rhol1;

                        FF0=1+veczz(iterpo,:)*ImFl0;
                        spec0=vecsin(iterpo,:)*rhol0;

                        FF2=1+veczz(iterpo,:)*ImFl2;
                        spec2=vecsin(iterpo,:)*rhol2;


                        %positivesemidefinite l=1, I=1
                        [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1;1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1);cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1]  == hermitian_semidefinite( 3 );
                        %positivesemidefinite l=0, I=0
                        [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0;1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0);cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0]  == hermitian_semidefinite( 3 );
                        %positivesemidefinite l=0, I=2
                        [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2;1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2);cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2]   == hermitian_semidefinite( 3 );


                    end



                    norm([intrho1m1s-QCD1m1; intrho10s-QCD10; intrho11s-QCD11].' )<= MSVZ;
                    norm([intrho01s-QCD01; intrho00s-QCD00; intrho02s-QCD02].')<= MSVZ0;
                    norm([intrho2m2s-QCD2m2; intrho2m1s-QCD2m1; intrho20s-QCD20].')<= MSVZ2;


                    %High energy limit
                    real((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2+imag((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2<= MFF0;
                    real((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2+imag((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2<= MFF;
                    real((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2+imag((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2<= MFF2;


                    cvx_end





                    vecprovdown(2)=f2P1;
                    vecprovdown(1)=f2S0;
                    if norm(vecprovdown-P)<normadown
                        sddown(:)=vecprovdown(:);
                        sdsvfdown(:)=x;
                        normadown=norm(sddown-P);
                        for iterform=1:m
                            FF1=1+veczz(iterform,:)*ImFl1;
                            spec1=vecsin(iterform,:)*rhol1;

                            sFormfactordown1(iterform)=cons1(iterform)*FF1;
                            sspectraldown1(iterform)=spec1;


                            FF0=1+veczz(iterform,:)*ImFl0;
                            spec0=vecsin(iterform,:)*rhol0;

                            sFormfactordown0(iterform)=cons0(iterform)*FF0;
                            sspectraldown0(iterform)=spec0;


                            FF2=1+veczz(iterform,:)*ImFl2;
                            spec2=vecsin(iterform,:)*rhol2;

                            sFormfactordown2(iterform)=cons2(iterform)*FF2;
                            sspectraldown2(iterform)=spec2;
                        end
                        sFF0(1,:)=ImFl0;
                        sFF1(1,:)=ImFl1;
                        sFF2(1,:)=ImFl2;
                        sspec0(1,:)=rhol0;
                        sspec1(1,:)=rhol1;
                        sspec2(1,:)=rhol2;

                        nh=zeros(m,v);
                        nh(:,:)=H(1,1,:,:);


                        smatrixdown(1,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(2,1,:,:);

                        smatrixdown(2,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(1,2,:,:);

                        smatrixdown(3,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(3,1,:,:);

                        smatrixdown(4,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(2,2,:,:);

                        smatrixdown(5,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(3,2,:,:);

                        smatrixdown(6,:)=1+1i*nh*x;




                    end

                    



                    sdsvfdown(:)=x;

                    for iterform=1:m
                        FF1=1+veczz(iterform,:)*ImFl1;
                        spec1=vecsin(iterform,:)*rhol1;

                        sFormfactordown1(iterform)=cons1(iterform)*FF1;
                        sspectraldown1(iterform)=spec1;


                        FF0=1+veczz(iterform,:)*ImFl0;
                        spec0=vecsin(iterform,:)*rhol0;

                        sFormfactordown0(iterform)=cons0(iterform)*FF0;
                        sspectraldown0(iterform)=spec0;


                        FF2=1+veczz(iterform,:)*ImFl2;
                        spec2=vecsin(iterform,:)*rhol2;

                        sFormfactordown2(iterform)=cons2(iterform)*FF2;
                        sspectraldown2(iterform)=spec2;
                    end
                    sFF0(1,:)=ImFl0;
                    sFF1(1,:)=ImFl1;
                    sFF2(1,:)=ImFl2;
                    sspec0(1,:)=rhol0;
                    sspec1(1,:)=rhol1;
                    sspec2(1,:)=rhol2;

                    nh=zeros(m,v);
                    nh(:,:)=H(1,1,:,:);

                    smatrixdown(1,:)=1+1i*nh*x;
                    phasedown1=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown1=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,1,:,:);

                    smatrixdown(2,:)=1+1i*nh*x;
                    phasedown2=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown2=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(1,2,:,:);

                    smatrixdown(3,:)=1+1i*nh*x;
                    phasedown3=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown3=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,1,:,:);

                    smatrixdown(4,:)=1+1i*nh*x;
                    phasedown4=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown4=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,2,:,:);

                    smatrixdown(5,:)=1+1i*nh*x;
                    phasedown5=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown5=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,2,:,:);

                    smatrixdown(6,:)=1+1i*nh*x;
                    phasedown6=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticitydown6=abs(1+1i*nh*x);


                    sAA(ni,2) = f2P1;

                    cvx_clear
                    cvx_begin  quiet

                    cvx_solver mosek

                    variable f0
                    variable rho1(nmax,nmax)
                    variable rho2(nmax,nmax) symmetric
                    variable sigma(2*nmax)
                    variable ImFl1(multnmax)
                    variable rhol1(multnmax)
                    variable ImFl0(multnmax)
                    variable rhol0(multnmax)
                    variable ImFl2(multnmax)
                    variable rhol2(multnmax)



                    x1=reshape(rho1,nmax*nmax,1);
                    x2=reshape(rho2,nmax*nmax,1);
                    x=transpose(cat(2,f0 ,transpose(x1),transpose(x2),transpose(sigma)));

                    % functional f2
                    f2S0=F00*x;
                    f2P1=F11*x;
                    fS0unphys=Fun00*x;
                    fP1unphys=Fun11*x;
                    fS2unphys=Fun02*x;


                    fD0unphys=Fun20*x;
                    fF1unphys=Fun13*x;
                    fD2unphys=Fun22*x;


                    %SVZ coefficients l=1 I=1
                    intrho10s=intrho10.'*(vecsin(1:is0,:)*rhol1);
                    intrho1m1s=intrho1m1.'*(vecsin(1:is0,:)*rhol1);
                    intrho11s=intrho11.'*(vecsin(1:is0,:)*rhol1);


                    %SVZ coefficients l=0 I=0
                    intrho00s=intrho00.'*(vecsin(1:is0,:)*rhol0);
                    intrho01s=intrho01.'*(vecsin(1:is0,:)*rhol0);
                    intrho02s=intrho02.'*(vecsin(1:is0,:)*rhol0);

                    %SVZ coefficients l=0 I=2
                    intrho2m2s=intrho2m2.'*(vecsin(1:is0,:)*rhol2);
                    intrho2m1s=intrho2m1.'*(vecsin(1:is0,:)*rhol2);
                    intrho20s=intrho20.'*(vecsin(1:is0,:)*rhol2);


                    maximize(f2P1);



                    subject to

                    f2S0==f2S0max;



                    % regularization
                    norm(x1,4) <= Mreg*10^(power);
                    norm(x2,4) <= Mreg*10^(power);


                    % quadratic cone
                    tensor_mult(imH(: ,1:lmax,:,:),x,[4],[1]).^2+tensor_mult(reH(: ,1:lmax,:,:),x,[4],[1]).^2 <= 2*tensor_mult(imH2(:,1:lmax,:,:),x,[4],[1])

                    %chiSB

                    norm([fS0unphys-rS0P1chi'.*fP1unphys; fS2unphys-rS2P1chi'.*fP1unphys])<=MtolCSB;

                    %optional
                    %norm([fP1unphys])<=MtolCSB;
                    %norm([fP1unphys,fD0unphys,fD2unphys])<=MtolCSB;

                    %Positive definite matrix
                    nh2=zeros(m,v);
                    nh3=zeros(m,v);
                    nh4=zeros(m,v);
                    nh2(:,:)=1i*H(2,1,:,:);
                    nh3(:,:)=1i*H(1,1,:,:);
                    nh4(:,:)=1i*H(1,2,:,:);
                    for iterpo=1:is0
                        FF1=1+veczz(iterpo,:)*ImFl1;
                        spec1=vecsin(iterpo,:)*rhol1;

                        FF0=1+veczz(iterpo,:)*ImFl0;
                        spec0=vecsin(iterpo,:)*rhol0;

                        FF2=1+veczz(iterpo,:)*ImFl2;
                        spec2=vecsin(iterpo,:)*rhol2;


                        %positivesemidefinite l=1, I=1
                        [1,1+nh2(iterpo,:)*x,cons1(iterpo)*FF1;1+conj(nh2(iterpo,:))*x,1,cons1(iterpo)*conj(FF1);cons1(iterpo)*conj(FF1),cons1(iterpo)*FF1,spec1]  == hermitian_semidefinite( 3 );
                        %positivesemidefinite l=0, I=0
                        [1,1+nh3(iterpo,:)*x,cons0(iterpo)*FF0;1+conj(nh3(iterpo,:))*x,1,cons0(iterpo)*conj(FF0);cons0(iterpo)*conj(FF0),cons0(iterpo)*FF0,spec0]  == hermitian_semidefinite( 3 );
                        %positivesemidefinite l=0, I=2
                        [1,1+nh4(iterpo,:)*x,cons2(iterpo)*FF2;1+conj(nh4(iterpo,:))*x,1,cons2(iterpo)*conj(FF2);cons2(iterpo)*conj(FF2),cons2(iterpo)*FF2,spec2]   == hermitian_semidefinite( 3 );


                    end


                    norm([intrho1m1s-QCD1m1; intrho10s-QCD10; intrho11s-QCD11].' )<= MSVZ;
                    norm([intrho01s-QCD01; intrho00s-QCD00; intrho02s-QCD02].')<= MSVZ0;
                    norm([intrho2m2s-QCD2m2; intrho2m1s-QCD2m1; intrho20s-QCD20].')<= MSVZ2;


                    %High energy limit
                    real((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2+imag((cons0(is0:m).')'.*(1+veczz(is0:m,:)*ImFl0)).^2<= MFF0;
                    real((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2+imag((cons1(is0:m).')'.*(1+veczz(is0:m,:)*ImFl1)).^2<= MFF;
                    real((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2+imag((cons2(is0:m).')'.*(1+veczz(is0:m,:)*ImFl2)).^2<= MFF2;


                    cvx_end




                    vecprovup(2)=f2P1;
                    vecprovup(1)=f2S0;
                    if norm(vecprovup-P)<normaup
                        sdup(:)=vecprovup(:);
                        sdsvfup(:)=x;
                        normaup=norm(sdup-P);
                        for iterform=1:m
                            FF1=1+veczz(iterform,:)*ImFl1;
                            spec1=vecsin(iterform,:)*rhol1;

                            sFormfactorup1(iterform)=cons1(iterform)*FF1;
                            sspectralup1(iterform)=spec1;

                            FF0=1+veczz(iterform,:)*ImFl0;
                            spec0=vecsin(iterform,:)*rhol0;

                            sFormfactorup0(iterform)=cons0(iterform)*FF0;
                            sspectralup0(iterform)=spec0;


                            FF2=1+veczz(iterform,:)*ImFl2;
                            spec2=vecsin(iterform,:)*rhol2;

                            sFormfactorup2(iterform)=cons2(iterform)*FF2;
                            sspectralup2(iterform)=spec2;

                        end
                        sFF0(2,:)=ImFl0;
                        sFF1(2,:)=ImFl1;
                        sFF2(2,:)=ImFl2;
                        sspec0(2,:)=rhol0;
                        sspec1(2,:)=rhol1;
                        sspec2(2,:)=rhol2;

                        nh=zeros(m,v);
                        nh(:,:)=H(1,1,:,:);

                        smatrixup(1,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(2,1,:,:);

                        smatrixup(2,:)=1+1i*nh*x;



                        nh=zeros(m,v);
                        nh(:,:)=H(1,2,:,:);

                        smatrixup(3,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(3,1,:,:);

                        smatrixup(4,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(2,2,:,:);

                        smatrixup(5,:)=1+1i*nh*x;


                        nh=zeros(m,v);
                        nh(:,:)=H(3,2,:,:);

                        smatrixup(6,:)=1+1i*nh*x;




                    end





                    sAA(ni,3) = f2P1;
                    fprintf('%f %f %f %f \n', ni, sAA(ni,1), sAA(ni,2),sAA(ni,3));
                    sfrs(:,ni)=x;

















                    sdsvfup(:)=x;

                    for iterform=1:m
                        FF1=1+veczz(iterform,:)*ImFl1;
                        spec1=vecsin(iterform,:)*rhol1;

                        sFormfactorup1(iterform)=cons1(iterform)*FF1;
                        sspectralup1(iterform)=spec1;

                        FF0=1+veczz(iterform,:)*ImFl0;
                        spec0=vecsin(iterform,:)*rhol0;

                        sFormfactorup0(iterform)=cons0(iterform)*FF0;
                        sspectralup0(iterform)=spec0;


                        FF2=1+veczz(iterform,:)*ImFl2;
                        spec2=vecsin(iterform,:)*rhol2;

                        sFormfactorup2(iterform)=cons2(iterform)*FF2;
                        sspectralup2(iterform)=spec2;

                    end
                    sFF0(2,:)=ImFl0;
                    sFF1(2,:)=ImFl1;
                    sFF2(2,:)=ImFl2;
                    sspec0(2,:)=rhol0;
                    sspec1(2,:)=rhol1;
                    sspec2(2,:)=rhol2;

                    nh=zeros(m,v);
                    nh(:,:)=H(1,1,:,:);

                    smatrixup(1,:)=1+1i*nh*x;
                    phaseup1=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup1=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,1,:,:);

                    smatrixup(2,:)=1+1i*nh*x;
                    phaseup2=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup2=abs(1+1i*nh*x);


                    nh=zeros(m,v);
                    nh(:,:)=H(1,2,:,:);

                    smatrixup(3,:)=1+1i*nh*x;
                    phaseup3=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup3=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,1,:,:);

                    smatrixup(4,:)=1+1i*nh*x;
                    phaseup4=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup4=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(2,2,:,:);

                    smatrixup(5,:)=1+1i*nh*x;
                    phaseup5=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup5=abs(1+1i*nh*x);

                    nh=zeros(m,v);
                    nh(:,:)=H(3,2,:,:);

                    smatrixup(6,:)=1+1i*nh*x;
                    phaseup6=unwrap(angle(1+1i*nh*x))*180/(2*pi);
                    elasticityup6=abs(1+1i*nh*x);







                    sfrs(:,ni)=x;




                end

                toc
                test0up=conj(sFormfactorup0).*sFormfactorup0;
                test0down=conj(sFormfactordown0).*sFormfactordown0;

                test1up=conj(sFormfactorup1).*sFormfactorup1;
                test1down=conj(sFormfactordown1).*sFormfactordown1;

                test2up=conj(sFormfactorup2).*sFormfactorup2;
                test2down=conj(sFormfactordown2).*sFormfactordown2;



                t=tiledlayout(7,3)


                testline=ones(m,1)*MFF;
                testline0=ones(m,1)*MFF0;
                mpi=140/1000;




                ax2 = nexttile;
                hold on

                plot(ax2,sqrt(vnu(4:is0))*mpi,test0up(4:is0))

                plot(ax2,sqrt(vnu(4:is0))*mpi,sspectralup0(4:is0))

                legend('Form Factor Up I=0 l=0','Spectral Density Up','Location','best')
                box off


                ax4 = nexttile;
                hold on
                plot(ax4,sqrt(vnu(4:is0))*mpi,test1up(4:is0))

                plot(ax4,sqrt(vnu(4:is0))*mpi,sspectralup1(4:is0))

                legend('Form Factor Up I=1 l=1','Spectral Density Up','Location','best')
                box off




                ax6 = nexttile;
                hold on
                plot(ax6,sqrt(vnu(4:is0))*mpi,test2up(4:is0))

                plot(ax6,sqrt(vnu(4:is0))*mpi,sspectralup2(4:is0))

                legend('Form Factor Up I=2 l=0 ','Spectral Density Up','Location','best')
                box off




                % Top plot
                ax1 = nexttile;
                hold on

                plot(ax1,sqrt(vnu(4:is0))*mpi,test0down(4:is0))
                plot(ax1,sqrt(vnu(4:is0))*mpi,sspectraldown0(4:is0))
                legend('Form Factor Down I=0 l=0','Spectral Density Down','Location','best')
                box off
                % Top plot

                axis([ax1 ax2],[0.2 2 0 8*10^(-4)])

                ax3 = nexttile;
                hold on

                plot(ax3,sqrt(vnu(4:is0))*mpi,test1down(4:is0))
                plot(ax3,sqrt(vnu(4:is0))*mpi,sspectraldown1(4:is0))
                legend('Form Factor Down I=1 l=1','Spectral Density Down','Location','best')
                box off
                % Top plot
                axis([ax3 ax4],[0.2 2 0 0.1])
                ax5 = nexttile;
                hold on

                plot(ax5,sqrt(vnu(4:is0))*mpi,test2down(4:is0))
                plot(ax5,sqrt(vnu(4:is0))*mpi,sspectraldown2(4:is0))
                legend('Form Factor Down I=2 l=0','Spectral Density Down','Location','best')
                box off
                % Top plot




                axis([ax5 ax6],[0.2 2 0 0.1])















                ax9 = nexttile;
                hold on

                pheno = importdata("pheno.mat");
                finpheno=find(pheno(:,1)>sqrt(vnu(is0))*mpi);
                %finpheno=finpheno(1);
                finpheno=size(pheno(:,1));
                plot(ax9,sqrt(vnu(1:is0))*mpi,phaseup1(1:is0))
                plot(ax9,sqrt(vnu(1:is0))*mpi,phasedown1(1:is0))
                scatter(ax9,pheno(1:finpheno,1),pheno(1:finpheno,2))

                legend('Phase Shift Up S0','Phase Shift Down','Pheno','Location','best')
                box off
                axis([ax9],[0.2 2 -20 190])

                ax10 = nexttile;
                hold on
                mpi=140/1000;
                pheno = importdata("phenop1.mat");
                finpheno=find(pheno(:,1)>sqrt(vnu(is0))*mpi);
                %finpheno=finpheno(1);
                finpheno=size(pheno(:,1));
                plot(ax10,sqrt(vnu(1:is0))*mpi,phaseup2(1:is0))
                plot(ax10,sqrt(vnu(1:is0))*mpi,phasedown2(1:is0))
                scatter(ax10,pheno(1:finpheno,1),pheno(1:finpheno,2))

                legend('Phase Shift Up P1','Phase Shift Down','Pheno','Location','best')
                box off


                ax11 = nexttile;
                hold on
                mpi=140/1000;
                pheno = importdata("phenos2.mat");
                finpheno=find(pheno(:,1)>sqrt(vnu(is0))*mpi);
                %finpheno=finpheno(1);
                finpheno=size(pheno(:,1));
                plot(ax11,sqrt(vnu(1:is0))*mpi,phaseup4(1:is0))
                plot(ax11,sqrt(vnu(1:is0))*mpi,phasedown4(1:is0))
                scatter(ax11,pheno(1:finpheno,1),pheno(1:finpheno,2))

                legend('Phase Shift Up S2','Phase Shift Down','Pheno','Location','best')
                box off


                ax13 = nexttile;
                hold on
                plot(ax13,sqrt(vnu(1:is0))*mpi,elasticityup1(1:is0))
                plot(ax13,sqrt(vnu(1:is0))*mpi,elasticitydown1(1:is0))
                legend('Elasticity Up S0','Elasticity Down','Location','best')
                box off




                ax14= nexttile;
                hold on
                plot(ax14,sqrt(vnu(1:is0))*mpi,elasticityup2(1:is0))
                plot(ax14,sqrt(vnu(1:is0))*mpi,elasticitydown2(1:is0))
                legend('Elasticity Up P1','Elasticity Down','Location','best')
                box off



                ax16= nexttile;
                hold on
                plot(ax16,sqrt(vnu(1:is0))*mpi,elasticityup4(1:is0))
                plot(ax16,sqrt(vnu(1:is0))*mpi,elasticitydown4(1:is0))
                legend('Elasticity Up S2','Elasticity Down','Location','best')
                box off


                ax12 = nexttile;
                hold on
                mpi=140/1000;
                pheno = importdata("phenod0.mat");
                finpheno=find(pheno(:,1)>sqrt(vnu(is0))*mpi);
                %finpheno=finpheno(1);
                finpheno=size(pheno(:,1));
                plot(ax12,sqrt(vnu)*mpi,phaseup3)
                plot(ax12,sqrt(vnu)*mpi,phasedown3)
                scatter(ax12,pheno(1:finpheno,1),pheno(1:finpheno,2))

                legend('Phase Shift Up D0','Phase Shift Down','Pheno','Location','best')
                box off
                axis([ax12],[0.2 3 -60 150])

                ax20 = nexttile;
                hold on
                mpi=140/1000;
                pheno = importdata("phenof1.mat");
                finpheno=find(pheno(:,1)>sqrt(vnu(is0))*mpi);
                %finpheno=finpheno(1);
                finpheno=size(pheno(:,1));
                plot(ax20,sqrt(vnu)*mpi,phaseup5)
                plot(ax20,sqrt(vnu)*mpi,phasedown5)
                scatter(ax20,pheno(1:finpheno,1),pheno(1:finpheno,2))

                legend('Phase Shift Up F1','Phase Shift Down','Pheno','Location','best')
                box off
                axis([ax20],[0.2 3 -6 6])

                ax19 = nexttile;
                hold on
                mpi=140/1000;
                pheno = importdata("phenod2.mat");
                finpheno=find(pheno(:,1)>sqrt(vnu(is0))*mpi);
                %finpheno=finpheno(1);
                finpheno=size(pheno(:,1));
                plot(ax19,sqrt(vnu(1:is0))*mpi,phaseup6(1:is0))
                plot(ax19,sqrt(vnu(1:is0))*mpi,phasedown6(1:is0))
                scatter(ax19,pheno(1:finpheno,1),pheno(1:finpheno,2))

                legend('Phase Shift Up D2','Phase Shift Down','Pheno','Location','best')
                box off

                axis([ax19],[0.2 2 -5       3])

                %axis([ax20],[0.2 2 0 10])

                ax15= nexttile;
                hold on
                plot(ax15,sqrt(vnu(1:is0))*mpi,elasticityup3(1:is0))
                plot(ax15,sqrt(vnu(1:is0))*mpi,elasticitydown3(1:is0))
                legend('Elasticity Up D0','Elasticity Down','Location','best')
                box off









                ax18= nexttile;
                hold on
                plot(ax18,sqrt(vnu(1:is0))*mpi,elasticityup5(1:is0))
                plot(ax18,sqrt(vnu(1:is0))*mpi,elasticitydown5(1:is0))
                legend('Elasticity Up F1','Elasticity Down','Location','best')
                box off

                ax17= nexttile;
                hold on
                plot(ax17,sqrt(vnu(1:is0))*mpi,elasticityup6(1:is0))
                plot(ax17,sqrt(vnu(1:is0))*mpi,elasticitydown6(1:is0))
                legend('Elasticity Up D2','Elasticity Down','Location','best')
                box off





                %
                % Bottom plot
                ax8 = nexttile;
                indexes=find(sAA(:,2)<0);

                hold on

                scatter(ax8,sAA(indexes,1),sAA(indexes,3),'LineWidth',1)
                scatter(ax8,sAA(indexes,1),sAA(indexes,2),'LineWidth',1)
                plot(P(1),P(2),'-o')
                box off
                grid(ax8,'on')


                axis([ax10],[0.2 2 -100 450])

                axis([ax11],[0.2 2 -35 0])

                axis([ax13 ax14 ax15 ax16 ax17 ax18],[0.2 2 0 1])






                exportgraphics(t,join(['Experiments/Images/spec=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2) 'multnmax=' num2str(multnmax) '.PDF'],""),'ContentType','vector','Resolution',1000)

                %MtolCSB=100;
                if MtolCSB==100
                    fid = fopen(join(['Experiments/MAT/SVZ/testnumberM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) '.txt'],""),'wt');
                    fprintf(fid, '%f',fstore);

                    fclose(fid);
                end
                %save space of f_00 f_11
                fid = fopen(join(['Experiments/MAT/SVZ/testconvSPfspanM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2) '.txt'],""),'wt');

                for ni=1:Mline
                    fprintf(fid, '%f %f %f %f \n', ni, sAA(ni,1), sAA(ni,2),sAA(ni,3));
                end
                fclose(fid);
                %save coefficents rho sigma f
                fid = fopen(join(['Experiments/MAT/SVZ/valfrhosigmaspanM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                for ni=1:Mline
                    fprintf(fid, '%f \n', sfrs(:,ni));
                end
                fclose(fid);

                %save coefficents rho sigma f nearest point up
                fid = fopen(join(['Experiments/MAT/SVZ/valfrhosigmaUPspanM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sdsvfup(:));
                fclose(fid);
                %save nearest point
                fid = fopen(join(['Experiments/MAT/SVZ/nearestUPpointM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sdup(:));
                fclose(fid);



                %save coefficents rho sigma f nearest point down
                fid = fopen(join(['Experiments/MAT/SVZ/valfrhosigmaDOWNspanM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sdsvfdown(:));
                fclose(fid);
                %save nearest point
                fid = fopen(join(['Experiments/MAT/SVZ/nearestDOWNpointM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sddown(:));
                fclose(fid);




                %save spectral density down
                fid = fopen(join(['Experiments/MAT/SVZ/specdensityDOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sspectraldown1);
                fclose(fid);

                %save spectral density up
                fid = fopen(join(['Experiments/MAT/SVZ/specdensityUPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sspectralup1);
                fclose(fid);

                %save Form Factor down
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorDOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', conj(sFormfactordown1).*sFormfactordown1);
                fclose(fid);
                %save Form Factor up
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorUPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', conj(sFormfactorup1).*sFormfactorup1);
                fclose(fid);



                %save spectral density 0 coefficients down
                fid = fopen(join(['Experiments/MAT/SVZ/specdensitycoeffS0DOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sspec0(1,:));
                fclose(fid);

                %save spectral density  0 coefficients up
                fid = fopen(join(['Experiments/MAT/SVZ/specdensitycoeffS0UPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sspec0(2,:));
                fclose(fid);

                %save Form Factor  0 coefficients down
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorcoeffS0DOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sFF0(1,:));
                fclose(fid);
                %save Form Factor  0 coefficients up
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorcoeffS0UPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sFF0(2,:));
                fclose(fid);


                %save spectral density  1 coefficients down
                fid = fopen(join(['Experiments/MAT/SVZ/specdensitycoeffP1DOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sspec1(1,:));
                fclose(fid);

                %save spectral density  1 coefficients up
                fid = fopen(join(['Experiments/MAT/SVZ/specdensitycoeffP1UPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sspec1(2,:));
                fclose(fid);

                %save Form Factor 1 coefficients down
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorcoeffP1DOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sFF1(1,:));
                fclose(fid);
                %save Form Factor  1 coefficients up
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorcoeffP1UPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sFF1(2,:));
                fclose(fid);


                %save spectral density  2 coefficients down
                fid = fopen(join(['Experiments/MAT/SVZ/specdensitycoeffD0DOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sspec2(1,:));
                fclose(fid);

                %save spectral density  2 coefficients up
                fid = fopen(join(['Experiments/MAT/SVZ/specdensitycoeffD0UPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sspec2(2,:));
                fclose(fid);








                %save Form Factor  2 coefficients down
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorcoeffD0DOWNM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sFF2(1,:));
                fclose(fid);
                %save Form Factor  2 coefficients up
                fid = fopen(join(['Experiments/MAT/SVZ/formfactorcoeffD0UPM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');
                fprintf(fid, '%f \n', sFF2(2,:));
                fclose(fid);




                %save xplot axis
                fid = fopen(join(['Experiments/MAT/SVZ/XplotaxisM=' num2str(M) 'nmax=' num2str(nmax) 'lmax=' num2str(lmax) 'MCSB=' num2str(MtolCSB) 'MSVZ0=' num2str(MSVZ0) 'MFF0=' num2str(MFF0) 'MSVZ1=' num2str(MSVZ) 'MFF1=' num2str(MFF) 'MSVZ2=' num2str(MSVZ2) 'MFF2=' num2str(MFF2)  '.txt'],""),'wt');

                fprintf(fid, '%f \n', sqrt(vnu)*mpi);
                fclose(fid);

            end
        end
    end
end

beep on; beep;
clear;

function C = tensor_mult(A, B, sum_idx_A, sum_idx_B)

sum_idx_A = reshape(sum_idx_A, 1, numel(sum_idx_A));
sum_idx_B = reshape(sum_idx_B, 1, numel(sum_idx_B));

size_A = size(A);
size_B = size(B);

num_size_A = numel(size_A);
num_size_B = numel(size_B);

perm_A = 1:num_size_A;
perm_B = 1:num_size_B;

num_idx = size(sum_idx_A, 2);

assert(isequal(size(sum_idx_A) , size(sum_idx_B) ));
assert(isequal(size_A(sum_idx_A), size_B(sum_idx_B)));

sum_dim = prod(size_A(sum_idx_A));

for i = 1:num_idx
    perm_A = perm_A(perm_A~=sum_idx_A(i));
    perm_B = perm_B(perm_B~=sum_idx_B(i));
end

perm_A = [perm_A, sum_idx_A];
perm_B = [sum_idx_B, perm_B];

size_A(sum_idx_A) = [];
size_B(sum_idx_B) = [];

if isempty(size_A)
    size_A = 1;
end
if isempty(size_B)
    size_B = 1;
end

C = squeeze(reshape(...
    reshape(permute(A, perm_A), [prod(size_A), sum_dim]) * ...
    reshape(permute(B, perm_B), [sum_dim, prod(size_B)]), ...
    [size_A,size_B]));

end

