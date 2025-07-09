clear all
close all
clc

MCS=input('Enter number of Monte Carlo iterations:')
mcs=1;
while mcs<MCS+1
    bdata_not_per_unit=[%Bus P(Kw) Q(Kvar)
        1 0 0
        2 120 72
        3 108 48
        4 144 96
        5 72 36
        6 72 24
        7 240 120
        8 240 120
        9 72 24
        10 72 24
        11 54 36
        12 72 42
        13 72 42
        14 144 96
        15 72 12
        16 72 24
        17 72 24
        18 108 48
        19 108 48
        20 108 48
        21 108 48
        22 108 48
        23 108 60
        24 504 240
        25 504 240
        26 72 30
        27 72 30
        28 72 24
        29 144 84
        30 240 720
        31 180 84
        32 252 120
        33 72 48 ];
    for i=2:33
        u(i)=rand;

        bdata_not_per_unit(i,2)=icdf('normal',u(i),bdata_not_per_unit(i,2),35);
        bdata_not_per_unit(i,3)=icdf('normal',u(i),bdata_not_per_unit(i,3),15);
    end
    ldata_not_per_unit=[
        % Inbus Outbus Resistance(ohm) Reactance(ohm)
        1 2 0.0922 0.0470
        2 3 0.4930 0.2511
        3 4 0.3660 0.1864
        4 5 0.3811 0.1941
        5 6 0.8191 0.7070
        6 7 0.1872 0.6188
        7 8 0.7114 0.2351
        8 9 1.0300 0.7400
        9 10 1.0440 0.7400
        10 11 0.1966 0.0650
        11 12 0.3744 0.1238
        12 13 1.4680 1.1550
        13 14 0.5416 0.7129
        14 15 0.5910 0.5260
        15 16 0.7463 0.5450
        16 17 1.2890 1.7210
        17 18 0.7320 0.5740
        2 19 0.1640 0.1565
        19 20 1.5042 1.3554
        20 21 0.4095 0.4784
        21 22 0.7089 0.9373
        3 23 0.4512 0.3083
        23 24 0.8980 0.7091
        24 25 0.8960 0.7011
        6 26 0.2030 0.1034
        26 27 0.2842 0.1447
        27 28 1.0590 0.9377
        28 29 0.8042 0.7006
        29 30 0.5075 0.2585
        30 31 0.9744 0.9630
        31 32 0.3105 0.3619
        32 33 0.3410 0.5302];
    sizbdata=size(bdata_not_per_unit);
    busnum=sizbdata(1,1);
    sizldata=size(ldata_not_per_unit);
    branchnum=sizldata(1,1);
    %per unit calculation:
    Sbase=10^3;
    Vbase=12.66*10^3;
    Zbase=Vbase^2/Sbase;
    bdata=bdata_not_per_unit;
    ldata=ldata_not_per_unit;
    for n=1:busnum
        bdata(n,2)=(bdata(n,2)*1000)/Sbase;
        bdata(n,3)=(bdata(n,3)*1000)/Sbase;
    end
    for n=1:branchnum
        ldata(n,3)=ldata(n,3)/Zbase;
        ldata(n,4)=ldata(n,4)/Zbase;
    end
    %per unit calculation finished
    terminatebus=zeros(busnum,1);
    intermediatebus=zeros(busnum,1);
    junctionbus=zeros(busnum,1);
    junctionnum=zeros(busnum,1);
    refbus=0;
    busI=zeros(busnum,1);
    v=ones(1,busnum);
    I=zeros(busnum,busnum);
    for k=1:busnum
        co=0;
        l=0;
        for n=1:branchnum
            if ldata(n,1)==k
                co=co+1;
            end
        end
        if co==0
            terminatebus(k,1)=k;
        elseif co>=2
            junctionbus(k,1)=k;
            junctionnum(k,1)=co;
        elseif co==1
            for m=1:branchnum
                l=l+1;
                if ldata(m,2)==k
                    intermediatebus(k,1)=k;
                    break
                elseif l==branchnum
                    refbus=k;
                end
            end
        end
    end
    junctionbus;

    intermediatebus;
    terminatebus;
    refbus;
    junctionnum;
    tempterminatebus=terminatebus;
    controljunctionnum=zeros(busnum,1);
    k=0;
    c=0;
    itecount=0;
    for s=1:15 %iteration
        itecount=itecount+1;
        %backward sweep
        while c==0
            k=k+1;
            juncnum=0;
            n=0;
            stop=0;
            previousI=0;
            if tempterminatebus(k,1)==k
                while(n<branchnum)&&(stop==0)
                    n=n+1;
                    if ldata(n,2)==k
                        a=ldata(n,1);
                        if a==refbus
                            c=1;
                        end
                        I(a,k)=busI(k)+(bdata(k,2)-1i*bdata(k,3))/conj(v(k))+previousI;
                        previousI=I(a,k);
                        tempterminatebus(k,1)=0;
                        if junctionbus(a,1)==a
                            busI(a)=busI(a)+I(a,k);
                            controljunctionnum(a,1)=controljunctionnum(a,1)+1;
                            if controljunctionnum(a,1)==junctionnum(a,1)
                                tempterminatebus(a,1)=a;
                            end
                            break
                        end
                        k=a;
                        n=0;
                    end
                end
                k=0;
            end
        end
        %end of backward sweep
        %forward sweep
        count=0;
        beforev=v;
        newldata=ldata;
        stop1=0;
        forwardbus=zeros(busnum,1);
        c=0;
        stopif=0;
        while stop1==0
            for k=1:branchnum
                if(newldata(k,1)==refbus)&&(stopif==0)
                    c=refbus;
                    if junctionbus(refbus,1)==refbus
                        forwardbus(refbus,1)=1;
                    end
                    stop2=0;
                    while 1
                        a=newldata(k,2);
                        z=newldata(k,3)+1i*newldata(k,4);
                        v(a)=v(c)-z*I(c,a);
                        newldata(k,:)=0;
                        if junctionbus(a,1)==a
                            forwardbus(a,1)=1;
                        end
                        if terminatebus(a,1)==a
                            stopif=1 ;
                            stop1=1;
                            break
                        end
                        for n=1:branchnum
                            if ldata(n,1)==a
                                c=a;
                                k=n;
                                break
                            end
                        end
                    end
                end
            end
        end
        stop3=0;
        while stop3==0
            stopif1=0;
            for k=1:busnum
                if (forwardbus(k,1)==1)&& (stopif1==0)
                    stop4=0;
                    while stop4==0
                        stopif2=0;
                        hhh=0;
                        counter=0;
                        for n=1:branchnum
                            counter=counter+1;
                            if (newldata(n,1)==k)&& (stopif2==0)
                                hhh=1;
                                c=k;
                                a=newldata(n,2);
                                z=newldata(n,3)+1i*newldata(n,4);
                                v(a)=v(c)-z*I(c,a);
                                newldata(n,:)=0;
                                if junctionbus(a,1)==a
                                    forwardbus(a,1)=1;
                                end
                                if terminatebus(a,1)==a
                                    if newldata==0
                                        stopif1=1;
                                        stop3=1;
                                        stopif2=1;
                                        stop4=1;
                                        break
                                    end
                                    stopif2=1;
                                    stop4=1;
                                end
                                k=a;
                                stopif2=1;
                            end
                            if (counter==branchnum) && (hhh==0)
                                stop4=1;
                            end
                        end
                    end
                end
            end
        end
        %end of forward sweep
        for p=2:busnum
            if abs(v(p)-beforev(p))<=0.000001
                count=count+1;
            end
        end
        if count==(busnum-1)
            break
        end
    end

    %data display
    voldisp=zeros(busnum,3);
    curdisp=zeros(branchnum,3);
    co=0;
    for n=1:busnum
        for k=1:busnum
            if I(n,k)~=0
                co=co+1;
                curdisp(co,1)=n;
                curdisp(co,2)=k;
                curdisp(co,3)=I(n,k);
            end
        end
    end
    for n=1:busnum
        voldisp(n,1)=n;
        voldisp(n,2)=abs(v(n));
        voldisp(n,3)=angle(v(n))*180;
    end
    co=0;
    W=zeros(1,branchnum);
    a=zeros(1,branchnum);
    b=zeros(1,branchnum);
    for n=1:busnum
        for k=1:busnum
            if I(n,k)~=0
                co=co+1;
                a(co)=n;
                b(co)=k;
            end
        end
    end
    DG = sparse(a,b,true,busnum,busnum);
    Ploss=0; Qloss=0;
    for k=1:branchnum
        z=curdisp(k,1);
        y=curdisp(k,2);
        for g=1:branchnum
            if ldata(g,1)==z && ldata(g,2)==y
                Ploss=Ploss+abs(curdisp(k,3))^2*ldata(g,3);
                %ploss=r*|I|^2
                Qloss=Qloss+abs(curdisp(k,3))^2*ldata(g,4);
                %qloss=x*|I|^2
                break
            end
        end
    end
    Inew=zeros(busnum,1);
    for l=1:busnum
        Inew(l,1)=(bdata(l,2)-1i*bdata(l,3))/conj(v(l));
    end
    %Outputs
    P_loss=abs(Ploss);
    Q_loss=abs(Qloss);
    Loss(mcs)=P_loss;
    for i=1:33
        voltage(i,mcs)=abs(v(i));
    end
        mcs = mcs + 1;

end

% After all iterations, calculate and display average losses and voltage results
average_loss = mean(Loss);  % Mean of power losses across all iterations
disp(['Average Active Power Loss (P_loss): ', num2str(average_loss), ' kW']);

% Display voltage magnitudes for each bus (mean over iterations)
average_voltage = mean(voltage, 2);  % Mean voltage for each bus over iterations
disp('Average Voltage Magnitudes at Each Bus:');
disp(average_voltage);

% Optionally, you can plot the voltage magnitudes to visualize the results
figure;
plot(1:33, average_voltage, '-o');
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Average Voltage Magnitudes for Each Bus after Monte Carlo Simulation');
grid on;

