% PME3380 - Modelagem de sistemas dinâmicos
% Exercício 1
% 
% Grupo W
% Arthur Poubel Bonamigo - 12549391
% João Victor Ramos Müller - 10705975
% Leonardo Torres de Oliveira Groot - 12624482
% Pedro Vicente Marino - 12553952

% Definição dos parâmetros
    mi = 1;
    g = 9.81;
    NUSP = [1, 5, 2, 2];
    D = sum(NUSP);
    lambda = 1 + (1/20)*(D);
    wp = 1;
    l1 = g/(wp^2);
    l2= l1/lambda;
    m1 = l1*mi;
    m2 = m1/lambda;

% Condições iniciais C1 e C2 como as linhas 1 e dois do vetor C
    C = [pi/90,pi/120,0,0;pi/3,pi/6,0,0];
    C1 = C(1,:)';
    C2 = C(2,:)';

% intervalo de tempo de numero de passos
    tin = [0 40];
    t = linspace(tin(1),tin(2),1000);

% integração numérica x_método(M1 ou M2)_(linear ou não linear)_Condição(1 ou 2)

    % integração numérica das duas condições de contorno pelo M1 - ode23
            % modelo linear
        sM1lin1=ode23(@linear,tin,C1);
        sM1lin2=ode23(@linear,tin,C2);
        xM1lin1 = deval(sM1lin1,t);
        xM1lin2 = deval(sM1lin2,t);
            % modelo não linear
        sM1nlin1=ode23(@naolinear,tin,C1);
        sM1nlin2=ode23(@naolinear,tin,C2);
        xM1nlin1 = deval(sM1nlin1,t);
        xM1nlin2 = deval(sM1nlin2,t);
    
    % integração numérica das duas condições de contorno pelo M2 - ode89
            % modelo linear
        sM2lin1=ode89(@linear,tin,C1);
        sM2lin2=ode89(@linear,tin,C2);
        xM2lin1 = deval(sM2lin1,t);
        xM2lin2 = deval(sM2lin2,t);
            % modelo não linear
        sM2nlin1=ode89(@naolinear,tin,C1);
        sM2nlin2=ode89(@naolinear,tin,C2);
        xM2nlin1 = deval(sM2nlin1,t);
        xM2nlin2 = deval(sM2nlin2,t);


% energias potencial, cinética e mecânica para os modelos linear e não linear
    for k = 1:length(t)
        % cinética
        TM1lin1(k) = (l1^2)*((xM1lin1(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM1lin1(1,k) - xM1lin1(2,k))*xM1lin1(3,k)*xM1lin1(4,k)/2 + ((l2*xM1lin1(4,k)).^2)*m2/6 ;
        TM1lin2(k) = (l1^2)*((xM1lin2(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM1lin2(1,k) - xM1lin2(2,k))*xM1lin2(3,k)*xM1lin2(4,k)/2 + ((l2*xM1lin2(4,k)).^2)*m2/6 ; ;
        TM1nlin1(k) = (l1^2)*((xM1nlin1(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM1nlin1(1,k) - xM1nlin1(2,k))*xM1nlin1(3,k)*xM1nlin1(4,k)/2 + ((l2*xM1nlin1(4,k)).^2)*m2/6; 
        TM1nlin2(k) = (l1^2)*((xM1nlin2(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM1nlin2(1,k) - xM1nlin2(2,k))*xM1nlin2(3,k)*xM1nlin2(4,k)/2 + ((l2*xM1nlin2(4,k)).^2)*m2/6 ;
        TM2lin1(k) = (l1^2)*((xM2lin1(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM2lin1(1,k) - xM2lin1(2,k))*xM2lin1(3,k)*xM2lin1(4,k)/2 + ((l2*xM2lin1(4,k)).^2)*m2/6 ;
        TM2lin2(k) = (l1^2)*((xM2lin2(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM2lin2(1,k) - xM2lin2(2,k))*xM2lin2(3,k)*xM2lin2(4,k)/2 + ((l2*xM2lin2(4,k)).^2)*m2/6 ;
        TM2nlin1(k) = (l1^2)*((xM2nlin1(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM2nlin1(1,k) - xM2nlin1(2,k))*xM2nlin1(3,k)*xM2nlin1(4,k)/2 + ((l2*xM2nlin1(4,k)).^2)*m2/6; 
        TM2nlin2(k) = (l1^2)*((xM2nlin2(3,k)).^2)*(((m1)/6) + ((m2)/2)) + l2*l1*m2*cos(xM2nlin2(1,k) - xM2nlin2(2,k))*xM2nlin2(3,k)*xM2nlin2(4,k)/2 + ((l2*xM2nlin2(4,k)).^2)*m2/6  ;
        % Potencial
        VM1lin1(k) = g*l1*m1*(1 - cos(xM1lin1(1,k)))/2 + g*m2*( l1*(1 - cos(xM1lin1(1,k))) +  (l2/2)*(1 - cos(xM1lin1(2,k))));
        VM1lin2(k) = g*l1*m1*(1 - cos(xM1lin2(1,k)))/2 + g*m2*( l1*(1 - cos(xM1lin2(1,k))) +  (l2/2)*(1 - cos(xM1lin2(2,k))));
        VM1nlin1(k) = g*l1*m1*(1 - cos(xM1nlin1(1,k)))/2 + g*m2*( l1*(1 - cos(xM1nlin1(1,k))) +  (l2/2)*(1 - cos(xM1nlin1(2,k))));
        VM1nlin2(k) = g*l1*m1*(1 - cos(xM1nlin2(1,k)))/2 + g*m2*( l1*(1 - cos(xM1nlin2(1,k))) +  (l2/2)*(1 - cos(xM1nlin2(2,k))));
        VM2lin1(k) = g*l1*m1*(1 - cos(xM2lin1(1,k)))/2 + g*m2*( l1*(1 - cos(xM2lin1(1,k))) +  (l2/2)*(1 - cos(xM2lin1(2,k))));
        VM2lin2(k) = g*l1*m1*(1 - cos(xM2lin2(1,k)))/2 + g*m2*( l1*(1 - cos(xM2lin2(1,k))) +  (l2/2)*(1 - cos(xM2lin2(2,k))));
        VM2nlin1(k) = g*l1*m1*(1 - cos(xM2nlin1(1,k)))/2 + g*m2*( l1*(1 - cos(xM2nlin1(1,k))) +  (l2/2)*(1 - cos(xM2nlin1(2,k))));
        VM2nlin2(k) = g*l1*m1*(1 - cos(xM2nlin2(1,k)))/2 + g*m2*( l1*(1 - cos(xM2nlin2(1,k))) +  (l2/2)*(1 - cos(xM2nlin2(2,k))));
        % mecânica
        EM1lin1(k) = TM1lin1(k)+ VM1lin1(k);
        EM1lin2(k) =  TM1lin2(k) + VM1lin2(k);
        EM1nlin1(k) = TM1nlin1(k) + VM1nlin1(k);
        EM1nlin2(k) =  TM1nlin2(k)+  VM1nlin2(k);
        EM2lin1(k) = TM2lin1(k) +VM2lin1(k);
        EM2lin2(k) = TM2lin2(k)+ VM2lin2(k);
        EM2nlin1(k) =  TM2nlin1(k) + VM2nlin1(k);
        EM2nlin2(k) = TM2nlin2(k)+ VM2nlin2(k);
    end



% % % % %  Escolha o gráfico que quer plotar e aperte ctrl+t 
% % % % %  para descomentar a seção

% Plotagem dos gráficos
    
%     % Gráficos gerais
%     
%         % Cenário C1
%         
%             figure(1)
%             plot(t,(xM1lin1(1,:)*180/pi),'r',t,(xM2lin1(1,:)*180/pi),'m',t,(xM1nlin1(1,:)*180/pi),'b',t,(xM2nlin1(1,:)*180/pi),'c')
%             %title('Gráfico 1: \theta_1 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M1','\theta_1 linear M2','\theta_1 não linear M1','\theta_1 não linear M2')
%             
%             figure(2)
%             plot(t,(xM1lin1(2,:)*180/pi),'r',t,(xM2lin1(2,:)*180/pi),'m',t,(xM1nlin1(2,:)*180/pi),'b',t,(xM2nlin1(2,:)*180/pi),'c')
%             %title('Gráfico 2: \theta_2 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M1','\theta_2 linear M2','\theta_2 não linear M1','\theta_2 não linear M2')
%             
%             figure(3)
%             plot(t,(xM1lin1(3,:)*180/pi),'r',t,(xM2lin1(3,:)*180/pi),'m',t,(xM1nlin1(3,:)*180/pi),'b',t,(xM2nlin1(3,:)*180/pi),'c')
%             %title('Gráfico 3:  d\theta_1/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M1','d\theta_1/dt linear M2','d\theta_1/dt não linear M1','d\theta_1/dt não linear M2')
%             
%             figure(4)
%             plot(t,(xM1lin1(4,:)*180/pi),'r',t,(xM2lin1(4,:)*180/pi),'m',t,(xM1nlin1(4,:)*180/pi),'b',t,(xM2nlin1(4,:)*180/pi),'c')
%             %title('Gráfico 4:  d\theta_2/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M1','d\theta_2/dt linear M2','d\theta_2/dt não linear M1','d\theta_2/dt não linear M2')
%             
%             figure(5)
%             plot(t,EM1lin1,'r',t,EM2lin1,'m',t,EM1nlin1,'b',t,EM2nlin1,'c')
%             %title('Gráfico 5: Energia mecânica no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('Energia mecânica (J)')
%             legend('Energia mecânica linear M1','Energia mecânica linear M2','Energia mecânica não linear M1','Energia mecânica não linear M2')
%             axis([tin(1) tin(2) 0.5 1])
%         
%         % Cenário C2
%         
%             figure(6)
%             plot(t,(xM1lin2(1,:)*180/pi),'r',t,(xM2lin2(1,:)*180/pi),'m',t,(xM1nlin2(1,:)*180/pi),'b',t,(xM2nlin2(1,:)*180/pi),'c')
%             %title('Gráfico 6: \theta_1 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M1','\theta_1 linear M2','\theta_1 não linear M1','\theta_1 não linear M2')
%             
%             figure(7)
%             plot(t,(xM1lin2(2,:)*180/pi),'r',t,(xM2lin2(2,:)*180/pi),'m',t,(xM1nlin2(2,:)*180/pi),'b',t,(xM2nlin2(2,:)*180/pi),'c')
%             %title('Gráfico 7: \theta_2 no cenário C2')
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M1','\theta_2 linear M2','\theta_2 não linear M1','\theta_2 não linear M2')
%             
%             figure(8)
%             plot(t,(xM1lin2(3,:)*180/pi),'r',t,(xM2lin2(3,:)*180/pi),'m',t,(xM1nlin2(3,:)*180/pi),'b',t,(xM2nlin2(3,:)*180/pi),'c')
%             %title('Gráfico 8:  d\theta_1/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M1','d\theta_1/dt linear M2','d\theta_1/dt não linear M1','d\theta_1/dt não linear M2')
%             
%             figure(9)
%             plot(t,(xM1lin2(4,:)*180/pi),'r',t,(xM2lin2(4,:)*180/pi),'m',t,(xM1nlin2(4,:)*180/pi),'b',t,(xM2nlin2(4,:)*180/pi),'c')
%             %title('Gráfico 9:  d\theta_2/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M1','d\theta_2/dt linear M2','d\theta_2/dt não linear M1','d\theta_2/dt não linear M2')
%             
%             figure(10)
%             plot(t,EM1lin2,'r',t,EM2lin2,'m',t,EM1nlin2,'b',t,EM2nlin2,'c')
%             %title('Gráfico 10: Energia mecânica no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('Energia mecânica (J)')
%             legend('Energia mecânica linear M1','Energia mecânica linear M2','Energia mecânica não linear M1','Energia mecânica não linear M2')
%             axis([0 tin(2) 500 650])
        
    % Gráficos específicos - item g

%         % C1 - theta_1
% 
%             figure(11)
%             plot(t,(xM1lin1(1,:)*180/pi),'r',t,(xM1nlin1(1,:)*180/pi),'b')
%             %title('Gráfico 11: \theta_1 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M1','\theta_1 não linear M1')
% 
%             figure(12)
%             plot(t,(xM2lin1(1,:)*180/pi),'m',t,(xM2nlin1(1,:)*180/pi),'c')
%             %title('Gráfico 12: \theta_1 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M2','\theta_1 não linear M2')
% 
%             figure(13)
%             plot(t,(xM1lin1(1,:)*180/pi),'r',t,(xM2lin1(1,:)*180/pi),'m')
%             %title('Gráfico 13: \theta_1 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M1','\theta_1 linear M2')
% 
%             figure(14)
%             plot(t,(xM1nlin1(1,:)*180/pi),'b',t,(xM2nlin1(1,:)*180/pi),'c')
%             %title('Gráfico 14: \theta_1 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 não linear M1','\theta_1 não linear M2')
 
%         % C1 - theta_2
% 
%             figure(15)
%             plot(t,(xM1lin1(2,:)*180/pi),'r',t,(xM1nlin1(2,:)*180/pi),'b')
%             %title('Gráfico 15: \theta_2 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M1','\theta_2 não linear M1')
%            
%             figure(16)
%             plot(t,(xM2lin1(2,:)*180/pi),'m',t,(xM2nlin1(2,:)*180/pi),'c')
%             %title('Gráfico 16: \theta_2 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M2','\theta_2 não linear M2')
% 
%             figure(17)
%             plot(t,(xM1lin1(2,:)*180/pi),'r',t,(xM2lin1(2,:)*180/pi),'m')
%             %title('Gráfico 18: \theta_2 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M1','\theta_2 linear M2')
% 
%             figure(18)
%             plot(t,(xM1nlin1(2,:)*180/pi),'b',t,(xM2nlin1(2,:)*180/pi),'c')
%             %title('Gráfico 17: \theta_2 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 não linear M1','\theta_2 não linear M2')

%         % C1 - dtheta_1
% 
%             figure(19)
%             plot(t,(xM1lin1(3,:)*180/pi),'r',t,(xM1nlin1(3,:)*180/pi),'b')
%             %title('Gráfico 21:  d\theta_1/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M1','d\theta_1/dt não linear M1')
% 
%             figure(20)
%             plot(t,(xM2lin1(3,:)*180/pi),'m',t,(xM2nlin1(3,:)*180/pi),'c')
%             %title('Gráfico 21:  d\theta_1/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M2','d\theta_1/dt não linear M2')
% 
%             figure(21)
%             plot(t,(xM1lin1(3,:)*180/pi),'r',t,(xM2lin1(3,:)*180/pi),'m')
%             %title('Gráfico 21:  d\theta_1/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M1','d\theta_1/dt linear M2')
% 
%             figure(22)
%             plot(t,(xM1nlin1(3,:)*180/pi),'b',t,(xM2nlin1(3,:)*180/pi),'c')
%             %title('Gráfico 22:  d\theta_1/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt não linear M1','d\theta_1/dt não linear M2')

%         % C1 - dtheta_2
% 
%             figure(23)
%             plot(t,(xM1lin1(4,:)*180/pi),'r',t,(xM1nlin1(4,:)*180/pi),'b')
%             %title('Gráfico 23:  d\theta_2/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M1','d\theta_2/dt não linear M1')
% 
%             figure(24)
%             plot(t,(xM2lin1(4,:)*180/pi),'m',t,(xM2nlin1(4,:)*180/pi),'c')
%             %title('Gráfico 24:  d\theta_2/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M2','d\theta_2/dt não linear M2')
% 
%             figure(25)
%             plot(t,(xM1lin1(4,:)*180/pi),'r',t,(xM2lin1(4,:)*180/pi),'m')
%             %title('Gráfico 25:  d\theta_2/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M1','d\theta_2/dt linear M2')
% 
%             figure(26)
%             plot(t,(xM1nlin1(4,:)*180/pi),'b',t,(xM2nlin1(4,:)*180/pi),'c')
%             %title('Gráfico 26:  d\theta_2/dt  no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt não linear M1','d\theta_2/dt não linear M2')


%         % C2 - theta_1
% 
%             figure(27)
%             plot(t,(xM1lin2(1,:)*180/pi),'r',t,(xM1nlin2(1,:)*180/pi),'b')
%             %title('Gráfico 11: \theta_1 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M1','\theta_1 não linear M1')
% 
%             figure(28)
%             plot(t,(xM2lin2(1,:)*180/pi),'m',t,(xM2nlin2(1,:)*180/pi),'c')
%             %title('Gráfico 12: \theta_1 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M2','\theta_1 não linear M2')
% 
%             figure(29)
%             plot(t,(xM1lin2(1,:)*180/pi),'r',t,(xM2lin2(1,:)*180/pi),'m')
%             %title('Gráfico 13: \theta_1 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 linear M1','\theta_1 linear M2')
% 
%             figure(30)
%             plot(t,(xM1nlin2(1,:)*180/pi),'b',t,(xM2nlin2(1,:)*180/pi),'c')
%             %title('Gráfico 14: \theta_1 no cenário C1')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_1 (graus)')
%             legend('\theta_1 não linear M1','\theta_1 não linear M2')

%         % C2 - theta_2
% 
%             figure(31)
%             plot(t,(xM1lin2(2,:)*180/pi),'r',t,(xM1nlin2(2,:)*180/pi),'b')
%             %title('Gráfico 15: \theta_2 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M1','\theta_2 não linear M1')
%            
%             figure(32)
%             plot(t,(xM2lin2(2,:)*180/pi),'m',t,(xM2nlin2(2,:)*180/pi),'c')
%             %title('Gráfico 16: \theta_2 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M2','\theta_2 não linear M2')
% 
%             figure(33)
%             plot(t,(xM1lin2(2,:)*180/pi),'r',t,(xM2lin2(2,:)*180/pi),'m')
%             %title('Gráfico 18: \theta_2 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 linear M1','\theta_2 linear M2')
% 
%             figure(34)
%             plot(t,(xM1nlin2(2,:)*180/pi),'b',t,(xM2nlin2(2,:)*180/pi),'c')
%             %title('Gráfico 17: \theta_2 no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('\theta_2 (graus)')
%             legend('\theta_2 não linear M1','\theta_2 não linear M2')

%         % C2 - dtheta_1
% 
%             figure(35)
%             plot(t,(xM1lin2(3,:)*180/pi),'r',t,(xM1nlin2(3,:)*180/pi),'b')
%             %title('Gráfico 35:  d\theta_1/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M1','d\theta_1/dt não linear M1')
% 
%             figure(36)
%             plot(t,(xM2lin2(3,:)*180/pi),'m',t,(xM2nlin2(3,:)*180/pi),'c')
%             %title('Gráfico 36:  d\theta_1/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M2','d\theta_1/dt não linear M2')
% 
%             figure(37)
%             plot(t,(xM1lin2(3,:)*180/pi),'r',t,(xM2lin2(3,:)*180/pi),'m')
%             %title('Gráfico 37:  d\theta_1/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt linear M1','d\theta_1/dt linear M2')
% 
%             figure(38)
%             plot(t,(xM1nlin2(3,:)*180/pi),'b',t,(xM2nlin2(3,:)*180/pi),'c')
%             %title('Gráfico 38:  d\theta_1/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_1/dt (graus/s)')
%             legend('d\theta_1/dt não linear M1','d\theta_1/dt não linear M2')

%         % C2 - dtheta_2
% 
%             figure(39)
%             plot(t,(xM1lin2(4,:)*180/pi),'r',t,(xM1nlin2(4,:)*180/pi),'b')
%             %title('Gráfico 39:  d\theta_2/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M1','d\theta_2/dt não linear M1')
% 
%             figure(40)
%             plot(t,(xM2lin2(4,:)*180/pi),'m',t,(xM2nlin2(4,:)*180/pi),'c')
%             %title('Gráfico 40:  d\theta_2/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M2','d\theta_2/dt não linear M2')
% 
%             figure(41)
%             plot(t,(xM1lin2(4,:)*180/pi),'r',t,(xM2lin2(4,:)*180/pi),'m')
%             %title('Gráfico 41:  d\theta_2/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt linear M1','d\theta_2/dt linear M2')
% 
%             figure(42)
%             plot(t,(xM1nlin2(4,:)*180/pi),'b',t,(xM2nlin2(4,:)*180/pi),'c')
%             %title('Gráfico 42:  d\theta_2/dt  no cenário C2')
%             grid("on")
%             xlabel('t(s)')
%             ylabel('d\theta_2/dt (graus/s)')
%             legend('d\theta_2/dt não linear M1','d\theta_2/dt não linear M2')



% %     % Gráficos específicos - item h
%         
%         figure(43)
%         plot(t,EM1nlin2,'b',t,EM2nlin2,'c')
%         %title('Gráfico 11: Energia mecânica no cenário C2 com o modelo não-linear')
%         grid("on")
%         xlabel('t(s)')
%         ylabel('Energia mecânica (J)')
%         legend('Energia mecânica não linear M1','Energia mecânica não linear M2')
%         axis([0 tin(2) 500 650])


function dy = linear(t,y)
%parâmetros
mi = 1;
g = 9.81;
NUSP= [1, 2, 5, 2];
D = sum(NUSP);
lambda = 1 + (1/20)*(D);
wp = 1;
l1 = g/(wp^2);
l2= l1/lambda;
m1 = l1*mi;
m2 = m1/lambda;
%equações da dinâmica
dy = [y(3);y(4); -3*(wp^2)*((2*lambda*y(1)) + (4*y(1)) - (3*y(2)))/((4*lambda) + 3);   3*lambda*(wp^2)*( (y(1)*((3*lambda) + 6)) - (y(2)*((2*lambda) + 6)))/((4*lambda) + 3)];   
end

function dy = naolinear(t,y)
%parâmetros
mi = 1;
g = 9.81;
NUSP= [1, 2, 5, 2];
D = sum(NUSP);
lambda = 1 + (1/20)*(D);
wp = 1;
l1 = g/(wp^2);
l2= l1/lambda;
m1 = l1*mi;
m2 = m1/lambda;
%equações da dinâmica
dy = [y(3);y(4); 3*((2*((lambda*wp)^2)*sin(y(1)))+(4*lambda*(wp^2)*sin(y(1)))-(3*lambda*(wp^2)*sin(y(2))*cos(y(1) - y(2)))+(3*lambda*sin(y(1) - y(2))*cos(y(1) - y(2))*(y(3).^2))+(2*sin(y(1) - y(2))*(y(4).^2)))/(lambda*((-4*lambda)+(9*((cos(y(1)-y(2)).^2))-12))) ;  (-3)*(     (3*((lambda*wp)^2)*sin(y(1))*cos(y(1)-y(2)))-(2*((lambda*wp)^2)*sin(y(2)))+(2*(lambda^2)*sin(y(1)-y(2))*(y(3).^2)) + (6*lambda*(wp^2)*sin(y(1))*cos(y(1)-y(2)))  - (6*lambda*(wp^2)*sin(y(2))) + (6*lambda*sin(y(1)-y(2))*(y(3).^2)) + (3*sin(y(1)-y(2))*cos(y(1)-y(2))*(y(4).^2)))/((-4*lambda)+(9*((cos(y(1)-y(2)).^2))-12))];   
end

