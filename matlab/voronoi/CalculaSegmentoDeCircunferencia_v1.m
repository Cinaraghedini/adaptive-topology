function [x,y] = CalculaSegmentoDeCircunferencia(p1, p2, pc, raio, pcGeral)
% function [x,y] = CalculateArc(p1, p2, pc, raio, pcGeral)
% Calcula segmento de curva entre dois pontos
% Parameters: p1 - (x,y)
%             p2 - (x,y)
%             pc - (x,y) - centroide da celula voronoi
%             raio - raio da circunferencia com centro em pc
%             pcGeral - ponto central resultante do circleFitting sobre
%             todos os pontos "verdes"
% Output:   x,y - Pontos do segmento de circunferencia
% Technological Institute of Aeronautics
% Author: Nicolas Pereira Borges - nicolas@ita.br
% Date: 06/10/2016

% Calcula o valor de theta
theta1 = atan2(p1(2) - pc(2), p1(1) - pc(1));
theta2 = atan2(p2(2) - pc(2), p2(1) - pc(1));

% Verifica casos 1, 3 e 4
if (p1(1) < pc(1) && p2(1) < pc(1)) ||... % Caso 1 - Dois pontos a esquerda de PC
        (p1(2) > pc(2) && p2(2) > pc(2)) ||... % Caso 3 - Dois pontos acima de PC
        (p1(2) < pc(2) && p2(2) < pc(2))       % Caso 4 - Dois pontos abaixo de PC
    
    % Os dois thetas devem ser positivos
    % Eh uma "curva" para esquerda
    
    % Se o theta1 for negativo, converte para positivo
    if theta1 < 0
        theta1 = 2*pi + theta1;
    end
    
    % Se o theta2 for negativo, converte para positivo
    if theta2 < 0
        theta2 = 2*pi + theta2;
    end
    
    % Verifica caso 2
elseif p1(1) > pc(1) &&  p2(1) > pc(1) % Caso 2 - Dois pontos a direita de PC
    
    % Se p1 esta acima de p2
    if p1(2) > p2(2)
        % Theta1 deve ser positivo e Theta2 negativo
        
        % Se o theta1 for negativo, converte para positivo
        if theta1 < 0
            theta1 = 2*pi + theta1;
        end
        
        % Se o theta2 for positivo, converte para negativo
        if theta2 > 0
            theta2 = theta2 - 2*pi;
        end
        
        % Se p2 esta acima de p1
    else
        % Theta2 deve ser positivo e Theta1 negativo
        
        % Se o theta1 for positivo, converte para negativo
        if theta1 > 0
            theta1 = theta1 - 2*pi;
        end
        
        % Se o theta2 for negativo, converte para positivo
        if theta2 < 0
            theta2 = 2*pi + theta2;
        end
        
    end
    % Caso 5 ou 6
else
    % Verifica se cell esta a esquerda do pc geral
    if pc(1) < pcGeral(1)
        % theta1 e theta2 devem ser positivos
        % Curva "a esquerda"
        
        % Se o theta1 for negativo, converte para positivo
        if theta1 < 0
            theta1 = 2*pi + theta1;
        end
        
        % Se o theta2 for negativo, converte para positivo
        if theta2 < 0
            theta2 = 2*pi + theta2;
        end
        % Se cell estiver a direita do pc geral
    else
        % Existira 1 theta positivo e 1 negativo
        % Curva "a direita"
        
        % Se p1 estiver acima de p2
        if p1(2) > p2(2)
            % theta1 devera ser positivo e theta2 negativo
            % Se o theta1 for negativo, converte para positivo
            if theta1 < 0
                theta1 = 2*pi + theta1;
            end
            
            % Se o theta2 for positivo, converte para negativo
            if theta2 > 0
                theta2 = theta2 - 2*pi;
            end
        else
            % Theta2 deve ser positivo e Theta1 negativo
            
            % Se o theta1 for positivo, converte para negativo
            if theta1 > 0
                theta1 = theta1 - 2*pi;
            end
            
            % Se o theta2 for negativo, converte para positivo
            if theta2 < 0
                theta2 = 2*pi + theta2;
            end
        end
    end
end

% Calcula range do theta
if theta1 < theta2
    th = theta1:0.1:theta2;
else
    th = theta1:-0.1:theta2;
end

% Calcula pontos do segmento da circunferencia
x = raio * cos(th) + pc(1);
y = raio * sin(th) + pc(2);

end