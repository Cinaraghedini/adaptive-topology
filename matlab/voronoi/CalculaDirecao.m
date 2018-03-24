function [ direcao ] = CalculaDirecao(xA, yA, xB, yB, xC, yC)
            
            %Situa??o 1
            alpha1 = abs(atan((xB - xA) / (yB - yA)));
            quadrante1 = 0;
            if xB >= xA && yB >= yA
                quadrante1 = 1;
            end
            if xB <= xA && yB >= yA
                alpha1 = 2.0 * pi - alpha1;
                quadrante1 = 2;
            end
            if xB <= xA&& yB <= yA
                alpha1 = pi + alpha1;
                quadrante1 = 3;
            end
            if xB >= xA && yB <= yA
                alpha1 = pi - alpha1;
                quadrante1 = 4;
            end
            
            %Situa??o 2
            alpha2 = abs(atan((xC - xB) / (yC - yB)));
            quadrante2 = 0;
            if xC >= xB && yC >= yB
                quadrante2 = 1;
            end
            if xC <= xB && yC >= yB
                alpha2 = 2.0 * pi - alpha2;
                quadrante2 = 2;
            end
            if xC <= xB && yC <= yB
                alpha2 = pi + alpha2;
                quadrante2 = 3;
            end
            if xC >= xB && yC <= yB
                alpha2 = pi - alpha2;
                quadrante2 = 4;
            end
            if alpha2 >= alpha1
                direcao = 'direita';
            else
                direcao = 'esquerda';
            end
            if (quadrante1 == 1 && quadrante2 == 2) || (quadrante1 == 2 && quadrante2 == 1)
                if strcmp(direcao, 'esquerda')
                    direcao = 'direita';
                else
                    direcao = 'esquerda';
                end
            end
            
        end