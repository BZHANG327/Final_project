function [ L,U ] = LU_ND( A,sign )
%This is the function compute the LU factorization based on  the ND 
%dissection technique
%sign = 0 means A is the matrix corresponding to square multigrids.
%sign = 1 means A is the matrix corresponding to more row multigrids.
%sign = -1 means A is the matrix corresponding to more col multigrids.
[N,M] = size(A);
if N~=M
    fprintf("error,A should be a square matrix");
    return
end
if N<=9
    [L,U] = lu(A);
else
    if sign == 0
        nrow = sqrt(N);
        %ncol = nrow;
        subrow = round(ceil(nrow/2)-mod(nrow,2));
        cutting1 = subrow^2;
        [L111,U111] = LU_ND(A(1:cutting1,1:cutting1),0);
        if mod(nrow,2)% if nrow is odd number
            cutting2 = cutting1 + subrow^2;
            [L122,U122] = LU_ND(A(cutting1+1:cutting2,cutting1+1:cutting2),0);
            cutting3 = cutting2+subrow;
            A131 = A(cutting2+1:cutting3,1:cutting1);
            A132 = A(cutting2+1:cutting3,cutting1+1:cutting2);
            A113 = A131';
            A123 = A132';
            A133 = A(cutting2+1:cutting3,cutting2+1:cutting3);
            S133 = A133-A131/U111/L111*A113-A132/U122/L122*A123;
            [L133,U133] = lu(S133);
            clear S133
            L11 = sparse([L111,zeros(cutting2-cutting1,cutting3-cutting1);...
                zeros(cutting2-cutting1,cutting1),L122,...
                zeros(cutting2-cutting1,cutting3-cutting2);...
                A131/U111,A132/U122,L133]);
            U11 = sparse([U111,zeros(cutting1,cutting2-cutting1),L111\A113;...
                zeros(cutting2-cutting1,cutting1),U122,L122\A123;...
                zeros(cutting3-cutting2,cutting2),U133]);
            clear L111 L122 U111 U122 A131 A132 A113 A123 A123 A133 ;
            cutting4 = cutting3 + subrow^2;
            cutting5 = cutting4 + subrow^2;
            cutting6 = cutting5 + subrow;
            cutting7 = cutting6 + nrow;
            [L211,U211] = LU_ND(A(cutting3+1:cutting4,cutting3+1:cutting4),0);
            [L222,U222] = LU_ND(A(cutting4+1:cutting5,cutting4+1:cutting5),0);
            A231 = A(cutting5+1:cutting6,cutting3+1:cutting4);
            A232 = A(cutting5+1:cutting6,cutting4+1:cutting5);
            A213 = A231';
            A223 = A232';
            A233 = A(cutting5+1:cutting6,cutting5+1:cutting6);
            S233 = A233-A231/U211/L211*A213-A232/U222/L222*A223;
            [L233,U233] = lu(S233);
            clear S233
            L22 = sparse([L211,zeros(cutting5-cutting4,cutting6-cutting4);...
                zeros(cutting5-cutting4,cutting4-cutting3),L222,...
                zeros(cutting5-cutting4,cutting6-cutting5);...
                A231/U211,A232/U222,L233]);
            U22 = sparse([U211,zeros(cutting4-cutting3,cutting5-cutting4),L211\A213;...
                zeros(cutting5-cutting4,cutting4-cutting3),U222,L222\A223;...
                zeros(cutting6-cutting5,cutting5-cutting3),U233]);
            clear L211 L222 U211 U222 A231 A232 A213 A223 A223 A233 ;
            A31 = A(cutting6+1:cutting7,1:cutting3);
            A32 = A(cutting6+1:cutting7,cutting3+1:cutting6);
            A33 = A(cutting6+1:cutting7,cutting6+1:cutting7);
            A13 = A31';
            A23 = A32';
            S33 = A33-A31/U11/L11*A13-A32/U22/L22*A23;
            [L33,U33] = lu(S33);
            clear S33;
            
            L = sparse([L11,zeros(cutting3,cutting7-cutting3);...
                zeros(cutting6-cutting3,cutting3),L22,...
                zeros(cutting6-cutting3,cutting7-cutting6);...
                A31/U11,A32/U22,L33]);
            U = sparse([U11,zeros(cutting3,cutting6-cutting3),L11\A13;...
                zeros(cutting6-cutting3,cutting3),U22,L22\A23;...
                zeros(cutting7-cutting6,cutting6),U33]);
            
        else
            cutting2 = cutting1 + subrow*(subrow-1);
            cutting3 = cutting2 + subrow;
            cutting4 =cutting3 + subrow*(subrow-1);
            cutting5 = cutting4 + (subrow-1)^2;
            cutting6 = cutting5 + subrow-1;
            cutting7 = cutting6 + nrow;
            [L122,U122] = LU_ND(A(cutting1+1:cutting2,cutting1+1:cutting2),-1);
            A131 = A(cutting2+1:cutting3,1:cutting1);
            A132 = A(cutting2+1:cutting3,cutting1+1:cutting2);
            A113 = A131';
            A123 = A132';
            A133 = A(cutting2+1:cutting3,cutting2+1:cutting3);
            S133 = A133-A131/U111/L111*A113-A132/U122/L122*A123;
            [L133,U133] = lu(S133);
            clear S133
            L11 = sparse([L111,zeros(cutting1,cutting3-cutting1);...
                zeros(cutting2-cutting1,cutting1),L122,...
                zeros(cutting2-cutting1,cutting3-cutting2);...
                A131/U111,A132/U122,L133]);
            U11 = sparse([U111,zeros(cutting1,cutting2-cutting1),L111\A113;...
                zeros(cutting2-cutting1,cutting1),U122,L122\A123;...
                zeros(cutting3-cutting2,cutting2),U133]);
            clear L111 L122 U111 U122 A131 A132 A113 A123 A123 A133 ;
            [L211,U211] = LU_ND(A(cutting3+1:cutting4,cutting3+1:cutting4),1);
            [L222,U222] = LU_ND(A(cutting4+1:cutting5,cutting4+1:cutting5),0);
            A231 = A(cutting5+1:cutting6,cutting3+1:cutting4);
            A232 = A(cutting5+1:cutting6,cutting4+1:cutting5);
            A213 = A231';
            A223 = A232';
            A233 = A(cutting5+1:cutting6,cutting5+1:cutting6);
            S233 = A233-A231/U211/L211*A213-A232/U222/L222*A223;
            [L233,U233] = lu(S233);
            clear S233
            L22 = sparse([L211,zeros(cutting4-cutting3,cutting6-cutting4);...
                zeros(cutting5-cutting4,cutting4-cutting3),L222,...
                zeros(cutting5-cutting4,cutting6-cutting5);...
                A231/U211,A232/U222,L233]);
            U22 = sparse([U211,zeros(cutting4-cutting3,cutting5-cutting4),L211\A213;...
                zeros(cutting5-cutting4,cutting4-cutting3),U222,L222\A223;...
                zeros(cutting6-cutting5,cutting5-cutting3),U233]);
            clear L211 L222 U211 U222 A231 A232 A213 A223 A223 A233 ;
            A31 = A(cutting6+1:cutting7,1:cutting3);
            A32 = A(cutting6+1:cutting7,cutting3+1:cutting6);
            A33 = A(cutting6+1:cutting7,cutting6+1:cutting7);
            A13 = A31';
            A23 = A32';
            S33 = A33-A31/U11/L11*A13-A32/U22/L22*A23;
            [L33,U33] = lu(S33);
            clear S33;
            L = sparse([L11,zeros(cutting3,cutting7-cutting3);...
                zeros(cutting6-cutting3,cutting3),L22,...
                zeros(cutting6-cutting3,cutting7-cutting6);...
                A31/U11,A32/U22,L33]);
            U = sparse([U11,zeros(cutting3,cutting6-cutting3),L11\A13;...
                zeros(cutting6-cutting3,cutting3),U22,L22\A23;...
                zeros(cutting7-cutting6,cutting6),U33]);
        end
    elseif sign == 1
        nrow = ceil(sqrt(N));
        %ncol = nrow-1;
        subrow = round(ceil(nrow/2)-mod(nrow,2));
        if mod(nrow,2) %odd rows
            cutting1 = subrow^2;
            cutting2 = cutting1 + subrow^2;
            cutting3 = cutting2 + subrow;
            cutting4 = cutting3 + subrow*(subrow-1);
            cutting5 = cutting4 + subrow*(subrow-1);
            cutting6 = cutting5+subrow-1;
            cutting7 = cutting6+ nrow;
            [L111,U111] = LU_ND(A(1:cutting1,1:cutting1),0);
            [L122,U122] = LU_ND(A(cutting1+1:cutting2,cutting1+1:cutting2),0);
            A131 = A(cutting2+1:cutting3,1:cutting1);
            A132 = A(cutting2+1:cutting3,cutting1+1:cutting2);
            A113 = A131';
            A123 = A132';
            A133 = A(cutting2+1:cutting3,cutting2+1:cutting3);
            S133 = A133-A131/U111/L111*A113-A132/U122/L122*A123;
            [L133,U133] = lu(S133);
            clear S133
            L11 = sparse([L111,zeros(cutting1,cutting3-cutting1);...
                zeros(cutting2-cutting1,cutting1),L122,...
                zeros(cutting2-cutting1,cutting3-cutting2);...
                A131/U111,A132/U122,L133]);
            U11 = sparse([U111,zeros(cutting1,cutting2-cutting1),L111\A113;...
                zeros(cutting2-cutting1,cutting1),U122,L122\A123;...
                zeros(cutting3-cutting2,cutting2),U133]);
            clear L111 L122 U111 U122 A131 A132 A113 A123 A123 A133;
            [L211,U211] = LU_ND(A(cutting3+1:cutting4,cutting3+1:cutting4),1);
            [L222,U222] = LU_ND(A(cutting4+1:cutting5,cutting4+1:cutting5),1);
            A231 = A(cutting5+1:cutting6,cutting3+1:cutting4);
            A232 = A(cutting5+1:cutting6,cutting4+1:cutting5);
            A213 = A231';
            A223 = A232';
            A233 = A(cutting5+1:cutting6,cutting5+1:cutting6);
            S233 = A233-A231/U211/L211*A213-A232/U222/L222*A223;
            [L233,U233] = lu(S233);
            clear S233
            L22 = sparse([L211,zeros(cutting4-cutting3,cutting6-cutting4);...
                zeros(cutting5-cutting4,cutting4-cutting3),L222,...
                zeros(cutting5-cutting4,cutting6-cutting5);...
                A231/U211,A232/U222,L233]);
            U22 = sparse([U211,zeros(cutting4-cutting3,cutting5-cutting4),L211\A213;...
                zeros(cutting5-cutting4,cutting4-cutting3),U222,L222\A223;...
                zeros(cutting6-cutting5,cutting5-cutting3),U233]);
            clear L211 L222 U211 U222 A231 A232 A213 A223 A223 A233 ;
            A31 = A(cutting6+1:cutting7,1:cutting3);
            A32 = A(cutting6+1:cutting7,cutting3+1:cutting6);
            A33 = A(cutting6+1:cutting7,cutting6+1:cutting7);
            A13 = A31';
            A23 = A32';
            S33 = A33-A31/U11/L11*A13-A32/U22/L22*A23;
            [L33,U33] = lu(S33);
            clear S33;
            L = sparse([L11,zeros(cutting3,cutting7-cutting3);...
                zeros(cutting6-cutting3,cutting3),L22,...
                zeros(cutting6-cutting3,cutting7-cutting6);...
                A31/U11,A32/U22,L33]);
            U = sparse([U11,zeros(cutting3,cutting6-cutting3),L11\A13;...
                zeros(cutting6-cutting3,cutting3),U22,L22\A23;...
                zeros(cutting7-cutting6,cutting6),U33]);    
        else
            cutting1 = subrow*(subrow-1);
            cutting2 = (subrow-1)^2 + cutting1;
            cutting3 = cutting2 + subrow-1;
            cutting4 = cutting3 + subrow*(subrow-1);
            cutting5 = cutting4 + (subrow-1)^2;
            cutting6 = cutting5+subrow-1;
            cutting7 = cutting6 + nrow;
            [L111,U111] = LU_ND(A(1:cutting1,1:cutting1),1);
            [L122,U122] = LU_ND(A(cutting1+1:cutting2,cutting1+1:cutting2),0);
            A131 = A(cutting2+1:cutting3,1:cutting1);
            A132 = A(cutting2+1:cutting3,cutting1+1:cutting2);
            A113 = A131';
            A123 = A132';
            A133 = A(cutting2+1:cutting3,cutting2+1:cutting3);
            S133 = A133-A131/U111/L111*A113-A132/U122/L122*A123;
            [L133,U133] = lu(S133);
            clear S133
            L11 = sparse([L111,zeros(cutting1,cutting3-cutting1);...
                zeros(cutting2-cutting1,cutting1),L122,...
                zeros(cutting2-cutting1,cutting3-cutting2);...
                A131/U111,A132/U122,L133]);
            U11 = sparse([U111,zeros(cutting1,cutting2-cutting1),L111\A113;...
                zeros(cutting2-cutting1,cutting1),U122,L122\A123;...
                zeros(cutting3-cutting2,cutting2),U133]);
            clear L111 L122 U111 U122 A131 A132 A113 A123 A123 A133;
            [L211,U211] = LU_ND(A(cutting3+1:cutting4,cutting3+1:cutting4),1);
            [L222,U222] = LU_ND(A(cutting4+1:cutting5,cutting4+1:cutting5),0);
            A231 = A(cutting5+1:cutting6,cutting3+1:cutting4);
            A232 = A(cutting5+1:cutting6,cutting4+1:cutting5);
            A213 = A231';
            A223 = A232';
            A233 = A(cutting5+1:cutting6,cutting5+1:cutting6);
            S233 = A233-A231/U211/L211*A213-A232/U222/L222*A223;
            [L233,U233] = lu(S233);
            clear S233
            L22 = sparse([L211,zeros(cutting4-cutting3,cutting6-cutting4);...
                zeros(cutting5-cutting4,cutting4-cutting3),L222,...
                zeros(cutting5-cutting4,cutting6-cutting5);...
                A231/U211,A232/U222,L233]);
            U22 = sparse([U211,zeros(cutting4-cutting3,cutting5-cutting4),L211\A213;...
                zeros(cutting5-cutting4,cutting4-cutting3),U222,L222\A223;...
                zeros(cutting6-cutting5,cutting5-cutting3),U233]);
            clear L211 L222 U211 U222 A231 A232 A213 A223 A223 A233 ;
            A31 = A(cutting6+1:cutting7,1:cutting3);
            A32 = A(cutting6+1:cutting7,cutting3+1:cutting6);
            A33 = A(cutting6+1:cutting7,cutting6+1:cutting7);
            A13 = A31';
            A23 = A32';
            S33 = A33-A31/U11/L11*A13-A32/U22/L22*A23;
            [L33,U33] = lu(S33);
            clear S33;
            L = sparse([L11,zeros(cutting3,cutting7-cutting3);...
                zeros(cutting6-cutting3,cutting3),L22,...
                zeros(cutting6-cutting3,cutting7-cutting6);...
                A31/U11,A32/U22,L33]);
            U = sparse([U11,zeros(cutting3,cutting6-cutting3),L11\A13;...
                zeros(cutting6-cutting3,cutting3),U22,L22\A23;...
                zeros(cutting7-cutting6,cutting6),U33]);    
        end
    elseif sign==-1
        ncol = ceil(sqrt(N));
        nrow = ncol-1;
        subrow = round(ceil(nrow/2)-mod(nrow,2));
        if mod(nrow,2)
            cutting1 = subrow*(subrow+1);
            cutting2 = cutting1 + subrow*(subrow+1);
            cutting3 = cutting2 + subrow+1;
            cutting4 = cutting3 + subrow^2;
            cutting5 = cutting4 + subrow^2;
            cutting6 = cutting5 + subrow;
            cutting7 = cutting6 + nrow;
            [L111,U111] = LU_ND(A(1:cutting1,1:cutting1),-1);
            [L122,U122] = LU_ND(A(cutting1+1:cutting2,cutting1+1:cutting2),-1);
            A131 = A(cutting2+1:cutting3,1:cutting1);
            A132 = A(cutting2+1:cutting3,cutting1+1:cutting2);
            A113 = A131';
            A123 = A132';
            A133 = A(cutting2+1:cutting3,cutting2+1:cutting3);
            S133 = A133-A131/U111/L111*A113-A132/U122/L122*A123;
            [L133,U133] = lu(S133);
            clear S133
            L11 = sparse([L111,zeros(cutting1,cutting3-cutting1);...
                zeros(cutting2-cutting1,cutting1),L122,...
                zeros(cutting2-cutting1,cutting3-cutting2);...
                A131/U111,A132/U122,L133]);
            U11 = sparse([U111,zeros(cutting1,cutting2-cutting1),L111\A113;...
                zeros(cutting2-cutting1,cutting1),U122,L122\A123;...
                zeros(cutting3-cutting2,cutting2),U133]);
            clear L111 L122 U111 U122 A131 A132 A113 A123 A123 A133;
            [L211,U211] = LU_ND(A(cutting3+1:cutting4,cutting3+1:cutting4),0);
            [L222,U222] = LU_ND(A(cutting4+1:cutting5,cutting4+1:cutting5),0);
            A231 = A(cutting5+1:cutting6,cutting3+1:cutting4);
            A232 = A(cutting5+1:cutting6,cutting4+1:cutting5);
            A213 = A231';
            A223 = A232';
            A233 = A(cutting5+1:cutting6,cutting5+1:cutting6);
            S233 = A233-A231/U211/L211*A213-A232/U222/L222*A223;
            [L233,U233] = lu(S233);
            clear S233
            L22 = sparse([L211,zeros(cutting4-cutting3,cutting6-cutting4);...
                zeros(cutting5-cutting4,cutting4-cutting3),L222,...
                zeros(cutting5-cutting4,cutting6-cutting5);...
                A231/U211,A232/U222,L233]);
            U22 = sparse([U211,zeros(cutting4-cutting3,cutting5-cutting4),L211\A213;...
                zeros(cutting5-cutting4,cutting4-cutting3),U222,L222\A223;...
                zeros(cutting6-cutting5,cutting5-cutting3),U233]);
            clear L211 L222 U211 U222 A231 A232 A213 A223 A223 A233 ;
            A31 = A(cutting6+1:cutting7,1:cutting3);
            A32 = A(cutting6+1:cutting7,cutting3+1:cutting6);
            A33 = A(cutting6+1:cutting7,cutting6+1:cutting7);
            A13 = A31';
            A23 = A32';
            S33 = A33-A31/U11/L11*A13-A32/U22/L22*A23;
            [L33,U33] = lu(S33);
            clear S33;
            L = sparse([L11,zeros(cutting3,cutting7-cutting3);...
                zeros(cutting6-cutting3,cutting3),L22,...
                zeros(cutting6-cutting3,cutting7-cutting6);...
                A31/U11,A32/U22,L33]);
            U = sparse([U11,zeros(cutting3,cutting6-cutting3),L11\A13;...
                zeros(cutting6-cutting3,cutting3),U22,L22\A23;...
                zeros(cutting7-cutting6,cutting6),U33]);    
            
        else
            cutting1 = subrow^2;
            cutting2 = cutting1+(subrow-1)*subrow;
            cutting3 = cutting2 + subrow;
            cutting4 = cutting3+subrow^2;
            cutting5 = cutting4+(subrow-1)*subrow;
            cutting6 = cutting5+subrow;
            cutting7 = cutting6+nrow;
            [L111,U111] = LU_ND(A(1:cutting1,1:cutting1),0);
            [L122,U122] = LU_ND(A(cutting1+1:cutting2,cutting1+1:cutting2),-1);
            A131 = A(cutting2+1:cutting3,1:cutting1);
            A132 = A(cutting2+1:cutting3,cutting1+1:cutting2);
            A113 = A131';
            A123 = A132';
            A133 = A(cutting2+1:cutting3,cutting2+1:cutting3);
            S133 = A133-A131/U111/L111*A113-A132/U122/L122*A123;
            [L133,U133] = lu(S133);
            clear S133
            L11 = sparse([L111,zeros(cutting1,cutting3-cutting1);...
                zeros(cutting2-cutting1,cutting1),L122,...
                zeros(cutting2-cutting1,cutting3-cutting2);...
                A131/U111,A132/U122,L133]);
            U11 = sparse([U111,zeros(cutting1,cutting2-cutting1),L111\A113;...
                zeros(cutting2-cutting1,cutting1),U122,L122\A123;...
                zeros(cutting3-cutting2,cutting2),U133]);
            clear L111 L122 U111 U122 A131 A132 A113 A123 A123 A133;
            [L211,U211] = LU_ND(A(cutting3+1:cutting4,cutting3+1:cutting4),0);
            [L222,U222] = LU_ND(A(cutting4+1:cutting5,cutting4+1:cutting5),-1);
            A231 = A(cutting5+1:cutting6,cutting3+1:cutting4);
            A232 = A(cutting5+1:cutting6,cutting4+1:cutting5);
            A213 = A231';
            A223 = A232';
            A233 = A(cutting5+1:cutting6,cutting5+1:cutting6);
            S233 = A233-A231/U211/L211*A213-A232/U222/L222*A223;
            [L233,U233] = lu(S233);
            clear S233
            L22 = sparse([L211,zeros(cutting4-cutting3,cutting6-cutting4);...
                zeros(cutting5-cutting4,cutting4-cutting3),L222,...
                zeros(cutting5-cutting4,cutting6-cutting5);...
                A231/U211,A232/U222,L233]);
            U22 = sparse([U211,zeros(cutting4-cutting3,cutting5-cutting4),L211\A213;...
                zeros(cutting5-cutting4,cutting4-cutting3),U222,L222\A223;...
                zeros(cutting6-cutting5,cutting5-cutting3),U233]);
            clear L211 L222 U211 U222 A231 A232 A213 A223 A223 A233 ;
            A31 = A(cutting6+1:cutting7,1:cutting3);
            A32 = A(cutting6+1:cutting7,cutting3+1:cutting6);
            A33 = A(cutting6+1:cutting7,cutting6+1:cutting7);
            A13 = A31';
            A23 = A32';
            S33 = A33-A31/U11/L11*A13-A32/U22/L22*A23;
            [L33,U33] = lu(S33);
            clear S33;
            L = sparse([L11,zeros(cutting3,cutting7-cutting3);...
                zeros(cutting6-cutting3,cutting3),L22,...
                zeros(cutting6-cutting3,cutting7-cutting6);...
                A31/U11,A32/U22,L33]);
            U = sparse([U11,zeros(cutting3,cutting6-cutting3),L11\A13;...
                zeros(cutting6-cutting3,cutting3),U22,L22\A23;...
                zeros(cutting7-cutting6,cutting6),U33]);    
        end
        
    end
            
end

