function [LF_EDT,LF_T20,LF_T30,LF_T60] = GetLinearFits(L,vectorTiempo)
    [~, Minus5DbIdx] = min(abs(L+5));
    [~, Minus15DbIdx] = min(abs(L+15));
    [~, Minus25DbIdx] = min(abs(L+25));
    [~, Minus35DbIdx] = min(abs(L+35));
    [~, Minus65DbIdx] = min(abs(L+65));
    
    X_EDT = vectorTiempo(Minus5DbIdx:Minus15DbIdx);
    X_T20 = vectorTiempo(Minus5DbIdx:Minus25DbIdx);
    X_T30 = vectorTiempo(Minus5DbIdx:Minus35DbIdx);
    X_T60 = vectorTiempo(Minus5DbIdx:Minus65DbIdx);
    
    L_EDT = L(Minus5DbIdx:Minus15DbIdx); 
    L_T20 = L(Minus5DbIdx:Minus25DbIdx); 
    L_T30 = L(Minus5DbIdx:Minus35DbIdx); 
    L_T60 = L(Minus5DbIdx:Minus65DbIdx); 
    
    P_EDT = polyfit(X_EDT,L_EDT,1);
    P_T20 = polyfit(X_T20,L_T20,1);
    P_T30 = polyfit(X_T30,L_T30,1);
    P_T60 = polyfit(X_T60,L_T60,1);
    
    LF_EDT = P_EDT(1)*vectorTiempo+P_EDT(2);
    LF_T20 = P_T20(1)*vectorTiempo+P_T20(2);
    LF_T30 = P_T30(1)*vectorTiempo+P_T30(2);
    LF_T60 = P_T60(1)*vectorTiempo+P_T60(2);
end