function [p, error_stroke,a,below,above] = brents(p,a,error_brents)
%function help here

if p.brent == 1
    p.x_stroke_brent(a) = p.x_stroke;
    if a >= 2   %a greater than 1
        
        if p.x_stroke_brent(a)< p.x_stroke_ref
            p.below = a;
        elseif p.x_stroke_brent(a) > p.x_stroke_ref
            p.above = a;
        end            
        
        if p.x_stroke_brent(a) < p.x_stroke_ref && p.x_stroke_brent(a-1) < p.x_stroke_ref
            %Secant Step
            p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
            
        elseif p.x_stroke_brent(a) < p.x_stroke_ref || p.x_stroke_brent(a-1) < p.x_stroke_ref
            %Bi-section Step
            p.P_brent(a+1) = (p.P_brent(a)+p.P_brent(a-1))/2;
            
        elseif p.x_stroke_brent(a) > p.x_stroke_ref && p.x_stroke_brent(a-1) > p.x_stroke_ref
            if abs(p.x_stroke_brent(a) - p.x_stroke_brent(a-1)) > 500*error_brents
                %Secant Step
                p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));  
                
            else
                %Bi-section Step using 'below'
                if p.below == 0
                    p.P_brent(a+1) = 0.5*p.P_brent(1);
                else
                    p.P_brent(a+1) = (p.P_brent(a)+p.P_brent(p.below))/2;
                end
            end
        end
        
        %p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) -
        %p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
        error_stroke = abs(p.x_stroke_ref - p.x_stroke_brent(a));
        
        a=a+1;
         
    else  %a == 1
  
        %p.P_electric = p.P_guess_2;
        error_stroke = 1;
        a = a+1;
        p.above = 0;
        p.below = 0;
        %keyboard
    end  % end entry logic
    
else  %if Brent's method is turned off
    
    error_stroke = 0;
    
end




end


  