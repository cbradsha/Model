function [p, error_stroke,a,below,above] = brents(p,a,error_brents)
%function help here

if p.brent == 1
    p.x_stroke_brent(a) = p.x_stroke;
    if a >= 2
        
        if p.x_stroke_brent(a)< p.x_max
            p.below = a;
        elseif p.x_stroke_brent(a) > p.x_max
            p.above = a;
        end            
        
        if p.x_stroke_brent(a) < p.x_max && p.x_stroke_brent(a-1) < p.x_max
            %Secant Step
            p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_max)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
            
        elseif p.x_stroke_brent(a) < p.x_max || p.x_stroke_brent(a-1) < p.x_max
            %Bi-section Step
            p.P_brent(a+1) = (p.P_brent(a)+p.P_brent(a-1))/2;
            
        elseif p.x_stroke_brent(a) > p.x_max && p.x_stroke_brent(a-1) > p.x_max
            if abs(p.x_stroke_brent(a) - p.x_stroke_brent(a-1)) > 500*error_brents
                %Secant Step
                p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_max)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));  
                
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
        %p.x_max)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
        error_stroke = abs(p.x_max - p.x_stroke_brent(a));

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


        %% old brents method
%         %For Craigs method, solutions must be below desired solution
%         if (p.x_stroke_brent(a-1)<p.x_max && p.x_stroke_brent(a)<p.x_max)
%             
%             if a == 2
%             
%                 if p.x_stroke_brent(a) > 0.5*p.x_max
%                     ln_P_brent = log(p.P_brent);
%                     coeff = polyfit(ln_P_brent,p.x_stroke_brent,1);            %fit linear fit
%                     p.P_brent(a+1) = exp((p.x_max - coeff(2))/coeff(1));
%                     p.P_brent(a+1) = 0.95*p.P_brent(a+1);
%                 else
%                     %if the stroke from the first two steps are very small,
%                     %just take a secant step
%                     p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_max)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
%                     
%                 end
%             
%             elseif a >= 3
%                 
%                 diff_stroke = abs(p.x_stroke_brent(a-1) - p.x_stroke_brent(a))/p.x_stroke_brent(a);
%                 
%                 if p.x_stroke_brent(a) > 0.8*p.x_max
%                     ln_P_brent = log(p.P_brent);
%                     coeff = polyfit(ln_P_brent,p.x_stroke_brent,1);            %fit linear fit
%                     p.P_brent(a+1) = exp((p.x_max - coeff(2))/coeff(1));
%                     
%                 else
%                     %p.P_brent(a+1) = 1.3*p.P_brent(a);
%                     
%                     p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_max)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
%                     
%                 end
%                 
%                 %error_stroke = 0;           
%                 
% %             elseif p.x_stroke_brent(a) >= p.x_max
% %                 
% %                 p.P_brent(a+1) = p.P_brent(a)*0.95;                        %reduce by 5%
% %                 
% %             elseif p.x_stroke_brent(a) < p.x_max
% %                 
% %                 if p.x_stroke_brent(a) > p.x_stroke_brent(a-1)
% %                     
% %                     ln_P_brent = ln(p.P_brent);
% %                     coeff = polyfit(ln_P_brent,p.x_stroke_brent,1);            %fit linear fit
% %                     p.P_brent(a+1) = exp((p.x_max - coeff(2))/coeff(1));
% %                     
% %                 elseif 
% 
%             end
%            
%             
%             error_stroke = abs(p.x_max - p.x_stroke_brent(a))/p.x_max;
%             a=a+1; 
%             
%         %if solutions are not both below desired solution   
%         else
%             
%             p.P_brent(a-1) = p.P_brent(a-1)*0.75;
%             p.P_brent(a) = p.P_brent(a)*0.75;
%             error_stroke = 1;
%             a=1;
%                         
%         end    %end straddle logic
%         %keyboard
