function [p, error_stroke,a,below,above] = brents(p,a,error_brents)
%function help here

tol = 1e-5;

if p.brent == 1
    p.x_stroke_brent(a) = p.x_stroke;
    if a >= 2   %a greater than 1
        
        if p.x_stroke_brent(a)< p.x_stroke_ref
            p.below = a;
        elseif p.x_stroke_brent(a) > p.x_stroke_ref
            p.above = a;
        end            
        
        %both previous pts are below desired stroke, secant only
        if p.x_stroke_brent(a) < p.x_stroke_ref && p.x_stroke_brent(a-1) < p.x_stroke_ref
            %Secant Step
            p.P_brent(a+1) = p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));
            p.method_n(a)=1;  %1 is secant, 2 is bisection
            
            
        %One below, one above (straddled), employ brents method    
        elseif (p.x_stroke_brent(a) < p.x_stroke_ref && p.x_stroke_brent(a-1) > p.x_stroke_ref) || (p.x_stroke_brent(a) > p.x_stroke_ref && p.x_stroke_brent(a-1) < p.x_stroke_ref)
            %If there is only one point below, then...
            if a==2
                %Bi-section Step
                p.P_brent(a+1) = (p.P_brent(a)+p.P_brent(a-1))/2;
                p.method_n(a)=2;
                
            else   %Brents method start
                %if the current iterate is above the reference
                if p.x_stroke_brent(a) >= p.x_stroke_ref
                    b_k = p.P_brent(a);
                    a_k = p.P_brent(p.below);
                    b_k_1 = p.P_brent(a-1);
                    %bisection step
                    m=(p.P_brent(a)+a_k)/2;
                    %Secant step
                    s=p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));

                    if s < b_k && s > m
                        p.P_brent(a+1)=s;
                        p.method_n(a)=1;
                        
                    else
                        p.P_brent(a+1)=m;
                        p.method_n(a)=2;
                        
                    end
                %if the current iterate is below the reference
                else
                    b_k = p.P_brent(a);
                    a_k = p.P_brent(p.above);
                    b_k_1 = p.P_brent(a-1);
                    %bisection step
                    m=(p.P_brent(a)+a_k)/2;
                    %Secant step
                    s=p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));

                    if s < b_k && s > m
                        p.P_brent(a+1)=s;
                        p.method_n(a)=1;
                        
                    else
                        p.P_brent(a+1)=m;
                        p.method_n(a)=2;
                        
                    end
                    

                    
                end
                
                %Brent's addition to the algorithmn
                if a>=3
                    if abs(p.P_brent(a)-p.P_brent(a-1)) < tol
                        if p.method_n(a-1) == 1  %used secant method
                            if abs(p.P_brent(a-1)-p.P_brent(a-2)) < tol
                                p.P_brent(a+1)=(p.P_brent(a)+a_k)/2; 
                                p.method_n(a)=2;
                                
                            end
                        else
                            p.P_brent(a+1)=(p.P_brent(a)+a_k)/2;
                            p.method_n(a)=2;
                            
                        end
                    elseif abs(s-b_k) > 0.5*abs(b_k-b_k_1)
                        if p.method_n(a-1) == 1  %used secant method
                            if abs(s-b_k) > 0.5*abs(b_k_1-p.P_brent(a-2))
                                p.P_brent(a+1)=(p.P_brent(a)+a_k)/2; 
                                p.method_n(a)=2;
                                
                            end
                        else
                            p.P_brent(a+1)=(p.P_brent(a)+a_k)/2;
                            p.method_n(a)=2;
                            
                        end
                    end
                end                

                
            end %brents method end

        
        %both above the desired stroke, if p.below does exist, employ brents method    
        elseif p.x_stroke_brent(a) > p.x_stroke_ref && p.x_stroke_brent(a-1) > p.x_stroke_ref
            if isfield(p,'below')==0
                p.P_brent(a+1)=p.P_brent(a-1)*0.4;
                
            else
                %if the current iterate is above the reference
                if p.x_stroke_brent(a) >= p.x_stroke_ref
                    b_k = p.P_brent(a);
                    a_k = p.P_brent(p.below);
                    b_k_1 = p.P_brent(a-1);
                    %bisection step
                    m=(p.P_brent(a)+a_k)/2;
                    %Secant step
                    s=p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));

                    if s < b_k && s > m
                        p.P_brent(a+1)=s;
                        p.method_n(a)=1;
                        
                    else
                        p.P_brent(a+1)=m;
                        p.method_n(a)=2;
                        
                    end
                %if the current iterate is below the reference
                else
                    b_k = p.P_brent(a);
                    a_k = p.P_brent(p.above);
                    b_k_1 = p.P_brent(a-1);
                    %bisection step
                    m=(p.P_brent(a)+a_k)/2;
                    %Secant step
                    s=p.P_brent(a) - (p.x_stroke_brent(a) - p.x_stroke_ref)*(p.P_brent(a) - p.P_brent(a-1))/(p.x_stroke_brent(a) - p.x_stroke_brent(a-1));

                    if s < b_k && s > m
                        p.P_brent(a+1)=s;
                        p.method_n(a)=1;
                        
                    else
                        p.P_brent(a+1)=m;
                        p.method_n(a)=2;
                        
                    end
                    
                end
                
                %Brent's addition to the algorithmn
                
                if a>=3
                    if abs(p.P_brent(a)-p.P_brent(a-1)) < tol
                        if p.method_n(a-1) == 1  %used secant method
                            if abs(p.P_brent(a-1)-p.P_brent(a-2)) < tol
                                p.P_brent(a+1)=(p.P_brent(a)+a_k)/2; 
                                p.method_n(a)=2;
                                
                            end
                        else
                            p.P_brent(a+1)=(p.P_brent(a)+a_k)/2;
                            p.method_n(a)=2;
                            
                        end
                    elseif abs(s-b_k) > 0.5*abs(b_k-b_k_1)
                        if p.method_n(a-1) == 1  %used secant method
                            if abs(s-b_k) > 0.5*abs(b_k_1-p.P_brent(a-2))
                                p.P_brent(a+1)=(p.P_brent(a)+a_k)/2;
                                p.method_n(a)=2;
                                
                            end
                        else
                            p.P_brent(a+1)=(p.P_brent(a)+a_k)/2;
                            p.method_n(a)=2;
                            
                        end
                    end
                end
                
            end %brents method end

        end
        
        error_stroke = abs(p.x_stroke_ref - p.x_stroke_brent(a));
        
        a=a+1;
        % 
        
    else  %a == 1
  
        %p.P_electric = p.P_guess_2;
        error_stroke = 1;
        a = a+1;
        %p.above = 0;
        %p.below = 0;
        p.method_n=0;
        %
    end  % end entry logic
    
else  %if Brent's method is turned off
    
    error_stroke = 0;
    
end




end


  