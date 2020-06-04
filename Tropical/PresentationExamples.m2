--Code for the presentation of the tropical group at the Warwick workshop
--4 June 2020





--isBalanced 

R = QQ[x, y];
I = ideal (x+y+1);
T = tropicalVariety(I);
assert (isBalancedCurves T == true)                                                     

U = tropicalCycle(fan T, {1, 2, 3});                                                    
assert (isBalancedCurves U == false)                                                    
                                                                      
								      