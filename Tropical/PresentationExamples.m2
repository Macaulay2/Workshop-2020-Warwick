--Code for the presentation of the tropical group at the Warwick workshop
--4 June 2020

needsPackage "Tropical"




--isBalanced 

R = QQ[x, y];
I = ideal (x+y+1);
T = tropicalVariety(I);
isBalancedCurves T == true

U = tropicalCycle(fan T, {1, 2, 3});
isBalancedCurves U
                                                                      
								      