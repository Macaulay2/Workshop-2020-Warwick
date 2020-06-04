-- -*- coding: utf-8 -*-
newPackage(
	"SubalgebraBases",
	AuxiliaryFiles => true,
    	Version => "0.1", 
    	Date => "November 24, 2006",
    	Authors => {{Name => "Mike Stillman", 
		  Email => "mike@math.cornell.edu", 
		  HomePage => "http://www.math.cornell.edu/~mike/"}},
    	Headline => "Canonical subalgebra bases (aka SAGBI/Khovanskii bases)",
	AuxiliaryFiles => true, -- set to true if package comes with auxiliary files
--  	DebuggingMode => false,		
  	DebuggingMode => true		 -- set to true only during development
    	)

export {}
exportMutable {}

-*
-- old code!
--load "SubalgebraBases/sagbi-common.m2"
--load "SubalgebraBases/sagbitop.m2"
--load "SubalgebraBases/sagbieng.m2"
--load "SubalgebraBases/sagbi-tests.m2"
*-
needs "./SubalgebraBases/classes.m2"
needs "./SubalgebraBases/service-functions.m2"
needs "./SubalgebraBases/main.m2"

beginDocumentation()
needs "./SubalgebraBases/documentation.m2"

end--