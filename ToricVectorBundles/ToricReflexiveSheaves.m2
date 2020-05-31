-- -*- coding: utf-8 -*-
------------------------------------------------------------------------------
-- Copyright 2017-20 Gregory G. Smith
--
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
--
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
-- FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
-- more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------------
newPackage(
  "ToricReflexiveSheaves",
  AuxiliaryFiles => true,
  Version => "0.3",
  Date => "31 May 2020",
  Authors => {{
      Name => "Gregory G. Smith", 
      Email => "ggsmith@mast.queensu.ca", 
      HomePage => "http://www.mast.queensu.ca/~ggsmith"}},
  Headline => "routines for working with toric reflexive sheaves",
  PackageExports => {"NormalToricVarieties", "HyperplaneArrangements"},
  PackageImports => {},
  DebuggingMode => true
  )


export {
    -- classes
    "ToricReflexiveSheaf",
    "ToricReflexiveSheafMap",
    -- methods   
    "toricReflexiveSheaf",
    "subspace",     
    "mukaiLazarsfeldBundle",
    "groundSet",
    "isArithmeticallyFree", 
    "isLocallyFree", 
    "associatedCharacters",
    "isGloballyGenerated",
    "subspacePoset",  -- is this necessary
    "toricTangentBundle",
    "toricCotangentBundle"
  }
------------------------------------------------------------------------------
-- CODE
------------------------------------------------------------------------------
load "ToricReflexiveSheaves/Code.m2"


------------------------------------------------------------------------------
-- DOCUMENTATION
------------------------------------------------------------------------------
beginDocumentation()
load "ToricReflexiveSheaves/Documentation.m2"


------------------------------------------------------------------------------
-- TESTS
------------------------------------------------------------------------------
load "ToricReflexiveSheaves/Tests.m2"

end---------------------------------------------------------------------------


------------------------------------------------------------------------------
-- methods to be added to package

dual ToricReflexiveSheaf := ToricReflexiveSheaf => E -> (  
ToricReflexiveSheaf ** ToricReflexiveSheaf := ToricReflexiveSheaf => (E,F) -> (
sheafHom (ToricReflexiveSheaf, ToricReflexiveSheaf) := ToricReflexiveSheaf => (E,F) -> (
euler ToricReflexiveSheaf := ZZ => E -> (
toricEuler ToricReflexiveSheaf := RingElement => E -> (

isAmple ToricReflexiveSheaf := Boolean => E -> (
isVeryAmple ToricReflexiveSheaf := Boolean => E -> (
cohomology ToricReflexiveSheaf := 

cokernel ToricReflexiveSheafMap
image ToricReflexiveSheafMap
coimage ToricReflexiveSheafMap
ToricReflexiveSheaf _ Array := ToricReflexiveSheafMap => (E,v) -> (
ToricReflexiveSheaf ^ Array := ToricReflexiveSheafMap => (E,v) -> (  
  
chowRing, chern classes, etc.

------------------------------------------------------------------------------

-- XXX
uninstallPackage "ToricReflexiveSheaves"
restart
installPackage "ToricReflexiveSheaves"
check ToricReflexiveSheaves

