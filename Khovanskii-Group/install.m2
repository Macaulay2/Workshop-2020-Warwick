-- workflow for installing this version of the package
-- run this as load "install.m2" from within a running instance of M2
uninstallPackage "SubalgebraBases"
restart
-- !!! assumes your M2 session starts inside of ".../-2020-Warwick/Khovanskii-Group/"
pathToPackage = "./SubalgebraBases.m2"
installPackage(
    "SubalgebraBases",
    FileName=>pathToPackage,
    RerunExamples => true
    )
check "SubalgebraBases"
help SubalgebraBases
-- viewHelp SubalgebraBases -- help in browser
