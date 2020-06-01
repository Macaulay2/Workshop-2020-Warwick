-- workflow for installing this version of the package
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