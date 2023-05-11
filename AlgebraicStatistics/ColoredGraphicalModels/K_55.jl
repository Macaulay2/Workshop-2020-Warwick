using LearningAlgebraicVarieties
using DelimitedFiles
using CSV
using DataFrames
using StatsPlots

SymL5=CSV.File("/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/SymL5.csv";header=false)

SymL5=SymL5|>DataFrame|>Matrix

Sym5Diag=DimensionDiagrams(SymL5, true)
savefig(Sym5Diag,"/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/Sym5Diag.png")

PSDL5=CSV.File("/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/PSDL5.csv";header=false)

PSDL5=PSDL5|>DataFrame|>Matrix
PSD5Diag=DimensionDiagrams(PSDL5, true)
savefig(PSD5Diag,"/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/PSD5Diag.png")

SymL4=CSV.File("/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/SymL4.csv";header=false)

SymL4=SymL4|>DataFrame|>Matrix
Sym4Diag=DimensionDiagrams(SymL4, true)
savefig(Sym4Diag,"/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/Sym4Diag.png")

PSDL4=CSV.File("/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/PSDL4.csv";header=false)

PSDL4=PSDL4|>DataFrame|>Matrix
PSD4Diag=DimensionDiagrams(PSDL4, true)
savefig(PSD4Diag,"/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/PSD4Diag.png")

SymL3=CSV.File("/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/SymL3.csv";header=false)


SymL3=SymL3|>DataFrame|>Matrix
Sym3Diag=DimensionDiagrams(SymL3, true)
savefig(Sym3Diag,"/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/Sym3Diag.png")

PSDL3=CSV.File("/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/PSDL3.csv";header=false)

PSDL3=PSDL3|>DataFrame|>Matrix
PSD3Diag=DimensionDiagrams(PSDL3, true)
savefig(PSD3Diag,"/Users/k87125/Dropbox/Math \
conferences/Workshop-2020-Warwick/AlgebraicStatistics/\
ColoredGraphicalModels/PSD3Diag.png")

FindEquations(PSDL4, :with_svd, 5, true)