library(OCNet)
OCN <- create_OCN(300,200)
OCN <- landscape_OCN(OCN)
OCN <- aggregate_OCN(OCN)
OCN <- rivergeometry_OCN(OCN)

library(R.matlab)

writeMat("subcatchments.mat", CTC = OCN$FD$toSC, X = OCN$FD$X, Y = OCN$FD$Y)

draw_subcatchments_OCN(OCN)

writeMat("FD.mat", FD.A = OCN$FD$A, FD.X = OCN$FD$X, FD.Y = OCN$FD$Y, FD.downNode = OCN$FD$downNode)
writeMat("SC.mat", SC = OCN$SC, SCX = OCN$SC$X, SCY = OCN$SC$Y)

writeMat("AG.mat", A = OCN$AG$A)

draw_contour_OCN(OCN)

