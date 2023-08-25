library(OCNet)
OCN <- create_OCN(80,100,initialNoCoolingPhase = 0.4,displayUpdates = 2,nUpdates = 100,nIter = 30000*50)
OCN <- landscape_OCN(OCN)
OCN <- aggregate_OCN(OCN,thrA = 55)
OCN <- rivergeometry_OCN(OCN)
draw_subcatchments_OCN(OCN)




library(R.matlab)

writeMat("subcatchments.mat", CTC = OCN$FD$toSC, X = OCN$FD$X, Y = OCN$FD$Y)



writeMat("FD.mat", FD.A = OCN$FD$A, FD.X = OCN$FD$X, FD.Y = OCN$FD$Y, FD.downNode = OCN$FD$downNode, FD.Z = OCN$FD$Z)
writeMat("SC.mat", SC = OCN$SC, SCX = OCN$SC$X, SCY = OCN$SC$Y)

writeMat("AG.mat", A = OCN$AG$A, downNode = OCN$AG$downNode)

writeMat("outlet.mat", outlet = OCN$FD$outlet)

draw_contour_OCN(OCN)

