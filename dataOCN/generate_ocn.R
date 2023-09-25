library(OCNet)

dimX = 80
dimY = 100
nIter = 500*dimX*dimY
OCN <- create_OCN(dimX,dimY,initialNoCoolingPhase = 0.3,displayUpdates = 2,
                  nUpdates = 100,nIter = nIter,
                  cellsize = 10000)
OCN <- landscape_OCN(OCN,slope0 = 0.000667)
OCN <- aggregate_OCN(OCN,thrA = 55*10000*10000)
OCN <- rivergeometry_OCN(OCN,widthMax = 400, depthMax = 7.50)
draw_subcatchments_OCN(OCN)


cellsize=10000


library(R.matlab)



writeMat("OCN_C.mat", CTC = OCN$FD$toSC, X = OCN$FD$X, Y = OCN$FD$Y, 
         FD_A = OCN$FD$A, FD_X = OCN$FD$X, FD_Y = OCN$FD$Y, 
         FD_downNode = OCN$FD$downNode, FD_Z = OCN$FD$Z,
         SC = OCN$SC, SCX = OCN$SC$X, SCY = OCN$SC$Y,
         A = OCN$AG$A, downNode = OCN$AG$downNode,
         outlet = OCN$FD$outlet,
         R_length = OCN$AG$leng, R_width = OCN$AG$width, R_depth = OCN$AG$depth,
         cellsize = cellsize,
         FD_width = OCN$RN$width,FD_depth = OCN$RN$depth)


#save("OCN_C.rda")

