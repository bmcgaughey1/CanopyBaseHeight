# code to demonstrate canopy base height methodology
#
library(lidR)
library(fusionwrapr)  # https://github.com/bmcgaughey1/fusionwrapr
library(rgl)

# Convert from LDA to LAS...only once -------------------------------------
#
# The resulting LAS files are not perfect. The number of returns for a pulse will be
# wrong (always 1) and some other fields in the LAS point records will not have valid values.
# This happens because FUSION's LDA format, developed before LAS format, does not include all
# of the information needed to create a "complete" LAS file.
#
# the lidR package will throw some warning when reading the converted files but I don't think
# the issues cause any problems.
#
# read point cloud for tree...manually change the file names for read and write
#treeLAS <- readLDA("extras/trees_clip_0000183.lda", type = "LAS", LASTemplate = "H:/T3_DroneLidar/Ba/Plot37/Ba_plot_37_001.laz")

# write LAS file
#writeLAS(treeLAS, "extras/tree_0000183.las")

# function to implement cone clipping
# basic idea is to specify the XYZ for the cone point and then the slope of the sides (or angle of the cone).
# putting the point below the ground simulates a conic frustum.
#
# below works with baseDia...if you have a baseDia and below=FALSE, points below the base of the frustum are removed
coneClip <- function(
    pts,
    coneX,
    coneY,
    coneZ,
    angle = 30,
    baseDia = 0.0,
    inside = TRUE,
    below = TRUE,
    xLabel = "x",
    yLabel = "y",
    zLabel = "z"
) {
  ipts <- pts[, c(xLabel, yLabel, zLabel)]
  colnames(ipts) <- c("x", "y", "z")

  cat("incoming: ", coneX, coneY, coneZ, "\n")

  # precompute some things
  tanAngle <- tan(angle / 2 * pi / 180)
  apexZ <- coneZ

  # modify coneZ to produce baseDia at the original coneZ
  if (baseDia > 0.0)
    apexZ <- coneZ - (baseDia / 2 / tanAngle)

  cat("modified: ", coneX, coneY, apexZ, "\n")

  # compute horizontal and vertical distance from apex XY to each point
  ipts$hdist <- sqrt((ipts$x - coneX)^2 + (ipts$y - coneY)^2)
  ipts$vdist <- ipts$z - apexZ

  # do the clip using the full cone
  if (inside) {
    cpts <- ipts[which(ipts$hdist <= (tanAngle * ipts$vdist)), ]
  } else {
    cpts <- ipts[which(ipts$hdist > (tanAngle * ipts$vdist)), ]
  }

  # do clip for points below the cone or frustum
  if (!below & baseDia > 0.0)
    cpts <- cpts[which(cpts$z >= coneZ), ]

  return(cpts)
}

computeCbh <- function(
    z,
    baseZ,
    showPlot = TRUE
    )
{
  # normalize...for tree 23, 272.602333 is the ground elevation under the high point XY
  #nz <- treeLAS$Z - minZ
  nz <- z - baseZ

  # compute percentiles
  p <- quantile(nz, probs = seq(0.0, 1, 0.01), names = FALSE)

  # normalize percentiles using P99 (stored in element 100)
  np <- p / p[100]

  # plot
 if (showPlot)
   plot(seq(1, 100, 1), np[-1], "l", xlab = "Percentile", ylab = "Normalized height")

  x <- c(1:99)
  slopes <- vector()
  # compute slopes
  for (j in 2:98) {
    x1 <- x[j - 1]
    x2 <- x[j]
    #x2 <- x[j + 1]    # compute slope, skipping center point
    y1 <- np[j - 1]
    y2 <- np[j]
    #y2 <- np[j + 1]    # compute slope, skipping center point
    slope_i <- (y2 - y1) / (x2 - x1)
    slopes <- append(slopes, slope_i)
  }

  # get max slope
  maxSlopeIndex <- which.max(slopes)
  cbhIndex <- maxSlopeIndex # + 1
  cbh <- p[cbhIndex + 1]

  if (showPlot) {
    abline(v = cbhIndex + 1, col = "red")
    abline(v = cbhIndex, col = "green")
  }

  return(cbh)
}

# Read tree ---------------------------------------------------------------
# throws a warning about the 'return number' greater than the 'number of returns'
treeLAS <- readLAS("extras/tree_0000121.las")

# get minimum Z to normalize...this could include ground points below tree base
# when tree is on steep slope
#minZ <- min(treeLAS$Z)

# normalize...for tree 23, 272.602333 is the ground elevation under the high point XY
#nz <- treeLAS$Z - minZ
nz <- treeLAS$Z - 221.4572

# compute percentiles
p <- quantile(nz, probs = seq(0.0, 1, 0.01), names = FALSE)

# normalize percentiles using P99 (stored in element 100)
np <- p / p[100]

# plot
plot(seq(0, 100, 1), np, "l", xlab = "Percentile", ylab = "Normalized height")

x <- c(1:99)
slopes <- vector()
# compute slopes
for (j in 2:99) {
  x1 <- x[j - 1]
  x2 <- x[j]
  y1 <- np[j -1]
  y2 <- np[j]
  slope_i <- (y2-y1)/(x2-x1)
  slopes <- append(slopes, slope_i)
}

# get max slope
maxSlopeIndex <- which.max(slopes[-1])
cbhIndex <- maxSlopeIndex + 1
cbh <- p[cbhIndex + 1]




# code for testing cone clip function -------------------------------------

# read high point XY for tree 23
treeNum <- 121

t <- read.csv("extras/TAO_metrics_highpts.csv", stringsAsFactors = FALSE)
highx <- t[t$Identifier == treeNum, 'High.point.X']
highy <- t[t$Identifier == treeNum, 'High.point.Y']
highz <- t[t$Identifier == treeNum, 'High.point.elevation']
df <- data.frame("X" = highx, "Y" = highy)
write.csv(df, "extras/t.csv", row.names = FALSE)

# get ground elevation under high point
SurfaceSample("extras/ground.dtm", "extras/t.csv", "extras/t_gnd.csv")

# read output file to get ground elevation
t <- read.csv("extras/t_gnd.csv", stringsAsFactors = FALSE)
groundz <- t[1, 'Value']

# get XYZ from LAS data
pts <- data.frame("x" = treeLAS@data$X,
                  "y" = treeLAS@data$Y,
                  "z" = treeLAS@data$Z)

cbh <- computeCbh(pts$z, groundz)

angle <- 10
baseDia <- 2
groundOffset <- 4

ptsClip <- coneClip(pts,
                    highx,
                    highy,
                    groundz + groundOffset,
                    angle = angle,
                    baseDia = baseDia,
                    below = FALSE)

cbh <- computeCbh(ptsClip$z, groundz)

#plot3d(treeLAS@data$X - highx, treeLAS@data$Y - highy, treeLAS@data$Z - groundz, size = 2, col = "cyan", aspect = c(1, 1, 10), decorate = FALSE)
plot3d(ptsClip$x - highx, ptsClip$y - highy, ptsClip$z - groundz, add = TRUE, size = 4, col = "black", aspect = c(1,1,10))


# rotate 45 degrees around x
source("Rcode/rotateXYZ.R")
dfr <- rotateXYZ(ptsClip$x, ptsClip$y, ptsClip$z, 0, 45, 45, highx, highy, groundz)
#dfr <- rotateXYZ(dfr$x, dfr$y, dfr$z, 0, 45, 0, highx, highy, groundz)


plot3d(dfr$x - highx, dfr$y - highy, dfr$z - groundz, add = TRUE, size = 4, col = "orange")

# display the cone with the points...need to compute coordinates
baseRadius <- baseDia / 2
topRadius <- (tan(angle / 2 * pi / 180) * (highz - (groundz + groundOffset) + (baseRadius / tan(angle / 2 * pi / 180))))
cpts <- data.frame("angle" = seq(0, 360, 10), "bX"=0, "bY"=0, "bZ"=0, "tX"=0, "tY"=0, "tZ"=0)
cpts$bX <- (baseRadius * cos(cpts$angle * pi / 180))
cpts$bY <- (baseRadius * sin(cpts$angle * pi / 180))
cpts$bZ <- groundOffset
cpts$tX <- topRadius * cos(cpts$angle * pi / 180)
cpts$tY <- topRadius * sin(cpts$angle * pi / 180)
cpts$tZ <- highz - groundz

# create circle pts
quads <- matrix(nrow = nrow(cpts) - 1, ncol = 4)
coords <- matrix(nrow = nrow(cpts) * 2, ncol = 3)
for (i in 1:nrow(cpts)) {
  coords[i, 1] <- cpts$bX[i]
  coords[i, 2] <- cpts$bY[i]
  coords[i, 3] <- cpts$bZ[i]
  coords[i + nrow(cpts), 1] <- cpts$tX[i]
  coords[i + nrow(cpts), 2] <- cpts$tY[i]
  coords[i + nrow(cpts), 3] <- cpts$tZ[i]
  if (i < nrow(cpts)) {
    quads[i, 1] <- i
    quads[i, 2] <- i + 1
    quads[i, 3] <- i + nrow(cpts) + 1
    quads[i, 4] <- i + nrow(cpts)
  }
}
wire3d(mesh3d(x = coords, quads = t(quads)), col = "red")

# disk at cbh
cpts$tZ <- cbh

tris <- matrix(nrow = nrow(cpts) - 1, ncol = 3)
coords <- matrix(nrow = nrow(cpts) + 1, ncol = 3)
coords[1, 1] <- 0
coords[1, 2] <- 0
coords[1, 3] <- cbh
for (i in 2:nrow(cpts)) {
  coords[i, 1] <- cpts$tX[i]
  coords[i, 2] <- cpts$tY[i]
  coords[i, 3] <- cpts$tZ[i]
#  if (i < nrow(cpts)) {
    tris[i - 1, 1] <- 1
    tris[i - 1, 2] <- i
    tris[i - 1, 3] <- i + 1
#  }
}
tris[nrow(cpts) - 1, 3] <- 2
wire3d(mesh3d(x = coords, triangles = t(tris)), col = "blue")
