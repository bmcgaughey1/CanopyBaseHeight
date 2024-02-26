# code to rotate XYZ data
#
# rotation order is z-y-x and each rotation is independent. This means that
# rotating around the y- and x-axes is not affected by rotation around the z-axis.
# For things like and azimuth and lean angle, you have to call the function twice.
rotateXYZ <- function(
    x,
    y,
    z,
    degx = 0.0,
    degy = 0.0,
    degz = 0.0,
    pivotx = 0.0,
    pivoty = 0.0,
    pivotz = 0.0,
    xlab = "x",
    ylab = "y",
    zlab = "z"
    )
{
  # if x is dataframe, ignore y and z
  if (is.data.frame(x)) {
    df <- x[, c(xlab, ylab, zlab)]
    colnames(df) <- c("x", "y", "z")
  }
  else {
    df <- data.frame("x" = x, "y" = y, "z" = z)
  }

  # convert angles to radians
  radx <- degx * pi / 180
  rady <- degy * pi / 180
  radz <- degz * pi / 180

  # do translation
  df$x <- df$x - pivotx
  df$y <- df$y - pivoty
  df$z <- df$z - pivotz

  dfi <- df
  colnames(dfi) <- c("xi", "yi", "zi")

  # do rotation around z
  dfi$xi <- df$x * cos(radz) - df$y * sin(radz)
  dfi$yi <- df$x * sin(radz) + df$y * cos(radz)
  dfi$zi <- df$z
  df <- dfi

  # do rotation around y
  dfi$xi <- df$x * cos(rady) - df$z * sin(rady)
  dfi$yi <- df$y
  dfi$zi <- df$x * sin(rady) + df$z * cos(rady)
  df <- dfi

  # do rotation around x
  dfi$xi <- df$x
  dfi$yi <- df$y * cos(radx) - df$z * sin(radx)
  dfi$zi <- df$y * sin(radx) + df$z * cos(radx)
  df <- dfi

  # do translation back to original location
  df$x <- df$x + pivotx
  df$y <- df$y + pivoty
  df$z <- df$z + pivotz

  invisible(return(df))
}
