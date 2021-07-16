

javaOptions += "-Xmx8G"
libraryDependencies ++= Seq(
  "org.locationtech"% "sfcurve" % "0.2.1-SNAPSHOT",
  uzaygezen,
  scalaTest % "test"
)
fork := true
