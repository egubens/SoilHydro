#' @export
SoilHydro_App = function(){
  appDir = system.file("SHApp", package = "SoilHydro")
  shiny::runApp(appDir, display.mode = "normal")
}
