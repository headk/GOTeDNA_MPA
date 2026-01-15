# GOTeDNA_V2
Most recent version of the Shiny app with modules

## Updates:

Currently, running the app takes a while and you cannot see the code it is running through. 

To ensure the code is running without error, run each file in the following order:
1) R/packages.R
2) R/dev_pipeline.R
3) ui.R
4) server.R

To run app:
shinyApp(ui, server)

The next steps are to precompute functions so that grids and calculations are done prior to loading the app (quicker running time), then fix the visual problems of the UI.R and www/style.css confusion (visual overlapping of menus), followed by creating diversity visuals and linking to OBIS data (once stored in OBIS).
