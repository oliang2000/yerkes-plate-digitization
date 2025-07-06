library(shiny)
library(reticulate)
library(png)

# <--- UI Definition --->
ui <- fluidPage(
  titlePanel("SDSS Quasar Explorer"),
  p("This application allows you to identify quasars 
  from the Sloan Digital Sky Survey (SDSS) within selected view. 
  Use the sliders below to define a sky region by 
  Right Ascension (RA) and Declination (Dec), 
  the maximum range is constrained to 30 to ensure query speed. 
  The app will then query the SDSS database, 
  display the SQL query being executed, 
  and plot the selected QSOs in an RA vs. Dec scatter plot. 
  You can also download the queried data."),
  tags$hr(), # horizontal rule for visual separation
  sidebarLayout(
    sidebarPanel(
      sliderInput("ra", "Right Ascension (deg)", min = 0, max = 360,
                  value = c(257.5, 260.0), step = 0.01),
      sliderInput("dec", "Declination (deg)", min = -90, max = 90,
                  value = c(42.75, 44.75), step = 0.01),
      br(),

      tags$h4("SQL Query"),
      verbatimTextOutput("sqltext"),
      br(),
      # download button for the queried data
      downloadButton("downloadData", "Download Data (CSV)")
    ),
    mainPanel(
      tags$h4("Location of quasars within the field"),
      plotOutput("scatter")
    )
  )
)

# <--- Server Logic --->
server <- function(input, output, session) {

  # initialize Python
  py_run_string("
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
import pandas as pd
")

  # max range for ra and dec
  max_range <- 30

  # Observe slider changes to enforce max range
  observeEvent(input$ra, {
    current_ra_range <- input$ra[[2]] - input$ra[[1]]
    if (current_ra_range > max_range) {
      new_ra_max <- input$ra[[1]] + max_range
      updateSliderInput(session, "ra", value = c(input$ra[[1]], new_ra_max))
    }
  }, ignoreNULL = TRUE)
  observeEvent(input$dec, {
    current_dec_range <- input$dec[[2]] - input$dec[[1]]
    if (current_dec_range > max_range) {
      new_dec_max <- input$dec[[1]] + max_range
      updateSliderInput(session, "dec", value = c(input$dec[[1]], new_dec_max))
    }
  }, ignoreNULL = TRUE)

  # Dynamically generate SQL query string
  sql_query_string <- reactive({
    ra1_fmt <- format(input$ra[[1]], scientific = FALSE, trim = TRUE)
    ra2_fmt <- format(input$ra[[2]], scientific = FALSE, trim = TRUE)
    dec1_fmt <- format(input$dec[[1]], scientific = FALSE, trim = TRUE)
    dec2_fmt <- format(input$dec[[2]], scientific = FALSE, trim = TRUE)

    paste0(
      "SELECT p.ra, p.dec, s.z, snMedian, \n",
      "p.u, p.g, p.r, p.i, p.z, \n",
      "p.modelmagerr_u, p.modelmagerr_g,\n",
      "p.modelmagerr_r, p.modelmagerr_i, p.modelmagerr_z\n",
      "FROM PhotoObj AS p, SpecObj AS s\n",
      "WHERE p.specObjID = s.specObjID\n",
      "  AND snMedian > 7.0\n",
      "  AND s.class = 'QSO'\n",
      "  AND (p.ra BETWEEN ", ra1_fmt, " AND ", ra2_fmt, ")\n",
      "  AND (p.dec BETWEEN ", dec1_fmt, " AND ", dec2_fmt, ")"
    )
  })

  # store fetched data
  qso_data <- reactiveVal(NULL)
  # changes in slider inputs and fetch data
  observeEvent(c(input$ra, input$dec), {
    req(input$ra, input$dec)
    query_to_execute <- sql_query_string()
    py$query_py <- query_to_execute
    # Python api call to fetch data
    py_run_string("
try:
    results_py = SDSS.query_sql(query_py)
    if results_py is not None:
        df_py = results_py.to_pandas()
    else:
        df_py = pd.DataFrame({'ra': [], 'dec': []})
except Exception as e:
    print(f'Error executing SDSS query for data download: {e}')
    df_py = pd.DataFrame({'ra': [], 'dec': []})
")
    # store the fetched data
    qso_data(reticulate::py$df_py)
  })

  # render SQL query text in the UI
  output$sqltext <- renderText({
    sql_query_string()
  })

  # Call Python to fetch data, plot it, save it, and then display it in R
  output$scatter <- renderPlot({
    req(input$ra, input$dec)
    query_to_execute <- sql_query_string()
    # pass the query string and plot parameters to Python's global scope
    py$query_py <- query_to_execute
    py$ra_min_py <- input$ra[[1]]
    py$ra_max_py <- input$ra[[2]]
    py$dec_min_py <- input$dec[[1]]
    py$dec_max_py <- input$dec[[2]]

    # temporary file path for the plot image
    plot_filepath <- tempfile(pattern = "sdss_qso_scatter", fileext = ".png")
    py$plot_filepath_py <- plot_filepath # Pass this path to Python
    # remove temp plot file
    on.exit({
      if (file.exists(plot_filepath)) {
        file.remove(plot_filepath)
      }
    })

    # <--- Python logic for querying data and generating the plot --->
    py_run_string("
try:
    results = SDSS.query_sql(query_py)
    if results is not None:
        df = results.to_pandas()
    else:
        print('SDSS query returned no results (None).')
        df = pd.DataFrame({'ra': [], 'dec': []})

except Exception as e:
    # catch any errors during the query execution
    print(f'Error executing SDSS query: {e}')
    df = pd.DataFrame({'ra': [], 'dec': []})

# plot
plt.figure(figsize=(15, 5))

if not df.empty:
    plt.scatter(df['ra'], df['dec'], s=50, c='blue', marker='*', alpha=0.8)
else:
    plt.text(0.5, 0.5, 'No data to display or query failed.', 
    horizontalalignment='center',
    verticalalignment='center', transform=plt.gca().transAxes, 
    fontsize=14, color='gray')

plt.xlabel('Right Ascension (deg)')
plt.ylabel('Declination (deg)')
plt.title('Quasar locations from SDSS')

# plot limits based on slider inputs
plt.xlim(ra_min_py, ra_max_py)
plt.ylim(dec_min_py, dec_max_py)
# correct aspect ratio
plt.gca().set_aspect('equal')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(plot_filepath_py, dpi = 300)
plt.close()
")

    # <--- Display the plot in R --->
    # if plot file was successfully created by Python
    if (file.exists(plot_filepath)) {
      img <- png::readPNG(plot_filepath)
      par(mar = c(0, 0, 0, 0))
      plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
      rasterImage(img, 0, 0, 1, 1)
    } else {
      plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes=FALSE)
      text(0.5, 0.5, "Plot could not be generated.", col = "red", cex = 1.2)
    }
  },
  res = 96, #resoultion
  )

  # Server logic for the download button
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sdss_qso_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(qso_data())) {
        write.csv(qso_data(), file, row.names = FALSE)
      } else {
        showNotification("No data to download. 
        Please adjust sliders and allow time for query.", type = "warning")
      }
    }
  )
}


shinyApp(ui, server)