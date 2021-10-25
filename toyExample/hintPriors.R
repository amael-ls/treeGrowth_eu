## A hint on how the priors look like
# Define UI for sliders
library(shiny)
ui = fluidPage(
	# ---- App title
	titlePanel("Priors preview"),

	# ---- Sidebar layout with input and output definitions
	sidebarLayout(
		# ---- Sidebar to demonstrate various slider options
		sidebarPanel(
			# ---- Input: Mean associated to normal distribution plot
			sliderInput("mean_dnorm", "mean (normal distrib):",
				min = -50, max = 50,
				value = 0, step = 0.1),

			# ---- Input: Std deviation associated to normal distribution plot
			sliderInput("sd_dnorm", "standard deviation (normal distrib):",
				min = 0, max = 20,
				value = 10, step = 0.01),

			# ---- Input: Range to plot normal distribution
			sliderInput("range_dnorm", "Range normal:",
				min = -80, max = 80,
				value = c(-40, 40)),
			
			# ---- Input: Mean associated to normal distribution plot
			sliderInput("mean_dgamma", "mean (gamma distrib):",
				min = 0, max = 200,
				value = 5, step = 0.1),

			# ---- Input: Std deviation associated to normal distribution plot
			sliderInput("var_dgamma", "variance (gamma distrib):",
				min = 0, max = 200,
				value = 5, step = 0.1),

			# ---- Input: Range to plot gamma distribution
			sliderInput("range_dgamma", "Range gamma:",
				min = 0, max = 200,
				value = c(0, 20)),

			# ---- Input: Integration range
			sliderInput("integration_range", "a: P(X < a):",
				min = 0, max = 200,
				value = 1, step = 0.1)
		),

		# ---- Main panel for displaying outputs
		mainPanel(
			# ---- Output: Plot of the distributions and display shape/scale for gamma distrib
			tableOutput("values"),
			plotOutput(outputId = "normal_pdf"),
			plotOutput(outputId = "gamma_pdf")
		)
	)
)

## Define server logic for slider examples
server = function(input, output) {
	# --- Reactive expression to display the value of shape and rate for gamma distrib
	sliderValues = reactive({
	data.frame(
	Name = c("Gamma distribution, shape:", "Gamma distribution, rate:", "P(x < a), gamma distrib:"),
	Value = as.character(c(
		paste(input$mean_dgamma, "^2/", input$var_dgamma, sep = ""),
		paste(input$mean_dgamma, "/", input$var_dgamma, sep = ""),
		integrate(
			function(x){return (dgamma(x, shape = input$mean_dgamma^2/input$var_dgamma, rate = input$mean_dgamma/input$var_dgamma))},
			lower = 0, upper = input$integration_range)$value)),
	stringsAsFactors = FALSE)
	})

	# ---- Show the values in an HTML table
	output$values <- renderTable({sliderValues()})

	output$normal_pdf = renderPlot({curve(dnorm(x, input$mean_dnorm, input$sd_dnorm),
		input$range_dnorm[1], input$range_dnorm[2], lwd = 5, col = "#0072B5",
		ylab = "Normal pdf")})

	output$gamma_pdf = renderPlot({
		curve(dgamma(x,
			shape = input$mean_dgamma^2/input$var_dgamma, rate = input$mean_dgamma/input$var_dgamma),
			input$range_dgamma[1], input$range_dgamma[2], lwd = 5, col = "#363945",
			ylab = "Gamma pdf")

		sh <<- input$mean_dgamma^2/input$var_dgamma # https://stackoverflow.com/questions/2628621/how-do-you-use-scoping-assignment-in-r
		ra <<- input$mean_dgamma/input$var_dgamma

		DescTools::Shade(
			dgamma(x, shape = sh, rate = ra),
			breaks = c(0, input$integration_range), col = "#A1A1A166", density = NA)
		})
}

shinyApp(ui, server)