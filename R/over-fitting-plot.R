library(ggplot2)
library(ggforce)
library(bezier)
library(tibble)
library(dplyr)

# save plot to tikz
setwd("~/Desktop/model-selection-uncertainty")
source("./R/aux/aux_plotting.R")

t <- seq(0, 1, length=100)
p_theory <- matrix(
  c(0,-100, 
    5, 87,
    10, -25,
    15, -80,
    20, 0), 
  nrow=5, 
  ncol=2, 
  byrow=TRUE
)
p_emp <- matrix(
  c(0,-90, 
    5, 110,
    10, 50,
    15, 30,
    20, 0),
  nrow=5, 
  ncol=2, 
  byrow=TRUE
)

bezier_points_theory <- bezier(t=t, p=p_theory)
bezier_points_emp <- bezier(t=t, p=p_emp)

colnames( bezier_points_theory ) <- c( "size", "elpd_out" )
colnames( bezier_points_emp ) <- c( "size", "elpd_in" )

data <- as_tibble( bezier_points_theory ) %>%
  add_column(bezier_points_emp[, "elpd_in"])
colnames( data ) <- c( "size", "elpd_out", "elpd_in" )

label_data <- data.frame(
  x = c(18, 3, 2, 10, 9),
  y = c(30, -50, 25, -35, 12),
  label = c("LOO-CV", "test", "selection-induced\n bias", "over-fitting", 
            "point of predictive\nsaturation"),
  colour = c("red", "grey", "black", "black", "black")
)
line_data_1 <- data |>
  filter(size >= 3.83 & size < 4) |>
  mutate(x = size, xend = size, y = elpd_out, yend = elpd_in)
line_data_2 <- data |>
  filter(size > 11 & size < 11.2) |>
  mutate(x = size, xend = size, y = elpd_out, yend = 0)

( p <- ggplot(data, aes(x=size, y=elpd_out)) + 
  geom_line(aes(y = elpd_out), colour = "grey") + 
  geom_line(aes(y = elpd_in), colour = "red") +
  geom_hline(aes(yintercept = 0), size = 0.2) +
  geom_segment(data = line_data_1, aes(x = x, y = y, xend = xend, yend = yend),
              arrow = arrow(length = unit(0.2, "cm"),
                            ends = "both", type = "closed")) +
  geom_curve(
    aes(x = 7, y = -5, xend = 11, yend = -25),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    curvature = -0.1
  ) +
  geom_mark_ellipse(aes(x = 6.5, y = 0), expand = unit(0.25, "cm")) +
    #geom_circle(aes(x0 = 6, y0 = 0, r = 2), inherit.aes = F) +
  geom_text(data = label_data, aes(x, y, label = label, colour = colour)) + 
  scale_fill_manual(values=c("white", "blue", "yellow")) +
  scale_colour_manual(values=c("black", "darkgrey", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none")  +
  labs(x = "Model size", y = "$Delta elpdPlain$") )

save_tikz_plot(p, width = 5, filename = "./tex/over-fit.tex")

