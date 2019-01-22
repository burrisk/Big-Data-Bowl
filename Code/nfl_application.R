library(tidyverse)
library(gganimate)
library(cowplot)
library(ggmap)
library(ggvoronoi)
library(nleqslv)
library(nloptr)
library(reshape)

file.tracking <- "Data/tracking_gameId_2017100112.csv"
tracking.example <- read_csv(file.tracking)

file.game <- "Data/games.csv"
games.sum <- read_csv(file.game) 

file.plays <- "Data/plays.csv"
plays.sum <- read_csv(file.plays) 

tracking.example.merged <- tracking.example %>% inner_join(games.sum) %>% inner_join(plays.sum) 


plot_play <- function(play, data){
  ## General field boundaries
  xmin <- 0
  xmax <- 160/3
  hash.right <- 38.35
  hash.left <- 12
  hash.width <- 3.3

  play_data <- data %>%
    filter(playId == play)

  ## Specific boundaries for a given play
  ymin <- max(round(min(play_data$x, na.rm = TRUE) - 10, -1), 0)
  ymax <- min(round(max(play_data$x, na.rm = TRUE) + 10, -1), 120)
  df.hash <- expand.grid(x = c(0, 23.36667, 29.96667, xmax), y = (10:110))
  df.hash <- df.hash %>% filter(!(floor(y %% 5) == 0))
  df.hash <- df.hash %>% filter(y < ymax, y > ymin)

  animate.play <- ggplot() +
    geom_point(data = play_data, aes(x = (xmax-y), y = x,
                                        colour = team, group = nflId, pch = team, size = team)) +
    geom_text(data = play_data, aes(x = (xmax-y), y = x, label = jerseyNumber), colour = "white",
              vjust = 0.36, size = 3.5) +
    scale_size_manual(values = c(6, 4, 6), guide = FALSE) +
    scale_shape_manual(values = c(19, 16, 19), guide = FALSE) +
    scale_colour_manual(values = c("black", "#654321", "darkorange1"), guide = FALSE) +
    annotate("text", x = df.hash$x[df.hash$x < 55/2],
             y = df.hash$y[df.hash$x < 55/2], label = "_", hjust = 0, vjust = -0.2) +
    annotate("text", x = df.hash$x[df.hash$x > 55/2],
             y = df.hash$y[df.hash$x > 55/2], label = "_", hjust = 1, vjust = -0.2) +
    annotate("segment", x = xmin,
             y = seq(max(10, ymin), min(ymax, 110), by = 5),
             xend =  xmax,
             yend = seq(max(10, ymin), min(ymax, 110), by = 5)) +
    annotate("text", x = rep(hash.left, 11), y = seq(10, 110, by = 10),
             label = c("G   ", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "   G"),
             angle = 270, size = 4) +
    annotate("text", x = rep((xmax - hash.left), 11), y = seq(10, 110, by = 10),
             label = c("   G", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "G   "),
             angle = 90, size = 4) +
    annotate("segment", x = c(xmin, xmin, xmax, xmax),
             y = c(ymin, ymax, ymax, ymin),
             xend = c(xmin, xmax, xmax, xmin),
             yend = c(ymax, ymax, ymin, ymin), colour = "black") +
    ylim(ymin, ymax) +
    coord_fixed() +
    theme_nothing() +
    transition_time(frame.id)  +
    ease_aes('linear') +
    NULL

  ## Ensure timing of play matches 10 frames-per-second
  play.length.ex <- length(unique(play_data$frame.id))
  animate(animate.play, fps = 10, nframe = play.length.ex)

}

td_plot <- plot_play(1540, tracking.example.merged)
td_plot

anim_save("figures/td_plot.gif", td_plot)

model <- load_model_hdf5("time_model.h5")

example.play <- tracking.example.merged %>% filter(playId == 1540)

plot_ownership <- function(play, frame){
  beaten <- play %>%
    filter(frame.id == frame) %>%
    arrange(nflId) %>%
    dplyr::rename(x_pos = y, y_pos = x) %>%
    select(nflId, displayName, jerseyNumber, team, x_pos, y_pos, s, dir)
  
  find_angle <- function(x0, y0, x1, y1, dir){
    init_angle <- dir * pi / 180
    diff_x <- x1 - x0; diff_y <- y1 - y0
    atan2(cos(init_angle) * diff_y - diff_x * sin(init_angle),
          diff_x * cos(init_angle) + diff_y * sin(init_angle))
  }
  
  xmin <- 0
  xmax <- 160/3
  hash.right <- 38.35
  hash.left <- 12
  hash.width <- 3.3
  
  ## Specific boundaries for a given play
  ymin <- max(round(min(example.play$x, na.rm = TRUE) - 10, -1), 0)
  ymax <- min(round(max(example.play$x, na.rm = TRUE) + 10, -1), 120)
  
  discretized_field <- expand.grid(x1 = seq(xmin, xmax, by = 0.5), y1 = seq(ymin, ymax, by = 0.5)) 
  player_expanded <- expand.grid.df(beaten %>%
                                      filter(displayName != "football"), discretized_field)
  player_expanded <- player_expanded %>%
    mutate(distance = sqrt((x1 - x_pos)^2 + (y1 - y_pos) ^ 2),
           angle = find_angle(x_pos, y_pos, x1, y1, dir)) %>%
    dplyr::rename(speed = s)
  
  times <- model %>% 
    predict(player_expanded[, c("distance", "speed", "angle")] %>% as.matrix())
  
  player_expanded$time <- times[,1]
  
  control <- player_expanded %>%
    group_by(x1, y1) %>%
    arrange(time) %>%
    filter(row_number() == 1)
  
  
  plot_beaten <- ggplot() +
    geom_raster(data = control, aes(x = xmax - x1, y = y1, fill = team), alpha = 0.5) +
    geom_point(data = beaten, aes(x = (xmax-x_pos), y = y_pos,
                                  colour = team, group = nflId, pch = team, size = team)) +
    geom_text(data = beaten, aes(x = (xmax-x_pos), y = y_pos, label = jerseyNumber), colour = "white",
              vjust = 0.36, size = 7) +
    scale_size_manual(values = c(12, 8, 12), guide = FALSE) +
    scale_shape_manual(values = c(19, 16, 19), guide = FALSE) +
    scale_colour_manual(values = c("black", "#654321", "darkorange"), guide = FALSE) +
    scale_fill_manual(values = c("grey", "orange"), guide = FALSE) +
    annotate("text", x = df.hash$x[df.hash$x < 55/2],
             y = df.hash$y[df.hash$x < 55/2], label = "_", hjust = 0, vjust = -0.2) +
    annotate("text", x = df.hash$x[df.hash$x > 55/2],
             y = df.hash$y[df.hash$x > 55/2], label = "_", hjust = 1, vjust = -0.2) +
    annotate("segment", x = xmin,
             y = seq(max(10, ymin), min(ymax, 110), by = 5),
             xend =  xmax,
             yend = seq(max(10, ymin), min(ymax, 110), by = 5)) +
    annotate("text", x = rep(hash.left, 11), y = seq(10, 110, by = 10),
             label = c("G   ", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "   G"),
             angle = 270, size = 8) +
    annotate("text", x = rep((xmax - hash.left), 11), y = seq(10, 110, by = 10),
             label = c("   G", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "G   "),
             angle = 90, size = 8) +
    annotate("segment", x = c(xmin, xmin, xmax, xmax),
             y = c(ymin, ymax, ymax, ymin),
             xend = c(xmin, xmax, xmax, xmin),
             yend = c(ymax, ymax, ymin, ymin), colour = "black") +
    ylim(ymin, ymax) +
    coord_fixed() +
    theme_nothing()
  plot_beaten
}

p15 <- plot_ownership(example.play, 15)
p38 <- plot_ownership(example.play, 38)
p44 <- plot_ownership(example.play, 44)
p50 <- plot_ownership(example.play, 50)

p15
p38
p44
p50

ggsave("figures/frame15.pdf", p15, device = "pdf")
ggsave("figures/frame38.pdf", p38, device = "pdf")
ggsave("figures/frame44.pdf", p44, device = "pdf")
ggsave("figures/frame50.pdf", p50, device = "pdf")


