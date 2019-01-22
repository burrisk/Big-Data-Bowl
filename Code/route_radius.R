library(tidyverse)
library(nloptr)
library(mice)
library(truncnorm)

time <- function(route_pars, start_pos, end_pos, bounds, inits, gradient = T){
  decel_constant <- bounds[1]; accel_constant <- bounds[2]; max_speed <- bounds[3]
  turn_speed <- route_pars[1]; turn_radius <- route_pars[2]; total_turn <- route_pars[3];
  accel_time <- route_pars[4]
  init_speed <- inits[1]; init_angle <- inits[2]
  decel_time <- (log(init_speed) - log(turn_speed)) / decel_constant
  turn_time <- abs(total_turn) * turn_radius / turn_speed
  objective <- decel_time + turn_time + accel_time
  if (gradient){
    grad <- c(-1 / (decel_constant * turn_speed) - abs(total_turn) * turn_radius / turn_speed ^ 2,
              abs(total_turn) / turn_speed,
              sign(total_turn) * turn_radius / turn_speed,
              1)
    return(list(
      "objective" = objective,
      "gradient" = grad
    ))
  }
  list("objective" = objective)
}

# time(route, start_pos, end_pos, bounds, inits)

feasible_route_BFGS <- function(start_pos, end_pos, bounds, inits, clockwise_turn = 1){
  x0 <- start_pos[1]; y0 <- start_pos[2]; x1 <- end_pos[1]; y1 <- end_pos[2]
  decel_constant <- bounds[1]; accel_constant <- bounds[2]; max_speed <- bounds[3]; max_ca <- bounds[4]
  init_speed <- inits[1]; init_angle <- inits[2]
  
  turn_speed <- init_speed / 10; turn_radius <- 2 * turn_speed ^ 2 / max_ca
  f <- function(pars){
    sum((displacement(c(turn_speed, turn_radius, pars[1], pars[2]), bounds, inits)$objective -
           (end_pos - start_pos))^2)
  }
  g <- function(pars){
    2 *  matrix((displacement(c(turn_speed, turn_radius, pars[1], pars[2]), bounds, inits)$objective
     - (end_pos - start_pos)), ncol = 2) %*% 
      displacement_gradient(c(turn_speed, turn_radius, pars[1], pars[2]), bounds, inits)[, 3:4]
  }
  if (clockwise_turn != 1){
    k <- optim(c(pi / 2, 1), f, g, method = "L-BFGS-B", lower = c(0, 0), upper = c(2 * pi, Inf),
               control = list("maxit" = 500))
  } else{
    k <- optim(c(-pi / 2, 1), f, g, method = "L-BFGS-B", lower = c(-2 * pi, 0), upper = c(0, Inf),
               control = list("maxit" = 500)) 
  }
  route_pars <- c(turn_speed, turn_radius, k$par[1], k$par[2])
  convergence <- ifelse(k$convergence == 0, 1, 0)
  list("route" = route_pars, "convergence" = convergence)
}

displacement <- function(route_pars, bounds, inits, gradient = F){
  # Extract Parameters
  decel_constant <- bounds[1]; accel_constant <- bounds[2]; max_speed <- bounds[3];
  turn_speed <- route_pars[1]; turn_radius <- route_pars[2]; total_turn <- route_pars[3];
  accel_time <- route_pars[4]
  init_speed <- inits[1]; init_angle <- inits[2]
  rotation_angle <- init_angle + total_turn
  
  x_decel <- cos(init_angle) / decel_constant * (init_speed - turn_speed)
  x_turn <- turn_radius * sign(total_turn) * (sin(rotation_angle) - sin(init_angle))
  x_accel <- cos(rotation_angle) * (max_speed * accel_time + (turn_speed - max_speed) / accel_constant
                                    * (1 - exp(-accel_constant * accel_time)))
  y_decel <- sin(init_angle) / decel_constant * (init_speed - turn_speed)
  y_turn <- -turn_radius * sign(total_turn) * (cos(rotation_angle) - cos(init_angle))
  y_accel <- sin(rotation_angle) * (max_speed * accel_time + (turn_speed - max_speed) / accel_constant
                                    * (1 - exp(-accel_constant * accel_time)))
  objective <- c(x_decel + x_turn + x_accel, y_decel + y_turn + y_accel)
  
  if (gradient){
    grad <- matrix(nrow = 2, ncol = 4)
    
    grad[1, 1] <- -cos(init_angle) / decel_constant + 
      cos(rotation_angle) / accel_constant * (1 - exp(-accel_constant * accel_time))
    grad[1, 2] <- sign(total_turn) * (sin(rotation_angle) - sin(init_angle))
    grad[1, 3] <- turn_radius * sign(total_turn) * cos(rotation_angle)  - sin(rotation_angle) * 
      (max_speed * accel_time + (turn_speed - max_speed) / accel_constant
                                                         * (1 - exp(-accel_constant * accel_time)))
    grad[1, 4] <- cos(rotation_angle) * ((turn_speed - max_speed) * exp(-accel_constant * accel_time)
                                                           + max_speed)
    grad[2, 1] <- -sin(init_angle) / decel_constant + 
      sin(rotation_angle) / accel_constant * (1 - exp(-accel_constant * accel_time))
    grad[2, 2] <-  -sign(total_turn) * (cos(rotation_angle) - cos(init_angle))
    grad[2, 3] <-  turn_radius * sign(total_turn) * sin(rotation_angle) +
      cos(rotation_angle) * (max_speed * accel_time + (turn_speed - max_speed) / accel_constant
                                          * (1 - exp(-accel_constant * accel_time)))
    grad[2, 4] <- sin(rotation_angle) * ((turn_speed - max_speed) * exp(-accel_constant * accel_time)
                                         + max_speed)
    return(list("objective" = objective, "gradient" = grad))
  } else{
    return(list("objective" = objective))
  }
}

displacement_gradient <- function(route_pars, bounds, inits){
  # Extract Parameters
  decel_constant <- bounds[1]; accel_constant <- bounds[2]; max_speed <- bounds[3];
  turn_speed <- route_pars[1]; turn_radius <- route_pars[2]; total_turn <- route_pars[3];
  accel_time <- route_pars[4]
  init_speed <- inits[1]; init_angle <- inits[2]
  rotation_angle <- init_angle + total_turn
  grad <- matrix(nrow = 2, ncol = 4)
  grad[1, 1] <- -cos(init_angle) / decel_constant + 
    cos(rotation_angle) / accel_constant * (1 - exp(-accel_constant * accel_time))
  grad[1, 2] <- sign(total_turn) * (sin(rotation_angle) - sin(init_angle))
  grad[1, 3] <- turn_radius * sign(total_turn) * cos(rotation_angle)  - sin(rotation_angle) * 
    (max_speed * accel_time + (turn_speed - max_speed) / accel_constant
     * (1 - exp(-accel_constant * accel_time)))
  grad[1, 4] <- cos(rotation_angle) * ((turn_speed - max_speed) * exp(-accel_constant * accel_time)
                                       + max_speed)
  grad[2, 1] <- -sin(init_angle) / decel_constant + 
    sin(rotation_angle) / accel_constant * (1 - exp(-accel_constant * accel_time))
  grad[2, 2] <-  -sign(total_turn) * (cos(rotation_angle) - cos(init_angle))
  grad[2, 3] <-  turn_radius * sign(total_turn) * sin(rotation_angle) +
    cos(rotation_angle) * (max_speed * accel_time + (turn_speed - max_speed) / accel_constant
                           * (1 - exp(-accel_constant * accel_time)))
  grad[2, 4] <- sin(rotation_angle) * ((turn_speed - max_speed) * exp(-accel_constant * accel_time)
                                       + max_speed)
  grad
}

movement_constraints <- function(route_pars, bounds, inits, gradient = F){
  # Unpack parameters
  turn_speed <- route_pars[1]; turn_radius <- route_pars[2]
  ca <- turn_speed ^ 2 / turn_radius
  
  constraints <- c(ca, -ca)
  if (gradient){
    grad <- rbind(c(2 * turn_speed / turn_radius, -(turn_speed / turn_radius) ^ 2, 0, 0),
                  c(-2 * turn_speed / turn_radius, (turn_speed / turn_radius) ^ 2, 0, 0))
    return(list("constraints" = constraints,
                "jacobian" = grad))
  } else{
    return(list("constraints" = constraints))
  }
}

equality_function <- function(route_pars, start_pos, end_pos, bounds, inits, gradient = T){
  d <- displacement(route_pars, bounds, inits, gradient = T)
  objective <- d$objective + start_pos - end_pos
  grad <- d$gradient
  list("constraints" = objective, "jacobian" = grad)
}

inequality_function <- function(route_pars, start_pos, end_pos, bounds, inits, gradient = T){
  mov <- movement_constraints(route_pars, bounds, inits, gradient = T)
  list("constraints" = c(mov$constraints[1] - bounds[4], mov$constraints[2] - bounds[4]),
       "jacobian" = mov$jacobian)
}

feasible_route_NLEQ <- function(start_pos, end_pos, bounds, inits, clockwise_turn = 1){
  x0 <- start_pos[1]; y0 <- start_pos[2]; x1 <- end_pos[1]; y1 <- end_pos[2]
  decel_constant <- bounds[1]; accel_constant <- bounds[2]; max_speed <- bounds[3]; max_ca <- bounds[4]
  init_speed <- inits[1]; init_angle <- inits[2]
  
  turn_speed <- 0.01; turn_radius <- 2 * turn_speed ^ 2 / max_ca
  diff_x <- x1 - x0; diff_y <- y1 - y0

  ending_difference <- function(pars){
    total_turn <- pars[1]
    accel_time <- exp(pars[2])
    route_pars <- c(turn_speed, turn_radius, total_turn, accel_time)
    movement <- displacement(route_pars, bounds, inits)$objective
    c(diff_x - movement[1], diff_y - movement[2])
  }
  initial_solns <- nleqslv(c(0, 0), ending_difference,
                           control = list("maxit" = 1000, "allowSingular" = TRUE),
                           method = "Newton")
  initial_total_turn <- initial_solns$x[1]
  while (abs(initial_total_turn) >= 2 * pi){
    initial_total_turn <- initial_total_turn + floor(abs(initial_total_turn) / (2 * pi)) *
      - sign(initial_total_turn) * 2 * pi 
  }
  if (-sign(initial_total_turn) != clockwise_turn){
    initial_total_turn <- initial_total_turn - clockwise_turn * 2 * pi
  }
  if (abs(initial_total_turn) > round(2 * pi, 6)){
    initial_total_turn <- 0
  }
  convergence <- 1
  if (initial_solns$termcd != 1){
    convergence <- 0
  }
  list("route" = c(turn_speed, turn_radius, initial_total_turn,
    exp(initial_solns$x[2])), "convergence" = convergence)
}

min_time <- function(start_pos, end_pos, bounds, inits){
  x0 <- start_pos[1]; y0 <- start_pos[2]; x1 <- end_pos[1]; y1 <- end_pos[2]
  init_speed <- inits[1]; init_angle <- inits[2]
  diff_x <- x1 - x0; diff_y <- y1 - y0
  find_angle <- function(diff_x, diff_y, init_angle){
    distance <- sqrt(diff_x^2 + diff_y^2)
    atan2(cos(init_angle) * diff_y - diff_x * sin(init_angle),
          diff_x * cos(init_angle) + diff_y * sin(init_angle))
  }
  angle_with_velocity <- find_angle(diff_x, diff_y, init_angle)
  clockwise_turn <- 2 * (angle_with_velocity < 0) - 1
  route_init <- feasible_route_BFGS(start_pos, end_pos, bounds, inits, clockwise_turn)
  if (route_init$convergence == 0){
    return(list("route" = rep(NA, 4), "time" = NA, "init_fail" = T))
  }
  solution <- nloptr(route_init$route, eval_f = time, lb = c(0, 0, min(-2 * pi * clockwise_turn, 0), 0),
                     ub = c(inits[1], Inf, max(-2 * pi * (clockwise_turn), 0), Inf),
                     eval_g_ineq = inequality_function,
                     eval_g_eq = equality_function,
                     start_pos = start_pos, end_pos = end_pos, bounds = bounds, inits = inits,
                     gradient = T,
                     opts = list("algorithm" = "NLOPT_LD_SLSQP", "local_opts" =
                                   list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1e-7),
                                 "xtol_rel" = 1e-7, "maxeval" = 1000))
  
  if (sum(abs(displacement(solution$solution, bounds, inits)$objective -
              (end_pos - start_pos))) > 0.1){
    return(list("route" = rep(NA, 4), "time" = NA, "init_fail" = F))
  }
  if (sum(solution$solution == route_init$route) == 4){
    return(list("route" = rep(NA, 4), "time" = NA, "init_fail" = F))
  }
  
  return(list("route" = solution$solution, "time" = solution$objective,
              "init_fail" = F))
}

create_timeplot <- function(distance, bounds){
  start_pos <- c(10.3, 27.2)
  end_pos <- c(10.3 + distance, 27.2)
  speeds <- seq(0.2, bounds[3], length = 20)
  angle <- seq(-pi, pi, length = 20)
  df_inits <- as.matrix(expand.grid(speeds, angle))
  df_times <- sapply(1:nrow(df_inits), function(i){
    inits <- df_inits[i, ]
    min_time(start_pos, end_pos, bounds, inits)
  })
  t <- matrix(unlist(df_times), ncol = 6, byrow = T)
  df_inits <- cbind(df_inits, t)
  df_plot <- tibble(Velo = df_inits[,1], Angle = df_inits[,2], Time = df_inits[, 7]) %>%
    mutate(Angle = Angle / pi * 180)
  data_impute <- mice(df_plot, m = 1, method = "rf")
  df_plot <- complete(data_impute)
  g <- ggplot(df_plot, aes(x = Angle, y = Velo, fill = Time)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", direction = 1, name = "Time (s)") +
    coord_polar(start = pi) +
    xlab("Initial Angle (degrees)") +
    ylab("Initial Velocity (m/s)")  +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"))
  list("g" = g, "data" = df_plot)
}

bounds <- c(2, 2, 9, 2 * pi)
plot2 <- create_timeplot(2, bounds)
plot5 <- create_timeplot(5, bounds)
plot15 <- create_timeplot(15, bounds)
plot40 <- create_timeplot(40, bounds)

ggsave("figures/2yards.pdf", plot2$g, device = "pdf")
ggsave("figures/5yards.pdf", plot5$g, device = "pdf")
ggsave("figures/15yards.pdf", plot15$g, device = "pdf")
ggsave("figures/40yards.pdf", plot40$g, device = "pdf")


min_time_dist <- function(distance, velocity, angle, bounds){
  start_pos <- c(10.3, 27.2)
  end_pos <- c(10.3 + distance, 27.2)
  inits <- c(velocity, angle)
  time <- min_time(start_pos, end_pos, bounds, inits)$time
  time
}


# Fit Neural Network
N <- 15000
x <- rtruncnorm(N, a = 0, b = Inf, mean = 20, sd = 20)
y <- runif(N, -pi, pi)
v <- runif(N, 0.05, bounds[3])
z <- sapply(1:N, function(i){
  min_time_dist(x[i], v[i], y[i], bounds)
})
mat <- tibble(distance = x, speed = v, angle = y, time = z)
write_csv(mat, "speed_est.csv")
