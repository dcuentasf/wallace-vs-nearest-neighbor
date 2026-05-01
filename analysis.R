library(ggplot2)

set.seed(123)

# Wallace rule function
melting_temp <- function(seq){
  seq <- unlist(strsplit(seq, split=""))
  counts <- table(seq)
  a <- sum(counts["A"],na.rm=T)
  c <- sum(counts["C"],na.rm=T)
  g <- sum(counts["G"],na.rm=T)
  t <- sum(counts["T"],na.rm=T)
  t_melt <- 4*(g+c)+2*(a+t)
  return(t_melt)
}


# Nearest-Neighbor model function
tm_nn <- function(primer,
                  conc_primer = 0.25e-6,
                  conc_na = 50e-3,
                  self_comp = FALSE){
  
  primer <- toupper(primer)
  bases <- strsplit(primer, "")[[1]]
  n <- length(bases)
  
  dinuc <- character(n - 1)
  for(i in 1:(n - 1)){
    dinuc[i] <- paste0(bases[i], bases[i+1])
  }
  
  nn_map <- c(
    "AA"="AA/TT", "TT"="AA/TT",
    "AT"="AT/TA",
    "TA"="TA/AT",
    "CA"="CA/GT", "TG"="CA/GT",
    "GT"="GT/CA", "AC"="GT/CA",
    "CT"="CT/GA", "AG"="CT/GA",
    "GA"="GA/CT", "TC"="GA/CT",
    "CG"="CG/GC",
    "GC"="GC/CG",
    "GG"="GG/CC", "CC"="GG/CC"
  )
  
  keys <- nn_map[dinuc]
  
  nn_params <- data.frame(
    pair = c("AA/TT","AT/TA","TA/AT","CA/GT","GT/CA",
             "CT/GA","GA/CT","CG/GC","GC/CG","GG/CC"),
    dH = c(-7.9, -7.2, -7.2, -8.5, -8.4,
           -7.8, -8.2, -10.6, -9.8, -8.0),
    dS = c(-22.2, -20.4, -21.3, -22.7, -22.4,
           -21.0, -22.2, -27.2, -24.4, -19.9)
  )
  
  idx <- match(keys, nn_params$pair)
  
  dH_total <- sum(nn_params$dH[idx])
  dS_total <- sum(nn_params$dS[idx])
  
  for (ext in c(bases[1], bases[n])) {
    if (ext %in% c("G","C")) {
      dH_total <- dH_total + 0.1
      dS_total <- dS_total - 2.8
    } else if (ext %in% c("A","T")) {
      dH_total <- dH_total + 2.3
      dS_total <- dS_total + 4.1
    } else {
      stop("Invalid base")
    }
  }
  
  if(self_comp){
    dS_total <- dS_total - 1.4
  }
  
  R <- 1.987
  
  Ct <- if(self_comp) conc_primer/2 else conc_primer/4
  
  Tm <- (1000 * dH_total) / (dS_total + R * log(Ct)) - 273.15 + 16.6 * log10(conc_na)
  
  return(Tm)
}

# Simulation of three groups of primers
small <- c(12:17)
medium <- c(18:24)
large <- c(25:30)

bases <- c("A","T","G")

s_primers <- vector()
m_primers <- vector()
l_primers <- vector()

for(i in 1:1000){
  primer_length <- sample(small,1)
  primer <- sample(bases,primer_length,replace=T)
  primer <- paste(primer,collapse="")
  s_primers[i] <- primer
}

for(i in 1:1000){
  primer_length <- sample(medium,1)
  primer <- sample(bases,primer_length,replace=T)
  primer <- paste(primer,collapse="")
  m_primers[i] <- primer
}

for(i in 1:1000){
  primer_length <- sample(large,1)
  primer <- sample(bases,primer_length,replace=T)
  primer <- paste(primer,collapse="")
  l_primers[i] <- primer
}

# Calculation of Tm with Wallace rule
wallace_s <- vector()
wallace_m <- vector()
wallace_l <- vector()

for(i in 1:1000){
  p <- s_primers[i]
  tm <- melting_temp(p)
  wallace_s[i] <- tm
}

for(i in 1:1000){
  p <- m_primers[i]
  tm <- melting_temp(p)
  wallace_m[i] <- tm
}

for(i in 1:1000){
  p <- l_primers[i]
  tm <- melting_temp(p)
  wallace_l[i] <- tm
}

# Calculation of Tm with Nearest-Neighbor
nn_s <- vector()
nn_m <- vector()
nn_l <- vector()

for(i in 1:1000){
  p <- s_primers[i]
  tm <- tm_nn(p)
  nn_s[i] <- tm
}

for(i in 1:1000){
  p <- m_primers[i]
  tm <- tm_nn(p)
  nn_m[i] <- tm
}

for(i in 1:1000){
  p <- l_primers[i]
  tm <- tm_nn(p)
  nn_l[i] <- tm
}

#Differences between Wallace and NN Tm
delta_s <- nn_s - wallace_s
delta_m <- nn_m - wallace_m
delta_l <- nn_l - wallace_l

# Absolute differences
abs_delta_s <- abs(delta_s)
abs_delta_m <- abs(delta_m)
abs_delta_l <- abs(delta_l)

df <- data.frame(
  group = rep(c("small","medium","large"), each = 1000),
  wallace = c(wallace_s, wallace_m, wallace_l),
  nn = c(nn_s, nn_m, nn_l),
  delta = c(delta_s, delta_m, delta_l),
  abs_delta = c(abs_delta_s, abs_delta_m, abs_delta_l)
)


# Boxplot
boxplot_plot <- ggplot(df, aes(x = group, y = abs_delta, fill = group)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  labs(
    title = expression("|"*Delta*T[m]*"|"~"by primer length group"),
    y = expression("|"*Delta*T[m]*"|"~"(°C)"),
    x = "Group"
  )

boxplot_plot
ggsave("boxplot.pdf", plot = boxplot_plot, width = 6, height = 5)


# Scatter plot
lims <- range(c(df$wallace, df$nn))

scatter_plot <- ggplot(df, aes(x = wallace, y = nn, color = group)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed(xlim = lims, ylim = lims) +
  labs(
    title = "Wallace vs NN",
    x = expression(T[m]~"Wallace"),
    y = expression(T[m]~"NN")
  )

scatter_plot
ggsave("scatter.pdf", plot = scatter_plot, width = 6, height = 5)


# Histogram
histogram_plot <- ggplot(df, aes(x = delta, fill = group)) +
  geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
  labs(
    title = expression("Distribution of "*Delta*T[m]),
    x = expression(Delta*T[m]~"(NN - Wallace)")
  )

histogram_plot
ggsave("histogram.pdf", plot = histogram_plot, width = 6, height = 5)


# Summary table with mean and confidence interval of |dTm|
summary_df <- aggregate(abs_delta ~ group, df, function(x){
  m <- mean(x)
  se <- sd(x) / sqrt(length(x))
  c(
    mean = m,
    ci_low = m - 1.96 * se,
    ci_high = m + 1.96 * se
  )
})

# Convert to proper data frame
summary_df <- do.call(data.frame, summary_df)

# Rename columns
colnames(summary_df) <- c("group", "mean_abs_delta", "ci_low", "ci_high")

summary_df

# Save as CSV
write.csv(summary_df, "summary_table.csv", row.names = FALSE)



