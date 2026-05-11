############################################################
# Mixed Models - Alzheimer Biomarker (Clean Version)
############################################################

library(nlme)
library(rgl)

set.seed(123)

############################################################
# Simulated dataset
############################################################

MRI <- factor(rep(1:9, each = 5))

Age <- round(runif(45, 55, 85), 1)
Tracer <- round(runif(45, 0.5, 2), 2)

# Random effects
b0 <- rnorm(9, 0, 6)   # intercepts
b1 <- rnorm(9, 0, 2)   # slopes (sensitivity)

Biomarker <- numeric(45)

for(i in 1:45){
  
  scanner <- as.numeric(MRI[i])
  
  Biomarker[i] <-
    120 -
    1.1 * Age[i] +
    6 * Tracer[i] +
    b0[scanner] +
    b1[scanner] * Tracer[i] +
    rnorm(1, 0, 4)
}

ALZ <- data.frame(MRI, Age, Tracer, Biomarker)

############################################################
# 1) RANDOM INTERCEPT MODEL
############################################################

M1 <- lme(Biomarker ~ Age + Tracer,
          random = ~1 | MRI,
          data = ALZ)

summary(M1)

############################################################
# 2) RANDOM INTERCEPT + SLOPE MODEL
############################################################

M2 <- lme(Biomarker ~ Age + Tracer,
          random = ~1 + Tracer | MRI,
          data = ALZ)

summary(M2)

############################################################
# 3D GRID FOR SURFACES
############################################################

Age_seq <- seq(min(ALZ$Age), max(ALZ$Age), length = 20)
Tracer_seq <- seq(min(ALZ$Tracer), max(ALZ$Tracer), length = 20)

grid <- expand.grid(Age = Age_seq,
                    Tracer = Tracer_seq,
                    MRI = levels(ALZ$MRI))

############################################################
# 3D PLOT - RANDOM INTERCEPT MODEL
############################################################

grid$Pred1 <- predict(M1, newdata = grid, level = 1)

cols <- rainbow(9)

open3d()
plot3d(ALZ$Age, ALZ$Tracer, ALZ$Biomarker,
       col = cols[ALZ$MRI],
       size = 6,
       xlab = "Age",
       ylab = "Tracer",
       zlab = "Biomarker",
       main = "Random Intercept Model")

for(i in 1:9){
  
  g <- grid[grid$MRI == i,]
  
  z <- matrix(g$Pred1, nrow = length(Age_seq))
  
  surface3d(Age_seq,
            Tracer_seq,
            z,
            color = cols[i],
            alpha = 0.4)
}

############################################################
# 3D PLOT - RANDOM INTERCEPT + SLOPE MODEL
############################################################

grid$Pred2 <- predict(M2, newdata = grid, level = 1)

open3d()
plot3d(ALZ$Age, ALZ$Tracer, ALZ$Biomarker,
       col = cols[ALZ$MRI],
       size = 6,
       xlab = "Age",
       ylab = "Tracer",
       zlab = "Biomarker",
       main = "Random Intercept + Slope Model")

for(i in 1:9){
  
  g <- grid[grid$MRI == i,]
  
  z <- matrix(g$Pred2, nrow = length(Age_seq))
  
  surface3d(Age_seq,
            Tracer_seq,
            z,
            color = cols[i],
            alpha = 0.4)
}


