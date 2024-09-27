library(ggplot2)

df.raw <- read.csv("U:/home/bnplab-admin/TMS_localization/Figures/contrasts.csv")
df <- data.frame(df.raw)
df$Response[startsWith(df$Receiver, "+")] <- sqrt(sqrt(df.raw$Response[startsWith(df.raw$Receiver, "+")]))

# Renaming:
df$Coil <- gsub("\\+", "h", df$Coil)
df$Coil <- gsub("\\-", "c", df$Coil)
df$Receiver <- gsub("\\+", "h", df$Receiver)
df$Receiver <- gsub("\\-", "c", df$Receiver)
df$Contrast <- gsub("\\+", "h", df$Contrast)
df$Contrast <- gsub("\\-", "c", df$Contrast)
df$Dominant <- ifelse(df$Hemisphere == "L", "dominant", "nondominant")
df$Dominant[df$Subject == "sub-006"] <- ifelse(df$Hemisphere[df$Subject == "sub-006"] == "L", "nondominant", "dominant")

# Example plots:
# hot vs. cold on Left:
(p <- ggplot(data=df[df$Receiver == "cFDI" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())
(p <- ggplot(data=df[df$Contrast == "hFDI_vs_cFDI" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())

(p <- ggplot(data=df[df$Contrast == "hAPB_vs_hFDI" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())


(p <- ggplot(data=df[df$Contrast == "hADM_vs_hFDI" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())



#library(lme4)
library(car)
library(performance)
library(nlme)

fit.contrast.model <- function(data, c1, c2, dominance) {
  df <- data[data$Contrast == sprintf("%s_vs_%s", c1, c2) & data$Dominant==dominance,]
  df$Coil <- relevel(factor(df$Coil), ref=c1)
  mdl <- lme(fixed = Response ~ Coil,
             random = ~ 1|Subject,
             weights = varIdent(form= ~ (1|Subject/Coil)), # maybe also across subjects? (1|Coil/Subject)?
             data=df, method="ML") # TODO: There is a problem when using REML --- check this again!
  return(mdl)
}



compare.spots <- function(data, c1, c2, dominance) {
  mdl <- fit.contrast.model(data, c1, c2, dominance)
  summary(mdl)
  c1.value <- fixef(mdl)["(Intercept)"]
  difference <- fixef(mdl)[sprintf("Coil%s", c2)]
  c2.value <- c1.value + difference
  r <- c(c1.value, c2.value, difference)
  names(r) <- c(c1, c2, "difference")
  return(r)
}


spot.difference.confInt <- function(data, c1, c2, dominance) {
  mdl <- fit.contrast.model(data, c1, c2, dominance)
  c.in <- intervals(mdl)
  name.c2 <- sprintf("Coil%s", c2)
  df <- data.frame("Contrast"  =c(sprintf("%s_vs_%s", c1, c2),sprintf("%s_vs_%s", c1, c2)),
                   "Hemisphere"=c(dominance, dominance),
                   "Coil"      =c(c1, c2),
                   "Receiver"  =c(c1, c1),
                   "Estimate"  =c(c.in$fixed["(Intercept)", "est."],  c.in$fixed[name.c2, "est."]), 
                   "Lower"     =c(c.in$fixed["(Intercept)", "lower"], c.in$fixed[name.c2, "lower"]), 
                   "Upper"     =c(c.in$fixed["(Intercept)", "upper"], c.in$fixed[name.c2, "upper"]),
                   "p.value"   =c(summary(mdl)$tTable[name.c2,"p-value"], summary(mdl)$tTable[name.c2,"p-value"]),
                   "t.value"   =c(summary(mdl)$tTable[name.c2,"t-value"], summary(mdl)$tTable[name.c2,"t-value"]))
  return(df)
}



mdl <- fit.contrast.model(df, "cFDI", "hFDI", "dominant")
mdl$apVar


library(stringr)
all.results <- data.frame()
for(c in unique(df$Contrast)) {
  contrastants <- str_split(c, "_vs_")[[1]]
  t.L <- spot.difference.confInt(df, contrastants[1], contrastants[2], "dominant")
  t.R <- spot.difference.confInt(df, contrastants[1], contrastants[2], "nondominant")
  all.results <- rbind(all.results, t.L, t.R)
}

write.csv(all.results,"B:/Projects/2023-01 HighDef/Results/Evaluation/group_level_spot_cross_relevance_dominant.csv", row.names = FALSE)
