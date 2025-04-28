library(ggplot2)

df.raw <- read.csv("U:/home/bnplab-admin/TMS_localization/Figures/contrasts.csv")
df <- data.frame(df)
df$Response[startsWith(df$Receiver, "+")] <- sqrt(sqrt(df.raw$Response[startsWith(df.raw$Receiver, "+")]))


# Example plots:
# hot vs. cold on Left:
(p <- ggplot(data=df[df$Receiver == "-FDI" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())
(p <- ggplot(data=df[df$Contrast == "CsE_FDI_in_uV_vs_SIHIscore_FDI" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())

(p <- ggplot(data=df[df$Contrast == "CsE_ADM_in_uV_vs_CsE_FDI_in_uV" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())
(p <- ggplot(data=df[df$Contrast == "CsE_FDI_in_uV_vs_CsE_ADM_in_uV" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())
(p <- ggplot(data=df[df$Contrast == "CsE_APB_in_uV_vs_CsE_FDI_in_uV" & df$Hemisphere == "L",], aes(y=Response, fill=Coil, x=Subject)) + geom_boxplot())

# Edited Contrast to instead be "+APB_vs_+FDI"