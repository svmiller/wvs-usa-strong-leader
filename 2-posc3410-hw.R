library(RCurl)

data <- getURL("https://raw.githubusercontent.com/svmiller/wvs-usa-strong-leader/master/wvs-usa-sl.csv")
Data <- read.csv(text = data)

M1 <- glm(sldummy ~ z_age + I(z_age^2) + female +
              hsedorless + z_ideo + z_incscale + gop + unemployed +
              + z_lemanc + gop*hsedorless,
            data=subset(Data), family=binomial(link = "logit"))

summary(M1)

M2 <- glm(sldummy ~ z_age + I(z_age^2) + female +
            hsedorless + z_ideo + z_incscale + gop + unemployed +
            + z_lemanc + gop*hsedorless,
          data=subset(Data, year == 2011), family=binomial(link = "logit"))

summary(M2)

M3 <- glm(hddummy ~ z_age + I(z_age^2) + female +
            hsedorless + z_ideo + z_incscale + gop + unemployed +
            + z_lemanc + gop*hsedorless,
          data=subset(Data), family=binomial(link = "logit"))

summary(M3)
