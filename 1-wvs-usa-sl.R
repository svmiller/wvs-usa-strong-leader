setwd("~/Dropbox/projects/blog-posts/wvs-usa-sl")
library(foreign)
library(RCurl)
library(countrycode)
library(plyr)
library(mirt)
library(car)
library(broom)
library(arm)
library(DataCombine)
library(reshape2)
library(ggplot2)
library(knitr)
library(stargazer)

augment.ranef.mer <- function(x,
                              ci.level=0.9,
                              reorder=TRUE,
                              order.var=1) {
  tmpf <- function(z) {
    if (is.character(order.var) && !order.var %in% names(z)) {
      order.var <- 1
      warning("order.var not found, resetting to 1")
    }
    ## would use plyr::name_rows, but want levels first
    zz <- data.frame(level=rownames(z),z,check.names=FALSE)
    if (reorder) {
      ## if numeric order var, add 1 to account for level column
      ov <- if (is.numeric(order.var)) order.var+1 else order.var
      zz$level <- reorder(zz$level, zz[,order.var+1], FUN=identity)
    }
    ## Q-Q values, for each column separately
    qq <- c(apply(z,2,function(y) {
      qnorm(ppoints(nrow(z)))[order(order(y))]
    }))
    rownames(zz) <- NULL
    pv   <- attr(z, "postVar")
    cols <- 1:(dim(pv)[1])
    se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    ## n.b.: depends on explicit column-major ordering of se/melt
    zzz <- cbind(melt(zz,id.vars="level",value.name="estimate"),
                 qq=qq,std.error=se)
    ## reorder columns:
    subset(zzz,select=c(variable, level, estimate, qq, std.error))
  }
  dd <- ldply(x,tmpf,.id="grp")
  ci.val <- -qnorm((1-ci.level)/2)
  transform(dd,
            p=2*pnorm(-abs(estimate/std.error)), ## 2-tailed p-val
            lb=estimate-ci.val*std.error,
            ub=estimate+ci.val*std.error)
}


WVS <- read.dta("~/Dropbox/data/wvs/WVS_Longitudinal_1981_2014_stata_v2015_04_18.dta", convert.factors = FALSE)
EVS <- read.csv("~/Dropbox/data/evs/evs-longitudinal.csv", na.strings = c(".a",".b",".c",".d",".e"))

data <- getURL("https://raw.githubusercontent.com/svmiller/wvsccodes/master/wvs-cow-ccodes-table.csv")
wvsccodes <- read.csv(text = data)

colnames(WVS) <- tolower(names(WVS))

ncolwvs <- dim(WVS)[2]+1

WVS$wave <- WVS$s002
WVS$wvsccode <- WVS$s003
WVS$year <- WVS$s020

USA <- subset(WVS, wvsccode == 840 & wave >= 3)

# Look at democratic peers first.

WVS <- join(WVS, wvsccodes, by=c("wvsccode"), type="left", match="first")

West <- subset(WVS, ccode == 2 | ccode == 20 | (ccode >= 200 & ccode < 300) 
               | (ccode > 300 & ccode <= 338) | ccode == 350 | ccode == 352 | (ccode >= 375 & ccode < 400) | ccode >= 900)
West <- with(West, data.frame(ccode, country, year, wave, e114, e117))
West$e114 <- with(West, recode(e114, "-5:-1=NA"))
West$e117 <- with(West, recode(e117, "-5:-1=NA"))
West$study <- "WVS"
levels(West$country)[levels(West$country)=="Great Britain"] <- "United Kingdom"

EVS$evsccode <- EVS$s003
EVS <- EVS[order(EVS$evsccode), ]

EVS$country <- countrycode(EVS$evsccode, "un", "country.name")

# Some light cleaning
EVS$country[EVS$evsccode == 197] <- "Northern Cyprus"
EVS$country[EVS$evsccode == 909] <- "Northern Ireland"
EVS$country[EVS$evsccode == 915] <- "Kosovo"

EVS$country[EVS$evsccode == 807] <- "Macedonia"
EVS$country[EVS$evsccode == 643] <- "Russia"
EVS$country[EVS$evsccode == 498] <- "Moldova"


EVS$ccode <- countrycode(EVS$country, "country.name", "cown")

# More light cleaning for ccodes
EVS$ccode[EVS$evsccode == 197] <- NA
EVS$ccode[EVS$evsccode == 909] <- NA
EVS$ccode[EVS$evsccode == 688] <- 345

EVS$wave <- EVS$s002evs # 1 = 1981-1984, 2 = 1990-1993, 3 = 1999-2001, 4 = 2008-2010
EVS$year <- EVS$s020

EVS <-  subset(EVS, ccode == 2 | ccode == 20 | (ccode >= 200 & ccode < 300) 
               | (ccode > 300 & ccode <= 338) | ccode == 350 | ccode == 352 | (ccode >= 375 & ccode < 400) | ccode >= 900)
EVS <- with(EVS, data.frame(ccode, country, year, wave, e114, e117))
EVS$study <- "EVS"

Peers <- rbind(West, EVS)
Peers$sldummy <- with(Peers, recode(e114, "-5:-1=NA; 1:2=1; 3:4=0"))
Peers$slsdummy <- with(Peers, recode(e114, "-5:-1=NA; 1=1; 2:4=0"))
Peers$hddummy <- with(Peers, recode(e117, "-5:-1=NA; 1:2=0; 3:4=1"))
Peers$hdsdummy <- with(Peers, recode(e117, "-5:-1=NA; 1:3=0; 4=1"))
Peers$wave[Peers$year == 1999 & Peers$study == "EVS"] <- 4
Peers$wave[Peers$year == 2000 & Peers$study == "EVS"] <- 4
Peers$wave[Peers$year == 2008 & Peers$study == "EVS"] <- 5
Peers$wave[Peers$year == 2009 & Peers$study == "EVS"] <- 5

peers <- aggregate(cbind(sldummy, slsdummy, hddummy, hdsdummy) ~ country + ccode + year + wave + study, Peers, mean)

write.table(peers,file="peers.csv",sep=",",row.names=F,na="")


usasums <- subset(peers, ccode == 2)
usasums$ccode <- usasums$year <- usasums$study <- NULL
peers <- subset(peers, ccode != 2)

peermedians <- aggregate(cbind(sldummy, slsdummy, hddummy, hdsdummy) ~ wave, peers, median)
peermedians$country <- "Democratic peers"

usapeers <- rbind(usasums, peermedians)

# Now look at just the U.S.
# -------------------------

USA$uid <- seq(1, nrow(USA))
USA$region <- USA$x048wvs

# In order: New England, Mid-Atlantic, South Atlantic, ESC, WSC, ENC, WNC, "Rocky Mountain State", Northwest, CA, AK, HI, Pacific
# These are almost full Census regions. My hunch is they only broke apart the Pacific and that Pacific might include US territories.
USA$region <- with(USA, ifelse(region == -2, NA, region))
USA$region <- with(USA, region - 840000)
USA$censusd <- USA$region # census division
USA$censusd <- with(USA, ifelse(censusd >= 9, 9, censusd))
USA$censusd[USA$censusd == 1] <- "New England"
USA$censusd[USA$censusd == 2] <- "Middle Atlantic"
USA$censusd[USA$censusd == 3] <- "South Atlantic"
USA$censusd[USA$censusd == 4] <- "East South Central"
USA$censusd[USA$censusd == 5] <- "West South Central"
USA$censusd[USA$censusd == 6] <- "East North Central"
USA$censusd[USA$censusd == 7] <- "West North Central"
USA$censusd[USA$censusd == 8] <- "Mountain"
USA$censusd[USA$censusd == 9] <- "Pacific"

USA$censusr <- USA$region
USA$censusr <- with(USA, ifelse(censusr >= 9, 9, censusr))
USA$censusr[USA$censusr == 1 | USA$censusr == 2] <- "Northeast"
USA$censusr[USA$censusr == 6 | USA$censusr == 7] <- "Midwest"
USA$censusr[USA$censusr == 3 | USA$censusr == 4 | USA$censusr == 5] <- "South"
USA$censusr[USA$censusr == 8 | USA$censusr == 9] <- "West"

USA$raceethnic <- USA$x051
USA$raceethnic <- with(USA, ifelse(raceethnic == -2, NA, raceethnic))
# WVS *really* made a meal of this variable.
USA$raceethnic[USA$x051 == 200] <- "Black" # "Black African"
USA$raceethnic[USA$x051 == 1250] <- "Hispanic" # Spanish/Hispanic
USA$raceethnic[USA$x051 == 1400] <- "White"
USA$raceethnic[USA$x051 == 8000] <- "Other"
USA$raceethnic[USA$x051 == 8001] <- "Mixed Race"
USA$raceethnic[USA$x051 == 840002] <- "White" # White, not-Hispanic
USA$raceethnic[USA$x051 == 840003] <- "Black" # Black, not-Hispanic
USA$raceethnic[USA$x051 == 840004] <- "Other" # Other, not-Hispanic
USA$raceethnic[USA$x051 == 840005] <- "Hispanic" # Hispanic
USA$raceethnic[USA$x051 == 840006] <- "Mixed Race" # 2+ Races, Non-Hispanic
USA$raceethnic[USA$x051 == 840007] <- "Other" # South Asian
USA$raceethnic[USA$x051 == 840008] <- "Other" # East Asian
USA$raceethnic[USA$x051 == 840009] <- "Other" # Arabic
USA$raceethnic[USA$x051 == 840010] <- "White"
USA$raceethnic[USA$x051 == 840011] <- "Black"
USA$raceethnic[USA$x051 == 840012] <- "Hispanic"

USA$white <- with(USA, ifelse(raceethnic == "White", 1, 0))
USA$black <- with(USA, ifelse(raceethnic == "Black", 1, 0))
USA$hispanic <- with(USA, ifelse(raceethnic == "Hispanic", 1, 0))

USA$party <- NA
USA$party[USA$e179wvs == 5] <- "Other"
USA$party[USA$e179wvs == 840001] <- "Republican"
USA$party[USA$e179wvs == 840002] <- "Democrat"
USA$party[USA$e179wvs == 840003] <- "Independent"
USA$party[USA$e179wvs == 840004] <- "Libertarian"
USA$party[USA$e179wvs == 840005] <- "Reform"

USA$gop <- with(USA, ifelse(party == "Republican", 1, 0))



USA$strongleader <- with(USA, recode(e114, "-5:-1=NA; 1=4; 2=3; 3=2; 4=1"))
USA$armyrule <- with(USA, recode(e115, "-5:-1=NA; 1=4; 2=3; 3=2; 4=1"))
USA$expertdecision <- with(USA, recode(e116, "-5:-1=NA; 1=4; 2=3; 3=2; 4=1"))
USA$havedem <- with(USA, recode(e117, "-5:-1=NA; 1=4; 2=3; 3=2; 4=1"))

USA$demtaxrich <- with(USA, recode(e224, "-5:-1=NA"))
USA$demrelig <- with(USA, recode(e225, "-5:-1=NA"))
USA$demchooseleaders <- with(USA, recode(e226, "-5:-1=NA"))
USA$demaidunemp <- with(USA, recode(e227, "-5:-1=NA"))
USA$demarmycoup <- with(USA, recode(e228, "-5:-1=NA"))
USA$demcivilrights <- with(USA, recode(e229, "-5:-1=NA"))
USA$demeconomy <- with(USA, recode(e230, "-5:-1=NA"))
USA$demcriminals <- with(USA, recode(e231, "-5:-1=NA"))
USA$dempeoplelaws <- with(USA, recode(e232, "-5:-1=NA"))
USA$demwomenrights <- with(USA, recode(e233, "-5:-1=NA"))
USA$demincome <- with(USA, recode(e234, "-5:-1=NA"))
USA$demimp <- with(USA, recode(e235, "-5:-1=NA"))
USA$democraticness <- with(USA, recode(e236, "-5:-1=NA"))


USA$age <- with(USA, recode(x003, "-5:-1=NA"))
USA$female <- with(USA, recode(x001, "-5:-1=NA; 1=0; 2=1"))
USA$unemployed <- with(USA, recode(x028, "-5:-1=NA; 1:6=0; 7=1; 8=0"))
USA$satisfin <- with(USA, ifelse(c006 < 0, NA, c006-1))
USA$ideo <- with(USA, ifelse(e033 < 0, NA, e033-1))
USA$socialclass <- with(USA, recode(x045, "-5:-1=NA; 1=4; 2=3; 3=2; 4=1; 5=0"))
USA$incscale <- with(USA, ifelse(x047 < 0, NA, x047-1))

USA$educat <- with(USA, recode(x025, "-5:-1=NA"))
USA$hsedorless <- with(USA, recode(x025, "-5:-1=NA; 1:6=1; 7:8=0"))
USA$lessthanhs <- with(USA, recode(x025, "-5:-1=NA; 1:3=1; 4:8=0"))
USA$collegeed <- with(USA, recode(x025, "-5:-1=NA; 1:7=0; 8=1"))

USA$educatr <- NA
USA$educatr <- with(USA, ifelse((educat >= 1 & educat <= 3) | educat == 5, "Did not finish HS", educatr))
USA$educatr <- with(USA, ifelse(educat == 4 | educat == 6, "HS grad", educatr))
USA$educatr <- with(USA, ifelse(educat == 7, "Some college", educatr))
USA$educatr <- with(USA, ifelse(educat == 8, "College educated", educatr))

# Get emancipative values.

USA$emancvalues <- USA$y020
USA$autonomy <- with(USA, recode(y021, "-5:-1=NA"))
USA$equality <- with(USA, recode(y022, "-5:-1=NA"))
USA$choice <- with(USA, recode(y023, "-5:-1=NA"))
USA$voice <- with(USA, recode(y024, "-5:-1=NA"))

Autonomy <- with(USA, data.frame(uid, autonomy, a029, a034, a042))

Autonomy[,3:ncol(Autonomy)] <- sapply(Autonomy[,3:ncol(Autonomy)],
                                      function(x)ifelse(x<=-1,NA,x))

Autonomy$removeme <- with(Autonomy, ifelse(is.na(a029) & is.na(a034) & is.na(a042), 1, 0))
Autonomy <- subset(Autonomy, removeme == 0)
Autonomy$removeme <- NULL

colnames(Autonomy) <- c("uid", "autonomy", "kid_ind", "kid_imag", "kid_obed")
Autonomy$kid_obed <- with(Autonomy, recode(kid_obed, "1=0;0=1"))

AutM <- mirt(Autonomy[ ,  3:ncol(Autonomy)], model = 1,
             itemtype = "graded", SE = TRUE, verbose = FALSE)

autscores <- fscores(AutM, full.scores = TRUE, full.scores.SE = TRUE)
Autonomy <- cbind(Autonomy, autscores)
Autonomy <- rename(Autonomy, c("F1" = "laut", "SE_F1" = "se_laut"))
with(Autonomy, cor(laut, autonomy, use="complete.obs"))

USA <- join(USA, Autonomy, by=c("uid", "autonomy"), type="left", match="first")

Equality <- with(USA, data.frame(uid, equality, c001, d059, d060))

Equality[,3:ncol(Equality)] <- sapply(Equality[,3:ncol(Equality)],
                                      function(x)ifelse(x<=-1,NA,x))

Equality$removeme <- with(Equality, ifelse(is.na(c001) & is.na(d059) & is.na(d060), 1, 0))
Equality <- subset(Equality, removeme == 0)
Equality$removeme <- NULL

colnames(Equality) <- c("uid", "equality", "menjob", "menleaders", "boycollege")

EquM <- mirt(Equality[ ,  3:ncol(Equality)], model = 1,
             itemtype = "graded", SE = TRUE, verbose = FALSE)

equscores <- fscores(EquM, full.scores = TRUE, full.scores.SE = TRUE)
Equality <- cbind(Equality, equscores)
Equality <- rename(Equality, c("F1" = "lequ", "SE_F1" = "se_lequ"))
with(Equality, cor(lequ, equality, use="complete.obs"))

USA <- join(USA, Equality, by=c("uid", "equality"), type="left", match="first")


Choice <- with(USA, data.frame(uid, choice, f118, f120, f121))

Choice[,3:ncol(Choice)] <- sapply(Choice[,3:ncol(Choice)],
                                  function(x)ifelse(x<=-1,NA,x))

Choice$removeme <- with(Choice, ifelse(is.na(f118) & is.na(f120) & is.na(f121), 1, 0))
Choice <- subset(Choice, removeme == 0)
Choice$removeme <- NULL

colnames(Choice) <- c("uid", "choice", "hj", "aj", "dj")


ChoM <- mirt(Choice[ ,  3:ncol(Choice)], model = 1,
             itemtype = "graded", SE = TRUE, verbose = FALSE)

choscores <- fscores(ChoM, full.scores = TRUE, full.scores.SE = TRUE)
Choice <- cbind(Choice, choscores)
Choice <- rename(Choice, c("F1" = "lcho", "SE_F1" = "se_lcho"))
with(Choice, cor(lcho, choice, use="complete.obs"))

USA <- join(USA, Choice, by=c("uid", "choice"), type="left", match="first")


Voice <- with(USA, data.frame(uid, voice, e001, e002, e003, e004))

Voice[,3:ncol(Voice)] <- sapply(Voice[,3:ncol(Voice)],
                                function(x)ifelse(x<=-1,NA,x))

Voice$acsay <- NA
Voice$acsay <- with(Voice, ifelse(e001 == 3, 2, acsay))
Voice$acsay <- with(Voice, ifelse(e002 == 3, 1, acsay))
Voice$acsay <- with(Voice, ifelse(e001 != 3 & e002 != 3 & !is.na(e001), 0, acsay))

Voice$apsay <- NA
Voice$apsay <- with(Voice, ifelse((e003 == 2  & e004 == 4) | (e003 == 4  & e004 == 2),
                                  3, apsay))
Voice$apsay <- with(Voice, ifelse((e003 == 2  & e004 != 4) | (e003 == 4  & e004 != 2),
                                  2, apsay))
Voice$apsay <- with(Voice, ifelse((e003 != 2  & e004 == 4) | (e003 != 4  & e004 == 2),
                                  1, apsay))
Voice$apsay <- with(Voice, ifelse((e003 != 2  & e004 != 4) & (e003 != 4  & e004 != 2),
                                  0, apsay))


Voice$removeme <- with(Voice, ifelse(is.na(acsay) & is.na(apsay), 1, 0))
Voice <- subset(Voice, removeme == 0)
Voice$removeme <- NULL

VoiM <- mirt(Voice[ ,  7:ncol(Voice)], model = 1,
             itemtype = "graded", SE = TRUE, verbose = FALSE)

voiscores <- fscores(VoiM, full.scores = TRUE, full.scores.SE = TRUE)
Voice <- cbind(Voice, voiscores)
Voice <- rename(Voice, c("F1" = "lvoi", "SE_F1" = "se_lvoi"))
with(Voice, cor(lvoi, voice, use="complete.obs"))

USA <- join(USA, Voice, by=c("uid", "voice"), type="left", match="first")

# duplicate emancvalues

Emanc <- with(USA, data.frame(uid, emancvalues, laut, lequ, lcho, lvoi))
Emanc$lemanc <- with(Emanc, (1/4)*(laut + lequ + lcho + lvoi))

with(Emanc, cor(emancvalues, lemanc, use="complete.obs"))

A1 <- lm(lemanc ~ lequ + lcho + lvoi, data=Emanc) # missing laut
A2 <- lm(lemanc ~ laut + lcho + lvoi, data=Emanc) # missing lequ
A3 <- lm(lemanc ~ laut + lequ + lvoi, data=Emanc) # missing lcho
A4 <- lm(lemanc ~ laut + lequ + lcho, data=Emanc) # missing lvoi
A1df <- tidy(A1)
A2df <- tidy(A2)
A3df <- tidy(A3)
A4df <- tidy(A4)

Emanc$lemanc <- with(Emanc, ifelse(is.na(laut) & is.na(lemanc),
                                   A1df[1,2] + A1df[2,2]*lequ +
                                     A1df[3,2]*lcho + A1df[4,2]*lvoi, lemanc))

Emanc$lemanc <- with(Emanc, ifelse(is.na(lequ) & is.na(lemanc),
                                   A2df[1,2] + A2df[2,2]*laut +
                                     A2df[3,2]*lcho + A2df[4,2]*lvoi, lemanc))

Emanc$lemanc <- with(Emanc, ifelse(is.na(lcho) & is.na(lemanc),
                                   A3df[1,2] + A3df[2,2]*laut +
                                     A3df[3,2]*lequ + A3df[4,2]*lvoi, lemanc))

Emanc$lemanc <- with(Emanc, ifelse(is.na(lvoi) & is.na(lemanc),
                                   A4df[1,2] + A4df[2,2]*laut +
                                     A4df[3,2]*lequ + A4df[4,2]*lcho, lemanc))

Emanc <- with(Emanc, data.frame(uid, lemanc))

USA <- join(USA, Emanc, by=c("uid"), type="left", match="first")

# Subset to just the columns we need. Also move uid around.
USA <- USA[ ,c(2:ncol(USA), 1)]
USA <- USA[ ,  1415:ncol(USA)]
USA <- MoveFront(USA, "uid")

USA$sldummy <- with(USA, recode(strongleader, "1:2=0; 3:4=1"))
USA$hddummy <- with(USA, recode(havedem, "1:2=1; 3:4=0"))

USA <- ddply(USA, c("wave"), transform, z_sl = arm::rescale(strongleader))

USA <- ddply(USA, c("wave"), transform, z_age = arm::rescale(age))
USA <- ddply(USA, c("wave"), transform, z_ideo = arm::rescale(ideo))
USA <- ddply(USA, c("wave"), transform, z_lemanc = arm::rescale(lemanc))
USA <- ddply(USA, c("wave"), transform, z_laut = arm::rescale(laut))
USA <- ddply(USA, c("wave"), transform, z_lequ = arm::rescale(lequ))
USA <- ddply(USA, c("wave"), transform, z_lcho = arm::rescale(lcho))
USA <- ddply(USA, c("wave"), transform, z_lvoi = arm::rescale(lvoi))
USA <- ddply(USA, c("wave"), transform, z_incscale = arm::rescale(incscale))


M1 <- glmer(sldummy ~ (1 | censusr) + (1 | year) + (1 | raceethnic) + (1 | educatr)  +
               (1 | censusr:year) + (1 | raceethnic:year) + (1 | educatr:year), data=USA,
             family=binomial(link = "logit"))

M1ranef <- augment(ranef(M1,condVar=TRUE))

show_ranef <- function(data, grp){
  require(ggplot2)
  ggplot(data[data$grp == grp,],aes(estimate,level,xmin=lb,xmax=ub))+
    geom_errorbarh(height=0)+
    geom_vline(xintercept=0,lty=2)+
    geom_point()+facet_wrap(~variable,scale="free_x") + 
    ylab("Levels of the Random Effect") +
    xlab("Estimated Intercept")
}

show_ranef(M1ranef, "year")
show_ranef(M1ranef, "educatr")

educatrsldum <- with(USA, table(sldummy, educatr))
kable(prop.table(educatrsldum, 2))

show_ranef(M1ranef, "raceethnic")

raceethnicsldum <- with(USA, table(sldummy, raceethnic))
kable(prop.table(raceethnicsldum, 2))

M2 <- glmer(sldummy ~ hsedorless + (1 | censusr) + (1 | year) + (1 + hsedorless | raceethnic) + (1 | educatr) +
              (1 | censusr:year) + (1 | raceethnic:year)
              , data=USA,
            family=binomial(link = "logit"),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e5)))

M2ranef <- augment(ranef(M2,condVar=TRUE))

levels(M2ranef$variable)[levels(M2ranef$variable)=="hsedorless"] <- "HS Education or Less"

show_ranef(M2ranef, "raceethnic")


M3 <- glmer(sldummy ~ z_age + I(z_age^2) + female +
              hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
            + z_laut + z_lcho + z_lequ + z_lvoi  + gop*hsedorless +
              (1  | censusr) + (1 | year) + (1 | raceethnic) +
              (1  | censusr:year) + (1 | raceethnic:year),
            data=subset(USA), family=binomial(link = "logit"),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e5)))

M4 <- glmer(sldummy ~ z_age + I(z_age^2) + female +
              hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
              + z_laut + z_lcho + z_lequ + z_lvoi  + gop*hsedorless +
              (1  | censusr) + (1 | raceethnic),
            data=subset(USA, year == 2011), family=binomial(link = "logit"),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e5)))

Table1 <- stargazer(M3, M4, type="html", style="ajps",
                   title="Mixed Effects Models of American Preferences for a Strong Leader",
                   covariate.labels=c("Age", "Age-square", "Female", 
                                      "High School Education or Less", "Ideology",
                                      "Ideology-square", "Income Scale", "Republican",
                                      "Unemployed", "Autonomy Values", "Choice Values",
                                      "Equality Values", "Voice Values",
                                      "Republican*HS Ed. or Less"),
                  dep.var.labels.include = FALSE, multicolumn = FALSE,
                   model.names=FALSE, omit.stat=c("aic","ll","bic"), omit="Constant",
                  add.lines = list(c("Years", "1995-2011", "2011"))
)

M5 <- glmer(hddummy ~ z_sl + z_age + I(z_age^2) + female +
              hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
              + z_lemanc  +
              (1  | censusr) + (1 | year) + (1 | raceethnic) +
              (1  | censusr:year) + (1 | raceethnic:year),
            data=subset(USA), family=binomial(link = "logit"),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e5)))

M6 <- lmer(demimp ~ z_sl + z_age + I(z_age^2) + female +
              hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
              + z_lemanc  +
              (1  | censusr) + (1 | raceethnic),
            data=subset(USA))

M7 <- lmer(dempeoplelaws ~ z_sl + z_age + I(z_age^2) + female +
             hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
             + z_lemanc  +
             (1  | censusr) + (1 | raceethnic),
           data=subset(USA))

M8 <- lmer(demwomenrights ~ z_sl + z_age + I(z_age^2) + female +
             hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
             + z_lemanc  +
             (1  | censusr) + (1 | raceethnic),
           data=subset(USA))

M9 <- lmer(demarmycoup ~ z_sl + z_age + I(z_age^2) + female +
             hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
             + z_lemanc  +
             (1  | censusr) + (1 | raceethnic),
           data=subset(USA))

M10 <- lmer(demchooseleaders ~ z_sl + z_age + I(z_age^2) + female +
             hsedorless + z_ideo + I(z_ideo^2) + z_incscale + gop + unemployed +
             + z_lemanc  +
             (1  | censusr) + (1 | raceethnic),
           data=subset(USA))


Table2 <- stargazer(M5, M6, M7, M8, M9, M10, type="html", style="ajps",
                    title="Mixed Effects Models of American Attitudes to Democracy",
                    covariate.labels=c("Strong Leader", "Age", "Age-square", "Female", 
                                       "High School Education or Less", "Ideology",
                                       "Ideology-square", "Income Scale", "Republican",
                                       "Unemployed", "Emancipative Values"),
                    dep.var.labels.include = TRUE, multicolumn = FALSE,
                    model.names=FALSE, omit.stat=c("aic","ll","bic"), omit="Constant",
                    dep.var.labels=c("Opposition to Democracy", "Democracy is Important",
                                     "People Can Change the Laws",
                                     "Women Have Same Rights as Men",
                                     "Army Takes Over Incompetent Government",
                                     "People Choose Leaders in Elections",
                                     "Religious Authorities Interpret Laws")
)

write.table(USA,file="wvs-usa-sl.csv",sep=",",row.names=F,na="")
