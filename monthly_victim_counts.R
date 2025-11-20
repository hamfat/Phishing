setwd("C:/Users/hamme/OneDrive - AUT University/Phishing_Study")


fraud.data<- read.csv("canada_fraud_edited.csv")
View(fraud.data)

fraud.data$date <- as.Date(fraud.data$Date_Received, "%m/%d/%Y")


#######Make a year column*********
fraud.data$year <- as.Date(fraud.data$date,'%m/%d/%Y')
fraud.data$year <- as.factor(format(fraud.data$year, '%Y'))


#######Make a month column*********
fraud.data$month<- as.Date(fraud.data$date, "%m/%d/%Y")
fraud.data$month <- as.factor(format(fraud.data$month, "%m"))


a<-table (fraud.data$year)
addmargins(a)

library(dplyr)


case_counts <- fraud.data %>%
  group_by(year, Number_of_Victims) %>%
  summarise(counts = n())
case_counts


case_counts_month <- fraud.data %>%
  group_by(month, Number_of_Victims) %>%
  summarise(counts = n())
case_counts_month
