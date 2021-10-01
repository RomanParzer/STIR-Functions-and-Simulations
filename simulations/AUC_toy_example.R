######################################################################
######################################################################
### AUC binary y
require(ROCR)
set.seed(1234)
x <- rnorm(100)
xnew <- rnorm(100)
y <- ifelse(-x+rnorm(100,0,0.5)>0,1,0)
ynew <- ifelse(-xnew+rnorm(100,0,0.5)>0,1,0)

fit_logit <- glm(y~.,family = binomial(link = "logit"),data = data.frame(X1=x,y=y))
fit_logit$coefficients[2]## negative
preds <- predict(fit_logit,data.frame(X1=xnew),type="link") ### default, same as:
preds2 <- fit_logit$coefficients[1] + fit_logit$coefficients[2]*xnew
sum((preds-preds2)^2)
probs <- predict(fit_logit,data.frame(X1=xnew),type="response")

### all give the same AUC value: (any strictly monotone (increasing) function of preds works)
performance(prediction(preds,ynew),measure="auc")@y.values[[1]] ## what I use
performance(prediction(probs,ynew),measure="auc")@y.values[[1]]
performance(prediction(5+3*preds,ynew),measure="auc")@y.values[[1]]
performance(prediction(1-exp(-preds),ynew),measure="auc")@y.values[[1]]
## does not work:
performance(prediction(xnew,ynew),measure="auc")@y.values[[1]]
## works:
performance(prediction(- xnew,ynew),measure="auc")@y.values[[1]]


######################################################################
######################################################################