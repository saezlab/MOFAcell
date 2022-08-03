facts <- MOFA2::get_factors(MOFAobject.trained) %>%
    pluck(1) %>%
    as.data.frame() %>%
    rownames_to_column("condition") %>%
    mutate(condition=as.factor(gsub(x = condition, pattern = "[|].*", replacement = "")))
facts

# Set random seed to make results reproducible:
set.seed(6)
# Calculate the size of each of the data sets:
data_set_size <- floor(nrow(facts)/2)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(facts), size = data_set_size)

# Assign the data to the correct sets
training <- facts[indexes,]
validation1 <- facts[-indexes,]

#import the package
library(randomForest)


# Perform training:
rf_classifier = randomForest(condition ~ .,
                             data=training,
                             ntree=100,
                             mtry=2,
                             importance=TRUE)


# Validation set assessment #2: ROC curves and AUC

# Needs to import ROCR package for ROC curve plotting:
library(ROCR)

# Calculate the probability of new observations belonging to each class
# prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
prediction_for_roc_curve <- stats::predict(rf_classifier,validation1[,-1],type="prob")

# Use pretty colours:
pretty_colours <- c("#F8766D","#00BA38","#619CFF")
# Specify the different classes
classes <- levels(validation1$condition)

# For each class
for (i in 1:3){
    # Define which observations belong to class[i]
    true_values <- ifelse(validation1[,1]==classes[i],1,0)
    # Assess the performance of classifier for class[i]
    pred <- prediction(prediction_for_roc_curve[,i],true_values)
    perf <- performance(pred, "tpr", "fpr")
    if (i==1)
    {
        plot(perf,main="ROC Curve",col=pretty_colours[i])
    }
    else
    {
        plot(perf,main="ROC Curve",col=pretty_colours[i],add=TRUE)
    }
    # Calculate the AUC and print it to screen
    auc.perf <- performance(pred, measure = "auc")
    print(auc.perf@y.values)
}

