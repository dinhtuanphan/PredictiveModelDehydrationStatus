# Predictive Model for Dehydration Status

In this project, a predictive model was built using the random forest algorithm in R to classify dehydration status based on physiological parameters. This is a research project on wearable sensors for physiological measurements that I participated in at Johns Hopkins Lab. We performed 300 clinical trials under different IRBs at Johns Hopkins Hospital.


## Prediction model for fluid loss during exertion

Weight change associated with fluid loss is an important metric for fatigue but is difficult to measure in real time. Since the sweat losses in this study were relatively small (<1.5% body weight), they would have no practical physiological impact. However, having demonstrated that the sweat chloride profiles are related to physiological parameters, we next created a prediction model to classify subjects according to their weight loss during exercise. For classification we set the weight loss cut-off to the median value: 0.74% bodyweight (group 1: ΔWgt <0.74%, group 2: ΔWgt >0.74%). We employed a random forest method using ΔC, ΔHR, ΔTcore, ΔRPE, onset time, exercise frequency, BMI, and gender as independent variables.

Two metrics commonly used to evaluate the performance of random forest models are the out-of-bag (OOB) error and the AUC (area under the ROC curve, 1: perfect classification, 0.5: no better than random classification). The OOB error, an estimation of the true prediction error, was 32% (group 1: 28% (7/25), group 2: 36% (9/25)) and the AUC value was 0.742. Interestingly, the sweat chloride concentration was the most effective variable contributing to the classification. Although a few previous studies have reported that hydration status can be classified or estimated by monitoring physiological parameters or by analyzing biofluids, these studies used post hoc laboratory-based analytical measurements. Here we show that chloride ion concentration and other physiological parameters recorded using wearable devices can be used to predict fluid loss due to sweating.

## Publication & Future Works
- The result of this research in our lab (Dr. Searson's Lab at JHU) that I work was published here: [Two Distinct Types of Sweat Profile in Healthy Subjects While Exercising at Constant Power Output Measured by a Wearable Sweat Sensor](https://www.nature.com/articles/s41598-019-54202-1#article-info)

- We are developing more comprehensive predictive models for healthy and CF subjects, based on unstructured datasets on physiological measurements, the GitHub repo here: [Machine Learning for Treatment Therapies](https://github.com/dinhtuanphan/MachineLearningForTherapies)

## File Structure
- MATLAB files used to first process raw data from multiple skin sensors (e.g., GSR, temp, heart rate, etc.)
- R files used to perform final data analysis and build machine learning models
