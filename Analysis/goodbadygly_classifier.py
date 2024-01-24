# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 09:15:05 2024

@author: Administrator
"""

from sklearn.svm import SVC
import pandas as pd
import sklearn.model_selection as model_selection
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn import svm
from sklearn.metrics import confusion_matrix 
from sklearn.preprocessing import scale
from pycaret.classification import *
import sklearn.preprocessing
from sklearn.preprocessing import StandardScaler


data = pd.read_csv("K:/shares\\di04_limet_bioinformatics\\PhD Pablo\\Publicaties\\WIP\\TARDIS\\data\\metabease_metabolomics\\positive\\feature_table.csv")

X = data[["AUC" , "MaxInt", "SNR", "peak_cor"]]

y = data[["rating"]]

X_train, X_test, y_train, y_test = model_selection.train_test_split(X, y, train_size=0.80, test_size=0.20, random_state=101)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


encoder = sklearn.preprocessing.LabelEncoder()

encoder.fit(y_train)
Y_train = encoder.transform(y_train)

# encoding test labels 
encoder.fit(y_test)
Y_test = encoder.transform(y_test)

svm_model_linear = SVC(kernel = 'linear', C = 1)
svm_model_linear.fit(X_train, Y_train) 



svm_predictions = svm_model_linear.predict(X_test) 
  
# model accuracy for X_test   
accuracy = svm_model_linear.score(X_test, y_test) 
  
# creating a confusion matrix 
cm = confusion_matrix(y_test, svm_predictions) 


