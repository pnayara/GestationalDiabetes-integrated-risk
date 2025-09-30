# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 16:14:47 2025

@author: Nayara
"""

#%% Pacotes

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, classification_report
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, recall_score
from statstests.process import stepwise
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import make_pipeline

#%% Upload de dados e exploração 

data = pd.read_excel('df_placenta.xlsx', sheet_name='data3')
data.info()

data_modif = data.copy()
data_modif['group'].unique()
mapping = {'Control' :0 , 'GDM' : 1}
data_modif['group_dummie'] = data_modif['group'].map(mapping)

data.groupby('group')['ogtt_1h'].describe().round(2).T
data.groupby('group')['ogtt_2h'].describe().round(2).T
data.groupby('group')['maternal_weight_adjusted'].describe().round(2).T
data.groupby('group')['BMI_adjust'].describe().round(2)

#%%  Criando classes de IMC

data_modif['BMI_adjust_cat'] = data_modif['BMI_adjust'].apply(
    lambda x : 'Normal' if x < 24 else 'Sobrepeso')

#%% Dummificar as categorias 

dummies = pd.get_dummies(data_modif['BMI_adjust_cat'], drop_first=True).astype('int')
data_category = pd.concat([data_modif, dummies], axis=1)

#data_category.to_csv('Dados_IMC_cat.csv', index = False)


#%% Corrigindo o nome das variaveis 

data_category.rename(columns ={'group_dummie': 'Grupos',
                              'BMI_adjust':'IMC',
                              'ogtt_1h':'oGTT_1h',
                              'ogtt_2h':'oGTT_2h',
                              'maternal_age': 'Idade',
                              'maternal_height':'Altura',
                              'maternal_weight_adjusted':'Peso'
                              }, inplace = True)

#%% Correlação e colinearidade entre variáveis 
#### Correlação

# Seleciona as variáveis e já renomeia
correlation_matrix = data_category[['oGTT_1h', 'oGTT_2h',
                                    'Idade','Altura',
                                    'Peso', 'IMC']].rename(
    columns={'oGTT_1h': 'Glicemia-1h',
             'oGTT_2h': 'Glicemia-2h'}
).corr()

# Plot
plt.figure(figsize=(8, 6), dpi =300)
sns.heatmap(correlation_matrix,
            annot=True,
            cmap='coolwarm',
            fmt=".2f",
            linewidths=0.5,
            square=True,
            cbar_kws={'shrink': 0.5, 'aspect': 15, 'pad': 0.02})
plt.yticks(rotation=0, fontsize=14 )
plt.xticks(rotation=45, fontsize=14)
plt.show()


####VIF - Variance Inflation Factor

X_VIF = data_category[['oGTT_1h', 'oGTT_2h',
                                    'Idade','Altura',
                                    'Peso', 'IMC']]
X_VIF_const = sm.add_constant(X_VIF)

VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
VIF["Tolerance"] = 1 / VIF["VIF"]

print(VIF)

#%% Modelos

## ogtt_1h, 'maternal_age', 'maternal_weight', 'maternal_height'

model_3 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Idade + Peso  + Altura',
                     data=data_category).fit()
#AIC / BIC
aic_m3= model_3.aic
bic_m3 = model_3.bic

print(model_3.summary())
print(f"Model 3 : AIC({aic_m3:.2f}) e BIC({bic_m3:.2f})")


#%% MODEL 4 

## ogtt_1h, 'maternal_age', 'BMI'

model_4 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Idade + IMC',
                     data=data_category).fit()
#AIC / BIC
aic_m4 = model_4.aic
bic_m4 = model_4.bic

print(model_4.summary())
print(f"Model 4 : AIC({aic_m4:.2f}) e BIC({bic_m4:.2f})")

#%% Model 7
## ogtt_1h, 'maternal_age', 'BMI_CAT'

model_7 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Idade + Sobrepeso',
                     data=data_category).fit()
#AIC / BIC
aic_m7 = model_7.aic
bic_m7 = model_7.bic

print(model_7.summary())
print(f"Model 7 : AIC({aic_m7:.2f}) e BIC({bic_m7:.2f})")

#%% Stepwise

# Estimação do modelo por meio do procedimento Stepwise
step_var_3 = stepwise(model_3, pvalue_limit=0.05)
step_var_4 = stepwise(model_4, pvalue_limit=0.05)
step_var_7 = stepwise(model_7, pvalue_limit=0.05)

#%% Modelos finais

#Modelo 1
model_1 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Peso',
                     data=data_category).fit()
#AIC / BIC
aic_m1= model_1.aic
bic_m1 = model_1.bic

print(model_1.summary())
print(f"Model 1 : AIC({aic_m1:.2f}) e BIC({bic_m1:.2f})")

#Modelo 2
model_2 = sm.Logit.from_formula('Grupos ~ oGTT_1h + IMC',
                     data=data_category).fit()
#AIC / BIC
aic_m2 = model_2.aic
bic_m2 = model_2.bic

print(model_2.summary())
print(f"Model 2 : AIC({aic_m2:.2f}) e BIC({bic_m2:.2f})")

#Modelo 3
model_3 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Sobrepeso',
                     data=data_category).fit()
#AIC / BIC
aic_m3 = model_3.aic
bic_m3 = model_3.bic

print(model_3.summary())
print(f"Model 3 : AIC({aic_m3:.2f}) e BIC({bic_m3:.2f})")

#%% Modelos de interação 
#### Modelo 4 : oGTT_1h * Peso

data_category["oGTT_1h_c"] = data_category["oGTT_1h"] - data_category["oGTT_1h"].mean()
data_category["Peso_c"] = data_category["Peso"] - data_category["Peso"].mean()
data_category['IMC_c'] = data_category['IMC'] - data_category['IMC'].mean()


modelo4 = sm.Logit.from_formula('Grupos ~ oGTT_1h_c + Peso_c + oGTT_1h_c:Peso_c',
                     data=data_category).fit()
print(modelo4.summary())

#### Modelo 5 : oGTT_1h * IMC

modelo5 = sm.Logit.from_formula('Grupos ~ oGTT_1h_c + IMC_c + oGTT_1h_c:IMC_c',
                     data=data_category).fit()
print(modelo5.summary())


modelo5 = sm.Logit.from_formula('Grupos ~ oGTT_1h_c + IMC_c',
                     data=data_category).fit()
print(modelo5.summary())

X_VIF = data_category[['oGTT_1h','IMC', 'oGTT_1h_c', 'IMC_c' ]]
X_VIF_const = sm.add_constant(X_VIF)
VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
print(VIF)

#%% VIF Modelos Finais 

### Modelo 1 

X_VIF = data_category[['oGTT_1h','Peso']]
X_VIF_const = sm.add_constant(X_VIF)

VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
VIF["Tolerance"] = 1 / VIF["VIF"]

print(VIF)


### Modelo 2

X_VIF = data_category[['oGTT_1h','IMC']]
X_VIF_const = sm.add_constant(X_VIF)

VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
VIF["Tolerance"] = 1 / VIF["VIF"]

print(VIF)

### Modelo 3

X_VIF = data_category[['oGTT_1h','Sobrepeso']]
X_VIF_const = sm.add_constant(X_VIF)

VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
VIF["Tolerance"] = 1 / VIF["VIF"]

print(VIF)

#%% ODDS Ratio

def calcular_or(model):
    params = model.params
    conf = model.conf_int()
    conf.columns = ['2.5%', '97.5%']

    or_values = np.exp(params)
    or_conf = np.exp(conf)
    
    results = pd.DataFrame({
        'Coeficiente': params,
        'Razão de Chaces': or_values,
        'IC 2.5%': or_conf['2.5%'],
        'IC 97.5%': or_conf['97.5%'],
        'p-value': model.pvalues})
    
    return results


or_m1 = calcular_or(model_1).round(3)
or_m2 = calcular_or(model_2).round(3)
or_m3 = calcular_or(model_3).round(3)

print("\n=== Modelo 1 ===")
print(or_m1.T)
print("\n=== Modelo 2 ===")
print(or_m2.T)
print("\n=== Modelo 3 ===")
print(or_m3.T)

#%% Analise de Previsão - Modelos 1 

#### MODELO 1
X_m1 = data_category[['oGTT_1h', 'Peso']]
y_m1 = data_category['Grupos']

#Traning and test data 70/30 (28/13)
X_train, X_test, y_train, y_test = train_test_split(
   X_m1, y_m1, test_size=0.3, random_state=42)

#Model adjust
X_train_const = sm.add_constant(X_train)
X_test_const = sm.add_constant(X_test)
model = sm.Logit(y_train, X_train_const).fit()

# Prediction
y_pred_proba = model.predict(X_test_const)
y_pred = (y_pred_proba > 0.5).astype(int)

# Matriz de confusão 
cm_m1 = confusion_matrix(y_test, y_pred)
print(cm_m1)

#Indicadores de Previsão 
accuracy_m1 = accuracy_score(y_test, y_pred)
sensitivity_m1 = cm_m1[1,1] / (cm_m1[1,1] + cm_m1[1,0])
specificity_m1 = cm_m1[0,0] / (cm_m1[0,0] + cm_m1[0,1])
auc_m1 = roc_auc_score(y_test, y_pred_proba)

print(model_1.summary())
print(f"Model 1: AIC({aic_m1:.2f}) e BIC({bic_m1:.2f})")
print(f"Model 1: Acurácia ({accuracy_m1:.2f})")
print(f"Model 1: Sensibilidade({sensitivity_m1:.2f})")
print(f"Model 1: Especificidade({specificity_m1:.2f})")
print(f"Model 1: AUC ROC({auc_m1:.2f})")


#### Cross Validation Modelo 1 

# Pipeline e validação cruzada
pipeline = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

# Métricas a calcular
specificity = make_scorer(recall_score, pos_label=0)
metrics = {
    "ROC AUC": "roc_auc",
    "Acurácia": "accuracy",
    "Sensibilidade": "recall",
    "Especificidade": specificity
}

# Loop para calcular médias e desvios
results = {}
for name, scorer in metrics.items():
    scores = cross_val_score(pipeline, X_m1, y_m1, cv=cv, scoring=scorer)
    results[name] = (scores.mean(), scores.std())

# Exibir resultados
print("\nResultados do Cross-Validation:")
for name, (mean, std) in results.items():
    print(f"{name}: {mean:.3f} ± {std:.3f}")

#### Curva ROC com predições de CV
y_pred_cv = cross_val_predict(pipeline, X_m1, y_m1, cv=cv, method="predict_proba")[:, 1]
fpr_cv, tpr_cv, _ = roc_curve(y_m1, y_pred_cv)
roc_auc_cv = auc(fpr_cv, tpr_cv)

# Plot da curva ROC
plt.figure(figsize=(6, 6))
plt.plot(fpr_cv, tpr_cv, color="green", lw=2, label=f'ROC CV (AUC = {roc_auc_cv:.2f})')
plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("1 - Especificidade (FPR)", fontsize = 14, fontweight ='bold')
plt.ylabel("Sensibilidade (TPR)",fontsize = 14, fontweight ='bold')
plt.title("Curva ROC com Validação Cruzada", fontsize = 16, fontweight ='bold' )
plt.legend(loc="lower right")
plt.show()


#%% Analise de Previsão - Modelos 2

#MODEL 2
X_m2 = data_category[['oGTT_1h','IMC']]
y_m2 = data_category['Grupos']

#Traning and test data
X_train, X_test, y_train, y_test = train_test_split(
   X_m2, y_m2, test_size=0.3, random_state=42)

#Model adjust
X_train_const = sm.add_constant(X_train)
X_test_const = sm.add_constant(X_test)
model = sm.Logit(y_train, X_train_const).fit()

# Prediction
y_pred_proba = model.predict(X_test_const)
y_pred = (y_pred_proba > 0.5).astype(int)

# Evaluation
cm_m2 = confusion_matrix(y_test, y_pred)
accuracy_m2 = accuracy_score(y_test, y_pred)
sensitivity_m2 = cm_m2[1,1] / (cm_m2[1,1] + cm_m2[1,0])
specificity_m2 = cm_m2[0,0] / (cm_m2[0,0] + cm_m2[0,1])
auc_m2 = roc_auc_score(y_test, y_pred_proba)

print(model_2.summary())
print(f"Model 2: AIC({aic_m2:.2f}) e BIC({bic_m2:.2f})")
print(f"Model 2: Acurácia ({accuracy_m2:.2f})")
print(f"Model 2: Sensibilidade({sensitivity_m2:.2f})")
print(f"Model 2: Especificidade({specificity_m2:.2f})")
print(f"Model 2: AUC ROC({auc_m2:.2f})")

#### Cross Validation do Modelo 2 

# Pipeline e validação cruzada
pipeline = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

# Métricas a calcular
specificity = make_scorer(recall_score, pos_label=0)
metrics = {
    "ROC AUC": "roc_auc",
    "Acurácia": "accuracy",
    "Sensibilidade": "recall",
    "Especificidade": specificity
}

# Loop para calcular médias e desvios
results = {}
for name, scorer in metrics.items():
    scores = cross_val_score(pipeline, X_m2, y_m2, cv=cv, scoring=scorer)
    results[name] = (scores.mean(), scores.std())

# Exibir resultados
print("\nResultados do Cross-Validation:")
for name, (mean, std) in results.items():
    print(f"{name}: {mean:.3f} ± {std:.3f}")

#### Curva ROC com predições de CV
y_pred_cv = cross_val_predict(pipeline, X_m2, y_m2, cv=cv, method="predict_proba")[:, 1]
fpr_cv, tpr_cv, _ = roc_curve(y_m2, y_pred_cv)
roc_auc_cv = auc(fpr_cv, tpr_cv)

# Plot da curva ROC
plt.figure(figsize=(6, 6))
plt.plot(fpr_cv, tpr_cv, color="green", lw=2, label=f'ROC CV (AUC = {roc_auc_cv:.2f})')
plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("1 - Especificidade (FPR)", fontsize = 14, fontweight ='bold')
plt.ylabel("Sensibilidade (TPR)",fontsize = 14, fontweight ='bold')
plt.title("Curva ROC com Validação Cruzada", fontsize = 16, fontweight ='bold' )
plt.legend(loc="lower right")
plt.show()

#%% Dois Modelos na mesma curva ROC 

pipeline = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

# Função para calcular FPR, TPR e AUC de um modelo/dataset
def get_roc_auc(X, y):
    y_pred_cv = cross_val_predict(pipeline, X, y, cv=cv, method="predict_proba")[:, 1]
    fpr, tpr, _ = roc_curve(y, y_pred_cv)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc

# Modelo 1
fpr1, tpr1, auc1 = get_roc_auc(X_m1, y_m1)

# Modelo 2
fpr2, tpr2, auc2 = get_roc_auc(X_m2, y_m2)

# Plot único
plt.figure(figsize=(6, 6))
plt.plot(fpr1, tpr1, color="green", lw=2, label=f'Modelo 1 (AUC = {auc1:.2f})')
plt.plot(fpr2, tpr2, color="blue", lw=2, label=f'Modelo 2 (AUC = {auc2:.2f})')
plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("1 - Especificidade (FPR)", fontsize=14, fontweight='bold')
plt.ylabel("Sensibilidade (TPR)", fontsize=14, fontweight='bold')
plt.title("Curvas ROC com Validação Cruzada", fontsize=16, fontweight='bold')
plt.legend(loc="lower right")
plt.show()

