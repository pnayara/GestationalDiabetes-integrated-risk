# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 11:09:51 2025

@author: Nayara
"""

#%% Pacotes 

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter


#%% Pacotes R 

ro.r('install.packages("hnp")')
base = importr('base')
r_version = base.R_Version()
hnp = importr('hnp')
stats = importr('stats')

print(f"Versão do R instalada: {r_version.rx2('version.string')[0]}")

#%% Configurações dos gráficos

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False 

#%% Dados clinicos 

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


data_modif['BMI_adjust_cat'] = data_modif['BMI_adjust'].apply(
    lambda x : 'Normal' if x < 24 else 'Sobrepeso')

##Dummificar as categorias 
dummies = pd.get_dummies(data_modif['BMI_adjust_cat'], drop_first=True).astype('int')
data_category = pd.concat([data_modif, dummies], axis=1)

#%% Conversor de df py2r

with localconverter(ro.default_converter + pandas2ri.converter):
    r_data = ro.conversion.py2rpy(data_category)
    ro.globalenv['dados'] = r_data

#%% Ajustar o modelo no R

ro.globalenv['dados'] = r_data
ro.r('''modelo1 <- glm(group_dummie ~ ogtt_1h + maternal_weight_adjusted,
              family = binomial(link = "logit"),
              data = dados) ''')


ro.r('''modelo2 <- glm(group_dummie ~ ogtt_1h + BMI_adjust,
              family = binomial(link = "logit"),
              data = dados) ''')


ro.r('''modelo3 <- glm(group_dummie ~ ogtt_1h + Sobrepeso,
              family = binomial(link = "logit"),
              data = dados) ''')


#%% Rodar o envelope simulado com hnp

#### MODELO 1
ro.r('''
hnp_resultado_modelo1 <- hnp(modelo1, resid.type = "deviance", sim = 1000)
plot(hnp_resultado_modelo1)
dev.off()
''')

ro.r('''
resultado1 <- hnp(modelo1, resid.type = "deviance", sim = 1000)
''')

#Extração
resultado1 = ro.r('resultado1')
observados = resultado1.rx2('x')
residuos_modelo1 = resultado1.rx2('residuals')
envelope_lower = resultado1.rx2('lower')
envelope_upper = resultado1.rx2('upper')
envelope_median = resultado1.rx2('median')

# Converte para pandas
with localconverter(default_converter + pandas2ri.converter):
    observados_py = ro.conversion.rpy2py(observados)
    residuos_modelo1_py = ro.conversion.rpy2py(residuos_modelo1)
    envelope_lower_py = ro.conversion.rpy2py(envelope_lower)
    envelope_upper_py = ro.conversion.rpy2py(envelope_upper)
    envelope_median_py = ro.conversion.rpy2py(envelope_median)

# Dataframe do Envelope
df_plot = pd.DataFrame({
    "observados": observados_py,
    "residuos": residuos_modelo1_py,
    "envelope_lower": envelope_lower_py,
    "envelope_upper": envelope_upper_py,
    "envelope_median": envelope_median_py
})

## PLOT
df_plot.sort_values(by="observados", inplace=True)
plt.figure(figsize=(8, 6), dpi=300)
plt.fill_between(df_plot["observados"],
                 df_plot["envelope_lower"],
                 df_plot["envelope_upper"],
                 color="lightgray", 
                 alpha=0.5,
                 label="Envelope 95%")

plt.plot(df_plot["observados"],
         df_plot["envelope_median"],
         color="red",
         linestyle="--",
         label="Mediana")
sns.scatterplot(x="observados",
                y="residuos",
                data=df_plot,
                color="blue",
                s=60,
                label="Resíduos")
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlabel("Valores teóricos",
           fontsize = 16)
plt.ylabel("Resíduos (deviance)",
           fontsize = 16)
plt.title("Envelope Simulado para Modelo 1",
          fontsize = 18,
          fontweight = 'bold')
sns.despine()
plt.legend(fontsize = 16)
plt.show()


#### MODELO 2
ro.r('''
hnp_resultado_modelo2 <- hnp(modelo2, resid.type = "deviance", sim = 1000)
plot(hnp_resultado_modelo2)
dev.off()
''')

ro.r('''
resultado2 <- hnp(modelo2, resid.type = "deviance", sim = 1000)
''')

#Extração
resultado2 = ro.r('resultado2')
observados = resultado2.rx2('x')
residuos_modelo2 = resultado2.rx2('residuals')
envelope_lower = resultado2.rx2('lower')
envelope_upper = resultado2.rx2('upper')
envelope_median = resultado2.rx2('median')

# Converte para pandas
with localconverter(default_converter + pandas2ri.converter):
    observados_py = ro.conversion.rpy2py(observados)
    residuos_modelo2_py = ro.conversion.rpy2py(residuos_modelo2)
    envelope_lower_py = ro.conversion.rpy2py(envelope_lower)
    envelope_upper_py = ro.conversion.rpy2py(envelope_upper)
    envelope_median_py = ro.conversion.rpy2py(envelope_median)

# Dataframe do Envelope
df_plot_2 = pd.DataFrame({
    "observados": observados_py,
    "residuos": residuos_modelo2_py,
    "envelope_lower": envelope_lower_py,
    "envelope_upper": envelope_upper_py,
    "envelope_median": envelope_median_py
})

## PLOT
df_plot_2.sort_values(by="observados", inplace=True)
plt.figure(figsize=(8, 6),dpi =300)
plt.fill_between(df_plot_2["observados"],
                 df_plot_2["envelope_lower"],
                 df_plot_2["envelope_upper"],
                 color="lightgray", 
                 alpha=0.5,
                 label="Envelope 95%")

plt.plot(df_plot_2["observados"],
         df_plot_2["envelope_median"],
         color="red",
         linestyle="--",
         label="Mediana")
sns.scatterplot(x="observados",
                y="residuos",
                data=df_plot_2,
                color="blue",
                s=60,
                label="Resíduos")
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlabel("Valores teóricos",
           fontsize = 16)
plt.ylabel("Resíduos (deviance)",
           fontsize = 16)
plt.title("Envelope Simulado para Modelo 2",
          fontsize = 18,
          fontweight = 'bold')
sns.despine()
plt.legend(fontsize = 16)
plt.show()



#### MODELO 3
ro.r('''
hnp_resultado_modelo3 <- hnp(modelo3, resid.type = "deviance", sim = 1000)
plot(hnp_resultado_modelo3)
dev.off()
''')

ro.r('''
resultado3 <- hnp(modelo3, resid.type = "deviance", sim = 1000)
''')

#Extração
resultado3 = ro.r('resultado3')
observados = resultado3.rx2('x')
residuos_modelo3 = resultado3.rx2('residuals')
envelope_lower = resultado3.rx2('lower')
envelope_upper = resultado3.rx2('upper')
envelope_median = resultado3.rx2('median')

# Converte para pandas
with localconverter(default_converter + pandas2ri.converter):
    observados_py = ro.conversion.rpy2py(observados)
    residuos_modelo3_py = ro.conversion.rpy2py(residuos_modelo3)
    envelope_lower_py = ro.conversion.rpy2py(envelope_lower)
    envelope_upper_py = ro.conversion.rpy2py(envelope_upper)
    envelope_median_py = ro.conversion.rpy2py(envelope_median)

# Dataframe do Envelope
df_plot_3 = pd.DataFrame({
    "observados": observados_py,
    "residuos": residuos_modelo3_py,
    "envelope_lower": envelope_lower_py,
    "envelope_upper": envelope_upper_py,
    "envelope_median": envelope_median_py
})

## PLOT
df_plot_3.sort_values(by="observados", inplace=True)
plt.figure(figsize=(8, 6))
plt.fill_between(df_plot_3["observados"],
                 df_plot_3["envelope_lower"],
                 df_plot_3["envelope_upper"],
                 color="lightgray", 
                 alpha=0.5,
                 label="Envelope 95%")

plt.plot(df_plot_3["observados"],
         df_plot_3["envelope_median"],
         color="red",
         linestyle="--",
         label="Mediana")
sns.scatterplot(x="observados",
                y="residuos",
                data=df_plot,
                color="blue",
                s=60,
                label="Resíduos")
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("Valores teóricos",
           fontsize = 14)
plt.ylabel("Resíduos (deviance)",
           fontsize = 14)
plt.title("Envelope Simulado para Modelo 3",
          fontsize = 16,
          fontweight = 'bold')
sns.despine()
plt.legend()
plt.show()

#%% Dados clinicos + Genes 

dados_genes = pd.read_excel('DadosClinicosGenes.xlsx')
dados_genes.info()

#### Conversor de df py2r
with localconverter(ro.default_converter + pandas2ri.converter):
    r_dados_genes = ro.conversion.py2rpy(dados_genes)


#%% Ajustar os modelo no R - Modelo 1 e Modelo 2 apenas 

ro.globalenv['dados_genes'] = r_dados_genes

ro.r('''modelo1.1 <- glm(Grupos  ~  oGTT_1h + Peso + CSH_group,
              family = binomial(link = "logit"),
              data = dados_genes) ''')

resumo = ro.r('summary(modelo1.1)')
print(resumo)


ro.r('''modelo2.1 <- glm(Grupos  ~  oGTT_1h + IMC + CSH_group,
              family = binomial(link = "logit"),
              data = dados_genes) ''')

resumo = ro.r('summary(modelo2.1)')
print(resumo)

#%% Rodar o envelope simulado com hnp - Modelos Com Genes 

#### MODELO 1
ro.r('''
hnp_resultado_modelo1.1 <- hnp(modelo1.1, resid.type = "deviance", sim = 1000)
plot(hnp_resultado_modelo1.1)
dev.off()
''')

ro.r('''
resultado1 <- hnp(modelo1.1, resid.type = "deviance", sim = 1000)
''')

#Extração
resultado1 = ro.r('resultado1')
observados = resultado1.rx2('x')
residuos_modelo1 = resultado1.rx2('residuals')
envelope_lower = resultado1.rx2('lower')
envelope_upper = resultado1.rx2('upper')
envelope_median = resultado1.rx2('median')

# Converte para pandas
with localconverter(default_converter + pandas2ri.converter):
    observados_py = ro.conversion.rpy2py(observados)
    residuos_modelo1_py = ro.conversion.rpy2py(residuos_modelo1)
    envelope_lower_py = ro.conversion.rpy2py(envelope_lower)
    envelope_upper_py = ro.conversion.rpy2py(envelope_upper)
    envelope_median_py = ro.conversion.rpy2py(envelope_median)

# Dataframe do Envelope
df_plot = pd.DataFrame({
    "observados": observados_py,
    "residuos": residuos_modelo1_py,
    "envelope_lower": envelope_lower_py,
    "envelope_upper": envelope_upper_py,
    "envelope_median": envelope_median_py
})

## PLOT
df_plot.sort_values(by="observados", inplace=True)
plt.figure(figsize=(8, 6))
plt.fill_between(df_plot["observados"],
                 df_plot["envelope_lower"],
                 df_plot["envelope_upper"],
                 color="lightgray", 
                 alpha=0.5,
                 label="Envelope 95%")

plt.plot(df_plot["observados"],
         df_plot["envelope_median"],
         color="red",
         linestyle="--",
         label="Mediana")
sns.scatterplot(x="observados",
                y="residuos",
                data=df_plot,
                color="blue",
                s=60,
                label="Resíduos")
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("Valores teóricos",
           fontsize = 14)
plt.ylabel("Resíduos (deviance)",
           fontsize = 14)
plt.title("Envelope Simulado para Modelo 1",
          fontsize = 16,
          fontweight = 'bold')
sns.despine()
plt.legend()
plt.show()


#### MODELO 2
ro.r('''
hnp_resultado_modelo2.1 <- hnp(modelo2.1, resid.type = "deviance", sim = 1000)
plot(hnp_resultado_modelo2.1)
dev.off()
''')

ro.r('''
resultado2 <- hnp(modelo2.1, resid.type = "deviance", sim = 1000)
''')

#Extração
resultado2 = ro.r('resultado2')
observados = resultado2.rx2('x')
residuos_modelo2 = resultado2.rx2('residuals')
envelope_lower = resultado2.rx2('lower')
envelope_upper = resultado2.rx2('upper')
envelope_median = resultado2.rx2('median')

# Converte para pandas
with localconverter(default_converter + pandas2ri.converter):
    observados_py = ro.conversion.rpy2py(observados)
    residuos_modelo2_py = ro.conversion.rpy2py(residuos_modelo2)
    envelope_lower_py = ro.conversion.rpy2py(envelope_lower)
    envelope_upper_py = ro.conversion.rpy2py(envelope_upper)
    envelope_median_py = ro.conversion.rpy2py(envelope_median)

# Dataframe do Envelope
df_plot_2 = pd.DataFrame({
    "observados": observados_py,
    "residuos": residuos_modelo2_py,
    "envelope_lower": envelope_lower_py,
    "envelope_upper": envelope_upper_py,
    "envelope_median": envelope_median_py
})

## PLOT
df_plot_2.sort_values(by="observados", inplace=True)
plt.figure(figsize=(8, 6), dpi =300)
plt.fill_between(df_plot_2["observados"],
                 df_plot_2["envelope_lower"],
                 df_plot_2["envelope_upper"],
                 color="lightgray", 
                 alpha=0.5,
                 label="Envelope 95%")

plt.plot(df_plot_2["observados"],
         df_plot_2["envelope_median"],
         color="red",
         linestyle="--",
         label="Mediana")
sns.scatterplot(x="observados",
                y="residuos",
                data=df_plot_2,
                color="blue",
                s=60,
                label="Resíduos")
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlabel("Valores teóricos",
           fontsize = 16)
plt.ylabel("Resíduos (deviance)",
           fontsize = 16)
plt.title("Envelope Simulado para Modelo Clínico-CSH",
          fontsize = 18,
          fontweight = 'bold')
sns.despine()
plt.legend(fontsize = 16)
plt.show()